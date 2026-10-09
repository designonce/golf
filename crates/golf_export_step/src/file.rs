//! A STEP file of parts, one per body.

use std::collections::HashMap;

use crate::error::StepError;
use crate::geometry::StepCurve;
use crate::geometry::StepSurface;
use crate::source::StepSource;
use crate::writer::Ref;
use crate::writer::Writer;
use crate::writer::logical;
use crate::writer::real;
use crate::writer::refs;
use crate::writer::string;

/// File-wide settings.
#[derive(Clone, Debug, PartialEq)]
pub struct StepOptions {
    /// The file's name, in its header.
    pub name: String,
    /// When the file was written, as ISO 8601 text, in its header. Empty if not
    /// given.
    pub time_stamp: String,
    /// How far apart two points may be and still be the same point, in
    /// millimetres.
    pub uncertainty: f64,
}

impl Default for StepOptions {
    fn default() -> Self {
        Self {
            name: String::new(),
            time_stamp: String::new(),
            uncertainty: 1e-6,
        }
    }
}

/// A STEP file being built: add bodies with [`Self::add_body`], then
/// [`Self::finish`] it.
pub struct StepFile {
    w: Writer,
    options: StepOptions,
    context: Ref,
    product_context: Ref,
    definition_context: Ref,
    products: Vec<Ref>,
}

impl StepFile {
    pub fn new(options: StepOptions) -> Result<Self, StepError> {
        let mut w = Writer::default();
        let context = representation_context(&mut w, options.uncertainty)?;
        let application = w.add("APPLICATION_CONTEXT('managed model based 3d engineering')");
        w.add(format!(
            "APPLICATION_PROTOCOL_DEFINITION('international standard','ap242_managed_model_based_3d_engineering',2014,{application})"
        ));
        let product_context = w.add(format!("PRODUCT_CONTEXT('',{application},'mechanical')"));
        let definition_context = w.add(format!(
            "PRODUCT_DEFINITION_CONTEXT('part definition',{application},'design')"
        ));
        Ok(Self {
            w,
            options,
            context,
            product_context,
            definition_context,
            products: Vec::new(),
        })
    }

    /// Adds `body` as a part named `name`.
    ///
    /// A body whose shells are all closed is written as a solid, its first shell
    /// outside and the rest voids; otherwise as a surface model.
    pub fn add_body<M: StepSource>(&mut self, name: &str, body: &M) -> Result<(), StepError> {
        let w = &mut self.w;
        let shape = write_body(w, name, body)?;
        let origin = w.placement(&golf_geom::Placement::<golf_manifold::World>::at(
            golf_manifold::Point::new(nalgebra::Vector3::zeros()),
        ))?;
        let kind = match shape.closed {
            true => "ADVANCED_BREP_SHAPE_REPRESENTATION",
            false => "MANIFOLD_SURFACE_SHAPE_REPRESENTATION",
        };
        let representation = w.add(format!(
            "{kind}({},({},{origin}),{})",
            string(name),
            shape.item,
            self.context
        ));

        let name = string(name);
        let product = w.add(format!(
            "PRODUCT({name},{name},'',({}))",
            self.product_context
        ));
        let formation = w.add(format!("PRODUCT_DEFINITION_FORMATION('','',{product})"));
        let definition = w.add(format!(
            "PRODUCT_DEFINITION('design','',{formation},{})",
            self.definition_context
        ));
        let definition_shape = w.add(format!("PRODUCT_DEFINITION_SHAPE('','',{definition})"));
        w.add(format!(
            "SHAPE_DEFINITION_REPRESENTATION({definition_shape},{representation})"
        ));
        self.products.push(product);
        Ok(())
    }

    /// The file's text.
    pub fn finish(mut self) -> String {
        if !self.products.is_empty() {
            self.w.add(format!(
                "PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,{})",
                refs(self.products)
            ));
        }
        self.w.finish(&self.options.name, &self.options.time_stamp)
    }
}

/// Millimetres, radians and steradians, to `uncertainty`.
fn representation_context(w: &mut Writer, uncertainty: f64) -> Result<Ref, StepError> {
    let length = w.add("(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.))");
    let angle = w.add("(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.))");
    let solid_angle = w.add("(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT())");
    let uncertainty = w.add(format!(
        "UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE({}),{length},'distance_accuracy_value','')",
        real(uncertainty)?
    ));
    Ok(w.add(format!(
        "(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT(({uncertainty}))\
         GLOBAL_UNIT_ASSIGNED_CONTEXT(({length},{angle},{solid_angle}))REPRESENTATION_CONTEXT('',''))"
    )))
}

/// A body's shape: the solid or surface model, and whether it's closed.
struct Shape {
    item: Ref,
    closed: bool,
}

fn write_body<M: StepSource>(w: &mut Writer, name: &str, body: &M) -> Result<Shape, StepError> {
    let points: HashMap<M::Vertex, _> = body.vertices().collect();
    let mut vertices: HashMap<M::Vertex, Ref> = HashMap::new();
    let mut edges: HashMap<M::Edge, Ref> = HashMap::new();
    for (id, edge) in body.edges() {
        let mut vertex = |v: M::Vertex, w: &mut Writer| -> Result<Ref, StepError> {
            if let Some(&r) = vertices.get(&v) {
                return Ok(r);
            }
            let point = w.point(points[&v].coords)?;
            let r = w.add(format!("VERTEX_POINT('',{point})"));
            vertices.insert(v, r);
            Ok(r)
        };
        let (start, end) = (vertex(edge.start, w)?, vertex(edge.end, w)?);
        let curve = edge.curve.write(w)?;
        edges.insert(
            id,
            w.add(format!("EDGE_CURVE('',{start},{end},{curve},.T.)")),
        );
    }

    // Faces are written per shell: a void's faces point into the void, but a
    // STEP void is a closed shell pointing out of it, used reversed, so they're
    // written flipped.
    let faces: HashMap<M::Face, _> = body.faces().collect();
    let shells: Vec<Vec<M::Face>> = body.shells().collect();
    if shells.iter().all(Vec::is_empty) {
        return Err(StepError::EmptyBody(name.to_string()));
    }
    let mut closed = true;
    let mut written_shells = Vec::new();
    for (index, shell) in shells.iter().enumerate() {
        let flip = index > 0;
        let mut uses: HashMap<M::Edge, usize> = HashMap::new();
        let mut face_refs = Vec::new();
        for face_id in shell {
            let face = &faces[face_id];
            if face.loops.is_empty() {
                return Err(StepError::UnboundedFace {
                    face: format!("{face_id:?}"),
                });
            }
            let surface = face.surface.write(w)?;
            let mut bounds = Vec::new();
            for l in &face.loops {
                let oriented: Vec<Ref> = l
                    .iter()
                    .map(|&(edge, reversed)| {
                        *uses.entry(edge).or_default() += 1;
                        w.add(format!(
                            "ORIENTED_EDGE('',*,*,{},{})",
                            edges[&edge],
                            logical(!reversed)
                        ))
                    })
                    .collect();
                let edge_loop = w.add(format!("EDGE_LOOP('',{})", refs(oriented)));
                bounds.push(w.add(format!("FACE_BOUND('',{edge_loop},{})", logical(!flip))));
            }
            let same_sense = face.same_sense != flip;
            face_refs.push(w.add(format!(
                "ADVANCED_FACE('',{},{surface},{})",
                refs(bounds),
                logical(same_sense)
            )));
        }
        closed &= uses.values().all(|&n| n == 2);
        written_shells.push(face_refs);
    }

    let name = string(name);
    if !closed {
        let shells: Vec<Ref> = written_shells
            .into_iter()
            .map(|faces| w.add(format!("OPEN_SHELL('',{})", refs(faces))))
            .collect();
        let item = w.add(format!(
            "SHELL_BASED_SURFACE_MODEL({name},{})",
            refs(shells)
        ));
        return Ok(Shape { item, closed });
    }
    let shells: Vec<Ref> = written_shells
        .into_iter()
        .map(|faces| w.add(format!("CLOSED_SHELL('',{})", refs(faces))))
        .collect();
    let outer = shells[0];
    let voids: Vec<Ref> = shells[1..]
        .iter()
        .map(|shell| w.add(format!("ORIENTED_CLOSED_SHELL('',*,{shell},.F.)")))
        .collect();
    let item = match voids.is_empty() {
        true => w.add(format!("MANIFOLD_SOLID_BREP({name},{outer})")),
        false => w.add(format!("BREP_WITH_VOIDS({name},{outer},{})", refs(voids))),
    };
    Ok(Shape { item, closed })
}
