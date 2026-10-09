//! A STEP file of parts, one per body.

use std::collections::HashMap;

use golf_color::Color;

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

/// A part or assembly written to a [`StepFile`], for placing in assemblies.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Product {
    definition: Ref,
    representation: Ref,
    /// The identity placement among the representation's items, which
    /// placements in assemblies map from.
    origin: Ref,
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
    /// One presentation style per colour.
    styles: HashMap<Color, Ref>,
    /// Every styled item, gathered into one presentation at the end.
    styled_items: Vec<Ref>,
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
            styles: HashMap::new(),
            styled_items: Vec::new(),
        })
    }

    /// Adds `body` as a part named `name`, returning the product for placing it
    /// in assemblies.
    ///
    /// A body whose shells are all closed is written as a solid, its first shell
    /// outside and the rest voids; otherwise as a surface model.
    pub fn add_body<M: StepSource>(&mut self, name: &str, body: &M) -> Result<Product, StepError> {
        let shape = write_body(&mut self.w, name, body)?;
        self.style_body(&shape, body.color())?;
        let origin = self.w.motion(&nalgebra::Isometry3::identity())?;
        let kind = match shape.closed {
            true => "ADVANCED_BREP_SHAPE_REPRESENTATION",
            false => "MANIFOLD_SURFACE_SHAPE_REPRESENTATION",
        };
        let representation = self.w.add(format!(
            "{kind}({},({},{origin}),{})",
            string(name),
            shape.item,
            self.context
        ));
        Ok(self.product(name, representation, origin))
    }

    /// Adds an assembly product named `name` placing each child product (a
    /// part or another assembly) by its motion from the assembly's origin. A
    /// product may be placed many times, in one assembly or several.
    pub fn add_group(
        &mut self,
        name: &str,
        children: &[(Product, nalgebra::Isometry3<f64>)],
    ) -> Result<Product, StepError> {
        let w = &mut self.w;
        let origin = w.motion(&nalgebra::Isometry3::identity())?;
        let axes = children
            .iter()
            .map(|(_, motion)| w.motion(motion))
            .collect::<Result<Vec<_>, _>>()?;
        let mut items = vec![origin];
        items.extend(&axes);
        let representation = w.add(format!(
            "SHAPE_REPRESENTATION({},{},{})",
            string(name),
            refs(items),
            self.context
        ));
        let assembly = self.product(name, representation, origin);
        // Each child placed: a usage of its product in this one, and its
        // representation related to this one's by the transformation taking its
        // origin to the placement.
        for (i, ((child, _), axis)) in children.iter().zip(axes).enumerate() {
            let w = &mut self.w;
            let usage = w.add(format!(
                "NEXT_ASSEMBLY_USAGE_OCCURRENCE({},'','',{},{},$)",
                string(&(i + 1).to_string()),
                assembly.definition,
                child.definition
            ));
            let usage_shape = w.add(format!("PRODUCT_DEFINITION_SHAPE('','',{usage})"));
            let transformation = w.add(format!(
                "ITEM_DEFINED_TRANSFORMATION('','',{},{axis})",
                child.origin
            ));
            let relationship = w.add(format!(
                "(REPRESENTATION_RELATIONSHIP('','',{},{})REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION({transformation})SHAPE_REPRESENTATION_RELATIONSHIP())",
                child.representation, assembly.representation
            ));
            w.add(format!(
                "CONTEXT_DEPENDENT_SHAPE_REPRESENTATION({relationship},{usage_shape})"
            ));
        }
        Ok(assembly)
    }

    /// Styles a written body: its colour on the whole, and faces' own colours
    /// over it.
    fn style_body(&mut self, shape: &Shape, color: Option<Color>) -> Result<(), StepError> {
        let body_item = match color {
            Some(color) => {
                let style = self.style(color)?;
                let item = self
                    .w
                    .add(format!("STYLED_ITEM('color',({style}),{})", shape.styled));
                self.styled_items.push(item);
                Some(item)
            }
            None => None,
        };
        for &(color, face) in &shape.face_colors {
            let style = self.style(color)?;
            let item = match body_item {
                Some(body_item) => self.w.add(format!(
                    "OVER_RIDING_STYLED_ITEM('color',({style}),{face},{body_item})"
                )),
                None => self.w.add(format!("STYLED_ITEM('color',({style}),{face})")),
            };
            self.styled_items.push(item);
        }
        Ok(())
    }

    /// A `PRESENTATION_STYLE_ASSIGNMENT` filling both sides of a surface with
    /// `color`, made once per colour.
    fn style(&mut self, color: Color) -> Result<Ref, StepError> {
        if let Some(&style) = self.styles.get(&color) {
            return Ok(style);
        }
        let w = &mut self.w;
        let [r, g, b] = color.unit_rgb();
        let rgb = w.add(format!(
            "COLOUR_RGB('',{},{},{})",
            real(r)?,
            real(g)?,
            real(b)?
        ));
        let fill_colour = w.add(format!("FILL_AREA_STYLE_COLOUR('',{rgb})"));
        let fill = w.add(format!("FILL_AREA_STYLE('',({fill_colour}))"));
        let mut sides = vec![w.add(format!("SURFACE_STYLE_FILL_AREA({fill})"))];
        if color.a.is_some() {
            let transparency = w.add(format!(
                "SURFACE_STYLE_TRANSPARENT({})",
                real(1.0 - color.opacity())?
            ));
            sides.push(w.add(format!(
                "SURFACE_STYLE_RENDERING_WITH_PROPERTIES(.NORMAL_SHADING.,{rgb},({transparency}))"
            )));
        }
        let side = w.add(format!("SURFACE_SIDE_STYLE('',{})", refs(sides)));
        let usage = w.add(format!("SURFACE_STYLE_USAGE(.BOTH.,{side})"));
        let style = w.add(format!("PRESENTATION_STYLE_ASSIGNMENT(({usage}))"));
        self.styles.insert(color, style);
        Ok(style)
    }

    /// A product named `name` whose shape is `representation`, with the
    /// identity placement `origin` among its items.
    fn product(&mut self, name: &str, representation: Ref, origin: Ref) -> Product {
        let w = &mut self.w;
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
        Product {
            definition,
            representation,
            origin,
        }
    }

    /// The file's text.
    pub fn finish(mut self) -> String {
        if !self.styled_items.is_empty() {
            self.w.add(format!(
                "MECHANICAL_DESIGN_GEOMETRIC_PRESENTATION_REPRESENTATION('',{},{})",
                refs(core::mem::take(&mut self.styled_items)),
                self.context
            ));
        }
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
    /// What the body's colour is styled on: the solid, or the first open shell.
    styled: Ref,
    /// Faces with their own colour, and their entities.
    face_colors: Vec<(Color, Ref)>,
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
    let mut face_colors = Vec::new();
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
            let face_ref = w.add(format!(
                "ADVANCED_FACE('',{},{surface},{})",
                refs(bounds),
                logical(same_sense)
            ));
            if let Some(color) = body.face_color(*face_id) {
                face_colors.push((color, face_ref));
            }
            face_refs.push(face_ref);
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
            refs(shells.iter().copied())
        ));
        return Ok(Shape {
            item,
            closed,
            styled: shells[0],
            face_colors,
        });
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
    Ok(Shape {
        item,
        closed,
        styled: item,
        face_colors,
    })
}
