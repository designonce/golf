//! Shells of faces into bodies.

use std::collections::HashMap;

use golf_brep::Body;
use golf_brep::Coedge;
use golf_brep::EdgeId;
use golf_brep::FaceId;
use golf_brep::FaceUv;
use golf_brep::Loop;
use golf_brep::VertexId;
use golf_color::Color;
use golf_geom::AnyCurve;
use golf_geom::AnyCurve2;
use golf_geom::Line;
use golf_geom::NurbsCurve;
use golf_manifold::Embedding;
use golf_manifold::Mapping;
use golf_manifold::Point;
use golf_manifold::Vector;
use golf_manifold::World;
use golf_step::Id;
use nalgebra::Vector2;
use nalgebra::Vector3;

use crate::error::At;
use crate::error::ImportWarning;
use crate::error::Read;
use crate::error::invalid;
use crate::error::unsupported;
use crate::geometry::Curve2;
use crate::geometry::Geometry;

/// Builds one body from shells, keeping what it can: a face that can't be
/// read is left out (the body is then open there) with a warning.
pub(crate) struct BodyBuilder<'a> {
    geometry: Geometry<'a>,
    colors: &'a HashMap<Id, Color>,
    body: Body<World>,
    vertices: HashMap<Id, VertexId>,
    /// Each edge read; `None` for one with no extent.
    edges: HashMap<Id, Option<ReadEdge>>,
    pub warnings: Vec<ImportWarning>,
}

impl<'a> BodyBuilder<'a> {
    pub fn new(geometry: Geometry<'a>, colors: &'a HashMap<Id, Color>) -> Self {
        Self {
            geometry,
            colors,
            body: Body::new(),
            vertices: HashMap::new(),
            edges: HashMap::new(),
            warnings: Vec::new(),
        }
    }

    pub fn finish(mut self, color_of: &[Id]) -> (Body<World>, Vec<ImportWarning>) {
        if let Some(color) = color_of.iter().find_map(|id| self.colors.get(id)) {
            self.body.set_color(Some(*color));
        }
        (self.body, self.warnings)
    }

    /// Adds a shell: `CLOSED_SHELL`, `OPEN_SHELL`, or an oriented one.
    pub fn add_shell(&mut self, id: Id) -> Read<()> {
        let data = self.geometry.data;
        let entity = data.entity(id).at(id)?;
        let record = entity
            .simple()
            .ok_or_else(|| invalid(id, "a complex shell"))?;
        let (faces, flip) = match record.name.as_str() {
            "CLOSED_SHELL" | "OPEN_SHELL" | "CONNECTED_FACE_SET" => {
                (record.references(1).at(id)?, false)
            }
            "ORIENTED_CLOSED_SHELL" | "ORIENTED_OPEN_SHELL" => {
                let element = record.reference(2).at(id)?;
                let faces = data.entity(element).at(element)?.records[0]
                    .references(1)
                    .at(element)?;
                (faces, !record.boolean(3).at(id)?)
            }
            other => return Err(unsupported(id, format!("the shell {other}"))),
        };
        let mut added = Vec::new();
        for face in faces {
            match self.add_face(face, flip) {
                Ok(face) => added.push(face),
                Err(warning) => self.warnings.push(warning),
            }
        }
        self.body
            .add_shell(added)
            .map_err(|e| invalid(id, e.to_string()))?;
        Ok(())
    }

    fn add_face(&mut self, id: Id, flip: bool) -> Read<FaceId> {
        let data = self.geometry.data;
        let entity = data.entity(id).at(id)?;
        let record = entity
            .record("ADVANCED_FACE")
            .or_else(|| entity.record("FACE_SURFACE"))
            .ok_or_else(|| unsupported(id, format!("the face {}", entity.name())))?;
        let bounds = record.references(1).at(id)?;
        let surface_id = record.reference(2).at(id)?;
        let same_sense = record.boolean(3).at(id)? != flip;

        // Outer bound first; each bound an EDGE_LOOP of ORIENTED_EDGEs.
        let mut bound_list = Vec::new();
        for bound in bounds {
            let b = data.entity(bound).at(bound)?;
            let r = b
                .records
                .iter()
                .find(|r| r.name.contains("BOUND"))
                .ok_or_else(|| unsupported(bound, format!("the bound {}", b.name())))?;
            let outer = r.name == "FACE_OUTER_BOUND";
            let loop_id = r.reference(1).at(bound)?;
            let orientation = r.boolean(2).at(bound)? != flip;
            bound_list.push((outer, loop_id, orientation));
        }
        bound_list.sort_by_key(|&(outer, ..)| !outer);

        // The vertices first, to size swept surfaces.
        let mut uses = Vec::new();
        for &(_, loop_id, orientation) in &bound_list {
            let entity = data.entity(loop_id).at(loop_id)?;
            if entity.is("VERTEX_LOOP") {
                // A loop of one vertex (a cone's apex): golf leaves poles out.
                continue;
            }
            let record = entity.expect("EDGE_LOOP").at(loop_id)?;
            let mut oriented = Vec::new();
            for oe in record.references(1).at(loop_id)? {
                let r = data.record(oe, "ORIENTED_EDGE").at(oe)?;
                oriented.push((r.reference(3).at(oe)?, r.boolean(4).at(oe)?));
            }
            if !orientation {
                oriented.reverse();
                for o in &mut oriented {
                    o.1 = !o.1;
                }
            }
            uses.push((loop_id, oriented));
        }
        let mut points = Vec::new();
        for (_, oriented) in &uses {
            for &(edge, _) in oriented {
                let r = data.record(edge, "EDGE_CURVE").at(edge)?;
                for v in [r.reference(1).at(edge)?, r.reference(2).at(edge)?] {
                    points.push(self.vertex_point(v)?);
                }
            }
        }
        let surface = self.geometry.surface(surface_id, &points)?;

        let mut loops = Vec::new();
        for (loop_id, oriented) in uses {
            let mut coedges = Vec::new();
            let mut candidates = Vec::new();
            for (edge, orientation) in oriented {
                let Some(ReadEdge {
                    id: edge_id,
                    against,
                    pcurves: curve,
                }) = self.edge(edge)?
                else {
                    // Degenerate: no extent to bound anything.
                    continue;
                };
                // The coedge runs golf's edge forwards when it runs the file's
                // edge the way golf's does.
                let reversed = orientation == against;
                let pcurves: Vec<AnyCurve2<FaceUv<World>>> = match surface.uv_scale {
                    Some(scale) => curve
                        .iter()
                        .filter(|&&(on, _)| on == surface_id)
                        .filter_map(|&(_, c)| {
                            self.geometry
                                .curve_2d(c)
                                .ok()
                                .and_then(|c| scaled(c, scale))
                        })
                        .collect(),
                    None => Vec::new(),
                };
                coedges.push(Coedge::new(edge_id, reversed));
                candidates.push(pcurves);
            }
            if coedges.is_empty() {
                continue;
            }
            choose_pcurves(&self.body, &mut coedges, candidates);
            let _ = loop_id;
            loops.push(Loop::new(coedges));
        }
        let face = self
            .body
            .add_face(surface.surface, same_sense, loops)
            .map_err(|e| invalid(id, e.to_string()))?;
        if let Some(color) = self.colors.get(&id) {
            self.body.set_face_color(face, Some(*color));
        }
        Ok(face)
    }

    fn vertex_point(&self, id: Id) -> Read<Vector3<f64>> {
        let r = self.geometry.data.record(id, "VERTEX_POINT").at(id)?;
        self.geometry.point::<3>(r.reference(1).at(id)?)
    }

    fn vertex(&mut self, id: Id) -> Read<VertexId> {
        if let Some(&v) = self.vertices.get(&id) {
            return Ok(v);
        }
        let point = self.vertex_point(id)?;
        let v = self
            .body
            .add_vertex(Point::new(point))
            .map_err(|e| invalid(id, e.to_string()))?;
        self.vertices.insert(id, v);
        Ok(v)
    }

    /// The golf edge for an `EDGE_CURVE`, made on first use; `None` for one
    /// with no extent.
    fn edge(&mut self, id: Id) -> Read<Option<ReadEdge>> {
        if !self.edges.contains_key(&id) {
            let made = self.make_edge(id)?;
            self.edges.insert(id, made);
        }
        Ok(self.edges[&id].clone())
    }

    fn make_edge(&mut self, id: Id) -> Read<Option<ReadEdge>> {
        let data = self.geometry.data;
        let r = data.record(id, "EDGE_CURVE").at(id)?;
        let (start, end) = (r.reference(1).at(id)?, r.reference(2).at(id)?);
        let same_sense = r.boolean(4).at(id)?;
        let curve = self.geometry.edge_curve(r.reference(3).at(id)?)?;
        // golf's edges run with their curve's parameter.
        let (from, to) = if same_sense {
            (start, end)
        } else {
            (end, start)
        };
        let (from_v, to_v) = (self.vertex(from)?, self.vertex(to)?);
        let p = |v: VertexId| self.body.vertex(v).point;
        let Some(range) =
            edge_range(&curve.curve, p(from_v), p(to_v), from == to).map_err(|m| invalid(id, m))?
        else {
            return Ok(None);
        };
        let edge = self
            .body
            .add_edge(curve.curve.clone(), range, from_v, to_v)
            .map_err(|e| invalid(id, e.to_string()))?;
        Ok(Some(ReadEdge {
            id: edge,
            against: !same_sense,
            pcurves: curve.pcurves,
        }))
    }
}

/// An edge read from a file.
#[derive(Clone)]
struct ReadEdge {
    id: EdgeId,
    /// Whether golf's edge runs against the file's.
    against: bool,
    /// The file's pcurves for it: the surface each is on, and its 2D curve.
    pcurves: Vec<(Id, Id)>,
}

/// The parameter range of `curve` from `start` to `end`, the whole period for
/// a closed edge; `None` for an edge with no extent.
fn edge_range(
    curve: &AnyCurve<World>,
    start: Point<World, 3>,
    end: Point<World, 3>,
    closed: bool,
) -> Result<Option<(f64, f64)>, String> {
    let project = |p: Point<World, 3>| {
        curve
            .project(p, None)
            .map(|t| t.coords.x)
            .map_err(|e| format!("its vertex couldn't be found on its curve: {e}"))
    };
    let axis = curve.domain().axes[0];
    let t0 = project(start)?;
    let mut t1 = project(end)?;
    match (axis.period(), closed) {
        (Some(period), true) => t1 = t0 + period,
        (Some(period), false) => {
            while t1 <= t0 {
                t1 += period;
            }
        }
        (None, true) if axis.min.is_finite() && axis.max.is_finite() => {
            return Ok(Some((axis.min, axis.max)));
        }
        (None, _) => {}
    }
    if t1 - t0 <= 1e-12 * (1.0 + t0.abs()) {
        if closed {
            return Ok(None);
        }
        return Err(format!(
            "its vertices are in the wrong order along its curve ({t0} to {t1})"
        ));
    }
    Ok(Some((t0, t1)))
}

/// A file's 2D curve as a pcurve in golf's surface parameters.
fn scaled(curve: Curve2, [su, sv]: [f64; 2]) -> Option<AnyCurve2<FaceUv<World>>> {
    let scale = |v: Vector2<f64>| Vector2::new(v.x * su, v.y * sv);
    match curve {
        Curve2::Line(origin, direction) => {
            Some(Line::new(Point::new(scale(origin)), Vector::new(scale(direction))).into())
        }
        Curve2::Nurbs(nurbs) => {
            let count = nurbs.control_point_count();
            let points = (0..count).map(|i| scale(nurbs.control_point(i))).collect();
            let weights = (0..count).map(|i| nurbs.weight(i)).collect();
            let curve = golf_nurbs::NurbsCurve::new(
                nurbs.degree(),
                nurbs.knots().knots().to_vec(),
                points,
                Some(weights),
            )
            .ok()?;
            Some(NurbsCurve::new(curve).into())
        }
    }
}

/// Gives each coedge one of its candidate pcurves: the only one, or for a
/// seam's two, the one continuing from the coedge before (or into the one
/// after). Coedges without candidates are left for healing.
fn choose_pcurves(
    body: &Body<World>,
    coedges: &mut [Coedge<World>],
    candidates: Vec<Vec<AnyCurve2<FaceUv<World>>>>,
) {
    let ends = |c: &Coedge<World>| {
        let range = body.edge(c.edge).range;
        if c.reversed {
            (range.1, range.0)
        } else {
            range
        }
    };
    let at = |p: &AnyCurve2<FaceUv<World>>, t: f64| p.apply(Point::new([t].into())).coords;
    let n = coedges.len();
    for (i, list) in candidates.iter().enumerate() {
        if list.len() == 1 {
            coedges[i].pcurve = Some(list[0].clone());
        }
    }
    for (i, list) in candidates.into_iter().enumerate() {
        if list.len() < 2 {
            continue;
        }
        let (start, end) = ends(&coedges[i]);
        let before = &coedges[(i + n - 1) % n];
        let after = &coedges[(i + 1) % n];
        let target = match (&before.pcurve, &after.pcurve) {
            (Some(p), _) => Some((at(p, ends(before).1), start)),
            (None, Some(p)) => Some((at(p, ends(after).0), end)),
            _ => None,
        };
        let chosen = match target {
            Some((uv, t)) => list
                .into_iter()
                .min_by(|a, b| (at(a, t) - uv).norm().total_cmp(&(at(b, t) - uv).norm())),
            None => None,
        };
        coedges[i].pcurve = chosen;
    }
}
