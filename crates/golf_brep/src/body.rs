use std::collections::HashMap;

use golf_geom::AnyCurve;
use golf_geom::AnySurface;
use golf_manifold::Embedding;
use golf_manifold::Mapping;
use golf_manifold::Point;
use golf_manifold::Space;

use super::Coedge;
use super::Edge;
use super::EdgeId;
use super::Face;
use super::FaceId;
use super::Loop;
use super::Shell;
use super::ShellId;
use super::TopologyError;
use super::Vertex;
use super::VertexId;

/// A boundary-represented body in space `S`.
///
/// Everything in it lives in one instance of `S`, named by the body's tag.
#[derive(Clone, Debug)]
pub struct Body<S: Space<3>> {
    pub(super) tag: S::Tag,
    pub(super) vertices: Vec<Vertex<S>>,
    pub(super) edges: Vec<Edge<S>>,
    pub(super) faces: Vec<Face<S>>,
    pub(super) shells: Vec<Shell>,
}

impl<S: Space<3, Tag = ()>> Default for Body<S> {
    fn default() -> Self {
        Self::new()
    }
}

impl<S: Space<3, Tag = ()>> Body<S> {
    pub fn new() -> Self {
        Self::with_tag(())
    }
}

impl<S: Space<3>> Body<S> {
    /// An empty body in the instance of `S` that `tag` names.
    pub fn with_tag(tag: S::Tag) -> Self {
        Self {
            tag,
            vertices: Vec::new(),
            edges: Vec::new(),
            faces: Vec::new(),
            shells: Vec::new(),
        }
    }

    pub fn tag(&self) -> S::Tag {
        self.tag
    }

    pub fn add_vertex(&mut self, point: Point<S, 3>) -> Result<VertexId, TopologyError> {
        if point.tag() != self.tag {
            return Err(TopologyError::WrongSpace);
        }
        self.vertices.push(Vertex { point });
        Ok(VertexId::new(self.vertices.len() - 1))
    }

    /// Adds the part of `curve` over `range`, from `start` to `end`.
    ///
    /// Whether the vertices lie on the curve is left to [`Self::validate`].
    pub fn add_edge(
        &mut self,
        curve: impl Into<AnyCurve<S>>,
        range: (f64, f64),
        start: VertexId,
        end: VertexId,
    ) -> Result<EdgeId, TopologyError> {
        let curve = curve.into();
        self.try_vertex(start)?;
        self.try_vertex(end)?;
        // Also rejects NaN.
        if range.0.partial_cmp(&range.1) != Some(core::cmp::Ordering::Less) {
            return Err(TopologyError::BadRange(range.0, range.1));
        }
        if curve.apply(Point::new([range.0].into())).tag() != self.tag {
            return Err(TopologyError::WrongSpace);
        }
        self.edges.push(Edge {
            curve,
            range,
            start,
            end,
        });
        Ok(EdgeId::new(self.edges.len() - 1))
    }

    /// Adds a face, checking each loop is a closed cycle of existing edges.
    pub fn add_face(
        &mut self,
        surface: impl Into<AnySurface<S>>,
        same_sense: bool,
        loops: Vec<Loop<S>>,
    ) -> Result<FaceId, TopologyError> {
        let surface = surface.into();
        if surface.apply(Point::new([0.0, 0.0].into())).tag() != self.tag {
            return Err(TopologyError::WrongSpace);
        }
        for (loop_index, l) in loops.iter().enumerate() {
            if l.coedges.is_empty() {
                return Err(TopologyError::EmptyLoop);
            }
            for coedge in &l.coedges {
                self.try_edge(coedge.edge)?;
            }
            for (coedge_index, pair) in l
                .coedges
                .iter()
                .zip(l.coedges.iter().cycle().skip(1))
                .enumerate()
            {
                if self.coedge_end(pair.0) != self.coedge_start(pair.1) {
                    return Err(TopologyError::OpenLoop {
                        loop_index,
                        coedge_index,
                    });
                }
            }
        }
        self.faces.push(Face {
            surface,
            same_sense,
            loops,
        });
        Ok(FaceId::new(self.faces.len() - 1))
    }

    pub fn add_shell(&mut self, faces: Vec<FaceId>) -> Result<ShellId, TopologyError> {
        for &face in &faces {
            self.try_face(face)?;
        }
        self.shells.push(Shell { faces });
        Ok(ShellId::new(self.shells.len() - 1))
    }

    pub fn vertex(&self, id: VertexId) -> &Vertex<S> {
        &self.vertices[id.index()]
    }

    pub fn edge(&self, id: EdgeId) -> &Edge<S> {
        &self.edges[id.index()]
    }

    pub fn face(&self, id: FaceId) -> &Face<S> {
        &self.faces[id.index()]
    }

    pub fn shell(&self, id: ShellId) -> &Shell {
        &self.shells[id.index()]
    }

    pub fn vertices(&self) -> impl ExactSizeIterator<Item = (VertexId, &Vertex<S>)> {
        self.vertices
            .iter()
            .enumerate()
            .map(|(i, v)| (VertexId::new(i), v))
    }

    pub fn edges(&self) -> impl ExactSizeIterator<Item = (EdgeId, &Edge<S>)> {
        self.edges
            .iter()
            .enumerate()
            .map(|(i, e)| (EdgeId::new(i), e))
    }

    pub fn faces(&self) -> impl ExactSizeIterator<Item = (FaceId, &Face<S>)> {
        self.faces
            .iter()
            .enumerate()
            .map(|(i, f)| (FaceId::new(i), f))
    }

    pub fn shells(&self) -> impl ExactSizeIterator<Item = (ShellId, &Shell)> {
        self.shells
            .iter()
            .enumerate()
            .map(|(i, s)| (ShellId::new(i), s))
    }

    /// The vertex a coedge starts from, following its direction.
    pub fn coedge_start(&self, coedge: &Coedge<S>) -> VertexId {
        let edge = self.edge(coedge.edge);
        if coedge.reversed {
            edge.end
        } else {
            edge.start
        }
    }

    /// The vertex a coedge ends at, following its direction.
    pub fn coedge_end(&self, coedge: &Coedge<S>) -> VertexId {
        let edge = self.edge(coedge.edge);
        if coedge.reversed {
            edge.start
        } else {
            edge.end
        }
    }

    /// Every use of each edge: the face, and whether it runs the edge reversed.
    pub fn edge_uses(&self) -> HashMap<EdgeId, Vec<(FaceId, bool)>> {
        let mut uses: HashMap<EdgeId, Vec<(FaceId, bool)>> = HashMap::new();
        for (face_id, face) in self.faces() {
            for coedge in face.loops.iter().flat_map(|l| &l.coedges) {
                uses.entry(coedge.edge)
                    .or_default()
                    .push((face_id, coedge.reversed));
            }
        }
        uses
    }

    /// The genus from the Euler–Poincaré formula `V − E + F − R = 2(S − G)`,
    /// where `R` counts loops beyond one per face. `None` if some face has no
    /// loops (the formula needs every face bounded) or the counts don't give a
    /// non-negative whole genus.
    pub fn genus(&self) -> Option<usize> {
        if self.faces.iter().any(|f| f.loops.is_empty()) {
            return None;
        }
        let loops: usize = self.faces.iter().map(|f| f.loops.len()).sum();
        let chi = self.vertices.len() as isize - self.edges.len() as isize
            + self.faces.len() as isize
            - (loops as isize - self.faces.len() as isize);
        let twice_genus = 2 * self.shells.len() as isize - chi;
        (twice_genus >= 0 && twice_genus % 2 == 0).then_some((twice_genus / 2) as usize)
    }

    /// Checks the body is a closed manifold whose geometry agrees with its
    /// topology, to within `tolerance`:
    ///
    /// - every face is in exactly one shell;
    /// - every edge is used once in each direction;
    /// - each edge's vertices are at its curve's ends;
    /// - each edge lies on the surfaces of the faces using it;
    /// - each pcurve, mapped through its face's surface, follows its edge.
    ///
    /// Reports every problem found rather than stopping at the first.
    pub fn validate(&self, tolerance: f64) -> Result<(), Vec<TopologyError>> {
        let mut errors = Vec::new();

        let mut shells_of = vec![0; self.faces.len()];
        for shell in &self.shells {
            for face in &shell.faces {
                shells_of[face.index()] += 1;
            }
        }
        for (i, &shells) in shells_of.iter().enumerate() {
            if shells != 1 {
                errors.push(TopologyError::FaceShells {
                    face: FaceId::new(i),
                    shells,
                });
            }
        }

        let uses = self.edge_uses();
        for (edge_id, edge) in self.edges() {
            let edge_uses = uses.get(&edge_id).map_or(&[][..], Vec::as_slice);
            let reversed = edge_uses.iter().filter(|(_, r)| *r).count();
            let forward = edge_uses.len() - reversed;
            if (forward, reversed) != (1, 1) {
                errors.push(TopologyError::EdgeUses {
                    edge: edge_id,
                    forward,
                    reversed,
                });
            }

            for (vertex, t) in [(edge.start, edge.range.0), (edge.end, edge.range.1)] {
                let distance = (edge.curve.apply(Point::new([t].into()))
                    - self.vertex(vertex).point)
                    .coords
                    .norm();
                if distance > tolerance {
                    errors.push(TopologyError::VertexOffCurve {
                        edge: edge_id,
                        vertex,
                        distance,
                    });
                }
            }
        }

        for (face_id, face) in self.faces() {
            for coedge in face.loops.iter().flat_map(|l| &l.coedges) {
                if let Err(error) = self.check_coedge(face_id, face, coedge, tolerance) {
                    errors.push(error);
                }
            }
        }

        if errors.is_empty() {
            Ok(())
        } else {
            Err(errors)
        }
    }

    /// Checks one coedge's edge lies on the face, and its pcurve follows it.
    fn check_coedge(
        &self,
        face_id: FaceId,
        face: &Face<S>,
        coedge: &Coedge<S>,
        tolerance: f64,
    ) -> Result<(), TopologyError> {
        const SAMPLES: usize = 8;
        let edge = self.edge(coedge.edge);
        let (t0, t1) = edge.range;
        for i in 0..=SAMPLES {
            let t = Point::new([t0 + (t1 - t0) * i as f64 / SAMPLES as f64].into());
            let on_edge = edge.curve.apply(t);
            let uv = face.surface.project(on_edge, None).map_err(|source| {
                TopologyError::Projection {
                    edge: coedge.edge,
                    face: face_id,
                    source,
                }
            })?;
            let distance = (face.surface.apply(uv) - on_edge).coords.norm();
            if distance > tolerance {
                return Err(TopologyError::EdgeOffFace {
                    edge: coedge.edge,
                    face: face_id,
                    distance,
                });
            }
            if let Some(pcurve) = &coedge.pcurve {
                // A pcurve shares its edge curve's parameter, by convention.
                let shared = Point::new(t.coords);
                let distance = (face.surface.apply(pcurve.apply(shared)) - on_edge)
                    .coords
                    .norm();
                if distance > tolerance {
                    return Err(TopologyError::PcurveOffEdge {
                        edge: coedge.edge,
                        face: face_id,
                        distance,
                    });
                }
            }
        }
        Ok(())
    }

    fn try_vertex(&self, id: VertexId) -> Result<&Vertex<S>, TopologyError> {
        self.vertices
            .get(id.index())
            .ok_or(TopologyError::MissingVertex(id))
    }

    fn try_edge(&self, id: EdgeId) -> Result<&Edge<S>, TopologyError> {
        self.edges
            .get(id.index())
            .ok_or(TopologyError::MissingEdge(id))
    }

    fn try_face(&self, id: FaceId) -> Result<&Face<S>, TopologyError> {
        self.faces
            .get(id.index())
            .ok_or(TopologyError::MissingFace(id))
    }
}

#[cfg(test)]
mod tests {
    use golf_geom::Line;
    use golf_geom::Placement;
    use golf_geom::Sphere;
    use golf_manifold::Vector;
    use golf_manifold::World;
    use nalgebra::Vector2;
    use nalgebra::Vector3;

    use super::*;
    use crate::cuboid;

    #[test]
    fn add_face_rejects_open_loops() {
        let mut body = Body::<World>::new();
        let p = |x: f64| Point::new(Vector3::new(x, 0.0, 0.0));
        let [a, b, c] = [0.0, 1.0, 2.0].map(|x| body.add_vertex(p(x)).unwrap());
        let line = |from: f64| Line::new(p(from), Vector::new(Vector3::x()));
        let ab = body.add_edge(line(0.0), (0.0, 1.0), a, b).unwrap();
        let bc = body.add_edge(line(1.0), (0.0, 1.0), b, c).unwrap();
        let plane = golf_geom::Plane::new(Placement::at(p(0.0)));
        let open = Loop::new(vec![Coedge::new(ab, false), Coedge::new(bc, false)]);
        assert_eq!(
            body.add_face(plane.clone(), true, vec![open]),
            Err(TopologyError::OpenLoop {
                loop_index: 0,
                coedge_index: 1
            })
        );
        let missing = Loop::new(vec![Coedge::new(EdgeId(7), false)]);
        assert_eq!(
            body.add_face(plane.clone(), true, vec![missing]),
            Err(TopologyError::MissingEdge(EdgeId(7)))
        );
        assert_eq!(
            body.add_face(plane, true, vec![Loop::new(vec![])]),
            Err(TopologyError::EmptyLoop)
        );
        assert_eq!(
            body.add_edge(line(0.0), (1.0, 1.0), a, b),
            Err(TopologyError::BadRange(1.0, 1.0))
        );
    }

    #[test]
    fn a_sphere_needs_no_edges() {
        let mut body = Body::<World>::new();
        let sphere = Sphere::new(Point::new(Vector3::zeros()), 1.0);
        let face = body.add_face(sphere, true, vec![]).unwrap();
        body.add_shell(vec![face]).unwrap();
        assert_eq!(body.validate(1e-9), Ok(()));
        assert_eq!(body.genus(), None);
    }

    #[test]
    fn pcurve_mismatch_is_reported() {
        let mut body = cuboid(
            Placement::<World>::at(Point::new(Vector3::zeros())),
            Vector3::new(1.0, 2.0, 3.0),
        );
        let coedge = &mut body.faces[0].loops[0].coedges[0];
        coedge.pcurve = Some(
            Line::new(
                Point::new(Vector2::new(0.5, 0.5)),
                Vector::new(Vector2::new(1.0, 0.0)),
            )
            .into(),
        );
        let errors = body.validate(1e-9).unwrap_err();
        assert!(
            matches!(
                errors[..],
                [TopologyError::PcurveOffEdge {
                    face: FaceId(0),
                    ..
                }]
            ),
            "{errors:?}"
        );
    }
}
