//! The entities a body is made of.

use golf_geom::AnyCurve;
use golf_geom::AnyCurve2;
use golf_geom::AnySurface;
use golf_manifold::Point;
use golf_manifold::Space;
use golf_manifold::Uv;

use super::EdgeId;
use super::FaceId;
use super::VertexId;

/// A point where edges meet.
#[derive(Clone, Debug)]
pub struct Vertex<S: Space<3>> {
    pub point: Point<S, 3>,
}

/// A bounded piece of a curve between two vertices.
#[derive(Clone, Debug)]
pub struct Edge<S: Space<3>> {
    pub curve: AnyCurve<S>,
    /// The curve parameters at `start` and `end`, increasing. A closed edge
    /// (`start == end`) spans a whole period.
    pub range: (f64, f64),
    pub start: VertexId,
    pub end: VertexId,
}

/// The parameter space of every face's surface, where pcurves live.
pub type FaceUv<S> = Uv<AnySurface<S>>;

/// One face's use of an edge, as part of a [`Loop`].
#[derive(Clone, Debug)]
pub struct Coedge<S: Space<3>> {
    pub edge: EdgeId,
    /// Whether the loop runs the edge from `end` to `start`.
    pub reversed: bool,
    /// The edge in the face's parameter space, sharing the edge curve's
    /// parameter: `surface(pcurve(t)) == curve(t)` over the edge's range.
    pub pcurve: Option<AnyCurve2<FaceUv<S>>>,
}

impl<S: Space<3>> Coedge<S> {
    pub fn new(edge: EdgeId, reversed: bool) -> Self {
        Self {
            edge,
            reversed,
            pcurve: None,
        }
    }

    pub fn with_pcurve(mut self, pcurve: impl Into<AnyCurve2<FaceUv<S>>>) -> Self {
        self.pcurve = Some(pcurve.into());
        self
    }
}

/// A cycle of coedges bounding a face.
#[derive(Clone, Debug)]
pub struct Loop<S: Space<3>> {
    pub coedges: Vec<Coedge<S>>,
}

impl<S: Space<3>> Loop<S> {
    pub fn new(coedges: Vec<Coedge<S>>) -> Self {
        Self { coedges }
    }
}

/// A bounded piece of a surface.
#[derive(Clone, Debug)]
pub struct Face<S: Space<3>> {
    pub surface: AnySurface<S>,
    /// Whether the face normal agrees with the surface normal.
    pub same_sense: bool,
    /// The face's boundary. A face covering a whole closed surface, such as a
    /// sphere or torus, has none.
    pub loops: Vec<Loop<S>>,
    /// The face's own colour, overriding its body's.
    pub color: Option<golf_color::Color>,
}

/// A connected set of faces, closed for a solid.
#[derive(Clone, Debug)]
pub struct Shell {
    pub faces: Vec<FaceId>,
}
