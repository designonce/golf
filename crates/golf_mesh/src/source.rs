//! What the mesher needs from a boundary representation.

use core::fmt::Debug;
use core::hash::Hash;

use golf_manifold::Embedding;
use golf_manifold::Mapping;
use golf_manifold::Point;
use golf_manifold::Space;

use crate::error::MeshError;
use crate::mesh::Mesh;

/// A boundary representation the mesher can tessellate.
///
/// Implement it for any B-rep (or anything shaped like one) to mesh it with
/// [`mesh`](crate::mesh). The mesher needs:
///
/// - vertices, as points;
/// - edges, each a piece of a curve between two vertices;
/// - faces, each a piece of a surface bounded by loops of coedges (uses of
///   edges, possibly with pcurves).
///
/// Orientation follows the usual B-rep conventions: a face's normal is its
/// surface's, flipped unless `same_sense`; walking a loop with the face normal
/// up, the face is on the left. A face with no loops covers its surface's
/// whole (finite) domain. A face on a periodic surface must be cut open by a
/// seam edge used twice, with pcurves on both uses.
pub trait MeshSource {
    /// The space the geometry, and so the mesh, is in.
    type Space: Space<3>;
    /// The geometry of edges.
    type Curve: Embedding<1, 3, To = Self::Space, From: Space<1, Tag = ()>>;
    /// The geometry of faces.
    type Surface: Embedding<2, 3, To = Self::Space, From: Space<2, Tag = ()>>;
    /// Edges in a face's parameter space. A pcurve shares its edge curve's
    /// parameter: `surface(pcurve(t)) == curve(t)` over the edge's range.
    type Pcurve: Mapping<1, 2, To = <Self::Surface as Mapping<2, 3>>::From, From: Space<1, Tag = ()>>;
    type Vertex: Copy + Eq + Hash + Debug;
    type Edge: Copy + Eq + Hash + Debug;
    type Face: Copy + Eq + Hash + Debug;

    fn vertices(&self) -> impl Iterator<Item = (Self::Vertex, Point<Self::Space, 3>)> + '_;
    fn edges(&self) -> impl Iterator<Item = (Self::Edge, EdgeData<'_, Self>)> + '_;
    fn faces(&self) -> impl Iterator<Item = (Self::Face, FaceData<'_, Self>)> + '_;
}

/// An edge: the part of `curve` over `range` (increasing), from `start` to `end`.
pub struct EdgeData<'a, M: MeshSource + ?Sized> {
    pub curve: &'a M::Curve,
    pub range: (f64, f64),
    pub start: M::Vertex,
    pub end: M::Vertex,
}

/// A face: the part of `surface` its `loops` bound.
pub struct FaceData<'a, M: MeshSource + ?Sized> {
    pub surface: &'a M::Surface,
    /// Whether the face normal agrees with the surface normal.
    pub same_sense: bool,
    pub loops: Vec<Vec<CoedgeData<'a, M>>>,
}

/// One use of an edge in a face's loop.
pub struct CoedgeData<'a, M: MeshSource + ?Sized> {
    pub edge: M::Edge,
    /// Whether the loop runs the edge from `end` to `start`.
    pub reversed: bool,
    pub pcurve: Option<&'a M::Pcurve>,
}

/// The mesh of a source `M`.
pub type MeshOf<M> = Mesh<<M as MeshSource>::Space, <M as MeshSource>::Face>;

/// Why a source `M` couldn't be meshed.
pub type MeshErrorOf<M> = MeshError<<M as MeshSource>::Face, <M as MeshSource>::Edge>;
