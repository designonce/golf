//! What the STEP writer needs from a boundary representation.

use core::fmt::Debug;
use core::hash::Hash;

use golf_color::Color;
use golf_manifold::Point;
use golf_manifold::Space;

use crate::geometry::StepCurve;
use crate::geometry::StepSurface;

/// A boundary representation that can be written as STEP.
///
/// Orientation follows the usual conventions, which are STEP's too: an edge
/// runs along its curve from `start` to `end`; a face's normal is its surface's,
/// flipped unless `same_sense`; walking a loop with the face normal up, the face
/// is on the left. Every face needs at least one loop.
pub trait StepSource {
    type Space: Space<3>;
    type Curve: StepCurve;
    type Surface: StepSurface;
    type Vertex: Copy + Eq + Hash + Debug;
    type Edge: Copy + Eq + Hash + Debug;
    type Face: Copy + Eq + Hash + Debug;

    fn vertices(&self) -> impl Iterator<Item = (Self::Vertex, Point<Self::Space, 3>)> + '_;
    fn edges(&self) -> impl Iterator<Item = (Self::Edge, StepEdge<'_, Self>)> + '_;
    fn faces(&self) -> impl Iterator<Item = (Self::Face, StepFace<'_, Self>)> + '_;
    /// Groups of faces forming shells: the outer shell first, then any voids,
    /// whose faces point into the void.
    fn shells(&self) -> impl Iterator<Item = Vec<Self::Face>> + '_;

    /// The whole body's colour.
    fn color(&self) -> Option<Color> {
        None
    }

    /// A face's own colour, overriding the body's.
    fn face_color(&self, _face: Self::Face) -> Option<Color> {
        None
    }
}

/// An edge along `curve` (in the direction of increasing parameter) from
/// `start` to `end`.
pub struct StepEdge<'a, M: StepSource + ?Sized> {
    pub curve: &'a M::Curve,
    pub start: M::Vertex,
    pub end: M::Vertex,
}

/// A face: the part of `surface` its `loops` bound. Each loop is a cycle of
/// edges, each with whether the loop runs it backwards.
pub struct StepFace<'a, M: StepSource + ?Sized> {
    pub surface: &'a M::Surface,
    pub same_sense: bool,
    pub loops: Vec<Vec<(M::Edge, bool)>>,
}
