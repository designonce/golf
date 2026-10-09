use golf_manifold::ProjectError;

/// Why a source could not be meshed. `F` and `E` are the source's face and
/// edge keys.
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum MeshError<F: core::fmt::Debug, E: core::fmt::Debug> {
    /// A face with no loops must cover its whole surface, which needs a finite
    /// domain.
    #[error("face {face:?} has no loops and an unbounded surface")]
    UnboundedFace { face: F },
    /// A loop goes all the way round a periodic surface, so it doesn't bound a
    /// region of the parameter space; the face needs a seam edge.
    #[error("loop {loop_index} of face {face:?} wraps round its surface; it needs a seam")]
    LoopWrapsSurface { face: F, loop_index: usize },
    /// The face's loops cross each other in parameter space.
    #[error("the loops of face {face:?} cross themselves in parameter space")]
    SelfIntersectingBoundary { face: F },
    /// A point of an edge without a pcurve couldn't be found on the face.
    #[error("projecting edge {edge:?} onto face {face:?} failed: {source}")]
    Projection {
        face: F,
        edge: E,
        source: ProjectError,
    },
    /// The triangulator rejected a point (NaN, or too large).
    #[error("triangulating face {face:?} failed: {message}")]
    Triangulation { face: F, message: String },
}
