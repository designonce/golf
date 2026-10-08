use golf_brep::EdgeId;
use golf_brep::FaceId;
use golf_manifold::ProjectError;

/// Why a body could not be meshed.
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum MeshError {
    /// A face with no loops must cover its whole surface, which needs a finite
    /// domain.
    #[error("{face} has no loops and an unbounded surface")]
    UnboundedFace { face: FaceId },
    /// A loop goes all the way round a periodic surface, so it doesn't bound a
    /// region of the parameter space; the face needs a seam edge.
    #[error("loop {loop_index} of {face} wraps round its surface; it needs a seam")]
    LoopWrapsSurface { face: FaceId, loop_index: usize },
    /// The face's loops cross each other in parameter space.
    #[error("the loops of {face} cross themselves in parameter space")]
    SelfIntersectingBoundary { face: FaceId },
    /// A point of an edge without a pcurve couldn't be found on the face.
    #[error("projecting {edge} onto {face} failed: {source}")]
    Projection {
        face: FaceId,
        edge: EdgeId,
        source: ProjectError,
    },
    /// The triangulator rejected a point (NaN, or too large).
    #[error("triangulating {face} failed: {message}")]
    Triangulation { face: FaceId, message: String },
}
