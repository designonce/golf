use golf_manifold::ProjectError;

use super::EdgeId;
use super::FaceId;
use super::VertexId;

/// What's wrong with a body, from construction or [`Body::validate`].
///
/// [`Body::validate`]: super::Body::validate
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum TopologyError {
    #[error("no vertex {0}")]
    MissingVertex(VertexId),
    #[error("no edge {0}")]
    MissingEdge(EdgeId),
    #[error("no face {0}")]
    MissingFace(FaceId),
    #[error("geometry is in a different instance of the body's space")]
    WrongSpace,
    #[error("edge range ({0}, {1}) is not increasing")]
    BadRange(f64, f64),
    #[error("a loop has no coedges")]
    EmptyLoop,
    #[error("loop {loop_index} breaks after coedge {coedge_index}")]
    OpenLoop {
        loop_index: usize,
        coedge_index: usize,
    },
    #[error("{vertex} is {distance} from the end of {edge}'s curve")]
    VertexOffCurve {
        edge: EdgeId,
        vertex: VertexId,
        distance: f64,
    },
    #[error(
        "{edge} is used {forward} times forwards and {reversed} backwards; a closed body uses it once each way"
    )]
    EdgeUses {
        edge: EdgeId,
        forward: usize,
        reversed: usize,
    },
    #[error("{face} is in {shells} shells")]
    FaceShells { face: FaceId, shells: usize },
    #[error("{edge} strays {distance} from {face}'s surface")]
    EdgeOffFace {
        edge: EdgeId,
        face: FaceId,
        distance: f64,
    },
    #[error("projecting {edge} onto {face} failed: {source}")]
    Projection {
        edge: EdgeId,
        face: FaceId,
        source: ProjectError,
    },
    #[error("{edge}'s pcurve on {face} strays {distance} from the edge")]
    PcurveOffEdge {
        edge: EdgeId,
        face: FaceId,
        distance: f64,
    },
}
