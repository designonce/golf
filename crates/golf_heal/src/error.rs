use golf_brep::EdgeId;
use golf_brep::FaceId;

/// Something [`crate::heal`] found and couldn't fix.
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
pub enum HealIssue {
    #[error("{edge} couldn't be projected onto {face}'s surface")]
    Projection { face: FaceId, edge: EdgeId },
    #[error("{edge} is {distance} from {face}'s surface, beyond the tolerance")]
    OffSurface {
        face: FaceId,
        edge: EdgeId,
        distance: f64,
    },
    #[error("{face} wraps round its surface but couldn't be cut open: {reason}")]
    NoSeam { face: FaceId, reason: &'static str },
}
