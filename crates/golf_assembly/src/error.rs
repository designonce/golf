use golf_frame::FrameId;
use golf_geom::TransformError;

use crate::id::PartId;

/// Why an assembly operation failed.
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum AssemblyError {
    #[error("no part {0}")]
    MissingPart(PartId),
    #[error("no frame {0:?} in this assembly")]
    MissingFrame(FrameId),
    /// Only groups can hold other frames, or be copied.
    #[error("frame {0:?} is a part instance, not a group")]
    NotAGroup(FrameId),
    /// The root is the assembly's own coordinates; it isn't placed.
    #[error("the root frame has no placement")]
    RootPlacement,
    #[error(transparent)]
    Transform(#[from] TransformError),
}
