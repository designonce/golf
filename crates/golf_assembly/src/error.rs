use golf_geom::TransformError;

use crate::id::GroupId;
use crate::id::PartId;

/// Why an assembly operation failed.
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum AssemblyError {
    #[error("no part {0}")]
    MissingPart(PartId),
    #[error("no group {0}")]
    MissingGroup(GroupId),
    /// Placing the group would make it contain itself.
    #[error("placing {child} in {parent} would make {child} contain itself")]
    Cycle { parent: GroupId, child: GroupId },
    #[error("group {group} has no child {index}")]
    MissingChild { group: GroupId, index: usize },
    #[error(transparent)]
    Transform(#[from] TransformError),
}
