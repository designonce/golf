//! Assemblies: parts placed in a tree of frames.
//!
//! An [`Assembly`] defines [parts](Assembly::add_part) (a named body in its own
//! coordinates) and places them in a [`FrameTree`](golf_frame::FrameTree):
//! every [group](Assembly::add_group) (sub-assembly) and every
//! [instance](Assembly::add_instance) of a part is a frame, placed by a rigid
//! motion in its parent group. A part placed many times is stored once; a
//! group can be [copied](Assembly::copy_group) to place it again.
//!
//! Frames give the relations between any two instances
//! ([`Assembly::transform`]), and moving a group moves everything in it.
//! Colours set on a group or instance override the parts' own below them.

mod assembly;
mod error;
mod id;
pub mod prelude;

pub use assembly::Assembly;
pub use assembly::Node;
pub use assembly::NodeKind;
pub use assembly::Occurrence;
pub use assembly::Part;
pub use error::AssemblyError;
pub use id::PartId;
