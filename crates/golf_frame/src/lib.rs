//! A hierarchy of rigid coordinate frames, and the transforms between them.

pub mod prelude;
mod tree;

pub use tree::Frame;
pub use tree::FrameId;
pub use tree::FrameTree;
pub use tree::Rigid;
