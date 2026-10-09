//! A minimal viewer for golf meshes, in Bevy.
//!
//! Add [`golf_mesh::Mesh`]es to a [`Viewer`], each face coloured as you choose,
//! then [`Viewer::run`] opens a window showing them, lit from the camera.
//! Drag with the left button to orbit, with the right to pan, and scroll to
//! zoom.
//!
//! Every mesh's coordinates are shown as they are, with z up, so meshes in
//! different frames should be brought into one first.

mod orbit;
pub mod prelude;
mod viewer;

pub use orbit::Orbit;
pub use viewer::Viewer;
