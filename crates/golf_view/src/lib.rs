//! A minimal viewer for golf meshes, in Bevy.
//!
//! Add [`golf_mesh::Mesh`]es to a [`Viewer`], each face coloured as you choose,
//! then [`Viewer::run`] opens a window showing them, lit from the camera.
//! Drag with the left button to orbit, with the right to pan, and scroll to
//! zoom; F frames everything again.
//!
//! In an app of your own, add the [`ViewerPlugin`] and put a [`Scene`] in
//! [`ShowScene`] whenever there's something new to show.
//!
//! Every mesh's coordinates are shown as they are, with z up, so meshes in
//! different frames should be brought into one first.

mod orbit;
mod plugin;
pub mod prelude;
mod scene;
mod viewer;

pub use orbit::Orbit;
pub use plugin::ShowScene;
pub use plugin::ViewerPlugin;
pub use scene::Scene;
pub use viewer::Viewer;
