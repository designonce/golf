#![cfg_attr(docsrs, feature(doc_cfg))]

//! The crates making up golf, each as a module, and their common items in a
//! prelude. Use them through the `golf` crate rather than directly.

pub mod prelude;

#[cfg(feature = "assembly")]
pub use golf_assembly as assembly;
pub use golf_brep as brep;
pub use golf_color as color;
#[cfg(feature = "export_step")]
pub use golf_export_step as export_step;
pub use golf_frame as frame;
pub use golf_geom as geom;
pub use golf_manifold as manifold;
#[cfg(feature = "mesh")]
pub use golf_mesh as mesh;
#[cfg(feature = "model")]
pub use golf_model as model;
pub use golf_nurbs as nurbs;
pub use golf_sketch as sketch;
#[cfg(feature = "view")]
pub use golf_view as view;
/// The linear algebra golf is written in, re-exported so its version always
/// matches.
pub use nalgebra;
