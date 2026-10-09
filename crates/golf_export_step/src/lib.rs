//! Writing boundary-represented bodies as STEP files (ISO 10303-21, AP242).
//!
//! Add bodies to a [`StepFile`]; each becomes a part whose shape is a
//! `MANIFOLD_SOLID_BREP` (or `BREP_WITH_VOIDS`) of `ADVANCED_FACE`s, or a
//! surface model if a shell isn't closed. Anything implementing [`StepSource`]
//! can be written; its geometry must implement [`StepCurve`] and
//! [`StepSurface`], which `golf_geom`'s types do, and which other geometry can
//! implement with the [`Writer`].
//!
//! `golf_brep`'s [`Body`](golf_brep::Body) is a [`StepSource`], and a
//! `golf_assembly` [`Assembly`](golf_assembly::Assembly) is written whole with
//! [`StepFile::add_assembly`]. Products can also be grouped by hand with
//! [`StepFile::add_group`], which nest and share parts and sub-assemblies.
//!
//! Coordinates are written as millimetres and angles as radians. Every body's
//! coordinates are written as they are, so bodies in different frames should be
//! brought into one first.

mod assembly;
mod body;
mod error;
mod file;
mod geometry;
pub mod prelude;
mod source;
mod writer;

pub use error::StepError;
pub use file::Product;
pub use file::StepFile;
pub use file::StepOptions;
pub use geometry::StepCurve;
pub use geometry::StepSurface;
pub use source::StepEdge;
pub use source::StepFace;
pub use source::StepSource;
pub use writer::Ref;
pub use writer::Writer;
