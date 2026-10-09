//! Reading STEP files into golf bodies and assemblies.
//!
//! [`read_step`] parses a file with `golf_step`, builds a [`golf_brep::Body`]
//! for each solid (`MANIFOLD_SOLID_BREP`, `BREP_WITH_VOIDS`,
//! `SHELL_BASED_SURFACE_MODEL`), heals it with `golf_heal` so every coedge has
//! a pcurve, and places them as the file's product structure says, in a
//! [`golf_assembly::Assembly`]. Lengths come out in millimetres and angles in
//! radians, whatever units the file uses; colours and product names are kept.
//!
//! Geometry is read exactly: analytic surfaces and curves as themselves, and
//! swept surfaces as NURBS. A file's own pcurves are kept where they're
//! accurate. Anything that can't be read (an unsupported surface, say) is
//! left out with an [`ImportWarning`], rather than failing the file.

mod body;
mod error;
mod geometry;
pub mod prelude;
mod read;
mod structure;
mod style;
mod units;

pub use error::ImportError;
pub use error::ImportWarning;
pub use error::Problem;
pub use read::ImportOptions;
pub use read::StepImport;
pub use read::read_step;
pub use read::read_step_with;
pub use units::Units;
