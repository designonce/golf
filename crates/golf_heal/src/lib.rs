//! Making boundary-represented bodies from outside sources complete and
//! consistent, so everything downstream can rely on them.
//!
//! Importers build bodies as their files describe them, which is often less
//! than golf's other crates expect. [`heal`] fills in what's missing:
//!
//! - every coedge gets a pcurve that follows its edge on its face's surface,
//!   keeping a file's own where it's accurate and computing the rest by
//!   projecting the edge, carried continuously around each loop so loops
//!   close in parameter space across periodic seams, and with poles taking
//!   the parameter of their neighbours;
//! - a face whose loops wrap round a periodic surface (a cylinder's side
//!   bounded by two circles, a sphere's cap by one), as files often leave
//!   them, is cut open by a seam edge along the surface's own iso-parameter
//!   curve, splitting a loop's edge where the seam meets it between vertices.
//!
//! What it finds and can't fix is reported in the [`HealReport`].

mod error;
mod heal;
mod pcurve;
mod seam;

pub use error::HealIssue;
pub use heal::HealOptions;
pub use heal::HealReport;
pub use heal::heal;
