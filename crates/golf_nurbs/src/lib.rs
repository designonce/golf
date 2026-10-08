//! Non-uniform rational B-spline (NURBS) curves and surfaces.
//!
//! Plain geometry over `nalgebra` vectors, with no notion of coordinate spaces;
//! `golf_brep` wraps these in its typed mappings.
//!
//! Algorithm numbers in comments (e.g. "A2.2") refer to Piegl & Tiller,
//! *The NURBS Book*, 2nd edition.

mod curve;
mod error;
mod homogeneous;
mod knots;
mod surface;

pub use curve::NurbsCurve;
pub use error::NurbsError;
pub use knots::KnotVector;
pub use surface::NurbsSurface;
