//! 2D sketches: closed profiles of lines and arcs, and regions bounded by them.
//!
//! Coordinates are a sketch plane's local `(x, y)`; placing the plane in space,
//! and turning regions into solids, is `golf_model`'s business.

mod corner;
mod error;
mod locus;
mod offset;
mod path;
mod profile;
mod region;
mod segment;

pub use error::SketchError;
pub use path::Path;
pub use profile::Profile;
pub use profile::ProfileBuilder;
pub use region::Region;
pub use segment::Segment;
