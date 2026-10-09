//! Building bodies with CAD-style operations.
//!
//! Draw closed [`Profile`]s of lines and arcs in a sketch plane, group them into
//! a [`Region`] (an outer boundary and holes), then turn the region into a solid
//! [`Body`](golf_brep::Body) with [`extrude`] or [`revolve`]. The [`primitives`]
//! are built the same way.
//!
//! A sketch plane is a [`Placement`](golf_geom::Placement): sketch coordinates
//! `(x, y)` are along its local x and y axes, and its local z is the sketch's
//! normal.
//!
//! Operations that need booleans (pockets, bosses onto existing bodies,
//! fillets) aren't here yet.

mod error;
mod extrude;
mod loft;
mod pattern;
pub mod primitives;
mod revolve;

pub use error::ModelError;
pub use extrude::extrude;
// The sketch types, for convenience.
pub use golf_sketch::Profile;
pub use golf_sketch::ProfileBuilder;
pub use golf_sketch::Region;
pub use golf_sketch::Segment;
pub use golf_sketch::SketchError;
pub use loft::loft;
pub use pattern::circular_pattern;
pub use pattern::linear_pattern;
pub use revolve::revolve;
pub use revolve::revolve_by;
