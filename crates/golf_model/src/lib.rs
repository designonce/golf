//! Building bodies with CAD-style operations.
//!
//! Sketch closed [`Profile`]s (lines and arcs, with fillets, chamfers and
//! offsets from `golf_sketch`) in a sketch plane, group them into a [`Region`]
//! (an outer boundary and holes), then make a solid
//! [`Body`](golf_brep::Body):
//!
//! - [`extrude`]: straight up the sketch's normal;
//! - [`extrude_with`]: with the edges round either end chamfered or filleted;
//! - [`extrude_drafted`]: with tapering sides;
//! - [`extrude_hollow`]: as a tube or cup with walls of a given thickness;
//! - [`revolve`], [`revolve_by`], [`revolve_about`]: round an axis in the
//!   sketch, fully or partly;
//! - [`loft`]: ruled between two matching regions in their own planes;
//! - [`pipe`]: a circle swept along a planar [`Path`].
//!
//! [`primitives`] (boxes, rounded and chamfered boxes, cylinders, cones,
//! spheres, tori) are built from these, and [`linear_pattern`] and
//! [`circular_pattern`] copy bodies. Bodies can be moved with
//! [`Body::transformed`](golf_brep::Body::transformed).
//!
//! A sketch plane is a [`Placement`](golf_geom::Placement): sketch coordinates
//! `(x, y)` are along its local x and y axes, and its local z is the sketch's
//! normal.
//!
//! Edges are filleted and chamfered as solids are built (in the sketch, or at
//! an extrusion's ends), not afterwards on an existing body: that, like
//! combining bodies, needs booleans, which aren't here yet.

mod ends;
mod error;
mod extrude;
mod hollow;
mod loft;
mod pattern;
mod pipe;
pub mod prelude;
pub mod primitives;
mod revolve;
mod stack;

pub use ends::EdgeTreatment;
pub use ends::ExtrudeEnds;
pub use ends::extrude_drafted;
pub use ends::extrude_with;
pub use error::ModelError;
pub use extrude::extrude;
// The sketch types, for convenience.
pub use golf_sketch::Path;
pub use golf_sketch::Profile;
pub use golf_sketch::ProfileBuilder;
pub use golf_sketch::Region;
pub use golf_sketch::Segment;
pub use golf_sketch::SketchError;
pub use hollow::extrude_hollow;
pub use loft::loft;
pub use pattern::circular_pattern;
pub use pattern::linear_pattern;
pub use pipe::pipe;
pub use revolve::revolve;
pub use revolve::revolve_about;
pub use revolve::revolve_by;
