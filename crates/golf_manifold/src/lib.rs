//! Coordinate spaces, mappings between them, and manifolds.
//!
//! Geometry is built from typed coordinate [`Space`]s and the [`Mapping`]s between
//! them. Every [`Point`] and [`Vector`] carries the space it lives in, so mixing
//! e.g. a sphere's uv coordinates with world coordinates is a type error.
//!
//! Geometry is generic over the space it lives in: a `Sphere<World>` is placed
//! directly in world space, while a `Sphere<Frame>` lives in a `FrameTree` (see `golf_frame`)
//! node and can be composed with that frame's placement.
//!

mod domain;
mod manifold;
mod mapping;
pub mod prelude;
mod space;

pub use domain::Axis;
pub use domain::Domain;
pub use manifold::Curve;
pub use manifold::Embedding;
pub use manifold::ProjectError;
pub use manifold::Surface;
pub use manifold::newton_project;
pub use mapping::Compose;
pub use mapping::Mapping;
pub use space::HasT;
pub use space::HasUv;
pub use space::Point;
pub use space::Space;
pub use space::T;
pub use space::Uv;
pub use space::Vector;
pub use space::World;
