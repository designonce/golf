//! Boundary Representation (B-REP)
//!
//! Geometry is built from typed coordinate [`Space`]s and the [`Mapping`]s between
//! them. Every [`Point`] and [`Vector`] carries the space it lives in, so mixing
//! e.g. a sphere's uv coordinates with world coordinates is a type error.
//!
//! Geometry is generic over the space it lives in: a `Sphere<World>` is placed
//! directly in world space, while a `Sphere<Frame>` lives in a [`frame::FrameTree`]
//! node and can be composed with that frame's placement.
//!
//! A curve on a surface (a pcurve) is a `Mapping<1, 2>` into the surface's uv space;
//! composing it with the surface gives the edge in the surface's space:
//!
//! ```
//! use golf_brep::{Compose, Mapping, Point, Uv, Vector, World};
//! use golf_brep::geom::{Line2, Sphere};
//! use nalgebra::{Vector1, Vector2, Vector3};
//!
//! let sphere = Sphere::<World>::new(Point::new(Vector3::zeros()), 1.0);
//! let pcurve = Line2::<Uv<Sphere<World>>>::new(Point::new(Vector2::zeros()), Vector::new(Vector2::x()));
//! let edge = Compose::<_, _, 2>(pcurve, sphere);
//! let _world = edge.apply(Point::new(Vector1::new(0.5)));
//! ```
//!
//! A pcurve in some other surface's uv space does not compose with the sphere:
//!
//! ```compile_fail
//! use golf_brep::{Compose, Mapping, Point, Uv, Vector, World};
//! use golf_brep::geom::{Line2, Plane, Sphere};
//! use nalgebra::{Vector1, Vector2, Vector3};
//!
//! let sphere = Sphere::<World>::new(Point::new(Vector3::zeros()), 1.0);
//! let pcurve = Line2::<Uv<Plane<World>>>::new(Point::new(Vector2::zeros()), Vector::new(Vector2::x()));
//! let edge = Compose::<_, _, 2>(pcurve, sphere);
//! let _world = edge.apply(Point::new(Vector1::new(0.5)));
//! ```

mod domain;
mod erased;
pub mod frame;
pub mod geom;
mod manifold;
mod mapping;
mod space;

pub use domain::Axis;
pub use domain::Domain;
pub use erased::AnyCurve;
pub use erased::AnyCurve2;
pub use erased::AnySurface;
pub use erased::DynEmbedding;
pub use manifold::Curve;
pub use manifold::Embedding;
pub use manifold::ProjectError;
pub use manifold::Surface;
pub use manifold::newton_project;
pub use mapping::Compose;
pub use mapping::Mapping;
pub use space::Point;
pub use space::Space;
pub use space::T;
pub use space::Uv;
pub use space::Vector;
pub use space::World;
