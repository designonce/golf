//! Concrete curves and surfaces: analytic shapes, NURBS, and type-erased
//! enums over them for topology to store.
//!
//! A curve on a surface (a pcurve) is a `Mapping<1, 2>` into the surface's uv space;
//! composing it with the surface gives the edge in the surface's space:
//!
//! ```
//! use golf_manifold::{Compose, Mapping, Point, Uv, Vector, World};
//! use golf_geom::{Line2, Sphere};
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
//! use golf_manifold::{Compose, Mapping, Point, Uv, Vector, World};
//! use golf_geom::{Line2, Plane, Sphere};
//! use nalgebra::{Vector1, Vector2, Vector3};
//!
//! let sphere = Sphere::<World>::new(Point::new(Vector3::zeros()), 1.0);
//! let pcurve = Line2::<Uv<Plane<World>>>::new(Point::new(Vector2::zeros()), Vector::new(Vector2::x()));
//! let edge = Compose::<_, _, 2>(pcurve, sphere);
//! let _world = edge.apply(Point::new(Vector1::new(0.5)));
//! ```

mod circle;
mod cone;
mod cylinder;
mod ellipse;
mod erased;
mod line;
mod nurbs;
mod placement;
mod plane;
mod sphere;
#[cfg(test)]
mod testing;
mod torus;

pub use circle::Circle;
pub use cone::Cone;
pub use cylinder::Cylinder;
pub use ellipse::Ellipse;
pub use erased::AnyCurve;
pub use erased::AnyCurve2;
pub use erased::AnySurface;
pub use erased::DynEmbedding;
pub use line::Line;
pub use line::Line2;
pub use nurbs::NurbsCurve;
pub use nurbs::NurbsSurface;
pub use placement::Placement;
pub use plane::Plane;
pub use sphere::Sphere;
pub use torus::Torus;
