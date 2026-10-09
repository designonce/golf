//! Boundary representation: vertices, edges and faces glued into shells.
//!
//! A [`Body`] owns its entities in arenas and refers to them by index ([`VertexId`],
//! [`EdgeId`], [`FaceId`], [`ShellId`]). Geometry is type-erased ([`AnyCurve`](golf_geom::AnyCurve),
//! [`AnySurface`](golf_geom::AnySurface)) so one body can mix kinds freely.
//!
//! # Orientation
//!
//! - An [`Edge`] runs along its curve over `range`, from `start` to `end`.
//! - A [`Face`]'s normal is its surface's normal, flipped unless `same_sense`.
//!   For a solid, face normals point out of the material.
//! - Each [`Loop`] of a face is a cycle of [`Coedge`]s, one face's use of an
//!   edge. Walking a loop with the face normal pointing up, the face is on the
//!   left: outer boundaries run anticlockwise, holes clockwise.
//! - In a closed manifold body every edge is used exactly twice, once in each
//!   direction (possibly by the same face, as along a cylinder's seam).

mod body;
mod entity;
mod error;
#[cfg(feature = "export_step")]
mod export_step;
mod id;
#[cfg(feature = "mesh")]
mod mesh;
mod primitives;

pub use body::Body;
pub use entity::Coedge;
pub use entity::Edge;
pub use entity::Face;
pub use entity::FaceUv;
pub use entity::Loop;
pub use entity::Shell;
pub use entity::Vertex;
pub use error::TopologyError;
pub use id::EdgeId;
pub use id::FaceId;
pub use id::ShellId;
pub use id::VertexId;
pub use primitives::cuboid;
pub use primitives::cylinder;
