//! Triangle meshes of boundary-represented bodies.
//!
//! [`mesh`] tessellates every face of a [`MeshSource`] (any boundary
//! representation that provides vertices, edges and faces) to a [`Tolerance`],
//! producing one [`Mesh`] in which faces meeting along an edge share that
//! edge's vertices, so a closed source gives a watertight mesh.
//!
//! The pipeline:
//!
//! 1. Each edge is sampled once, adaptively, and both faces using it reuse the
//!    samples.
//! 2. Each face's loops are taken into its surface's parameter space, through
//!    their pcurves or by projection. Gaps along a singular bound (a pole) are
//!    bridged with points that all raise to the pole's vertex. A face with no
//!    loops is bounded by its whole (finite) domain.
//! 3. The loops are triangulated by constrained Delaunay triangulation in uv
//!    space, scaled so distances roughly match the surface's, keeping the
//!    triangles inside the loops.
//! 4. Triangles that sag from the surface, turn too sharply or are too long are
//!    split, then the triangulation is refined for shape without touching the
//!    boundary.
//! 5. Vertices are raised to the surface; triangles collapsed at poles are
//!    dropped, and winding and normals follow the face's orientation.

mod boundary;
mod error;
mod mesh;
pub mod prelude;
mod sample;
mod source;
mod tessellate;
mod tolerance;
mod triangulate;

pub use error::MeshError;
pub use mesh::Mesh;
pub use mesh::Triangle;
pub use source::CoedgeData;
pub use source::EdgeData;
pub use source::FaceData;
pub use source::MeshErrorOf;
pub use source::MeshOf;
pub use source::MeshSource;
pub use tessellate::mesh;
pub use tolerance::Tolerance;
