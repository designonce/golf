use golf_brep::TopologyError;
use golf_sketch::SketchError;

/// Why an operation couldn't produce a body.
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum ModelError {
    #[error(transparent)]
    Sketch(#[from] SketchError),
    #[error("extrusion distance must be finite and non-zero")]
    ZeroDistance,
    #[error("revolution angle {0} must be in (0, 2π]")]
    BadAngle(f64),
    /// Revolving needs the profile on one side of the axis (sketch x >= 0).
    #[error("segment {index} crosses the axis of revolution")]
    CrossesAxis { index: usize },
    /// An arc whose circle reaches the axis would revolve to a self-intersecting
    /// (spindle) torus.
    #[error("segment {index} would revolve to a spindle torus")]
    SpindleTorus { index: usize },
    /// Lofted regions need the same number of holes, and each profile the
    /// same number of segments as its counterpart.
    #[error("lofted regions don't match segment for segment")]
    LoftMismatch,
    /// A loft's sketch planes face opposite ways along it, or it doesn't leave
    /// the first plane.
    #[error("lofted sketch planes must face the same way along the loft")]
    LoftFacing,
    /// Building NURBS geometry failed.
    #[error("NURBS construction failed: {0}")]
    Nurbs(String),
    #[error(transparent)]
    Topology(#[from] TopologyError),
}
