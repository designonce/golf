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
    /// Revolving needs the profile on one side of the axis (sketch x >= 0).
    #[error("segment {index} crosses the axis of revolution")]
    CrossesAxis { index: usize },
    /// An arc whose circle reaches the axis would revolve to a self-intersecting
    /// (spindle) torus.
    #[error("segment {index} would revolve to a spindle torus")]
    SpindleTorus { index: usize },
    #[error(transparent)]
    Topology(#[from] TopologyError),
}
