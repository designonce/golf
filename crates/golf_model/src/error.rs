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
    /// Edge treatments and drafts follow the profile round, so a sharp corner
    /// must be between two lines.
    #[error("corner {index} is sharp and not between two lines; fillet it in the sketch first")]
    UnsupportedCorner { index: usize },
    /// Insetting the profile (for a draft, chamfer or fillet) shrinks a segment
    /// to nothing or turns it inside out.
    #[error("insetting the profile collapses segment {index}")]
    InsetCollapses { index: usize },
    /// A fillet round a convex arc must be no larger than the arc's radius, or
    /// exactly equal (a spherical corner), not in between.
    #[error("the fillet is too large for arc {index}")]
    FilletTooLarge { index: usize },
    /// A pipe's arcs must be wider than the pipe, and its lines long enough
    /// for the mitres at their ends.
    #[error("the pipe is too thick for path segment {index}")]
    PipeTooThick { index: usize },
    /// A chamfer or fillet needs positive, finite sizes; a draft must be less
    /// than a right angle.
    #[error("invalid edge treatment or draft")]
    BadTreatment,
    /// The end treatments reach further than the extrusion is tall.
    #[error("the end treatments are taller than the extrusion")]
    TreatmentsTooTall,
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
