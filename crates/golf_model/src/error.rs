use golf_brep::TopologyError;

/// Why a sketch or operation couldn't produce a body.
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum ModelError {
    #[error("a profile needs at least one segment")]
    EmptyProfile,
    /// A line of zero length, or an arc of zero radius or sweep.
    #[error("segment {index} is degenerate")]
    DegenerateSegment { index: usize },
    /// Segment `index` doesn't end where the next one starts.
    #[error("segment {index} ends {gap} from where the next begins")]
    OpenProfile { index: usize, gap: f64 },
    /// A full circle is a profile on its own; it can't join other segments.
    #[error("a full circle must be a profile's only segment")]
    CircleInChain,
    /// An arc's end point isn't on the circle its start and centre define.
    #[error("arc end is {distance} off its circle")]
    ArcEndOffCircle { distance: f64 },
    /// The profile encloses no area.
    #[error("profile encloses no area")]
    ZeroArea,
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
