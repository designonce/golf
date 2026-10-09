/// Why a sketch is invalid.
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum SketchError {
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
    /// A shape's dimensions don't describe it (a negative radius, too few
    /// sides, three collinear points for an arc...).
    #[error("invalid shape: {reason}")]
    InvalidShape { reason: &'static str },
    /// The profile has no corner with this index.
    #[error("profile has no corner {corner}")]
    NoSuchCorner { corner: usize },
    /// The corner is already smooth (its segments meet tangentially), so has
    /// nothing to fillet or chamfer.
    #[error("corner {corner} is smooth")]
    SmoothCorner { corner: usize },
    /// No fillet of the radius touches both segments at the corner.
    #[error("no fillet fits corner {corner}")]
    NoFillet { corner: usize },
    /// The fillet or chamfer at the corner would consume a whole neighbouring
    /// segment, or overlap the treatment of the next corner.
    #[error("the treatment at corner {corner} is too large for its segments")]
    TooLarge { corner: usize },
    /// Offsetting would shrink segment `index` to nothing or turn it inside out.
    #[error("offsetting collapses segment {index}")]
    OffsetCollapses { index: usize },
}
