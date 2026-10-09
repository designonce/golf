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
}
