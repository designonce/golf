/// Why something couldn't be written as STEP.
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum StepError {
    /// Geometry STEP has no entity for, or that this writer doesn't know.
    #[error("{0} can't be written as STEP")]
    Unsupported(String),
    #[error("{0} is not finite and can't be written as STEP")]
    NonFinite(f64),
    #[error("a direction of zero length can't be written as STEP")]
    ZeroDirection,
    /// STEP faces need at least one bound, so a face covering a whole closed
    /// surface must be cut open by a seam first.
    #[error("face {face} has no bounds; STEP faces need at least one")]
    UnboundedFace { face: String },
    #[error("body {0} has no faces")]
    EmptyBody(String),
}
