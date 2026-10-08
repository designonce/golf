/// Why a NURBS could not be built or modified.
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum NurbsError {
    #[error("degree must be at least 1")]
    ZeroDegree,
    /// A degree-`p` B-spline needs at least `2(p + 1)` knots.
    #[error("degree {degree} needs at least {} knots, found {knots}", 2 * (degree + 1))]
    TooFewKnots { degree: usize, knots: usize },
    #[error("knots must be finite and non-decreasing")]
    InvalidKnots,
    /// A knot inside the domain is repeated more than `degree` times, or one at
    /// an end more than `degree + 1` times.
    #[error("knot {knot} is repeated {multiplicity} times, too many for the degree")]
    MultiplicityTooHigh { knot: f64, multiplicity: usize },
    /// The knots leave no parameter range to evaluate over.
    #[error("knots leave an empty parameter domain")]
    EmptyDomain,
    /// The number of control points (or weights) doesn't match the knots.
    #[error("expected {expected} control points, found {found}")]
    ControlPointCount { expected: usize, found: usize },
    #[error("weight {index} is {weight}; weights must be finite and positive")]
    InvalidWeight { index: usize, weight: f64 },
    /// The parameter must lie strictly inside the domain.
    #[error("parameter {parameter} is not strictly inside the domain")]
    ParameterOutsideDomain { parameter: f64 },
}
