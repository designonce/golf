/// How closely a mesh must follow the geometry.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Tolerance {
    /// Largest distance allowed between a segment or facet and the true curve
    /// or surface (the sag).
    pub chord: f64,
    /// Largest angle, in radians, a curve's tangent or a surface's normal may
    /// turn across one segment or facet.
    pub angle: f64,
    /// Longest segment or facet edge allowed, if limited.
    pub max_length: Option<f64>,
}

impl Tolerance {
    pub fn new(chord: f64, angle: f64) -> Self {
        Self {
            chord,
            angle,
            max_length: None,
        }
    }

    pub fn with_max_length(self, max_length: f64) -> Self {
        Self {
            max_length: Some(max_length),
            ..self
        }
    }
}
