use core::f64::consts::TAU;

use nalgebra::Vector2;

/// A piece of a sketch profile, in sketch coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Segment {
    Line {
        start: Vector2<f64>,
        end: Vector2<f64>,
    },
    /// The arc of the circle about `center` from `start_angle`, turning
    /// `sweep` radians: anticlockwise if positive. A sweep of ±2π is a full
    /// circle.
    Arc {
        center: Vector2<f64>,
        radius: f64,
        start_angle: f64,
        sweep: f64,
    },
}

impl Segment {
    pub fn start(&self) -> Vector2<f64> {
        match *self {
            Self::Line { start, .. } => start,
            Self::Arc { .. } => self.at(0.0),
        }
    }

    pub fn end(&self) -> Vector2<f64> {
        match *self {
            Self::Line { end, .. } => end,
            Self::Arc { .. } => self.at(1.0),
        }
    }

    /// The point a fraction `s` of the way along.
    pub fn at(&self, s: f64) -> Vector2<f64> {
        match *self {
            Self::Line { start, end } => start + (end - start) * s,
            Self::Arc {
                center,
                radius,
                start_angle,
                sweep,
            } => {
                let angle = start_angle + sweep * s;
                center + Vector2::new(angle.cos(), angle.sin()) * radius
            }
        }
    }

    pub fn length(&self) -> f64 {
        match *self {
            Self::Line { start, end } => (end - start).norm(),
            Self::Arc { radius, sweep, .. } => radius * sweep.abs(),
        }
    }

    /// Whether this is a whole circle, which closes on itself.
    pub fn is_full_circle(&self) -> bool {
        matches!(*self, Self::Arc { sweep, .. } if (sweep.abs() - TAU).abs() <= 1e-12 * TAU)
    }

    /// The same segment traversed the other way.
    pub fn reversed(&self) -> Self {
        match *self {
            Self::Line { start, end } => Self::Line {
                start: end,
                end: start,
            },
            Self::Arc {
                center,
                radius,
                start_angle,
                sweep,
            } => Self::Arc {
                center,
                radius,
                start_angle: start_angle + sweep,
                sweep: -sweep,
            },
        }
    }

    /// This segment's term in the profile's signed area, `½∮(x dy − y dx)`.
    pub(crate) fn area_term(&self) -> f64 {
        match *self {
            Self::Line { start, end } => (start.x * end.y - end.x * start.y) / 2.0,
            Self::Arc {
                center,
                radius,
                start_angle: a,
                sweep: s,
            } => {
                let b = a + s;
                (radius * radius * s
                    + radius * (center.x * (b.sin() - a.sin()) - center.y * (b.cos() - a.cos())))
                    / 2.0
            }
        }
    }
}
