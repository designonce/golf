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

    /// The unit tangent a fraction `s` of the way along, in the direction of
    /// travel.
    pub fn tangent_at(&self, s: f64) -> Vector2<f64> {
        match *self {
            Self::Line { start, end } => (end - start).normalize(),
            Self::Arc {
                start_angle, sweep, ..
            } => {
                let angle = start_angle + sweep * s;
                Vector2::new(-angle.sin(), angle.cos()) * sweep.signum()
            }
        }
    }

    /// The point a distance `d` along (arc length for arcs).
    pub fn point_at_length(&self, d: f64) -> Vector2<f64> {
        self.at(d / self.length())
    }

    /// The part between fractions `s0` and `s1` of the way along.
    pub fn between(&self, s0: f64, s1: f64) -> Self {
        if s0 == 0.0 && s1 == 1.0 {
            return *self;
        }
        match *self {
            Self::Line { .. } => Self::Line {
                start: self.at(s0),
                end: self.at(s1),
            },
            Self::Arc {
                center,
                radius,
                start_angle,
                sweep,
            } => Self::Arc {
                center,
                radius,
                start_angle: start_angle + sweep * s0,
                sweep: sweep * (s1 - s0),
            },
        }
    }

    /// How far along (as a fraction) the point of the underlying line or
    /// circle nearest `p` is: outside `[0, 1]` if it's off the segment. For an
    /// arc, angles are measured forwards from the start, so a point just
    /// before the start reads as just below 0.
    pub(crate) fn fraction_of(&self, p: Vector2<f64>) -> f64 {
        match *self {
            Self::Line { start, end } => {
                (p - start).dot(&(end - start)) / (end - start).norm_squared()
            }
            Self::Arc {
                center,
                start_angle,
                sweep,
                ..
            } => {
                let angle = (p.y - center.y).atan2(p.x - center.x);
                let mut ahead = ((angle - start_angle) * sweep.signum()).rem_euclid(TAU);
                if ahead > TAU - 1e-9 {
                    ahead -= TAU;
                }
                ahead / sweep.abs()
            }
        }
    }

    /// The point of the underlying line or circle nearest `p`.
    pub(crate) fn foot(&self, p: Vector2<f64>) -> Vector2<f64> {
        match *self {
            Self::Line { .. } => self.at(self.fraction_of(p)),
            Self::Arc { center, radius, .. } => center + (p - center).normalize() * radius,
        }
    }

    /// The segment moved `k` to its left (right if negative): lines shift along
    /// their normal, arcs change radius. `None` if an arc's radius would vanish.
    pub(crate) fn left_offset(&self, k: f64) -> Option<Self> {
        match *self {
            Self::Line { start, end } => {
                let direction = (end - start).normalize();
                let shift = Vector2::new(-direction.y, direction.x) * k;
                Some(Self::Line {
                    start: start + shift,
                    end: end + shift,
                })
            }
            Self::Arc {
                center,
                radius,
                start_angle,
                sweep,
            } => {
                // An anticlockwise arc's left is towards its centre.
                let radius = radius - k * sweep.signum();
                (radius > 1e-12 * (1.0 + center.norm())).then_some(Self::Arc {
                    center,
                    radius,
                    start_angle,
                    sweep,
                })
            }
        }
    }

    pub fn translated(&self, v: Vector2<f64>) -> Self {
        match *self {
            Self::Line { start, end } => Self::Line {
                start: start + v,
                end: end + v,
            },
            Self::Arc {
                center,
                radius,
                start_angle,
                sweep,
            } => Self::Arc {
                center: center + v,
                radius,
                start_angle,
                sweep,
            },
        }
    }

    /// Rotated anticlockwise by `angle` about `about`.
    pub fn rotated(&self, angle: f64, about: Vector2<f64>) -> Self {
        let (sin, cos) = angle.sin_cos();
        let turn = |p: Vector2<f64>| {
            let d = p - about;
            about + Vector2::new(cos * d.x - sin * d.y, sin * d.x + cos * d.y)
        };
        match *self {
            Self::Line { start, end } => Self::Line {
                start: turn(start),
                end: turn(end),
            },
            Self::Arc {
                center,
                radius,
                start_angle,
                sweep,
            } => Self::Arc {
                center: turn(center),
                radius,
                start_angle: start_angle + angle,
                sweep,
            },
        }
    }

    /// Scaled by `factor` about `about`.
    ///
    /// # Panics
    /// If `factor` isn't positive.
    pub fn scaled(&self, factor: f64, about: Vector2<f64>) -> Self {
        assert!(factor > 0.0, "scale factor {factor} must be positive");
        let scale = |p: Vector2<f64>| about + (p - about) * factor;
        match *self {
            Self::Line { start, end } => Self::Line {
                start: scale(start),
                end: scale(end),
            },
            Self::Arc {
                center,
                radius,
                start_angle,
                sweep,
            } => Self::Arc {
                center: scale(center),
                radius: radius * factor,
                start_angle,
                sweep,
            },
        }
    }

    /// Reflected in the line through `point` along `direction`. The segment
    /// still runs from (the image of) its start to its end, so arcs turn the
    /// other way.
    pub fn mirrored(&self, point: Vector2<f64>, direction: Vector2<f64>) -> Self {
        let u = direction.normalize();
        let reflect = |p: Vector2<f64>| {
            let d = p - point;
            point + u * (2.0 * d.dot(&u)) - d
        };
        match *self {
            Self::Line { start, end } => Self::Line {
                start: reflect(start),
                end: reflect(end),
            },
            Self::Arc {
                center,
                radius,
                start_angle,
                sweep,
            } => Self::Arc {
                center: reflect(center),
                radius,
                start_angle: 2.0 * u.y.atan2(u.x) - start_angle,
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

#[cfg(test)]
mod tests {
    use core::f64::consts::FRAC_PI_2;
    use core::f64::consts::PI;

    use super::*;

    fn arc() -> Segment {
        Segment::Arc {
            center: Vector2::new(1.0, 2.0),
            radius: 2.0,
            start_angle: 0.0,
            sweep: -FRAC_PI_2,
        }
    }

    #[test]
    fn tangents_follow_the_direction_of_travel() {
        let line = Segment::Line {
            start: Vector2::zeros(),
            end: Vector2::new(3.0, 4.0),
        };
        assert_eq!(line.tangent_at(0.5), Vector2::new(0.6, 0.8));
        // Clockwise from angle 0: heading down.
        assert!((arc().tangent_at(0.0) - Vector2::new(0.0, -1.0)).norm() < 1e-15);
        assert!((arc().tangent_at(1.0) - Vector2::new(-1.0, 0.0)).norm() < 1e-15);
        assert!((arc().point_at_length(PI / 2.0) - arc().at(0.5)).norm() < 1e-15);
    }

    #[test]
    fn fractions_invert_at() {
        for s in [0.0, 0.3, 1.0] {
            assert!((arc().fraction_of(arc().at(s)) - s).abs() < 1e-12);
        }
        assert!(arc().fraction_of(Vector2::new(1.0, 4.0)) > 1.0);
        let piece = arc().between(0.25, 0.75);
        assert!((piece.start() - arc().at(0.25)).norm() < 1e-15);
        assert!((piece.end() - arc().at(0.75)).norm() < 1e-15);
    }

    #[test]
    fn transforms_move_the_ends() {
        let about = Vector2::new(-1.0, 0.5);
        let rotated = arc().rotated(0.7, about);
        let turn = |p: Vector2<f64>| {
            let d = p - about;
            about
                + Vector2::new(
                    0.7f64.cos() * d.x - 0.7f64.sin() * d.y,
                    0.7f64.sin() * d.x + 0.7f64.cos() * d.y,
                )
        };
        assert!((rotated.start() - turn(arc().start())).norm() < 1e-14);
        assert!((rotated.end() - turn(arc().end())).norm() < 1e-14);
        let scaled = arc().scaled(3.0, about);
        assert!((scaled.end() - (about + (arc().end() - about) * 3.0)).norm() < 1e-14);
        // Mirrored in the y axis: x negated, still from start to end.
        let mirrored = arc().mirrored(Vector2::zeros(), Vector2::y());
        let flip = |p: Vector2<f64>| Vector2::new(-p.x, p.y);
        assert!((mirrored.start() - flip(arc().start())).norm() < 1e-14);
        assert!((mirrored.end() - flip(arc().end())).norm() < 1e-14);
        assert!((mirrored.at(0.5) - flip(arc().at(0.5))).norm() < 1e-14);
    }
}
