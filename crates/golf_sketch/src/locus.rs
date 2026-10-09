//! The full lines and circles segments lie on, and where they cross.

use nalgebra::Vector2;

use crate::segment::Segment;

/// An unbounded line or a whole circle.
#[derive(Clone, Copy, Debug)]
pub(crate) enum Locus {
    Line {
        point: Vector2<f64>,
        direction: Vector2<f64>,
    },
    Circle {
        center: Vector2<f64>,
        radius: f64,
    },
}

impl Locus {
    /// The line or circle `segment` lies on.
    pub(crate) fn of(segment: &Segment) -> Self {
        match *segment {
            Segment::Line { start, end } => Self::Line {
                point: start,
                direction: (end - start).normalize(),
            },
            Segment::Arc { center, radius, .. } => Self::Circle { center, radius },
        }
    }

    /// Where the two cross; a touching pair counts once.
    pub(crate) fn intersect(&self, other: &Self) -> Vec<Vector2<f64>> {
        match (*self, *other) {
            (
                Self::Line {
                    point: p,
                    direction: d,
                },
                Self::Line {
                    point: q,
                    direction: e,
                },
            ) => {
                let cross = d.x * e.y - d.y * e.x;
                if cross.abs() <= 1e-12 {
                    return Vec::new();
                }
                let w = q - p;
                let t = (w.x * e.y - w.y * e.x) / cross;
                vec![p + d * t]
            }
            (Self::Line { point, direction }, Self::Circle { center, radius })
            | (Self::Circle { center, radius }, Self::Line { point, direction }) => {
                // The foot of the centre on the line, and either side of it.
                let foot = point + direction * (center - point).dot(&direction);
                let h2 = radius * radius - (center - foot).norm_squared();
                touching(h2, radius).map_or(Vec::new(), |h| {
                    if h == 0.0 {
                        vec![foot]
                    } else {
                        vec![foot - direction * h, foot + direction * h]
                    }
                })
            }
            (
                Self::Circle {
                    center: c1,
                    radius: r1,
                },
                Self::Circle {
                    center: c2,
                    radius: r2,
                },
            ) => {
                let d = (c2 - c1).norm();
                if d <= 1e-12 * (1.0 + r1 + r2) {
                    return Vec::new();
                }
                let u = (c2 - c1) / d;
                let a = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
                let base = c1 + u * a;
                let perp = Vector2::new(-u.y, u.x);
                touching(r1 * r1 - a * a, r1.max(r2)).map_or(Vec::new(), |h| {
                    if h == 0.0 {
                        vec![base]
                    } else {
                        vec![base - perp * h, base + perp * h]
                    }
                })
            }
        }
    }
}

/// The half-chord `√h2`, treating a hair below zero (rounding at a tangency)
/// as zero, or `None` if the two miss.
fn touching(h2: f64, scale: f64) -> Option<f64> {
    let tolerance = 1e-12 * (1.0 + scale * scale);
    match h2 {
        h2 if h2 > tolerance => Some(h2.sqrt()),
        h2 if h2 >= -tolerance => Some(0.0),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crossings() {
        let line = Locus::Line {
            point: Vector2::new(0.0, 1.0),
            direction: Vector2::x(),
        };
        let circle = Locus::Circle {
            center: Vector2::zeros(),
            radius: 2.0,
        };
        let x = 3f64.sqrt();
        assert_eq!(
            line.intersect(&circle),
            vec![Vector2::new(-x, 1.0), Vector2::new(x, 1.0)]
        );
        let tangent = Locus::Line {
            point: Vector2::new(5.0, 2.0),
            direction: Vector2::x(),
        };
        assert_eq!(tangent.intersect(&circle), vec![Vector2::new(0.0, 2.0)]);
        let other = Locus::Circle {
            center: Vector2::new(2.0, 0.0),
            radius: 2.0,
        };
        let points = circle.intersect(&other);
        assert_eq!(points.len(), 2);
        for p in points {
            assert!(
                (p.norm() - 2.0).abs() < 1e-12
                    && ((p - Vector2::new(2.0, 0.0)).norm() - 2.0).abs() < 1e-12
            );
        }
        let vertical = Locus::Line {
            point: Vector2::new(1.0, 0.0),
            direction: Vector2::y(),
        };
        assert_eq!(line.intersect(&vertical), vec![Vector2::new(1.0, 1.0)]);
    }
}
