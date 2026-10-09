use core::f64::consts::PI;
use core::f64::consts::TAU;

use nalgebra::Vector2;

use crate::error::ModelError;
use crate::segment::Segment;

/// A closed chain of segments in a sketch: each starts where the last ended,
/// and the last ends where the first starts.
#[derive(Clone, Debug, PartialEq)]
pub struct Profile {
    segments: Vec<Segment>,
}

/// How close consecutive segments' ends must be, relative to the profile's size.
const JOIN_TOLERANCE: f64 = 1e-9;

impl Profile {
    /// Checks `segments` form a closed, non-degenerate chain enclosing some area.
    pub fn new(segments: Vec<Segment>) -> Result<Self, ModelError> {
        if segments.is_empty() {
            return Err(ModelError::EmptyProfile);
        }
        let size = segments
            .iter()
            .map(|s| s.start().norm().max(s.end().norm()) + s.length())
            .fold(0.0, f64::max);
        let tolerance = JOIN_TOLERANCE * (1.0 + size);
        for (index, segment) in segments.iter().enumerate() {
            let degenerate = match *segment {
                Segment::Line { .. } => segment.length() <= tolerance,
                Segment::Arc { radius, sweep, .. } => {
                    radius <= tolerance
                        || sweep.abs() * radius <= tolerance
                        || sweep.abs() > TAU * (1.0 + 1e-12)
                }
            };
            if degenerate {
                return Err(ModelError::DegenerateSegment { index });
            }
        }
        if segments.len() > 1 && segments.iter().any(Segment::is_full_circle) {
            return Err(ModelError::CircleInChain);
        }
        for (index, segment) in segments.iter().enumerate() {
            let next = &segments[(index + 1) % segments.len()];
            let gap = (segment.end() - next.start()).norm();
            if gap > tolerance {
                return Err(ModelError::OpenProfile { index, gap });
            }
        }
        let profile = Self { segments };
        if profile.signed_area().abs() <= tolerance * tolerance {
            return Err(ModelError::ZeroArea);
        }
        Ok(profile)
    }

    /// The polygon through `points`, closed back to the first.
    pub fn polygon(points: &[Vector2<f64>]) -> Result<Self, ModelError> {
        let segments = (0..points.len())
            .map(|i| Segment::Line {
                start: points[i],
                end: points[(i + 1) % points.len()],
            })
            .collect();
        Self::new(segments)
    }

    /// The axis-aligned rectangle with opposite corners `a` and `b`.
    pub fn rectangle(a: Vector2<f64>, b: Vector2<f64>) -> Result<Self, ModelError> {
        Self::polygon(&[a, Vector2::new(b.x, a.y), b, Vector2::new(a.x, b.y)])
    }

    /// The circle about `center`, starting and ending at angle π (its point
    /// furthest along −x).
    pub fn circle(center: Vector2<f64>, radius: f64) -> Result<Self, ModelError> {
        Self::new(vec![Segment::Arc {
            center,
            radius,
            start_angle: PI,
            sweep: TAU,
        }])
    }

    /// Starts drawing a profile at `start`.
    pub fn builder(start: Vector2<f64>) -> ProfileBuilder {
        ProfileBuilder {
            start,
            cursor: start,
            segments: Vec::new(),
            error: None,
        }
    }

    pub fn segments(&self) -> &[Segment] {
        &self.segments
    }

    /// The enclosed area: positive if the profile runs anticlockwise.
    pub fn signed_area(&self) -> f64 {
        self.segments.iter().map(Segment::area_term).sum()
    }

    /// The same profile traversed the other way.
    pub fn reversed(&self) -> Self {
        Self {
            segments: self.segments.iter().rev().map(Segment::reversed).collect(),
        }
    }
}

/// Draws a profile segment by segment, from a pen position.
#[derive(Clone, Debug)]
pub struct ProfileBuilder {
    start: Vector2<f64>,
    cursor: Vector2<f64>,
    segments: Vec<Segment>,
    /// The first mistake, reported by [`Self::close`].
    error: Option<ModelError>,
}

impl ProfileBuilder {
    /// A line from the pen to `end`.
    pub fn line_to(mut self, end: Vector2<f64>) -> Self {
        self.segments.push(Segment::Line {
            start: self.cursor,
            end,
        });
        self.cursor = end;
        self
    }

    /// An arc from the pen to `end` about `center`, anticlockwise if
    /// `anticlockwise`, else clockwise. `end` must be as far from `center` as the
    /// pen is; [`Self::close`] reports it if not.
    pub fn arc_to(mut self, end: Vector2<f64>, center: Vector2<f64>, anticlockwise: bool) -> Self {
        let from = self.cursor - center;
        let to = end - center;
        let distance = (to.norm() - from.norm()).abs();
        if distance > JOIN_TOLERANCE * (1.0 + from.norm()) && self.error.is_none() {
            self.error = Some(ModelError::ArcEndOffCircle { distance });
        }
        let start_angle = from.y.atan2(from.x);
        let mut sweep = (to.y.atan2(to.x) - start_angle).rem_euclid(TAU);
        if !anticlockwise {
            sweep -= TAU;
        }
        self.segments.push(Segment::Arc {
            center,
            radius: from.norm(),
            start_angle,
            sweep,
        });
        self.cursor = end;
        self
    }

    /// Finishes the profile, adding a line back to the start if the pen isn't
    /// there already.
    pub fn close(mut self) -> Result<Profile, ModelError> {
        if let Some(error) = self.error {
            return Err(error);
        }
        if (self.cursor - self.start).norm() > JOIN_TOLERANCE * (1.0 + self.start.norm()) {
            let start = self.start;
            self = self.line_to(start);
        }
        Profile::new(self.segments)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rectangle_area_and_orientation() {
        let r = Profile::rectangle(Vector2::new(1.0, 1.0), Vector2::new(4.0, 3.0)).unwrap();
        assert!((r.signed_area() - 6.0).abs() < 1e-12);
        assert!((r.reversed().signed_area() + 6.0).abs() < 1e-12);
    }

    #[test]
    fn circle_area() {
        let c = Profile::circle(Vector2::new(2.0, -1.0), 1.5).unwrap();
        assert!((c.signed_area() - PI * 2.25).abs() < 1e-12);
        assert!((c.reversed().signed_area() + PI * 2.25).abs() < 1e-12);
    }

    #[test]
    fn stadium_from_the_builder() {
        // A 4 × 2 rectangle with semicircular ends of radius 1.
        let p = Profile::builder(Vector2::new(0.0, -1.0))
            .line_to(Vector2::new(4.0, -1.0))
            .arc_to(Vector2::new(4.0, 1.0), Vector2::new(4.0, 0.0), true)
            .line_to(Vector2::new(0.0, 1.0))
            .arc_to(Vector2::new(0.0, -1.0), Vector2::new(0.0, 0.0), true)
            .close()
            .unwrap();
        assert_eq!(p.segments().len(), 4);
        assert!((p.signed_area() - (8.0 + PI)).abs() < 1e-12);
    }

    #[test]
    fn concave_arc_subtracts_area() {
        // A 2 × 2 square with a semicircular bite of radius 0.5 from its top.
        let p = Profile::builder(Vector2::new(0.0, 0.0))
            .line_to(Vector2::new(2.0, 0.0))
            .line_to(Vector2::new(2.0, 2.0))
            .line_to(Vector2::new(1.5, 2.0))
            .arc_to(Vector2::new(0.5, 2.0), Vector2::new(1.0, 2.0), false)
            .line_to(Vector2::new(0.0, 2.0))
            .close()
            .unwrap();
        assert!((p.signed_area() - (4.0 - PI / 8.0)).abs() < 1e-12);
    }

    #[test]
    fn invalid_profiles_are_rejected() {
        let open = vec![
            Segment::Line {
                start: Vector2::new(0.0, 0.0),
                end: Vector2::new(1.0, 0.0),
            },
            Segment::Line {
                start: Vector2::new(1.0, 0.0),
                end: Vector2::new(1.0, 1.0),
            },
        ];
        assert!(matches!(
            Profile::new(open),
            Err(ModelError::OpenProfile { index: 1, .. })
        ));
        assert_eq!(Profile::new(vec![]), Err(ModelError::EmptyProfile));
        assert_eq!(
            Profile::polygon(&[
                Vector2::new(0.0, 0.0),
                Vector2::new(1.0, 0.0),
                Vector2::new(2.0, 0.0)
            ]),
            Err(ModelError::ZeroArea)
        );
        let off = Profile::builder(Vector2::new(1.0, 0.0))
            .arc_to(Vector2::new(0.0, 2.0), Vector2::zeros(), true)
            .close();
        assert!(matches!(off, Err(ModelError::ArcEndOffCircle { .. })));
    }
}
