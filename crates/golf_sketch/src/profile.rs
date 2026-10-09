use core::f64::consts::PI;
use core::f64::consts::TAU;

use nalgebra::Vector2;

use crate::error::SketchError;
use crate::path::Path;
use crate::segment::Segment;

/// A closed chain of segments in a sketch: each starts where the last ended,
/// and the last ends where the first starts.
#[derive(Clone, Debug, PartialEq)]
pub struct Profile {
    segments: Vec<Segment>,
}

/// How close consecutive segments' ends must be, relative to the profile's size.
pub(crate) const JOIN_TOLERANCE: f64 = 1e-9;

impl Profile {
    /// Checks `segments` form a closed, non-degenerate chain enclosing some area.
    pub fn new(segments: Vec<Segment>) -> Result<Self, SketchError> {
        if segments.is_empty() {
            return Err(SketchError::EmptyProfile);
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
                return Err(SketchError::DegenerateSegment { index });
            }
        }
        if segments.len() > 1 && segments.iter().any(Segment::is_full_circle) {
            return Err(SketchError::CircleInChain);
        }
        for (index, segment) in segments.iter().enumerate() {
            let next = &segments[(index + 1) % segments.len()];
            let gap = (segment.end() - next.start()).norm();
            if gap > tolerance {
                return Err(SketchError::OpenProfile { index, gap });
            }
        }
        let profile = Self { segments };
        if profile.signed_area().abs() <= tolerance * tolerance {
            return Err(SketchError::ZeroArea);
        }
        Ok(profile)
    }

    /// The polygon through `points`, closed back to the first.
    pub fn polygon(points: &[Vector2<f64>]) -> Result<Self, SketchError> {
        let segments = (0..points.len())
            .map(|i| Segment::Line {
                start: points[i],
                end: points[(i + 1) % points.len()],
            })
            .collect();
        Self::new(segments)
    }

    /// The axis-aligned rectangle with opposite corners `a` and `b`.
    pub fn rectangle(a: Vector2<f64>, b: Vector2<f64>) -> Result<Self, SketchError> {
        Self::polygon(&[a, Vector2::new(b.x, a.y), b, Vector2::new(a.x, b.y)])
    }

    /// The circle about `center`, starting and ending at angle π (its point
    /// furthest along −x).
    pub fn circle(center: Vector2<f64>, radius: f64) -> Result<Self, SketchError> {
        Self::new(vec![Segment::Arc {
            center,
            radius,
            start_angle: PI,
            sweep: TAU,
        }])
    }

    /// The regular polygon with `sides` corners on the circle about `center`
    /// of `circumradius`, the first at angle `rotation`, anticlockwise.
    pub fn regular_polygon(
        center: Vector2<f64>,
        circumradius: f64,
        sides: usize,
        rotation: f64,
    ) -> Result<Self, SketchError> {
        if sides < 3 {
            return Err(SketchError::InvalidShape {
                reason: "a polygon needs at least 3 sides",
            });
        }
        if circumradius <= 0.0 {
            return Err(SketchError::InvalidShape {
                reason: "a polygon's radius must be positive",
            });
        }
        let points: Vec<_> = (0..sides)
            .map(|i| {
                let angle = rotation + TAU * i as f64 / sides as f64;
                center + Vector2::new(angle.cos(), angle.sin()) * circumradius
            })
            .collect();
        Self::polygon(&points)
    }

    /// The axis-aligned rectangle with opposite corners `a` and `b`, its
    /// corners rounded to `radius` (a plain rectangle if 0). The radius may be
    /// up to half the shorter side.
    pub fn rounded_rectangle(
        a: Vector2<f64>,
        b: Vector2<f64>,
        radius: f64,
    ) -> Result<Self, SketchError> {
        let (lo, hi) = (a.inf(&b), a.sup(&b));
        let size = hi - lo;
        if radius < 0.0 || radius > size.x.min(size.y) / 2.0 * (1.0 + 1e-12) {
            return Err(SketchError::InvalidShape {
                reason: "corner radius must be between 0 and half the shorter side",
            });
        }
        if radius == 0.0 {
            return Self::rectangle(lo, hi);
        }
        let r = radius;
        let corners = [
            (Vector2::new(hi.x - r, lo.y + r), -PI / 2.0),
            (Vector2::new(hi.x - r, hi.y - r), 0.0),
            (Vector2::new(lo.x + r, hi.y - r), PI / 2.0),
            (Vector2::new(lo.x + r, lo.y + r), PI),
        ];
        let tolerance = JOIN_TOLERANCE * (1.0 + hi.norm().max(lo.norm()));
        let mut segments = Vec::new();
        for (i, &(center, start_angle)) in corners.iter().enumerate() {
            let arc = Segment::Arc {
                center,
                radius: r,
                start_angle,
                sweep: PI / 2.0,
            };
            // The straight side leading into this corner, unless the corners
            // meet (radius half the side).
            let (previous, _) = corners[(i + 3) % 4];
            let previous_end = previous + Vector2::new(start_angle.cos(), start_angle.sin()) * r;
            if (arc.start() - previous_end).norm() > tolerance {
                segments.push(Segment::Line {
                    start: previous_end,
                    end: arc.start(),
                });
            }
            segments.push(arc);
        }
        Self::new(segments)
    }

    /// A slot (stadium): the points within `radius` of the segment from
    /// `center_a` to `center_b`. A circle if the two coincide.
    pub fn slot(
        center_a: Vector2<f64>,
        center_b: Vector2<f64>,
        radius: f64,
    ) -> Result<Self, SketchError> {
        if radius <= 0.0 {
            return Err(SketchError::InvalidShape {
                reason: "a slot's radius must be positive",
            });
        }
        let length = (center_b - center_a).norm();
        if length == 0.0 {
            return Self::circle(center_a, radius);
        }
        let d = (center_b - center_a) / length;
        let angle = d.y.atan2(d.x);
        let n = Vector2::new(-d.y, d.x) * radius;
        Self::new(vec![
            Segment::Line {
                start: center_a - n,
                end: center_b - n,
            },
            Segment::Arc {
                center: center_b,
                radius,
                start_angle: angle - PI / 2.0,
                sweep: PI,
            },
            Segment::Line {
                start: center_b + n,
                end: center_a + n,
            },
            Segment::Arc {
                center: center_a,
                radius,
                start_angle: angle + PI / 2.0,
                sweep: PI,
            },
        ])
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

    pub fn translated(&self, v: Vector2<f64>) -> Self {
        self.map(|s| s.translated(v))
    }

    /// Rotated anticlockwise by `angle` about `about`.
    pub fn rotated(&self, angle: f64, about: Vector2<f64>) -> Self {
        self.map(|s| s.rotated(angle, about))
    }

    /// Scaled by `factor` about `about`.
    ///
    /// # Panics
    /// If `factor` isn't positive.
    pub fn scaled(&self, factor: f64, about: Vector2<f64>) -> Self {
        self.map(|s| s.scaled(factor, about))
    }

    /// Reflected in the line through `point` along `direction`, which reverses
    /// its orientation.
    pub fn mirrored(&self, point: Vector2<f64>, direction: Vector2<f64>) -> Self {
        self.map(|s| s.mirrored(point, direction))
    }

    fn map(&self, f: impl Fn(&Segment) -> Segment) -> Self {
        Self {
            segments: self.segments.iter().map(f).collect(),
        }
    }

    /// The indices of corners where segments meet at an angle: corner `i` is
    /// where segment `i` starts.
    pub fn sharp_corners(&self) -> Vec<usize> {
        (0..self.segments.len())
            .filter(|&i| !self.is_smooth(i))
            .collect()
    }

    /// Whether the segments meeting at corner `i` are tangent there.
    pub fn is_smooth(&self, i: usize) -> bool {
        let (incoming, outgoing) = self.at_corner(i);
        let (a, b) = (incoming.tangent_at(1.0), outgoing.tangent_at(0.0));
        (a.x * b.y - a.y * b.x).abs() <= 1e-9 && a.dot(&b) > 0.0
    }

    /// The segments arriving at and leaving corner `i`.
    pub(crate) fn at_corner(&self, i: usize) -> (&Segment, &Segment) {
        let n = self.segments.len();
        (&self.segments[(i + n - 1) % n], &self.segments[i])
    }
}

/// Draws a profile segment by segment, from a pen position.
#[derive(Clone, Debug)]
pub struct ProfileBuilder {
    start: Vector2<f64>,
    cursor: Vector2<f64>,
    segments: Vec<Segment>,
    /// The first mistake, reported by [`Self::close`].
    error: Option<SketchError>,
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
            self.error = Some(SketchError::ArcEndOffCircle { distance });
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

    /// An arc from the pen through `through` to `end`.
    pub fn arc_through(self, through: Vector2<f64>, end: Vector2<f64>) -> Self {
        let (p, q, e) = (self.cursor, through, end);
        let cross = (q - p).x * (e - p).y - (q - p).y * (e - p).x;
        let size = (q - p).norm().max((e - p).norm());
        if cross.abs() <= 1e-12 * size * size {
            return self.fail(SketchError::InvalidShape {
                reason: "an arc through three collinear points",
            });
        }
        // The circumcentre.
        let (b, c) = (q - p, e - p);
        let d = 2.0 * cross;
        let center = p + Vector2::new(
            c.y * b.norm_squared() - b.y * c.norm_squared(),
            b.x * c.norm_squared() - c.x * b.norm_squared(),
        ) / d;
        // Points met anticlockwise round a circle turn anticlockwise.
        self.arc_to(end, center, cross > 0.0)
    }

    /// An arc from the pen to `end`, leaving tangent to the previous segment.
    pub fn tangent_arc_to(self, end: Vector2<f64>) -> Self {
        let Some(last) = self.segments.last() else {
            return self.fail(SketchError::InvalidShape {
                reason: "a tangent arc needs a segment before it",
            });
        };
        let t = last.tangent_at(1.0);
        let n = Vector2::new(-t.y, t.x);
        let w = end - self.cursor;
        let along = 2.0 * n.dot(&w);
        if along.abs() <= 1e-12 * w.norm() {
            return self.fail(SketchError::InvalidShape {
                reason: "a tangent arc to a point straight ahead is a line",
            });
        }
        // The centre is on the normal, as far from the pen as from `end`.
        let rho = w.norm_squared() / along;
        let center = self.cursor + n * rho;
        self.arc_to(end, center, rho > 0.0)
    }

    fn fail(mut self, error: SketchError) -> Self {
        if self.error.is_none() {
            self.error = Some(error);
        }
        self
    }

    /// Finishes the profile, adding a line back to the start if the pen isn't
    /// there already.
    pub fn close(mut self) -> Result<Profile, SketchError> {
        if let Some(error) = self.error {
            return Err(error);
        }
        if (self.cursor - self.start).norm() > JOIN_TOLERANCE * (1.0 + self.start.norm()) {
            let start = self.start;
            self = self.line_to(start);
        }
        Profile::new(self.segments)
    }

    /// Finishes an open path where the pen is, without closing it.
    pub fn path(self) -> Result<Path, SketchError> {
        if let Some(error) = self.error {
            return Err(error);
        }
        Path::new(self.segments)
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
    fn constructors() {
        let hexagon = Profile::regular_polygon(Vector2::new(1.0, 1.0), 2.0, 6, 0.3).unwrap();
        assert!((hexagon.signed_area() - 1.5 * 3f64.sqrt() * 4.0).abs() < 1e-12);
        assert!(Profile::regular_polygon(Vector2::zeros(), 1.0, 2, 0.0).is_err());

        let (w, h, r) = (4.0, 3.0, 0.5);
        let rounded = Profile::rounded_rectangle(Vector2::new(w, h), Vector2::zeros(), r).unwrap();
        assert_eq!(rounded.segments().len(), 8);
        assert!((rounded.signed_area() - (w * h - (4.0 - PI) * r * r)).abs() < 1e-12);
        // At half the shorter side the short sides vanish: a slot.
        let full = Profile::rounded_rectangle(Vector2::zeros(), Vector2::new(w, h), 1.5).unwrap();
        assert_eq!(full.segments().len(), 6);
        assert!((full.signed_area() - (w * h - (4.0 - PI) * 2.25)).abs() < 1e-12);
        assert!(Profile::rounded_rectangle(Vector2::zeros(), Vector2::new(w, h), 1.6).is_err());
        assert_eq!(
            Profile::rounded_rectangle(Vector2::zeros(), Vector2::new(w, h), 0.0)
                .unwrap()
                .segments()
                .len(),
            4
        );

        let slot = Profile::slot(Vector2::new(1.0, 1.0), Vector2::new(4.0, 5.0), 0.5).unwrap();
        assert!((slot.signed_area() - (2.0 * 0.5 * 5.0 + PI * 0.25)).abs() < 1e-12);
        assert!(slot.sharp_corners().is_empty());
    }

    #[test]
    fn three_point_and_tangent_arcs() {
        // A semicircle through the top of the unit circle, closed by a line.
        let half = Profile::builder(Vector2::new(1.0, 0.0))
            .arc_through(Vector2::new(0.0, 1.0), Vector2::new(-1.0, 0.0))
            .close()
            .unwrap();
        assert!((half.signed_area() - PI / 2.0).abs() < 1e-12);
        // Clockwise through the bottom: the same area, negative.
        let below = Profile::builder(Vector2::new(-1.0, 0.0))
            .arc_through(Vector2::new(0.0, -1.0), Vector2::new(1.0, 0.0))
            .close()
            .unwrap();
        assert!((below.signed_area() - PI / 2.0).abs() < 1e-12);

        // A line, then a tangent arc turning back: a D of a 2-long side and a
        // semicircle of radius 1.
        let d = Profile::builder(Vector2::new(0.0, -1.0))
            .line_to(Vector2::new(2.0, -1.0))
            .tangent_arc_to(Vector2::new(2.0, 1.0))
            .line_to(Vector2::new(0.0, 1.0))
            .close()
            .unwrap();
        assert!((d.signed_area() - (4.0 + PI / 2.0)).abs() < 1e-12);
        assert!(d.is_smooth(1) && d.is_smooth(2));

        let collinear = Profile::builder(Vector2::zeros())
            .arc_through(Vector2::new(1.0, 0.0), Vector2::new(2.0, 0.0))
            .close();
        assert!(matches!(collinear, Err(SketchError::InvalidShape { .. })));
        assert!(matches!(
            Profile::builder(Vector2::zeros())
                .tangent_arc_to(Vector2::x())
                .close(),
            Err(SketchError::InvalidShape { .. })
        ));
    }

    #[test]
    fn transforms() {
        let rounded =
            Profile::rounded_rectangle(Vector2::zeros(), Vector2::new(4.0, 3.0), 0.5).unwrap();
        let area = rounded.signed_area();
        let about = Vector2::new(1.0, -1.0);
        assert!((rounded.rotated(1.1, about).signed_area() - area).abs() < 1e-12);
        assert!((rounded.translated(about).signed_area() - area).abs() < 1e-12);
        assert!((rounded.scaled(3.0, about).signed_area() - 9.0 * area).abs() < 1e-11);
        let mirrored = rounded.mirrored(about, Vector2::new(1.0, 1.0));
        assert!((mirrored.signed_area() + area).abs() < 1e-12);
        for p in [rounded.rotated(1.1, about), mirrored] {
            assert_eq!(Profile::new(p.segments().to_vec()).unwrap(), p);
        }
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
            Err(SketchError::OpenProfile { index: 1, .. })
        ));
        assert_eq!(Profile::new(vec![]), Err(SketchError::EmptyProfile));
        assert_eq!(
            Profile::polygon(&[
                Vector2::new(0.0, 0.0),
                Vector2::new(1.0, 0.0),
                Vector2::new(2.0, 0.0)
            ]),
            Err(SketchError::ZeroArea)
        );
        let off = Profile::builder(Vector2::new(1.0, 0.0))
            .arc_to(Vector2::new(0.0, 2.0), Vector2::zeros(), true)
            .close();
        assert!(matches!(off, Err(SketchError::ArcEndOffCircle { .. })));
    }
}
