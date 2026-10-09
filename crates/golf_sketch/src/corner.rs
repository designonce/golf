//! Fillets and chamfers at a profile's corners.

use core::f64::consts::TAU;

use nalgebra::Vector2;

use crate::error::SketchError;
use crate::locus::Locus;
use crate::profile::Profile;
use crate::segment::Segment;

/// How close to its ends a trimmed segment may be cut before it counts as
/// consumed.
const TRIM_TOLERANCE: f64 = 1e-9;

/// What replaces a corner: the fractions its neighbours are cut back to, and
/// the segment joining the cuts.
struct Treatment {
    /// Where the incoming segment now ends.
    incoming_end: f64,
    /// Where the outgoing segment now starts.
    outgoing_start: f64,
    piece: Segment,
}

impl Profile {
    /// Rounds the named corners with tangent arcs of `radius`, trimming the
    /// segments either side. Corner `i` is where segment `i` starts.
    pub fn fillet(&self, corners: &[usize], radius: f64) -> Result<Self, SketchError> {
        // Also rejects NaN.
        if radius.is_nan() || radius <= 0.0 {
            return Err(SketchError::InvalidShape {
                reason: "a fillet's radius must be positive",
            });
        }
        self.check_corners(corners)?;
        self.treat(corners, |corner, incoming, outgoing| {
            fillet(corner, incoming, outgoing, radius)
        })
    }

    /// Rounds every sharp corner with a tangent arc of `radius`.
    pub fn fillet_all(&self, radius: f64) -> Result<Self, SketchError> {
        self.fillet(&self.sharp_corners(), radius)
    }

    /// Cuts the named corners with a line `distance` back along each segment
    /// (arc length for arcs).
    pub fn chamfer(&self, corners: &[usize], distance: f64) -> Result<Self, SketchError> {
        self.chamfer_asymmetric(corners, distance, distance)
    }

    /// Cuts every sharp corner with a line `distance` back along each segment.
    pub fn chamfer_all(&self, distance: f64) -> Result<Self, SketchError> {
        self.chamfer(&self.sharp_corners(), distance)
    }

    /// Cuts the named corners with a line from `incoming` back along the
    /// segment arriving at the corner to `outgoing` along the one leaving it.
    pub fn chamfer_asymmetric(
        &self,
        corners: &[usize],
        incoming: f64,
        outgoing: f64,
    ) -> Result<Self, SketchError> {
        if incoming.is_nan() || outgoing.is_nan() || incoming <= 0.0 || outgoing <= 0.0 {
            return Err(SketchError::InvalidShape {
                reason: "chamfer distances must be positive",
            });
        }
        self.check_corners(corners)?;
        self.treat(corners, |corner, a, b| {
            let incoming_end = 1.0 - incoming / a.length();
            let outgoing_start = outgoing / b.length();
            if incoming_end <= TRIM_TOLERANCE || outgoing_start >= 1.0 - TRIM_TOLERANCE {
                return Err(SketchError::TooLarge { corner });
            }
            Ok(Treatment {
                incoming_end,
                outgoing_start,
                piece: Segment::Line {
                    start: a.at(incoming_end),
                    end: b.at(outgoing_start),
                },
            })
        })
    }

    /// Every named corner exists and is sharp.
    fn check_corners(&self, corners: &[usize]) -> Result<(), SketchError> {
        for &corner in corners {
            if corner >= self.segments().len() {
                return Err(SketchError::NoSuchCorner { corner });
            }
            if self.is_smooth(corner) {
                return Err(SketchError::SmoothCorner { corner });
            }
        }
        Ok(())
    }

    /// Replaces each named corner by what `treatment` makes of it, trimming the
    /// segments either side.
    fn treat(
        &self,
        corners: &[usize],
        treatment: impl Fn(usize, &Segment, &Segment) -> Result<Treatment, SketchError>,
    ) -> Result<Self, SketchError> {
        let segments = self.segments();
        let n = segments.len();
        let mut starts = vec![0.0; n];
        let mut ends = vec![1.0; n];
        let mut pieces: Vec<Option<Segment>> = vec![None; n];
        // Which corner trimmed each end, to blame when a segment is consumed.
        let mut trimmed_by = vec![(None, None); n];
        for &corner in corners {
            let incoming = (corner + n - 1) % n;
            let t = treatment(corner, &segments[incoming], &segments[corner])?;
            ends[incoming] = t.incoming_end;
            starts[corner] = t.outgoing_start;
            trimmed_by[incoming].1 = Some(corner);
            trimmed_by[corner].0 = Some(corner);
            pieces[corner] = Some(t.piece);
        }
        let mut out = Vec::with_capacity(n + corners.len());
        for i in 0..n {
            if starts[i] >= ends[i] - TRIM_TOLERANCE {
                let corner = trimmed_by[i]
                    .1
                    .or(trimmed_by[i].0)
                    .expect("only trimmed segments shrink");
                return Err(SketchError::TooLarge { corner });
            }
            out.extend(pieces[i]);
            out.push(segments[i].between(starts[i], ends[i]));
        }
        Profile::new(out)
    }
}

/// The fillet of `radius` at `corner`, between `a` (arriving) and `b`
/// (leaving).
///
/// Its centre is `radius` from both segments, on the side the path turns
/// towards (inside the corner), so it lies on both segments' offset curves on
/// that side; of their crossings, the nearest the corner whose feet fall on
/// both segments is taken.
fn fillet(corner: usize, a: &Segment, b: &Segment, radius: f64) -> Result<Treatment, SketchError> {
    let (ta, tb) = (a.tangent_at(1.0), b.tangent_at(0.0));
    let turn = (ta.x * tb.y - ta.y * tb.x).signum();
    let no_fillet = SketchError::NoFillet { corner };
    let (Some(a_offset), Some(b_offset)) =
        (a.left_offset(turn * radius), b.left_offset(turn * radius))
    else {
        return Err(no_fillet);
    };
    let at = a.end();
    let mut centres = Locus::of(&a_offset).intersect(&Locus::of(&b_offset));
    centres.sort_by(|p, q| (p - at).norm().total_cmp(&(q - at).norm()));
    for centre in centres {
        let (foot_a, foot_b) = (a.foot(centre), b.foot(centre));
        let (incoming_end, outgoing_start) = (a.fraction_of(foot_a), b.fraction_of(foot_b));
        let on_a = incoming_end > -TRIM_TOLERANCE && incoming_end < 1.0 + TRIM_TOLERANCE;
        let on_b = outgoing_start > -TRIM_TOLERANCE && outgoing_start < 1.0 + TRIM_TOLERANCE;
        if !(on_a && on_b) {
            continue;
        }
        if incoming_end <= TRIM_TOLERANCE || outgoing_start >= 1.0 - TRIM_TOLERANCE {
            return Err(SketchError::TooLarge { corner });
        }
        let start_angle = angle(foot_a - centre);
        // Round the centre the way the path turns.
        let raw = angle(foot_b - centre) - start_angle;
        let sweep = match turn > 0.0 {
            true => raw.rem_euclid(TAU),
            false => -(-raw).rem_euclid(TAU),
        };
        return Ok(Treatment {
            incoming_end,
            outgoing_start,
            piece: Segment::Arc {
                center: centre,
                radius,
                start_angle,
                sweep,
            },
        });
    }
    Err(no_fillet)
}

fn angle(v: Vector2<f64>) -> f64 {
    v.y.atan2(v.x)
}

#[cfg(test)]
mod tests {
    use core::f64::consts::PI;

    use super::*;

    fn square(a: f64) -> Profile {
        Profile::rectangle(Vector2::zeros(), Vector2::new(a, a)).unwrap()
    }

    /// Consecutive segments meet, and do so tangentially at each fillet.
    fn assert_valid(p: &Profile) {
        let rebuilt = Profile::new(p.segments().to_vec()).unwrap();
        assert_eq!(&rebuilt, p);
    }

    fn assert_tangent_joins(p: &Profile, arc_radius: f64) {
        let segments = p.segments();
        let n = segments.len();
        for i in 0..n {
            let (s, next) = (&segments[i], &segments[(i + 1) % n]);
            assert!((s.end() - next.start()).norm() < 1e-12);
            let is_fillet =
                |seg: &Segment| matches!(*seg, Segment::Arc { radius, .. } if radius == arc_radius);
            if is_fillet(s) || is_fillet(next) {
                assert!(
                    (s.tangent_at(1.0) - next.tangent_at(0.0)).norm() < 1e-12,
                    "joint {i}"
                );
            }
        }
    }

    #[test]
    fn filleted_square() {
        let (a, r) = (3.0, 0.5);
        let p = square(a).fillet_all(r).unwrap();
        assert_eq!(p.segments().len(), 8);
        assert!((p.signed_area() - (a * a - (4.0 - PI) * r * r)).abs() < 1e-12);
        assert_valid(&p);
        assert_tangent_joins(&p, r);
        assert!(p.sharp_corners().is_empty());
        // Clockwise too.
        let q = square(a).reversed().fillet_all(r).unwrap();
        assert!((q.signed_area() + (a * a - (4.0 - PI) * r * r)).abs() < 1e-12);
    }

    #[test]
    fn concave_corner_fillet() {
        // An L: its inner corner (index 3, at (1, 1)) is concave.
        let l = Profile::polygon(&[
            Vector2::new(0.0, 0.0),
            Vector2::new(2.0, 0.0),
            Vector2::new(2.0, 1.0),
            Vector2::new(1.0, 1.0),
            Vector2::new(1.0, 2.0),
            Vector2::new(0.0, 2.0),
        ])
        .unwrap();
        let r = 0.25;
        let p = l.fillet(&[3], r).unwrap();
        // A concave fillet adds the corner's (1 − π/4)r².
        assert!((p.signed_area() - (3.0 + (1.0 - PI / 4.0) * r * r)).abs() < 1e-12);
        assert_valid(&p);
        assert_tangent_joins(&p, r);
    }

    #[test]
    fn line_arc_and_arc_arc_fillets() {
        // A semicircle on a line: two line–arc corners. Then a lens of two arcs.
        let half = Profile::builder(Vector2::new(-2.0, 0.0))
            .line_to(Vector2::new(2.0, 0.0))
            .arc_to(Vector2::new(-2.0, 0.0), Vector2::zeros(), true)
            .close()
            .unwrap();
        let p = half.fillet_all(0.3).unwrap();
        assert_eq!(p.segments().len(), 4);
        assert_valid(&p);
        assert_tangent_joins(&p, 0.3);
        assert!(p.signed_area() < half.signed_area());

        // A lens: two arcs of radius √2 meeting at (0, ±1).
        let lens = Profile::builder(Vector2::new(0.0, -1.0))
            .arc_to(Vector2::new(0.0, 1.0), Vector2::new(-1.0, 0.0), true)
            .arc_to(Vector2::new(0.0, -1.0), Vector2::new(1.0, 0.0), true)
            .close()
            .unwrap();
        let p = lens.fillet_all(0.2).unwrap();
        assert_eq!(p.segments().len(), 4);
        assert_valid(&p);
        assert_tangent_joins(&p, 0.2);
    }

    #[test]
    fn chamfered_square() {
        let (a, d) = (3.0, 0.5);
        let p = square(a).chamfer_all(d).unwrap();
        assert_eq!(p.segments().len(), 8);
        assert!((p.signed_area() - (a * a - 2.0 * d * d)).abs() < 1e-12);
        assert_valid(&p);
        // Asymmetric: each corner loses a d1 × d2 / 2 triangle.
        let q = square(a).chamfer_asymmetric(&[0, 2], 0.5, 1.0).unwrap();
        assert!((q.signed_area() - (a * a - 2.0 * 0.25)).abs() < 1e-12);
        assert_valid(&q);
    }

    #[test]
    fn chamfer_along_an_arc_uses_arc_length() {
        let half = Profile::builder(Vector2::new(-2.0, 0.0))
            .line_to(Vector2::new(2.0, 0.0))
            .arc_to(Vector2::new(-2.0, 0.0), Vector2::zeros(), true)
            .close()
            .unwrap();
        // Corner 1 is at (2, 0): 0.5 back along the line, π/4 · 2 along the arc.
        let p = half.chamfer_asymmetric(&[1], 0.5, PI / 2.0).unwrap();
        assert!((p.segments()[1].start() - Vector2::new(1.5, 0.0)).norm() < 1e-12);
        assert!((p.segments()[1].end() - Vector2::new(2f64.sqrt(), 2f64.sqrt())).norm() < 1e-12);
    }

    #[test]
    fn errors() {
        let s = square(1.0);
        assert_eq!(
            s.fillet(&[7], 0.1),
            Err(SketchError::NoSuchCorner { corner: 7 })
        );
        assert!(matches!(
            s.fillet_all(0.6),
            Err(SketchError::TooLarge { .. })
        ));
        assert!(matches!(
            s.chamfer_all(0.6),
            Err(SketchError::TooLarge { .. })
        ));
        let round = s.fillet_all(0.2).unwrap();
        assert!(matches!(
            round.fillet(&[1], 0.1),
            Err(SketchError::SmoothCorner { .. })
        ));
        assert!(matches!(
            s.fillet(&[0], -1.0),
            Err(SketchError::InvalidShape { .. })
        ));
        // A circle has no sharp corners: nothing to do.
        let c = Profile::circle(Vector2::zeros(), 1.0).unwrap();
        assert_eq!(c.fillet_all(0.1).unwrap(), c);
    }
}
