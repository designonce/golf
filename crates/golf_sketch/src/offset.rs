//! Offsetting a profile: moving its boundary a constant distance.

use crate::error::SketchError;
use crate::locus::Locus;
use crate::profile::Profile;
use crate::segment::Segment;

/// How close to its ends a trimmed segment may be cut before it counts as
/// collapsed.
const TRIM_TOLERANCE: f64 = 1e-9;

impl Profile {
    /// The profile moved `distance` outwards (inwards if negative), whichever
    /// way it runs: a positive distance grows the enclosed area.
    ///
    /// Lines move along their normals and arcs change radius about their
    /// centres. Where offset segments part at a corner a round join (an arc
    /// about the corner of radius `|distance|`) fills the gap; where they
    /// overlap both are trimmed back to where they cross; tangent joints still
    /// meet. Fails if a segment would shrink to nothing or turn inside out;
    /// self-intersections between segments that aren't neighbours aren't
    /// detected.
    pub fn offset(&self, distance: f64) -> Result<Self, SketchError> {
        if distance == 0.0 {
            return Ok(self.clone());
        }
        // The material is on the left of an anticlockwise profile, so outwards
        // is to the right.
        let left = if self.signed_area() > 0.0 {
            -distance
        } else {
            distance
        };
        let segments = self.segments();
        let n = segments.len();
        let shifted = segments
            .iter()
            .enumerate()
            .map(|(index, s)| {
                s.left_offset(left)
                    .ok_or(SketchError::OffsetCollapses { index })
            })
            .collect::<Result<Vec<_>, _>>()?;

        let mut starts = vec![0.0; n];
        let mut ends = vec![1.0; n];
        let mut joins: Vec<Option<Segment>> = vec![None; n];
        for i in 0..n {
            if self.is_smooth(i) {
                continue;
            }
            let incoming = (i + n - 1) % n;
            let (ta, tb) = (
                segments[incoming].tangent_at(1.0),
                segments[i].tangent_at(0.0),
            );
            let cross = ta.x * tb.y - ta.y * tb.x;
            if cross.abs() <= 1e-9 {
                // Doubling straight back: no sensible offset.
                return Err(SketchError::OffsetCollapses { index: i });
            }
            let corner = segments[i].start();
            if cross * left < 0.0 {
                // Offset to the outside of the turn: the ends part, so join them
                // round the corner, turning as the path does.
                let from = shifted[incoming].end() - corner;
                joins[i] = Some(Segment::Arc {
                    center: corner,
                    radius: left.abs(),
                    start_angle: from.y.atan2(from.x),
                    sweep: cross.atan2(ta.dot(&tb)),
                });
            } else {
                // To the inside: the ends overlap, so cut both where they cross.
                let mut crossings =
                    Locus::of(&shifted[incoming]).intersect(&Locus::of(&shifted[i]));
                crossings.sort_by(|p, q| (p - corner).norm().total_cmp(&(q - corner).norm()));
                let x = crossings
                    .first()
                    .ok_or(SketchError::OffsetCollapses { index: i })?;
                ends[incoming] = shifted[incoming].fraction_of(*x);
                starts[i] = shifted[i].fraction_of(*x);
                if ends[incoming] > 1.0 + TRIM_TOLERANCE || starts[i] < -TRIM_TOLERANCE {
                    return Err(SketchError::OffsetCollapses { index: i });
                }
            }
        }
        let mut out = Vec::with_capacity(2 * n);
        for i in 0..n {
            if starts[i] >= ends[i] - TRIM_TOLERANCE {
                return Err(SketchError::OffsetCollapses { index: i });
            }
            out.extend(joins[i]);
            out.push(shifted[i].between(starts[i], ends[i]));
        }
        Profile::new(out)
    }
}

#[cfg(test)]
mod tests {
    use core::f64::consts::PI;

    use nalgebra::Vector2;

    use super::*;

    fn square(a: f64) -> Profile {
        Profile::rectangle(Vector2::zeros(), Vector2::new(a, a)).unwrap()
    }

    #[test]
    fn square_outwards_gets_round_corners() {
        let (a, d) = (2.0, 0.5);
        for p in [square(a), square(a).reversed()] {
            let grown = p.offset(d).unwrap();
            assert_eq!(grown.segments().len(), 8);
            assert!((grown.signed_area().abs() - (a * a + 4.0 * a * d + PI * d * d)).abs() < 1e-12);
            assert_eq!(grown.signed_area().signum(), p.signed_area().signum());
        }
    }

    #[test]
    fn square_inwards_stays_square() {
        let (a, d) = (2.0, 0.5);
        let shrunk = square(a).offset(-d).unwrap();
        assert_eq!(shrunk.segments().len(), 4);
        assert!((shrunk.signed_area() - (a - 2.0 * d).powi(2)).abs() < 1e-12);
        assert!((shrunk.segments()[0].start() - Vector2::new(d, d)).norm() < 1e-12);
        assert_eq!(
            square(a).offset(-1.0),
            Err(SketchError::OffsetCollapses { index: 0 })
        );
    }

    #[test]
    fn rounded_shapes_offset_exactly() {
        // A circle and a slot: tangent joints, so only radii change.
        let c = Profile::circle(Vector2::new(1.0, 1.0), 2.0).unwrap();
        assert!((c.offset(0.5).unwrap().signed_area() - PI * 6.25).abs() < 1e-12);
        assert!((c.offset(-0.5).unwrap().signed_area() - PI * 2.25).abs() < 1e-12);
        assert_eq!(
            c.offset(-2.0),
            Err(SketchError::OffsetCollapses { index: 0 })
        );
        let slot = Profile::slot(Vector2::zeros(), Vector2::new(3.0, 0.0), 1.0).unwrap();
        let grown = slot.offset(0.5).unwrap();
        assert_eq!(grown.segments().len(), 4);
        assert!((grown.signed_area() - (2.0 * 1.5 * 3.0 + PI * 2.25)).abs() < 1e-12);
    }

    #[test]
    fn concave_corners_trim() {
        // An L grown by d: its concave corner trims; its five convex ones round.
        let l = Profile::polygon(&[
            Vector2::new(0.0, 0.0),
            Vector2::new(2.0, 0.0),
            Vector2::new(2.0, 1.0),
            Vector2::new(1.0, 1.0),
            Vector2::new(1.0, 2.0),
            Vector2::new(0.0, 2.0),
        ])
        .unwrap();
        let d = 0.25;
        let grown = l.offset(d).unwrap();
        assert_eq!(grown.segments().len(), 11);
        // Area 3, a d-wide strip along the perimeter (8), a quarter disc at
        // each of the five convex corners, less the d × d square where the
        // strips overlap at the concave one.
        assert!(
            (grown.signed_area() - (3.0 + 8.0 * d + 5.0 * PI / 4.0 * d * d - d * d)).abs() < 1e-12
        );
        // Shrunk: the strips overlap at the five convex corners (each squared
        // off), and a quarter disc goes round the concave one.
        let shrunk = l.offset(-d).unwrap();
        assert_eq!(shrunk.segments().len(), 7);
        assert!(
            (shrunk.signed_area() - (3.0 - 8.0 * d + 5.0 * d * d - PI / 4.0 * d * d)).abs() < 1e-12
        );
    }
}
