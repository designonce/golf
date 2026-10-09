//! Adaptive sampling of a parametrised curve to a tolerance.

use nalgebra::Vector3;

use crate::tolerance::Tolerance;

/// Parameters from `t0` to `t1` (both included) such that the polyline through
/// `f` at them meets `tolerance`.
///
/// Starts from `initial` equal segments, enough to tell a closed curve from a
/// point, then halves any segment that sags, turns or stretches too far. Sag is
/// probed at a quarter, half and three quarters, which catches S-bends a
/// midpoint alone would miss. `df` is the derivative of `f`.
pub(crate) fn sample(
    f: impl Fn(f64) -> Vector3<f64>,
    df: impl Fn(f64) -> Vector3<f64>,
    range: (f64, f64),
    initial: usize,
    tolerance: &Tolerance,
) -> Vec<f64> {
    sample_with(f, df, range, initial, tolerance, |_, _| true)
}

/// [`sample`], also halving any segment from `a` to `b` for which `also(a, b)`
/// is false.
pub(crate) fn sample_with(
    f: impl Fn(f64) -> Vector3<f64>,
    df: impl Fn(f64) -> Vector3<f64>,
    (t0, t1): (f64, f64),
    initial: usize,
    tolerance: &Tolerance,
    also: impl Fn(f64, f64) -> bool,
) -> Vec<f64> {
    const MAX_DEPTH: usize = 24;
    let settled = |a: f64, b: f64| {
        let (pa, pb) = (f(a), f(b));
        if tolerance
            .max_length
            .is_some_and(|max| (pb - pa).norm() > max)
        {
            return false;
        }
        let (ta, tb) = (df(a), df(b));
        if ta.norm() > 0.0 && tb.norm() > 0.0 && ta.angle(&tb) > tolerance.angle {
            return false;
        }
        [0.25, 0.5, 0.75]
            .iter()
            .all(|&s| distance_to_segment(f(a + (b - a) * s), pa, pb) <= tolerance.chord)
            && also(a, b)
    };

    fn split(
        a: f64,
        b: f64,
        depth: usize,
        settled: &impl Fn(f64, f64) -> bool,
        out: &mut Vec<f64>,
    ) {
        if depth < MAX_DEPTH && !settled(a, b) {
            let mid = (a + b) / 2.0;
            split(a, mid, depth + 1, settled, out);
            split(mid, b, depth + 1, settled, out);
        } else {
            out.push(b);
        }
    }

    let initial = initial.max(1);
    let mut out = vec![t0];
    for i in 0..initial {
        let a = t0 + (t1 - t0) * i as f64 / initial as f64;
        let b = if i + 1 == initial {
            t1
        } else {
            t0 + (t1 - t0) * (i + 1) as f64 / initial as f64
        };
        split(a, b, 0, &settled, &mut out);
    }
    out
}

pub(crate) fn distance_to_segment(p: Vector3<f64>, a: Vector3<f64>, b: Vector3<f64>) -> f64 {
    let ab = b - a;
    let length_squared = ab.norm_squared();
    if length_squared == 0.0 {
        return (p - a).norm();
    }
    let s = ((p - a).dot(&ab) / length_squared).clamp(0.0, 1.0);
    (p - (a + ab * s)).norm()
}

#[cfg(test)]
mod tests {
    use core::f64::consts::TAU;

    use super::*;

    #[test]
    fn circle_samples_meet_the_chord_tolerance() {
        let f = |t: f64| Vector3::new(t.cos(), t.sin(), 0.0) * 2.0;
        let df = |t: f64| Vector3::new(-t.sin(), t.cos(), 0.0) * 2.0;
        let tolerance = Tolerance::new(1e-3, 1.0);
        let ts = sample(f, df, (0.0, TAU), 3, &tolerance);
        assert_eq!((ts[0], *ts.last().unwrap()), (0.0, TAU));
        for w in ts.windows(2) {
            let sag = distance_to_segment(f((w[0] + w[1]) / 2.0), f(w[0]), f(w[1]));
            assert!(sag <= 1e-3);
        }
        // A sag of 1e-3 on radius 2 needs segments of about 0.063 rad: ~100.
        assert!((90..=260).contains(&ts.len()), "{}", ts.len());
    }

    #[test]
    fn straight_lines_need_one_segment_unless_too_long() {
        let f = |t: f64| Vector3::new(t, 2.0 * t, 0.0);
        let df = |_: f64| Vector3::new(1.0, 2.0, 0.0);
        let tolerance = Tolerance::new(1e-3, 0.1);
        assert_eq!(sample(f, df, (0.0, 1.0), 1, &tolerance), vec![0.0, 1.0]);
        let ts = sample(f, df, (0.0, 1.0), 1, &tolerance.with_max_length(0.6));
        assert_eq!(ts.len(), 5);
    }
}
