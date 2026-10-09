//! Fitting pcurves to projected samples, and moving them by whole periods.

use golf_brep::FaceUv;
use golf_geom::AnyCurve2;
use golf_geom::Line;
use golf_geom::NurbsCurve;
use golf_manifold::Point;
use golf_manifold::Space;
use golf_manifold::Vector;
use nalgebra::Vector2;

/// The pcurve through `samples` (edge parameter, uv), increasing in
/// parameter: a line if they lie on one within `uv_tolerance` (at the rate
/// the parameter runs), else the piecewise-linear curve through them, which
/// shares the edge's parameter exactly at each sample.
pub(crate) fn fit<S: Space<3>>(
    samples: &[(f64, Vector2<f64>)],
    uv_tolerance: f64,
) -> Option<AnyCurve2<FaceUv<S>>> {
    let (&(t0, a), &(t1, b)) = (samples.first()?, samples.last()?);
    // Also rejects NaN.
    if t1.partial_cmp(&t0) != Some(core::cmp::Ordering::Greater) {
        return None;
    }
    let rate = (b - a) / (t1 - t0);
    let origin = a - rate * t0;
    if samples
        .iter()
        .all(|&(t, uv)| (origin + rate * t - uv).norm() <= uv_tolerance)
    {
        return Some(Line::new(Point::new(origin), Vector::new(rate)).into());
    }
    let mut knots = vec![t0];
    knots.extend(samples.iter().map(|&(t, _)| t));
    knots.push(t1);
    let points = samples.iter().map(|&(_, uv)| uv).collect();
    let curve = golf_nurbs::NurbsCurve::new(1, knots, points, None).ok()?;
    Some(NurbsCurve::new(curve).into())
}

/// `pcurve` moved by `shift` in parameter space.
pub(crate) fn shifted<S: Space<3>>(
    pcurve: &AnyCurve2<FaceUv<S>>,
    shift: Vector2<f64>,
) -> Option<AnyCurve2<FaceUv<S>>> {
    match pcurve {
        AnyCurve2::Line(line) => {
            Some(Line::new(Point::new(line.origin.coords + shift), line.direction).into())
        }
        AnyCurve2::Nurbs(nurbs) => {
            let g = nurbs.geometry();
            let count = g.control_point_count();
            let points = (0..count).map(|i| g.control_point(i) + shift).collect();
            let weights = (0..count).map(|i| g.weight(i)).collect();
            let moved = golf_nurbs::NurbsCurve::new(
                g.degree(),
                g.knots().knots().to_vec(),
                points,
                Some(weights),
            )
            .ok()?;
            Some(NurbsCurve::new(moved).into())
        }
        _ => None,
    }
}
