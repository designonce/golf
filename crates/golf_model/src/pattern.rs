//! Copies of a body in a pattern.
//!
//! The copies are separate bodies: merging them into one needs booleans.

use golf_brep::Body;
use golf_geom::TransformError;
use golf_manifold::Point;
use golf_manifold::Space;
use golf_manifold::Vector;
use nalgebra::Isometry3;
use nalgebra::Translation3;
use nalgebra::Unit;
use nalgebra::UnitQuaternion;

/// `count` copies of `body` (the first in place), each `step` further along.
pub fn linear_pattern<S: Space<3>>(
    body: &Body<S>,
    step: Vector<S, 3>,
    count: usize,
) -> Result<Vec<Body<S>>, TransformError> {
    (0..count)
        .map(|i| {
            body.transformed(&Isometry3::from_parts(
                Translation3::from(step.coords * i as f64),
                UnitQuaternion::identity(),
            ))
        })
        .collect()
}

/// `count` copies of `body` (the first in place) turned about the axis through
/// `center` along `axis`, spread evenly over `angle` radians: a full turn puts
/// the last copy one step short of the first, anything less puts it at
/// `angle`.
pub fn circular_pattern<S: Space<3>>(
    body: &Body<S>,
    center: Point<S, 3>,
    axis: Vector<S, 3>,
    count: usize,
    angle: f64,
) -> Result<Vec<Body<S>>, TransformError> {
    let full = (angle.abs() - core::f64::consts::TAU).abs() <= 1e-12;
    let steps = if full || count < 2 {
        count.max(1)
    } else {
        count - 1
    };
    let axis = Unit::new_normalize(axis.coords);
    (0..count)
        .map(|i| {
            let turn = UnitQuaternion::from_axis_angle(&axis, angle * i as f64 / steps as f64);
            // About `center`: move it to the origin, turn, move it back.
            let about = Translation3::from(center.coords);
            let motion = Isometry3::from_parts(about, turn)
                * Isometry3::from_parts(about.inverse(), UnitQuaternion::identity());
            body.transformed(&motion)
        })
        .collect()
}
