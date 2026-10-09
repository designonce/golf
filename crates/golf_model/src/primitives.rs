//! Simple solids, built by sketching and extruding or revolving.

use core::f64::consts::FRAC_PI_2;
use core::f64::consts::PI;

use golf_brep::Body;
use golf_geom::Placement;
use golf_manifold::Space;
use golf_sketch::Profile;
use golf_sketch::Region;
use golf_sketch::Segment;
use nalgebra::Vector2;
use nalgebra::Vector3;

use crate::ends::EdgeTreatment;
use crate::ends::ExtrudeEnds;
use crate::ends::extrude_with;
use crate::error::ModelError;
use crate::extrude::extrude;
use crate::revolve::revolve;

/// A box with a corner at the placement's origin, extending `size` along its
/// local axes.
pub fn cuboid<S: Space<3>>(
    placement: &Placement<S>,
    size: Vector3<f64>,
) -> Result<Body<S>, ModelError> {
    let base = Profile::rectangle(Vector2::zeros(), size.xy())?;
    extrude(placement, &base.into(), size.z)
}

/// A box like [`cuboid`] with every edge rounded to `radius` (at most half the
/// smallest side): its corners are spherical.
pub fn rounded_cuboid<S: Space<3>>(
    placement: &Placement<S>,
    size: Vector3<f64>,
    radius: f64,
) -> Result<Body<S>, ModelError> {
    let base = Profile::rounded_rectangle(Vector2::zeros(), size.xy(), radius)?;
    extrude_with(
        placement,
        &base.into(),
        size.z,
        ExtrudeEnds::both(EdgeTreatment::Fillet(radius)),
    )
}

/// A box like [`cuboid`] with every edge chamfered by `distance`: the upright
/// edges in the sketch, then the ends, mitred where they meet.
pub fn chamfered_cuboid<S: Space<3>>(
    placement: &Placement<S>,
    size: Vector3<f64>,
    distance: f64,
) -> Result<Body<S>, ModelError> {
    let base = Profile::rectangle(Vector2::zeros(), size.xy())?.chamfer_all(distance)?;
    extrude_with(
        placement,
        &base.into(),
        size.z,
        ExtrudeEnds::both(EdgeTreatment::chamfer(distance)),
    )
}

/// A cylinder about the placement's local z axis, from its origin up to
/// `height`.
pub fn cylinder<S: Space<3>>(
    placement: &Placement<S>,
    radius: f64,
    height: f64,
) -> Result<Body<S>, ModelError> {
    extrude(
        placement,
        &Profile::circle(Vector2::zeros(), radius)?.into(),
        height,
    )
}

/// A cone (or frustum, if `top_radius` isn't zero) about the placement's local
/// z axis, from `bottom_radius` at its origin to `top_radius` at `height`.
pub fn cone<S: Space<3>>(
    placement: &Placement<S>,
    bottom_radius: f64,
    top_radius: f64,
    height: f64,
) -> Result<Body<S>, ModelError> {
    let mut points = vec![Vector2::zeros(), Vector2::new(bottom_radius, 0.0)];
    if top_radius > 0.0 {
        points.push(Vector2::new(top_radius, height));
    }
    points.push(Vector2::new(0.0, height));
    revolve(&axial_sketch(placement), &Profile::polygon(&points)?.into())
}

/// A sphere about the placement's origin, with its poles on the local z axis.
pub fn sphere<S: Space<3>>(placement: &Placement<S>, radius: f64) -> Result<Body<S>, ModelError> {
    let half_disc = Profile::new(vec![
        Segment::Arc {
            center: Vector2::zeros(),
            radius,
            start_angle: -FRAC_PI_2,
            sweep: PI,
        },
        Segment::Line {
            start: Vector2::new(0.0, radius),
            end: Vector2::new(0.0, -radius),
        },
    ])?;
    revolve(&axial_sketch(placement), &half_disc.into())
}

/// A torus about the placement's local z axis: a tube of `minor_radius`
/// round a circle of `major_radius` in the local xy plane.
pub fn torus<S: Space<3>>(
    placement: &Placement<S>,
    major_radius: f64,
    minor_radius: f64,
) -> Result<Body<S>, ModelError> {
    let tube = Profile::circle(Vector2::new(major_radius, 0.0), minor_radius)?;
    revolve(&axial_sketch(placement), &Region::from(tube))
}

/// The sketch plane containing `placement`'s local x and z axes, with sketch x
/// along local x and sketch y along local z, so revolving about the sketch's
/// y axis revolves about the placement's z.
pub fn axial_sketch<S: Space<3>>(placement: &Placement<S>) -> Placement<S> {
    Placement::from_axes(
        placement.origin,
        placement.vector(-Vector3::y()),
        placement.vector(Vector3::x()),
    )
}
