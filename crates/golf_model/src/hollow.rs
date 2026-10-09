use golf_brep::Body;
use golf_geom::Placement;
use golf_geom::Plane;
use golf_manifold::Space;
use golf_sketch::Profile;
use golf_sketch::Region;

use crate::error::ModelError;
use crate::extrude::extrude;
use crate::extrude::shifted;
use crate::extrude::walls;

/// Extrudes `profile`, drawn in the sketch plane `plane`, up its normal by
/// `height` as a hollow with walls `wall` thick: a tube open at both ends if
/// `floor` is `None`, else a cup open at the top with a floor that thick.
///
/// The hollow is the profile offset inwards by `wall`, with sharp convex
/// corners staying sharp and concave ones rounding (see
/// [`Profile::offset`](golf_sketch::Profile::offset)).
pub fn extrude_hollow<S: Space<3>>(
    plane: &Placement<S>,
    profile: &Profile,
    height: f64,
    wall: f64,
    floor: Option<f64>,
) -> Result<Body<S>, ModelError> {
    if !height.is_finite() || height <= 0.0 {
        return Err(ModelError::ZeroDistance);
    }
    if !wall.is_finite() || wall <= 0.0 {
        return Err(ModelError::BadTreatment);
    }
    let region = Region::from(profile.clone());
    let hollow = region.outer().offset(-wall)?;
    let Some(floor) = floor else {
        return extrude(
            plane,
            &Region::new(region.outer().clone(), vec![hollow]),
            height,
        );
    };
    if !floor.is_finite() || floor <= 0.0 || floor >= height {
        return Err(ModelError::TreatmentsTooTall);
    }
    // The hollow, as a hole, runs clockwise: its walls face into it, and its
    // own bottom loop runs anticlockwise round the floor, which faces up.
    let hollow = Region::new(region.outer().clone(), vec![hollow]).holes()[0].clone();
    let mut body = Body::with_tag(plane.origin.tag());
    let outside = walls(&mut body, plane, height, region.outer().segments())?;
    let inside = walls(
        &mut body,
        &shifted(plane, floor),
        height - floor,
        hollow.segments(),
    )?;
    let mut faces = outside.sides;
    faces.extend(inside.sides);
    faces.push(body.add_face(Plane::new(*plane), false, vec![outside.bottom])?);
    faces.push(body.add_face(
        Plane::new(shifted(plane, height)),
        true,
        vec![outside.top, inside.top],
    )?);
    faces.push(body.add_face(Plane::new(shifted(plane, floor)), true, vec![inside.bottom])?);
    body.add_shell(faces)?;
    Ok(body)
}
