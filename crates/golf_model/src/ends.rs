//! Extrusions with treated (chamfered or filleted) end edges, and drafted
//! extrusions.

use golf_brep::Body;
use golf_geom::Placement;
use golf_manifold::Space;
use golf_sketch::Region;

use crate::error::ModelError;
use crate::stack::Band;
use crate::stack::Ring;
use crate::stack::build;

/// What to do to the edges round one end of an extrusion.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum EdgeTreatment {
    /// Leave them sharp.
    #[default]
    None,
    /// Cut them off: `inset` in from the sides across the end, and `height`
    /// down the sides.
    Chamfer { inset: f64, height: f64 },
    /// Round them with this radius.
    Fillet(f64),
}

impl EdgeTreatment {
    /// A chamfer of equal inset and height.
    pub fn chamfer(distance: f64) -> Self {
        Self::Chamfer {
            inset: distance,
            height: distance,
        }
    }

    /// How far down the sides it reaches.
    fn height(&self) -> f64 {
        match *self {
            Self::None => 0.0,
            Self::Chamfer { height, .. } => height,
            Self::Fillet(radius) => radius,
        }
    }

    /// The cap's inset, and the band joining it to the sides.
    fn cap(&self) -> Option<(f64, Band)> {
        match *self {
            Self::None => None,
            Self::Chamfer { inset, .. } => Some((inset, Band::Ruled)),
            Self::Fillet(radius) => Some((radius, Band::Round)),
        }
    }

    fn check(&self) -> Result<(), ModelError> {
        let valid = |x: f64| x.is_finite() && x > 0.0;
        let ok = match *self {
            Self::None => true,
            Self::Chamfer { inset, height } => valid(inset) && valid(height),
            Self::Fillet(radius) => valid(radius),
        };
        match ok {
            true => Ok(()),
            false => Err(ModelError::BadTreatment),
        }
    }
}

/// The treatments of an extrusion's two ends.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ExtrudeEnds {
    /// The end in the sketch plane.
    pub bottom: EdgeTreatment,
    /// The far end.
    pub top: EdgeTreatment,
}

impl ExtrudeEnds {
    /// The same treatment at both ends.
    pub fn both(treatment: EdgeTreatment) -> Self {
        Self {
            bottom: treatment,
            top: treatment,
        }
    }
}

/// Extrudes `region` from the sketch plane `plane` up its normal by `height`,
/// chamfering or filleting the edges round either end.
///
/// The treatments follow the profiles round, so every corner of a profile must
/// be smooth (tangent), or a sharp corner between two lines (where the treated
/// faces meet in a mitre); fillet other sharp corners in the sketch first.
///
/// A fillet's radius round a convex arc must be less than half the arc's
/// radius (the fillet sweeps a torus), or exactly equal to it (a spherical
/// corner, as on a fully rounded box).
pub fn extrude_with<S: Space<3>>(
    plane: &Placement<S>,
    region: &Region,
    height: f64,
    ends: ExtrudeEnds,
) -> Result<Body<S>, ModelError> {
    if !height.is_finite() || height <= 0.0 {
        return Err(ModelError::ZeroDistance);
    }
    ends.bottom.check()?;
    ends.top.check()?;
    let (bottom, top) = (ends.bottom.height(), height - ends.top.height());
    if top < bottom {
        return Err(ModelError::TreatmentsTooTall);
    }
    let mut rings = Vec::new();
    let mut bands = Vec::new();
    if let Some((inset, band)) = ends.bottom.cap() {
        rings.push(Ring { z: 0.0, inset });
        bands.push(band);
    }
    rings.push(Ring {
        z: bottom,
        inset: 0.0,
    });
    if top > bottom {
        bands.push(Band::Ruled);
        rings.push(Ring { z: top, inset: 0.0 });
    }
    if let Some((inset, band)) = ends.top.cap() {
        bands.push(band);
        rings.push(Ring { z: height, inset });
    }
    if bands.is_empty() {
        return Err(ModelError::TreatmentsTooTall);
    }
    build(plane, region, &rings, &bands)
}

/// Extrudes `region` from the sketch plane `plane` up its normal by `height`,
/// with sides tapering in by `draft` radians from the normal (out, if
/// negative): a line sweeps a sloping plane, an arc a cone.
///
/// Corners must be smooth, or sharp between two lines (as for
/// [`extrude_with`]).
pub fn extrude_drafted<S: Space<3>>(
    plane: &Placement<S>,
    region: &Region,
    height: f64,
    draft: f64,
) -> Result<Body<S>, ModelError> {
    if !height.is_finite() || height <= 0.0 {
        return Err(ModelError::ZeroDistance);
    }
    if !draft.is_finite() || draft.abs() >= core::f64::consts::FRAC_PI_2 {
        return Err(ModelError::BadTreatment);
    }
    let rings = [
        Ring { z: 0.0, inset: 0.0 },
        Ring {
            z: height,
            inset: height * draft.tan(),
        },
    ];
    build(plane, region, &rings, &[Band::Ruled])
}
