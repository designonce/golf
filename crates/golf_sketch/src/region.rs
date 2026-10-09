use nalgebra::Vector2;

use crate::error::SketchError;
use crate::profile::Profile;

/// An area of a sketch: an outer profile less any holes.
///
/// Orientation is normalised: the outer profile runs anticlockwise and holes
/// clockwise, so the region is always on the left. Holes are assumed to lie
/// inside the outer profile and not to touch it or each other.
#[derive(Clone, Debug, PartialEq)]
pub struct Region {
    outer: Profile,
    holes: Vec<Profile>,
}

impl Region {
    pub fn new(outer: Profile, holes: Vec<Profile>) -> Self {
        let anticlockwise = |p: Profile| {
            if p.signed_area() < 0.0 {
                p.reversed()
            } else {
                p
            }
        };
        let clockwise = |p: Profile| {
            if p.signed_area() > 0.0 {
                p.reversed()
            } else {
                p
            }
        };
        Self {
            outer: anticlockwise(outer),
            holes: holes.into_iter().map(clockwise).collect(),
        }
    }

    pub fn outer(&self) -> &Profile {
        &self.outer
    }

    pub fn holes(&self) -> &[Profile] {
        &self.holes
    }

    /// The outer profile, then the holes.
    pub fn profiles(&self) -> impl Iterator<Item = &Profile> {
        core::iter::once(&self.outer).chain(&self.holes)
    }

    /// The enclosed area, holes subtracted.
    pub fn area(&self) -> f64 {
        self.profiles().map(Profile::signed_area).sum()
    }

    pub fn translated(&self, v: Vector2<f64>) -> Self {
        self.map(|p| p.translated(v))
    }

    /// Rotated anticlockwise by `angle` about `about`.
    pub fn rotated(&self, angle: f64, about: Vector2<f64>) -> Self {
        self.map(|p| p.rotated(angle, about))
    }

    /// Scaled by `factor` about `about`.
    ///
    /// # Panics
    /// If `factor` isn't positive.
    pub fn scaled(&self, factor: f64, about: Vector2<f64>) -> Self {
        self.map(|p| p.scaled(factor, about))
    }

    /// Reflected in the line through `point` along `direction`, orientations
    /// re-normalised.
    pub fn mirrored(&self, point: Vector2<f64>, direction: Vector2<f64>) -> Self {
        self.map(|p| p.mirrored(point, direction))
    }

    /// The region with its boundary moved `distance` into the surrounding
    /// space (out of it if negative): the outer profile grows and holes shrink.
    /// See [`Profile::offset`].
    pub fn offset(&self, distance: f64) -> Result<Self, SketchError> {
        let outer = self.outer.offset(distance)?;
        let holes = self
            .holes
            .iter()
            .map(|hole| hole.offset(-distance))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self::new(outer, holes))
    }

    fn map(&self, f: impl Fn(&Profile) -> Profile) -> Self {
        Self::new(f(&self.outer), self.holes.iter().map(f).collect())
    }
}

impl From<Profile> for Region {
    fn from(outer: Profile) -> Self {
        Self::new(outer, Vec::new())
    }
}

#[cfg(test)]
mod tests {
    use core::f64::consts::PI;

    use super::*;

    fn plate() -> Region {
        let outer = Profile::rectangle(Vector2::zeros(), Vector2::new(4.0, 3.0)).unwrap();
        let hole = Profile::circle(Vector2::new(2.0, 1.5), 0.5).unwrap();
        Region::new(outer, vec![hole])
    }

    #[test]
    fn transforms_keep_area_and_orientation() {
        let r = plate();
        let area = 12.0 - PI * 0.25;
        for moved in [
            r.translated(Vector2::new(1.0, -2.0)),
            r.rotated(0.4, Vector2::new(3.0, 1.0)),
            r.mirrored(Vector2::new(1.0, 1.0), Vector2::new(1.0, 2.0)),
        ] {
            assert!((moved.area() - area).abs() < 1e-12);
            assert!(moved.outer().signed_area() > 0.0);
            assert!(moved.holes()[0].signed_area() < 0.0);
        }
        assert!((r.scaled(2.0, Vector2::zeros()).area() - 4.0 * area).abs() < 1e-11);
    }

    #[test]
    fn offset_grows_the_material() {
        let d = 0.25;
        let grown = plate().offset(d).unwrap();
        let outer = 12.0 + 14.0 * d + PI * d * d;
        let hole = PI * (0.5 - d).powi(2);
        assert!((grown.area() - (outer - hole)).abs() < 1e-12);
        assert_eq!(
            plate().offset(0.5),
            Err(SketchError::OffsetCollapses { index: 0 })
        );
    }
}
