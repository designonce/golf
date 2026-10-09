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
}

impl From<Profile> for Region {
    fn from(outer: Profile) -> Self {
        Self::new(outer, Vec::new())
    }
}
