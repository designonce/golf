use nalgebra::Vector2;

use crate::error::SketchError;
use crate::profile::JOIN_TOLERANCE;
use crate::profile::Profile;
use crate::segment::Segment;

/// A chain of segments, each starting where the last ended: open, or closed
/// if the last ends where the first starts. Draw one with
/// [`Profile::builder`](crate::Profile::builder) and
/// [`ProfileBuilder::path`](crate::ProfileBuilder::path).
#[derive(Clone, Debug, PartialEq)]
pub struct Path {
    segments: Vec<Segment>,
}

impl Path {
    /// Checks `segments` form a chain of non-degenerate segments.
    pub fn new(segments: Vec<Segment>) -> Result<Self, SketchError> {
        if segments.is_empty() {
            return Err(SketchError::EmptyProfile);
        }
        let size = segments
            .iter()
            .map(|s| s.start().norm().max(s.end().norm()) + s.length())
            .fold(0.0, f64::max);
        let tolerance = JOIN_TOLERANCE * (1.0 + size);
        for (index, segment) in segments.iter().enumerate() {
            if segment.length() <= tolerance {
                return Err(SketchError::DegenerateSegment { index });
            }
        }
        for (index, pair) in segments.windows(2).enumerate() {
            let gap = (pair[0].end() - pair[1].start()).norm();
            if gap > tolerance {
                return Err(SketchError::OpenProfile { index, gap });
            }
        }
        Ok(Self { segments })
    }

    pub fn segments(&self) -> &[Segment] {
        &self.segments
    }

    pub fn start(&self) -> Vector2<f64> {
        self.segments[0].start()
    }

    pub fn end(&self) -> Vector2<f64> {
        self.segments[self.segments.len() - 1].end()
    }

    /// Whether it ends where it starts.
    pub fn is_closed(&self) -> bool {
        let size = self.start().norm().max(1.0);
        (self.end() - self.start()).norm() <= JOIN_TOLERANCE * size
    }

    pub fn length(&self) -> f64 {
        self.segments.iter().map(Segment::length).sum()
    }
}

impl From<Profile> for Path {
    /// A profile as a closed path.
    fn from(profile: Profile) -> Self {
        Self {
            segments: profile.segments().to_vec(),
        }
    }
}

#[cfg(test)]
mod tests {
    use core::f64::consts::PI;

    use super::*;

    #[test]
    fn open_and_closed_paths() {
        let path = Profile::builder(Vector2::zeros())
            .line_to(Vector2::new(2.0, 0.0))
            .arc_to(Vector2::new(3.0, 1.0), Vector2::new(2.0, 1.0), true)
            .path()
            .unwrap();
        assert!(!path.is_closed());
        assert!((path.length() - (2.0 + PI / 2.0)).abs() < 1e-12);
        let ring = Path::from(Profile::circle(Vector2::zeros(), 1.0).unwrap());
        assert!(ring.is_closed());
        let broken = vec![
            Segment::Line {
                start: Vector2::zeros(),
                end: Vector2::x(),
            },
            Segment::Line {
                start: Vector2::y(),
                end: Vector2::new(1.0, 1.0),
            },
        ];
        assert!(matches!(
            Path::new(broken),
            Err(SketchError::OpenProfile { index: 0, .. })
        ));
    }
}
