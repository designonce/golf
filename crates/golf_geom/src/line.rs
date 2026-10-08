use golf_manifold::Domain;
use golf_manifold::Embedding;
use golf_manifold::HasT;
use golf_manifold::Mapping;
use golf_manifold::Point;
use golf_manifold::ProjectError;
use golf_manifold::Space;
use golf_manifold::T;
use golf_manifold::Vector;
use nalgebra::SMatrix;
use nalgebra::Vector1;

/// The straight line `origin + t·direction` in an `N`-dimensional space `S`.
///
/// In a surface's uv space it is a pcurve; see [`Line2`].
#[derive(Clone, Debug)]
pub struct Line<S: Space<N>, const N: usize> {
    pub origin: Point<S, N>,
    pub direction: Vector<S, N>,
}

/// A line in a 2D space, e.g. a pcurve in a surface's uv space.
pub type Line2<S> = Line<S, 2>;

impl<S: Space<N>, const N: usize> Line<S, N> {
    pub fn new(origin: Point<S, N>, direction: Vector<S, N>) -> Self {
        Self { origin, direction }
    }
}

impl<S: Space<N>, const N: usize> HasT for Line<S, N> {}

impl<S: Space<N>, const N: usize> Mapping<1, N> for Line<S, N> {
    type From = T<Line<S, N>>;
    type To = S;

    fn apply(&self, t: Point<Self::From, 1>) -> Point<S, N> {
        self.origin + self.direction.scale(t.coords.x)
    }

    fn jacobian(&self, _t: Point<Self::From, 1>) -> SMatrix<f64, N, 1> {
        self.direction.coords
    }
}

impl<S: Space<N>, const N: usize> Embedding<1, N> for Line<S, N> {
    fn domain(&self) -> Domain<1> {
        Domain::unbounded()
    }

    /// A zero-length direction makes every `t` equally near, so the hint is returned.
    fn project(
        &self,
        p: Point<S, N>,
        hint: Option<Point<Self::From, 1>>,
    ) -> Result<Point<Self::From, 1>, ProjectError> {
        let length_squared = self.direction.coords.norm_squared();
        if length_squared == 0.0 {
            return hint.ok_or(ProjectError::Ambiguous);
        }
        let t = (p - self.origin).coords.dot(&self.direction.coords) / length_squared;
        Ok(Point::new(Vector1::new(t)))
    }
}

#[cfg(test)]
mod tests {
    use golf_manifold::Curve;
    use golf_manifold::Uv;
    use golf_manifold::World;
    use nalgebra::Vector2;
    use nalgebra::Vector3;

    use super::*;
    use crate::Plane;
    use crate::testing::assert_jacobian;
    use crate::testing::assert_matches_newton;
    use crate::testing::assert_round_trip;

    fn line() -> Line<World, 3> {
        Line::new(
            Point::new(Vector3::new(1.0, 0.0, -1.0)),
            Vector::new(Vector3::new(2.0, 1.0, 0.5)),
        )
    }

    #[test]
    fn round_trip_and_jacobian() {
        let l = line();
        for at in [[0.0], [2.5], [-40.0]] {
            assert_round_trip(&l, at);
            assert_jacobian(&l, at);
        }
        assert!(
            (l.tangent(Point::new(Vector1::new(0.0))).coords - l.direction.coords.normalize())
                .norm()
                < 1e-12
        );
    }

    #[test]
    fn off_line_points_match_newton() {
        assert_matches_newton(&line(), Point::new(Vector3::new(4.0, -3.0, 2.0)), [0.0]);
    }

    #[test]
    fn degenerate_line_needs_hint() {
        let l = Line2::<Uv<Plane<World>>>::new(
            Point::new(Vector2::zeros()),
            Vector::new(Vector2::zeros()),
        );
        let p = Point::new(Vector2::new(1.0, 1.0));
        assert_eq!(l.project(p, None), Err(ProjectError::Ambiguous));
        let hint = Point::new(Vector1::new(3.0));
        assert_eq!(l.project(p, Some(hint)), Ok(hint));
    }
}
