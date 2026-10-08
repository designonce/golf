use core::f64::consts::PI;

use nalgebra::Matrix3x1;
use nalgebra::Vector1;
use nalgebra::Vector3;

use super::Placement;
use crate::domain::Axis;
use crate::domain::Domain;
use crate::manifold::Embedding;
use crate::manifold::ProjectError;
use crate::mapping::Mapping;
use crate::space::Point;
use crate::space::Space;
use crate::space::T;

/// A circle in the local xy plane of its placement, about the local origin.
///
/// `t ∈ [-π, π)` is the angle from local x, and is periodic.
#[derive(Clone, Debug)]
pub struct Circle<S: Space<3>> {
    pub placement: Placement<S>,
    pub radius: f64,
}

impl<S: Space<3>> Circle<S> {
    pub fn new(placement: Placement<S>, radius: f64) -> Self {
        Self { placement, radius }
    }
}

impl<S: Space<3>> Space<1> for T<Circle<S>> {
    type Tag = ();
}

impl<S: Space<3>> Mapping<1, 3> for Circle<S> {
    type From = T<Circle<S>>;
    type To = S;

    fn apply(&self, t: Point<Self::From, 1>) -> Point<S, 3> {
        let (t, r) = (t.coords.x, self.radius);
        self.placement
            .point(Vector3::new(r * t.cos(), r * t.sin(), 0.0))
    }

    fn jacobian(&self, t: Point<Self::From, 1>) -> Matrix3x1<f64> {
        let (t, r) = (t.coords.x, self.radius);
        self.placement.rotation_matrix() * Vector3::new(-r * t.sin(), r * t.cos(), 0.0)
    }
}

impl<S: Space<3>> Embedding<1, 3> for Circle<S> {
    fn domain(&self) -> Domain<1> {
        Domain::new([Axis::periodic(-PI, PI)])
    }

    /// On the axis every `t` is equally near, so the hint is kept.
    fn project(
        &self,
        p: Point<S, 3>,
        hint: Option<Point<Self::From, 1>>,
    ) -> Result<Point<Self::From, 1>, ProjectError> {
        let d = self.placement.to_local(p);
        let t = match d.xy().norm() <= 1e-15 * (d.norm() + self.radius) {
            true => hint.ok_or(ProjectError::Ambiguous)?.coords.x,
            false => d.y.atan2(d.x),
        };
        Ok(self.domain().wrap(Point::new(Vector1::new(t))))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::testing::assert_jacobian;
    use crate::geom::testing::assert_matches_newton;
    use crate::geom::testing::assert_round_trip;
    use crate::manifold::Curve;
    use crate::space::Vector;
    use crate::space::World;

    fn circle() -> Circle<World> {
        Circle::new(
            Placement::from_axes(
                Point::new(Vector3::new(2.0, 0.0, 1.0)),
                Vector::new(Vector3::new(1.0, 1.0, 0.0)),
                Vector::new(Vector3::z()),
            ),
            1.5,
        )
    }

    #[test]
    fn round_trip_and_jacobian() {
        let c = circle();
        for at in [[0.0], [1.0], [-PI + 1e-9], [PI - 1e-9]] {
            assert_round_trip(&c, at);
            assert_jacobian(&c, at);
        }
    }

    #[test]
    fn tangent_runs_anticlockwise_about_local_z() {
        let c = circle();
        let local = c
            .placement
            .vector_to_local(c.tangent(Point::new(Vector1::new(0.0))));
        assert!((local - Vector3::y()).norm() < 1e-12);
    }

    #[test]
    fn axis_keeps_hint() {
        let c = circle();
        let on_axis = c.placement.point(Vector3::new(0.0, 0.0, 2.0));
        assert_eq!(c.project(on_axis, None), Err(ProjectError::Ambiguous));
        let hint = Point::new(Vector1::new(0.4));
        assert_eq!(c.project(on_axis, Some(hint)), Ok(hint));
    }

    #[test]
    fn off_circle_points_match_newton() {
        let c = circle();
        for (local, start) in [
            (Vector3::new(3.0, 1.0, 2.0), [0.0]),
            (Vector3::new(-0.2, 0.1, -1.0), [2.9]),
        ] {
            assert_matches_newton(&c, c.placement.point(local), start);
        }
    }
}
