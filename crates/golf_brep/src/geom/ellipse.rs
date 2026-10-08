use core::f64::consts::PI;

use nalgebra::Matrix3x1;
use nalgebra::Vector1;
use nalgebra::Vector3;

use super::Placement;
use crate::domain::Axis;
use crate::domain::Domain;
use crate::manifold::Embedding;
use crate::manifold::ProjectError;
use crate::manifold::newton_project;
use crate::mapping::Mapping;
use crate::space::Point;
use crate::space::Space;
use crate::space::T;

/// An ellipse in the local xy plane of its placement, major axis along local x.
///
/// `t ∈ [-π, π)` is the eccentric angle: the point is `(a cos t, b sin t, 0)`.
#[derive(Clone, Debug)]
pub struct Ellipse<S: Space<3>> {
    pub placement: Placement<S>,
    pub major_radius: f64,
    pub minor_radius: f64,
}

impl<S: Space<3>> Ellipse<S> {
    pub fn new(placement: Placement<S>, major_radius: f64, minor_radius: f64) -> Self {
        Self {
            placement,
            major_radius,
            minor_radius,
        }
    }
}

impl<S: Space<3>> Space<1> for T<Ellipse<S>> {
    type Tag = ();
}

impl<S: Space<3>> Mapping<1, 3> for Ellipse<S> {
    type From = T<Ellipse<S>>;
    type To = S;

    fn apply(&self, t: Point<Self::From, 1>) -> Point<S, 3> {
        let t = t.coords.x;
        self.placement.point(Vector3::new(
            self.major_radius * t.cos(),
            self.minor_radius * t.sin(),
            0.0,
        ))
    }

    fn jacobian(&self, t: Point<Self::From, 1>) -> Matrix3x1<f64> {
        let t = t.coords.x;
        self.placement.rotation_matrix()
            * Vector3::new(
                -self.major_radius * t.sin(),
                self.minor_radius * t.cos(),
                0.0,
            )
    }
}

impl<S: Space<3>> Embedding<1, 3> for Ellipse<S> {
    fn domain(&self) -> Domain<1> {
        Domain::new([Axis::periodic(-PI, PI)])
    }

    /// Newton iteration from the hint, or without one from the nearest of 16
    /// evenly spaced samples.
    fn project(
        &self,
        p: Point<S, 3>,
        hint: Option<Point<Self::From, 1>>,
    ) -> Result<Point<Self::From, 1>, ProjectError> {
        const SEEDS: usize = 16;
        let start = hint.unwrap_or_else(|| {
            let distance = |t| (self.apply(t) - p).coords.norm_squared();
            (0..SEEDS)
                .map(|i| Point::new(Vector1::new(-PI + 2.0 * PI * i as f64 / SEEDS as f64)))
                .min_by(|a, b| distance(*a).total_cmp(&distance(*b)))
                .unwrap()
        });
        newton_project(self, p, start)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::testing::assert_jacobian;
    use crate::geom::testing::assert_round_trip;
    use crate::space::Vector;
    use crate::space::World;

    fn ellipse() -> Ellipse<World> {
        Ellipse::new(
            Placement::from_axes(
                Point::new(Vector3::new(0.0, 1.0, -2.0)),
                Vector::new(Vector3::new(0.0, 0.0, 1.0)),
                Vector::new(Vector3::new(1.0, 1.0, 0.0)),
            ),
            3.0,
            1.0,
        )
    }

    #[test]
    fn round_trip_and_jacobian() {
        let e = ellipse();
        for at in [[0.0], [0.7], [-2.0], [-PI + 1e-9], [PI - 1e-9]] {
            assert_round_trip(&e, at);
            assert_jacobian(&e, at);
        }
    }

    #[test]
    fn off_ellipse_points_match_dense_sampling() {
        let e = ellipse();
        for local in [
            Vector3::new(4.0, 2.0, 0.0),
            Vector3::new(1.0, 0.3, 1.0),
            Vector3::new(-2.0, -5.0, -0.5),
            Vector3::new(2.9, 0.05, 0.0),
        ] {
            let p = e.placement.point(local);
            let q = e.project(p, None).unwrap();
            let best = (0..20_000)
                .map(|i| -PI + 2.0 * PI * i as f64 / 20_000.0)
                .map(|t| (e.apply(Point::new(Vector1::new(t))) - p).coords.norm())
                .fold(f64::INFINITY, f64::min);
            let found = (e.apply(q) - p).coords.norm();
            assert!(
                found <= best + 1e-9,
                "{local}: found {found}, sampled {best}"
            );
        }
    }

    #[test]
    fn hint_picks_the_local_minimum() {
        let e = ellipse();
        // From the centre the two ends of the minor axis are equally near.
        let centre = e.placement.origin;
        for hint in [1.4, -1.4] {
            let q = e
                .project(centre, Some(Point::new(Vector1::new(hint))))
                .unwrap();
            assert!((q.coords.x - hint.signum() * PI / 2.0).abs() < 1e-9);
        }
    }

    #[test]
    fn circular_ellipse_matches_circle() {
        let placement = ellipse().placement;
        let e = Ellipse::new(placement, 2.0, 2.0);
        let c = super::super::Circle::new(placement, 2.0);
        let p = placement.point(Vector3::new(1.0, -3.0, 0.7));
        let (qe, qc) = (e.project(p, None).unwrap(), c.project(p, None).unwrap());
        assert!((qe.coords - qc.coords).norm() < 1e-9);
    }
}
