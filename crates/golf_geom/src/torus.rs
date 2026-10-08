use core::f64::consts::PI;

use golf_manifold::Axis;
use golf_manifold::Domain;
use golf_manifold::Embedding;
use golf_manifold::HasUv;
use golf_manifold::Mapping;
use golf_manifold::Point;
use golf_manifold::ProjectError;
use golf_manifold::Space;
use golf_manifold::Uv;
use nalgebra::Matrix3x2;
use nalgebra::Vector2;
use nalgebra::Vector3;

use super::Placement;

/// A torus about local z in space `S`.
///
/// `u ∈ [-π, π)` is the angle about local z from local x, and `v ∈ [-π, π)` the
/// angle around the tube from its outermost point. Both are periodic. Assumes a
/// ring torus (`major_radius > minor_radius`).
#[derive(Clone, Debug)]
pub struct Torus<S: Space<3>> {
    pub placement: Placement<S>,
    /// Distance from the axis to the centre of the tube.
    pub major_radius: f64,
    /// Radius of the tube.
    pub minor_radius: f64,
}

impl<S: Space<3>> Torus<S> {
    pub fn new(placement: Placement<S>, major_radius: f64, minor_radius: f64) -> Self {
        Self {
            placement,
            major_radius,
            minor_radius,
        }
    }
}

impl<S: Space<3>> HasUv for Torus<S> {}

impl<S: Space<3>> Mapping<2, 3> for Torus<S> {
    type From = Uv<Torus<S>>;
    type To = S;

    fn apply(&self, uv: Point<Self::From, 2>) -> Point<S, 3> {
        let (u, v, r) = (uv.coords.x, uv.coords.y, self.minor_radius);
        let rho = self.major_radius + r * v.cos();
        self.placement
            .point(Vector3::new(rho * u.cos(), rho * u.sin(), r * v.sin()))
    }

    fn jacobian(&self, uv: Point<Self::From, 2>) -> Matrix3x2<f64> {
        let (u, v, r) = (uv.coords.x, uv.coords.y, self.minor_radius);
        let rho = self.major_radius + r * v.cos();
        let du = Vector3::new(-rho * u.sin(), rho * u.cos(), 0.0);
        let dv = Vector3::new(-r * v.sin() * u.cos(), -r * v.sin() * u.sin(), r * v.cos());
        self.placement.rotation_matrix() * Matrix3x2::from_columns(&[du, dv])
    }
}

impl<S: Space<3>> Embedding<2, 3> for Torus<S> {
    fn domain(&self) -> Domain<2> {
        Domain::new([Axis::periodic(-PI, PI), Axis::periodic(-PI, PI)])
    }

    /// On the axis every `u` is equally near, so the hint's `u` is kept; on the
    /// tube's core circle every `v` is, so the hint's `v` is kept.
    fn project(
        &self,
        p: Point<S, 3>,
        hint: Option<Point<Self::From, 2>>,
    ) -> Result<Point<Self::From, 2>, ProjectError> {
        let d = self.placement.to_local(p);
        let rho = d.xy().norm();
        let tolerance = 1e-15 * (d.norm() + self.major_radius);
        let u = match rho <= tolerance {
            true => hint.ok_or(ProjectError::Ambiguous)?.coords.x,
            false => d.y.atan2(d.x),
        };
        let (across, up) = (rho - self.major_radius, d.z);
        let v = match across.hypot(up) <= tolerance {
            true => hint.ok_or(ProjectError::Ambiguous)?.coords.y,
            false => up.atan2(across),
        };
        Ok(self.domain().wrap(Point::new(Vector2::new(u, v))))
    }
}

#[cfg(test)]
mod tests {
    use golf_manifold::Surface;
    use golf_manifold::Vector;
    use golf_manifold::World;

    use super::*;
    use crate::testing::assert_jacobian;
    use crate::testing::assert_matches_newton;
    use crate::testing::assert_round_trip;

    fn uv(u: f64, v: f64) -> Point<Uv<Torus<World>>, 2> {
        Point::new(Vector2::new(u, v))
    }

    fn torus() -> Torus<World> {
        Torus::new(
            Placement::from_axes(
                Point::new(Vector3::new(-1.0, 0.5, 2.0)),
                Vector::new(Vector3::new(1.0, 2.0, 2.0)),
                Vector::new(Vector3::z()),
            ),
            3.0,
            1.0,
        )
    }

    #[test]
    fn round_trip_and_jacobian() {
        let t = torus();
        for at in [[0.0, 0.0], [1.0, 2.0], [-PI + 1e-9, PI - 1e-9], [2.5, -1.2]] {
            assert_round_trip(&t, at);
            assert_jacobian(&t, at);
        }
    }

    #[test]
    fn normal_points_out_of_tube() {
        let t = torus();
        let at = uv(0.9, 2.2);
        let core = t
            .placement
            .point(Vector3::new(3.0 * 0.9f64.cos(), 3.0 * 0.9f64.sin(), 0.0));
        let out = (t.apply(at) - core).coords / t.minor_radius;
        assert!((t.normal(at).coords - out).norm() < 1e-12);
    }

    #[test]
    fn axis_keeps_hint_u() {
        let t = torus();
        let on_axis = t.placement.point(Vector3::new(0.0, 0.0, 0.5));
        assert_eq!(t.project(on_axis, None), Err(ProjectError::Ambiguous));
        let q = t.project(on_axis, Some(uv(1.0, 0.0))).unwrap();
        assert_eq!(q.coords.x, 1.0);
        assert!((q.coords.y - 0.5f64.atan2(-3.0)).abs() < 1e-12);
    }

    #[test]
    fn core_circle_keeps_hint_v() {
        let t = torus();
        let on_core = t.placement.point(Vector3::new(0.0, 3.0, 0.0));
        assert_eq!(t.project(on_core, None), Err(ProjectError::Ambiguous));
        let q = t.project(on_core, Some(uv(0.0, -2.0))).unwrap();
        assert!((q.coords - Vector2::new(PI / 2.0, -2.0)).norm() < 1e-12);
    }

    #[test]
    fn off_surface_points_match_newton() {
        let t = torus();
        for (local, start) in [
            (Vector3::new(6.0, 1.0, 2.0), [0.0, 0.5]),
            (Vector3::new(0.5, 2.6, 0.3), [1.5, 2.5]),
            (Vector3::new(-1.0, -1.0, -0.5), [-2.0, 2.8]),
        ] {
            assert_matches_newton(&t, t.placement.point(local), start);
        }
    }
}
