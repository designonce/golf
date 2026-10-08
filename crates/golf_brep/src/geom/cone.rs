use core::f64::consts::FRAC_PI_2;
use core::f64::consts::PI;

use nalgebra::Matrix3x2;
use nalgebra::Vector2;
use nalgebra::Vector3;

use super::Placement;
use crate::domain::Axis;
use crate::domain::Domain;
use crate::manifold::Embedding;
use crate::manifold::ProjectError;
use crate::mapping::Mapping;
use crate::space::Point;
use crate::space::Space;
use crate::space::Uv;

/// A cone about local z in space `S`.
///
/// `u ∈ [-π, π)` is the angle from local x and `v` the height along local z; the
/// radius at height `v` is `radius + v·tan(half_angle)`. The domain stops at the
/// apex, which is a singular bound: below it for a positive half angle, above it
/// for a negative one. A zero half angle is a cylinder.
#[derive(Clone, Debug)]
pub struct Cone<S: Space<3>> {
    pub placement: Placement<S>,
    /// Radius at `v = 0`.
    pub radius: f64,
    /// Angle between the axis and a generator line, in `(-π/2, π/2)`.
    pub half_angle: f64,
}

impl<S: Space<3>> Cone<S> {
    pub fn new(placement: Placement<S>, radius: f64, half_angle: f64) -> Self {
        debug_assert!(
            half_angle.abs() < FRAC_PI_2,
            "half angle {half_angle} out of range"
        );
        Self {
            placement,
            radius,
            half_angle,
        }
    }

    /// Height of the apex, or `None` for a zero half angle.
    pub fn apex_height(&self) -> Option<f64> {
        let slope = self.half_angle.tan();
        (slope != 0.0).then(|| -self.radius / slope)
    }

    fn radius_at(&self, v: f64) -> f64 {
        self.radius + v * self.half_angle.tan()
    }
}

impl<S: Space<3>> Space<2> for Uv<Cone<S>> {
    type Tag = ();
}

impl<S: Space<3>> Mapping<2, 3> for Cone<S> {
    type From = Uv<Cone<S>>;
    type To = S;

    fn apply(&self, uv: Point<Self::From, 2>) -> Point<S, 3> {
        let (u, v) = (uv.coords.x, uv.coords.y);
        let rho = self.radius_at(v);
        self.placement
            .point(Vector3::new(rho * u.cos(), rho * u.sin(), v))
    }

    fn jacobian(&self, uv: Point<Self::From, 2>) -> Matrix3x2<f64> {
        let (u, v) = (uv.coords.x, uv.coords.y);
        let (rho, slope) = (self.radius_at(v), self.half_angle.tan());
        let du = Vector3::new(-rho * u.sin(), rho * u.cos(), 0.0);
        let dv = Vector3::new(slope * u.cos(), slope * u.sin(), 1.0);
        self.placement.rotation_matrix() * Matrix3x2::from_columns(&[du, dv])
    }
}

impl<S: Space<3>> Embedding<2, 3> for Cone<S> {
    fn domain(&self) -> Domain<2> {
        let v = match self.apex_height() {
            None => Axis::UNBOUNDED,
            Some(apex) if self.half_angle > 0.0 => {
                Axis::bounded(apex, f64::INFINITY).singular_at_min()
            }
            Some(apex) => Axis::bounded(f64::NEG_INFINITY, apex).singular_at_max(),
        };
        Domain::new([Axis::periodic(-PI, PI), v])
    }

    /// Projects onto the generator line in the half-plane through `p`, clamped at
    /// the apex. On the axis every `u` is equally near, so the hint's `u` is kept;
    /// if the nearest point is the apex any `u` will do and `0` is used without one.
    fn project(
        &self,
        p: Point<S, 3>,
        hint: Option<Point<Self::From, 2>>,
    ) -> Result<Point<Self::From, 2>, ProjectError> {
        let d = self.placement.to_local(p);
        let (rho, slope) = (d.xy().norm(), self.half_angle.tan());
        let domain = self.domain();
        let v = domain.axes[1].clamp(((rho - self.radius) * slope + d.z) / (1.0 + slope * slope));
        let u = match rho <= 1e-15 * (d.norm() + self.radius) {
            true if Some(v) == self.apex_height() => hint.map_or(0.0, |h| h.coords.x),
            true => hint.ok_or(ProjectError::Ambiguous)?.coords.x,
            false => d.y.atan2(d.x),
        };
        Ok(domain.wrap(Point::new(Vector2::new(u, v))))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::testing::assert_jacobian;
    use crate::geom::testing::assert_matches_newton;
    use crate::geom::testing::assert_round_trip;
    use crate::manifold::Surface;
    use crate::space::Vector;
    use crate::space::World;

    fn uv(u: f64, v: f64) -> Point<Uv<Cone<World>>, 2> {
        Point::new(Vector2::new(u, v))
    }

    fn cone(half_angle: f64) -> Cone<World> {
        Cone::new(
            Placement::from_axes(
                Point::new(Vector3::new(1.0, 1.0, -1.0)),
                Vector::new(Vector3::new(0.0, -1.0, 2.0)),
                Vector::new(Vector3::x()),
            ),
            2.0,
            half_angle,
        )
    }

    #[test]
    fn round_trip_and_jacobian() {
        for half_angle in [0.4, -0.4, 0.0] {
            let c = cone(half_angle);
            for at in [[0.0, 0.0], [1.0, 1.5], [-PI + 1e-9, -1.0], [PI - 1e-9, 2.0]] {
                assert_round_trip(&c, at);
                assert_jacobian(&c, at);
            }
        }
    }

    #[test]
    fn domain_stops_at_apex() {
        let up = cone(0.4).domain().axes[1];
        let apex = -2.0 / 0.4f64.tan();
        assert!((up.min - apex).abs() < 1e-12 && up.max == f64::INFINITY && up.singular_at_min);
        let down = cone(-0.4).domain().axes[1];
        assert!((down.max + apex).abs() < 1e-12 && down.singular_at_max);
        assert_eq!(cone(0.0).domain().axes[1], Axis::UNBOUNDED);
    }

    #[test]
    fn normal_points_away_from_axis() {
        let c = cone(0.3);
        let at = uv(0.7, 1.0);
        let local = c.placement.vector_to_local(c.normal(at));
        let radial = Vector2::new(0.7f64.cos(), 0.7f64.sin());
        assert!(local.xy().dot(&radial) > 0.0);
        assert!((local.z + 0.3f64.sin()).abs() < 1e-12);
    }

    #[test]
    fn nearest_point_is_not_radial() {
        let c = cone(0.5);
        let p = c.placement.point(Vector3::new(5.0, 0.0, 0.0));
        let q = c.project(p, None).unwrap();
        // Radial projection would give v = 0; the generator is nearer higher up.
        assert!(q.coords.y > 0.1);
        assert_matches_newton(&c, p, [0.2, 1.0]);
    }

    #[test]
    fn points_past_the_apex_clamp_to_it() {
        let c = cone(0.5);
        let apex = c.apex_height().unwrap();
        let below = c.placement.point(Vector3::new(0.1, 0.2, apex - 3.0));
        let q = c.project(below, None).unwrap();
        assert_eq!(q.coords.y, apex);
        assert!(c.domain().is_singular(q, 0.0));
        let on_axis = c.placement.point(Vector3::new(0.0, 0.0, apex - 1.0));
        assert_eq!(
            c.project(on_axis, None).unwrap().coords,
            Vector2::new(0.0, apex)
        );
    }

    #[test]
    fn axis_inside_cone_keeps_hint_u() {
        let c = cone(0.5);
        let on_axis = c.placement.point(Vector3::new(0.0, 0.0, 1.0));
        assert_eq!(c.project(on_axis, None), Err(ProjectError::Ambiguous));
        assert_eq!(
            c.project(on_axis, Some(uv(1.5, 0.0))).unwrap().coords.x,
            1.5
        );
    }

    #[test]
    fn off_surface_points_match_newton() {
        for half_angle in [0.3, -0.6] {
            let c = cone(half_angle);
            for (local, start) in [
                (Vector3::new(4.0, 1.0, 0.5), [0.1, 0.5]),
                (Vector3::new(-1.0, 0.5, 0.2), [2.8, 0.0]),
            ] {
                assert_matches_newton(&c, c.placement.point(local), start);
            }
        }
    }
}
