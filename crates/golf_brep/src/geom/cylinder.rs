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

/// A cylinder about local z in space `S`: `u ∈ [-π, π)` is the angle from local x,
/// `v` the height along local z.
#[derive(Clone, Debug)]
pub struct Cylinder<S: Space<3>> {
    pub placement: Placement<S>,
    pub radius: f64,
}

impl<S: Space<3>> Cylinder<S> {
    pub fn new(placement: Placement<S>, radius: f64) -> Self {
        Self { placement, radius }
    }
}

impl<S: Space<3>> Space<2> for Uv<Cylinder<S>> {
    type Tag = ();
}

impl<S: Space<3>> Mapping<2, 3> for Cylinder<S> {
    type From = Uv<Cylinder<S>>;
    type To = S;

    fn apply(&self, uv: Point<Self::From, 2>) -> Point<S, 3> {
        let (u, v, r) = (uv.coords.x, uv.coords.y, self.radius);
        self.placement
            .point(Vector3::new(r * u.cos(), r * u.sin(), v))
    }

    fn jacobian(&self, uv: Point<Self::From, 2>) -> Matrix3x2<f64> {
        let (u, r) = (uv.coords.x, self.radius);
        let du = Vector3::new(-r * u.sin(), r * u.cos(), 0.0);
        self.placement.rotation_matrix() * Matrix3x2::from_columns(&[du, Vector3::z()])
    }
}

impl<S: Space<3>> Embedding<2, 3> for Cylinder<S> {
    fn domain(&self) -> Domain<2> {
        Domain::new([Axis::periodic(-PI, PI), Axis::UNBOUNDED])
    }

    /// On the axis every `u` is equally near, so the hint's `u` is kept.
    fn project(
        &self,
        p: Point<S, 3>,
        hint: Option<Point<Self::From, 2>>,
    ) -> Result<Point<Self::From, 2>, ProjectError> {
        let d = self.placement.to_local(p);
        let u = match d.xy().norm() <= 1e-15 * (d.norm() + self.radius) {
            true => hint.ok_or(ProjectError::Ambiguous)?.coords.x,
            false => d.y.atan2(d.x),
        };
        Ok(self.domain().wrap(Point::new(Vector2::new(u, d.z))))
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::Isometry3;
    use nalgebra::Vector3;

    use super::*;
    use crate::frame::FrameTree;
    use crate::geom::testing::assert_jacobian;
    use crate::geom::testing::assert_matches_newton;
    use crate::geom::testing::assert_round_trip;
    use crate::manifold::Surface;
    use crate::space::Vector;
    use crate::space::World;

    fn uv(u: f64, v: f64) -> Point<Uv<Cylinder<World>>, 2> {
        Point::new(Vector2::new(u, v))
    }

    fn cylinder() -> Cylinder<World> {
        Cylinder::new(
            Placement::from_axes(
                Point::new(Vector3::new(0.5, -1.0, 2.0)),
                Vector::new(Vector3::new(1.0, 0.0, 1.0)),
                Vector::new(Vector3::y()),
            ),
            1.25,
        )
    }

    #[test]
    fn round_trip_and_jacobian() {
        let c = cylinder();
        for at in [[0.0, 0.0], [1.0, -3.0], [-PI + 1e-9, 2.0], [PI - 1e-9, 0.5]] {
            assert_round_trip(&c, at);
            assert_jacobian(&c, at);
        }
    }

    #[test]
    fn normal_points_outward() {
        let c = cylinder();
        let at = uv(0.8, 1.5);
        let axis_point = c.placement.point(Vector3::new(0.0, 0.0, 1.5));
        let radial = (c.apply(at) - axis_point).coords / c.radius;
        assert!((c.normal(at).coords - radial).norm() < 1e-12);
    }

    #[test]
    fn axis_keeps_hint_u() {
        let c = cylinder();
        let on_axis = c.placement.point(Vector3::new(0.0, 0.0, 4.0));
        assert_eq!(c.project(on_axis, None), Err(ProjectError::Ambiguous));
        let q = c.project(on_axis, Some(uv(2.0, 0.0))).unwrap();
        assert!((q.coords - Vector2::new(2.0, 4.0)).norm() < 1e-12);
    }

    #[test]
    fn off_surface_points_match_newton() {
        let c = cylinder();
        for (local, start) in [
            (Vector3::new(3.0, 1.0, -2.0), [0.5, -1.5]),
            (Vector3::new(0.2, -0.3, 1.0), [-1.2, 1.2]),
            (Vector3::new(-20.0, 0.1, 5.0), [3.0, 4.0]),
        ] {
            assert_matches_newton(&c, c.placement.point(local), start);
        }
    }

    #[test]
    fn tags_flow_through_frames() {
        let mut tree = FrameTree::new();
        let part = tree.add(tree.root(), Isometry3::translation(1.0, 2.0, 3.0));
        let c = Cylinder::new(Placement::at(part.point(Vector3::zeros())), 2.0);
        let p = c.apply(Point::new(Vector2::new(0.0, 1.0)));
        assert_eq!(p.tag(), part);
        assert_eq!(c.normal(Point::new(Vector2::new(0.0, 1.0))).tag(), part);
        let q = c
            .project(part.point(Vector3::new(5.0, 0.0, 1.0)), None)
            .unwrap();
        assert!((q.coords - Vector2::new(0.0, 1.0)).norm() < 1e-12);
    }

    #[test]
    #[should_panic(expected = "different instances")]
    fn projecting_from_another_frame_panics() {
        let mut tree = FrameTree::new();
        let part = tree.add(tree.root(), Isometry3::identity());
        let c = Cylinder::new(Placement::at(part.point(Vector3::zeros())), 2.0);
        let _ = c.project(tree.root().point(Vector3::x()), None);
    }
}
