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

/// A sphere in space `S`, parametrised by longitude `u ∈ [-π, π)` about local z
/// from local x, and latitude `v ∈ [-π/2, π/2]`.
///
/// `u` is periodic; the poles `v = ±π/2` (on local z) are singular.
#[derive(Clone, Debug)]
pub struct Sphere<S: Space<3>> {
    pub placement: Placement<S>,
    pub radius: f64,
}

impl<S: Space<3>> Sphere<S> {
    /// A sphere about `center`, with its poles along `S`'s z axis.
    pub fn new(center: Point<S, 3>, radius: f64) -> Self {
        Self::from_placement(Placement::at(center), radius)
    }

    pub fn from_placement(placement: Placement<S>, radius: f64) -> Self {
        Self { placement, radius }
    }

    pub fn center(&self) -> Point<S, 3> {
        self.placement.origin
    }
}

impl<S: Space<3>> Space<2> for Uv<Sphere<S>> {
    type Tag = ();
}

impl<S: Space<3>> Mapping<2, 3> for Sphere<S> {
    type From = Uv<Sphere<S>>;
    type To = S;

    fn apply(&self, uv: Point<Self::From, 2>) -> Point<S, 3> {
        let (u, v) = (uv.coords.x, uv.coords.y);
        let dir = Vector3::new(v.cos() * u.cos(), v.cos() * u.sin(), v.sin());
        self.placement.point(dir * self.radius)
    }

    fn jacobian(&self, uv: Point<Self::From, 2>) -> Matrix3x2<f64> {
        let (u, v, r) = (uv.coords.x, uv.coords.y, self.radius);
        let du = Vector3::new(-v.cos() * u.sin(), v.cos() * u.cos(), 0.0) * r;
        let dv = Vector3::new(-v.sin() * u.cos(), -v.sin() * u.sin(), v.cos()) * r;
        self.placement.rotation_matrix() * Matrix3x2::from_columns(&[du, dv])
    }
}

impl<S: Space<3>> Embedding<2, 3> for Sphere<S> {
    fn domain(&self) -> Domain<2> {
        Domain::new([
            Axis::periodic(-PI, PI),
            Axis::bounded(-FRAC_PI_2, FRAC_PI_2)
                .singular_at_min()
                .singular_at_max(),
        ])
    }

    /// At the centre every point is equally near, so the hint is returned as is.
    /// On the axis through the poles the hint's longitude is kept.
    fn project(
        &self,
        p: Point<S, 3>,
        hint: Option<Point<Self::From, 2>>,
    ) -> Result<Point<Self::From, 2>, ProjectError> {
        let d = self.placement.to_local(p);
        let r = d.norm();
        if r == 0.0 {
            return hint
                .map(|h| self.domain().wrap(h))
                .ok_or(ProjectError::Ambiguous);
        }
        let v = (d.z / r).clamp(-1.0, 1.0).asin();
        let u = match d.xy().norm() <= 1e-15 * r {
            true => hint.map_or(0.0, |h| h.coords.x),
            false => d.y.atan2(d.x),
        };
        Ok(self.domain().wrap(Point::new(Vector2::new(u, v))))
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::Isometry3;
    use nalgebra::Translation3;
    use nalgebra::UnitQuaternion;

    use super::*;
    use crate::frame::FrameTree;
    use crate::geom::testing::assert_jacobian;
    use crate::geom::testing::assert_round_trip;
    use crate::manifold::Surface;
    use crate::manifold::newton_project;
    use crate::mapping::Compose;
    use crate::space::Vector;
    use crate::space::World;

    fn uv<S: Space<3>>(u: f64, v: f64) -> Point<Uv<Sphere<S>>, 2> {
        Point::new(Vector2::new(u, v))
    }

    fn sphere(center: [f64; 3], radius: f64) -> Sphere<World> {
        Sphere::new(Point::new(Vector3::from(center)), radius)
    }

    #[test]
    fn project_inverts_apply() {
        let s = sphere([1.0, 2.0, 3.0], 2.0);
        let p = uv(0.7, -0.3);
        assert!((s.project(s.apply(p), None).unwrap().coords - p.coords).norm() < 1e-12);
    }

    #[test]
    fn project_wraps_into_domain() {
        let s = sphere([0.0; 3], 1.0);
        let q = s.project(s.apply(uv(3.0 * PI - 0.2, 0.1)), None).unwrap();
        assert!((q.coords - Vector2::new(PI - 0.2, 0.1)).norm() < 1e-12);
    }

    #[test]
    fn project_keeps_hint_longitude_at_pole() {
        let s = sphere([1.0, 1.0, 1.0], 2.0);
        let pole = Point::new(Vector3::new(1.0, 1.0, 5.0));
        let q = s.project(pole, Some(uv(1.2, 0.0))).unwrap();
        assert!((q.coords - Vector2::new(1.2, FRAC_PI_2)).norm() < 1e-12);
        assert!(s.domain().is_singular(q, 1e-12));
    }

    #[test]
    fn project_from_center_needs_hint() {
        let s = sphere([1.0, 2.0, 3.0], 2.0);
        assert_eq!(s.project(s.center(), None), Err(ProjectError::Ambiguous));
        let hint = uv(0.5, 0.5);
        assert_eq!(s.project(s.center(), Some(hint)), Ok(hint));
    }

    #[test]
    fn project_within_rejects_distant_points() {
        let s = sphere([0.0; 3], 1.0);
        let p = Point::new(Vector3::new(3.0, 0.0, 0.0));
        assert!(s.project_within(p, None, 2.5).is_ok());
        assert!(matches!(
            s.project_within(p, None, 1.0),
            Err(ProjectError::TooFar { distance, .. }) if (distance - 2.0).abs() < 1e-12
        ));
    }

    #[test]
    fn newton_matches_closed_form() {
        let s = sphere([1.0, -2.0, 0.5], 1.5);
        let p = Point::new(Vector3::new(3.0, 1.0, 2.0));
        let exact = s.project(p, None).unwrap();
        let start = uv(exact.coords.x + 0.3, exact.coords.y - 0.2);
        let newton = newton_project(&s, p, start).unwrap();
        assert!((newton.coords - exact.coords).norm() < 1e-10);
    }

    #[test]
    fn newton_matches_closed_form_inside_and_far_away() {
        let s = sphere([0.5, 0.0, -1.0], 2.0);
        for distance in [0.3, 1.9, 2.1, 5.0, 40.0] {
            for (u, v) in [(0.1, 0.2), (2.5, -1.0), (-3.0, 1.3), (-1.2, -0.4)] {
                let dir = (s.apply(uv(u, v)) - s.center()).coords / s.radius;
                let p = Point::new(s.center().coords + dir * distance);
                let newton = newton_project(&s, p, uv(u + 0.25, v - 0.15)).unwrap();
                let exact = s.project(p, None).unwrap();
                assert!(
                    (newton.coords - exact.coords).norm() < 1e-9,
                    "distance {distance}, uv ({u}, {v}): {newton:?} != {exact:?}"
                );
            }
        }
    }

    #[test]
    fn newton_converges_at_pole() {
        let s = sphere([0.0; 3], 1.0);
        let p = Point::new(Vector3::new(0.0, 0.0, 4.0));
        let q = newton_project(&s, p, uv(0.7, 1.2)).unwrap();
        assert!((q.coords.y - FRAC_PI_2).abs() < 1e-10);
    }

    #[test]
    fn rotated_sphere_round_trips_and_matches_finite_difference() {
        let s = Sphere::from_placement(
            Placement::<World>::from_axes(
                Point::new(Vector3::new(1.0, -1.0, 2.0)),
                Vector::new(Vector3::new(1.0, 1.0, 1.0)),
                Vector::new(Vector3::x()),
            ),
            1.5,
        );
        for at in [
            [0.3, 0.2],
            [-PI + 1e-9, -1.0],
            [PI - 1e-9, 1.2],
            [2.0, -1.5],
        ] {
            assert_round_trip(&s, at);
            assert_jacobian(&s, at);
        }
        let pole = s.apply(uv(0.0, FRAC_PI_2));
        let q = s.project(pole, Some(uv(-2.0, 0.0))).unwrap();
        assert!((q.coords - Vector2::new(-2.0, FRAC_PI_2)).norm() < 1e-9);
    }

    #[test]
    fn jacobian_matches_finite_difference() {
        let s = sphere([0.0; 3], 1.5);
        let p = uv(0.4, 0.9);
        let h = 1e-6;
        let j = s.jacobian(p);
        let du = (s.apply(uv(0.4 + h, 0.9)) - s.apply(uv(0.4 - h, 0.9))).coords / (2.0 * h);
        let dv = (s.apply(uv(0.4, 0.9 + h)) - s.apply(uv(0.4, 0.9 - h))).coords / (2.0 * h);
        assert!((j.column(0) - du).norm() < 1e-8);
        assert!((j.column(1) - dv).norm() < 1e-8);
    }

    #[test]
    fn normal_points_outward() {
        let s = sphere([1.0, 0.0, 0.0], 3.0);
        let p = uv(1.1, 0.2);
        let radial = (s.apply(p) - s.center()).coords.normalize();
        assert!((s.normal(p).coords - radial).norm() < 1e-12);
    }

    #[test]
    fn sphere_in_frame_composes_with_placement() {
        let mut tree = FrameTree::new();
        let placement = Isometry3::from_parts(
            Translation3::new(5.0, 0.0, 0.0),
            UnitQuaternion::from_scaled_axis(Vector3::new(0.0, 0.0, 1.0)),
        );
        let part = tree.add(tree.root(), placement);
        let s = Sphere::new(part.point(Vector3::new(0.0, 1.0, 0.0)), 2.0);
        let placed = Compose::<_, _, 3>(s, tree.transform(part, tree.root()));

        let p = uv(0.3, 0.4);
        let world = placed.apply(p);
        assert_eq!(world.tag(), tree.root());
        let center = tree.transform(part, tree.root()).apply(placed.0.center());
        assert!(((world - center).coords.norm() - 2.0).abs() < 1e-12);
        // Normals rotate with the frame but stay radial.
        let normal = placed.0.normal(p);
        let rotated = placed.1.push_forward(placed.0.apply(p), normal);
        assert!((rotated.coords - (world - center).coords / 2.0).norm() < 1e-12);
    }
}
