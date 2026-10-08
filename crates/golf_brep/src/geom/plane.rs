use nalgebra::Matrix3x2;
use nalgebra::Vector2;
use nalgebra::Vector3;

use super::Placement;
use crate::domain::Domain;
use crate::manifold::Embedding;
use crate::manifold::ProjectError;
use crate::mapping::Mapping;
use crate::space::Point;
use crate::space::Space;
use crate::space::Uv;

/// A plane in space `S`: `u` along local x, `v` along local y, normal along local z.
#[derive(Clone, Debug)]
pub struct Plane<S: Space<3>> {
    pub placement: Placement<S>,
}

impl<S: Space<3>> Plane<S> {
    pub fn new(placement: Placement<S>) -> Self {
        Self { placement }
    }
}

impl<S: Space<3>> Space<2> for Uv<Plane<S>> {
    type Tag = ();
}

impl<S: Space<3>> Mapping<2, 3> for Plane<S> {
    type From = Uv<Plane<S>>;
    type To = S;

    fn apply(&self, uv: Point<Self::From, 2>) -> Point<S, 3> {
        self.placement.point(uv.coords.push(0.0))
    }

    fn jacobian(&self, _uv: Point<Self::From, 2>) -> Matrix3x2<f64> {
        self.placement
            .rotation_matrix()
            .fixed_columns::<2>(0)
            .into()
    }
}

impl<S: Space<3>> Embedding<2, 3> for Plane<S> {
    fn domain(&self) -> Domain<2> {
        Domain::unbounded()
    }

    fn project(
        &self,
        p: Point<S, 3>,
        _hint: Option<Point<Self::From, 2>>,
    ) -> Result<Point<Self::From, 2>, ProjectError> {
        let d: Vector3<f64> = self.placement.to_local(p);
        Ok(Point::new(Vector2::new(d.x, d.y)))
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

    fn plane() -> Plane<World> {
        Plane::new(Placement::from_axes(
            Point::new(Vector3::new(1.0, 2.0, 3.0)),
            Vector::new(Vector3::new(0.0, 1.0, 1.0)),
            Vector::new(Vector3::x()),
        ))
    }

    #[test]
    fn round_trip_and_jacobian() {
        let s = plane();
        for at in [[0.0, 0.0], [3.0, -2.0], [-1e3, 7.5]] {
            assert_round_trip(&s, at);
            assert_jacobian(&s, at);
        }
    }

    #[test]
    fn normal_is_local_z() {
        let s = plane();
        let n = s.normal(Point::new(Vector2::new(0.4, 0.1)));
        assert!((n.coords - Vector3::new(0.0, 1.0, 1.0).normalize()).norm() < 1e-12);
    }

    #[test]
    fn off_plane_points_match_newton() {
        let s = plane();
        let p = Point::new(Vector3::new(-4.0, 7.0, 0.5));
        assert_matches_newton(&s, p, [0.0, 0.0]);
    }
}
