//! Rigid motions of geometry within its own space.

use golf_manifold::Point;
use golf_manifold::Space;
use golf_manifold::Vector;
use nalgebra::Isometry3;

use crate::circle::Circle;
use crate::cone::Cone;
use crate::cylinder::Cylinder;
use crate::ellipse::Ellipse;
use crate::erased::AnyCurve;
use crate::erased::AnySurface;
use crate::line::Line;
use crate::nurbs::NurbsCurve;
use crate::nurbs::NurbsSurface;
use crate::placement::Placement;
use crate::plane::Plane;
use crate::sphere::Sphere;
use crate::torus::Torus;

/// Why geometry couldn't be moved.
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum TransformError {
    /// Geometry the transform doesn't know how to move, such as a `Custom`
    /// erased curve or surface.
    #[error("{0} can't be transformed")]
    Unsupported(String),
}

/// Geometry that can be moved rigidly.
///
/// The motion acts on coordinates within the geometry's own space, so tags
/// are unchanged. Parametrisations are kept exactly: the moved geometry at a
/// parameter is the motion applied to the original at that parameter, so
/// pcurves and parameter ranges stay valid.
pub trait Transform {
    fn transformed(&self, motion: &Isometry3<f64>) -> Result<Self, TransformError>
    where
        Self: Sized;
}

impl<S: Space<3>> Transform for Point<S, 3> {
    fn transformed(&self, motion: &Isometry3<f64>) -> Result<Self, TransformError> {
        Ok(Point::with_tag(
            motion.transform_point(&self.coords.into()).coords,
            self.tag(),
        ))
    }
}

impl<S: Space<3>> Transform for Vector<S, 3> {
    fn transformed(&self, motion: &Isometry3<f64>) -> Result<Self, TransformError> {
        Ok(Vector::with_tag(
            motion.transform_vector(&self.coords),
            self.tag(),
        ))
    }
}

impl<S: Space<3>> Transform for Placement<S> {
    fn transformed(&self, motion: &Isometry3<f64>) -> Result<Self, TransformError> {
        Ok(Placement::new(
            self.origin.transformed(motion)?,
            motion.rotation * self.rotation,
        ))
    }
}

impl<S: Space<3>> Transform for Line<S, 3> {
    fn transformed(&self, motion: &Isometry3<f64>) -> Result<Self, TransformError> {
        Ok(Line::new(
            self.origin.transformed(motion)?,
            self.direction.transformed(motion)?,
        ))
    }
}

/// Geometry that is a placement and some fixed sizes: moving it moves only
/// the placement.
macro_rules! placed {
    ($($ty:ident),*) => {$(
        impl<S: Space<3>> Transform for $ty<S> {
            fn transformed(&self, motion: &Isometry3<f64>) -> Result<Self, TransformError> {
                let mut moved = self.clone();
                moved.placement = self.placement.transformed(motion)?;
                Ok(moved)
            }
        }
    )*};
}

placed!(Circle, Ellipse, Plane, Sphere, Cylinder, Cone, Torus);

impl<S: Space<3>> Transform for NurbsCurve<S, 3> {
    /// NURBS are affinely invariant, so moving the control points moves the
    /// curve.
    fn transformed(&self, motion: &Isometry3<f64>) -> Result<Self, TransformError> {
        let moved = self
            .geometry()
            .map_points(|p| motion.transform_point(&p.into()).coords);
        Ok(NurbsCurve::with_tag(moved, self.tag()))
    }
}

impl<S: Space<3>> Transform for NurbsSurface<S> {
    fn transformed(&self, motion: &Isometry3<f64>) -> Result<Self, TransformError> {
        let moved = self
            .geometry()
            .map_points(|p| motion.transform_point(&p.into()).coords);
        Ok(NurbsSurface::with_tag(moved, self.tag()))
    }
}

impl<S: Space<3>> Transform for AnyCurve<S> {
    fn transformed(&self, motion: &Isometry3<f64>) -> Result<Self, TransformError> {
        Ok(match self {
            Self::Line(c) => c.transformed(motion)?.into(),
            Self::Circle(c) => c.transformed(motion)?.into(),
            Self::Ellipse(c) => c.transformed(motion)?.into(),
            Self::Nurbs(c) => c.transformed(motion)?.into(),
            other => return Err(TransformError::Unsupported(format!("{other:?}"))),
        })
    }
}

impl<S: Space<3>> Transform for AnySurface<S> {
    fn transformed(&self, motion: &Isometry3<f64>) -> Result<Self, TransformError> {
        Ok(match self {
            Self::Plane(s) => s.transformed(motion)?.into(),
            Self::Sphere(s) => s.transformed(motion)?.into(),
            Self::Cylinder(s) => s.transformed(motion)?.into(),
            Self::Cone(s) => s.transformed(motion)?.into(),
            Self::Torus(s) => s.transformed(motion)?.into(),
            Self::Nurbs(s) => s.transformed(motion)?.into(),
            other => return Err(TransformError::Unsupported(format!("{other:?}"))),
        })
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use golf_manifold::Mapping;
    use golf_manifold::World;
    use nalgebra::SVector;
    use nalgebra::Translation3;
    use nalgebra::UnitQuaternion;
    use nalgebra::Vector3;

    use super::*;

    fn motion() -> Isometry3<f64> {
        Isometry3::from_parts(
            Translation3::new(3.0, -1.0, 2.5),
            UnitQuaternion::from_scaled_axis(Vector3::new(0.4, -0.7, 1.1)),
        )
    }

    fn placement() -> Placement<World> {
        Placement::from_axes(
            Point::new(Vector3::new(1.0, 2.0, -0.5)),
            Vector::new(Vector3::new(0.2, -0.3, 1.0)),
            Vector::new(Vector3::x()),
        )
    }

    /// Moving then evaluating equals evaluating then moving, for points and
    /// for the jacobian's columns (which only rotate).
    fn assert_moves<G, const N: usize>(original: &G, params: &[[f64; N]])
    where
        G: Transform + Mapping<N, 3, To = World>,
        G::From: Space<N, Tag = ()>,
    {
        let m = motion();
        let moved = original.transformed(&m).unwrap();
        for &at in params {
            let p = Point::new(SVector::from(at));
            let expected = m.transform_point(&original.apply(p).coords.into()).coords;
            assert!((moved.apply(p).coords - expected).norm() < 1e-12, "{at:?}");
            let expected_jacobian = m.rotation.to_rotation_matrix() * original.jacobian(p);
            assert!(
                (moved.jacobian(p) - expected_jacobian).norm() < 1e-12,
                "{at:?}"
            );
        }
    }

    const UV: [[f64; 2]; 4] = [[0.0, 0.0], [0.7, -0.3], [-2.5, 0.9], [3.0, 1.4]];
    const T: [[f64; 1]; 4] = [[0.0], [0.4], [-1.7], [2.9]];

    #[test]
    fn curves_move_with_their_parametrisation() {
        assert_moves(
            &Line::new(placement().origin, Vector::new(Vector3::new(1.0, 2.0, 3.0))),
            &T,
        );
        assert_moves(&Circle::new(placement(), 1.5), &T);
        assert_moves(&Ellipse::new(placement(), 2.0, 0.5), &T);
        let arc = golf_nurbs::NurbsCurve::circular_arc(
            Vector3::zeros(),
            Vector3::x(),
            Vector3::y(),
            2.0,
            0.3,
            2.0,
        );
        assert_moves(
            &NurbsCurve::<World, 3>::new(arc),
            &[[0.0], [0.3], [0.8], [1.0]],
        );
        let any: AnyCurve<World> = Circle::new(placement(), 1.0).into();
        assert_moves(&any, &T);
    }

    #[test]
    fn surfaces_move_with_their_parametrisation() {
        assert_moves(&Plane::new(placement()), &UV);
        assert_moves(&Sphere::from_placement(placement(), 2.0), &UV);
        assert_moves(&Cylinder::new(placement(), 1.5), &UV);
        assert_moves(&Cone::new(placement(), 1.0, 0.3), &UV);
        assert_moves(&Torus::new(placement(), 3.0, 1.0), &UV);
        let patch = golf_nurbs::NurbsSurface::new(
            (1, vec![0.0, 0.0, 1.0, 1.0]),
            (1, vec![0.0, 0.0, 1.0, 1.0]),
            vec![
                vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.0, 1.0, 0.5)],
                vec![Vector3::new(1.0, 0.0, 0.2), Vector3::new(1.0, 1.0, 1.0)],
            ],
            Some(vec![vec![1.0, 2.0], vec![0.5, 1.0]]),
        )
        .unwrap();
        assert_moves(
            &NurbsSurface::<World>::new(patch),
            &[[0.0, 0.0], [0.3, 0.6], [1.0, 1.0]],
        );
        let any: AnySurface<World> = Torus::new(placement(), 3.0, 1.0).into();
        assert_moves(&any, &UV);
    }

    #[test]
    fn moved_nurbs_keeps_its_domain() {
        let full = golf_nurbs::NurbsCurve::circular_arc(
            Vector3::zeros(),
            Vector3::x(),
            Vector3::y(),
            1.0,
            0.0,
            0.0,
        );
        let curve = NurbsCurve::<World, 3>::new(full);
        let moved = curve.transformed(&motion()).unwrap();
        use golf_manifold::Embedding;
        assert_eq!(moved.domain(), curve.domain());
        assert!(moved.domain().axes[0].periodic);
    }

    #[test]
    fn custom_geometry_is_unsupported() {
        let custom = AnyCurve::<World>::Custom(Arc::new(Circle::new(placement(), 1.0)));
        assert!(matches!(
            custom.transformed(&motion()),
            Err(TransformError::Unsupported(_))
        ));
    }
}
