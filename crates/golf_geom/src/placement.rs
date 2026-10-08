use golf_manifold::Point;
use golf_manifold::Space;
use golf_manifold::Vector;
use nalgebra::Matrix3;
use nalgebra::Rotation3;
use nalgebra::UnitQuaternion;
use nalgebra::Vector3;

/// An orthonormal frame in space `S`: where a shape's canonical local coordinates sit.
///
/// `rotation` takes local axes to `S`; `origin` is the local origin.
pub struct Placement<S: Space<3>> {
    pub origin: Point<S, 3>,
    pub rotation: UnitQuaternion<f64>,
}

impl<S: Space<3>> Clone for Placement<S> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<S: Space<3>> Copy for Placement<S> {}

impl<S: Space<3>> core::fmt::Debug for Placement<S> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("Placement")
            .field("origin", &self.origin)
            .field("rotation", &self.rotation)
            .finish()
    }
}

impl<S: Space<3>> Placement<S> {
    pub fn new(origin: Point<S, 3>, rotation: UnitQuaternion<f64>) -> Self {
        Self { origin, rotation }
    }

    /// Local axes aligned with `S`'s.
    pub fn at(origin: Point<S, 3>) -> Self {
        Self::new(origin, UnitQuaternion::identity())
    }

    /// Local z along `z_axis`, and local x along the part of `x_axis` perpendicular to it.
    ///
    /// # Panics
    /// If `z_axis` is zero or `x_axis` is parallel to it.
    pub fn from_axes(origin: Point<S, 3>, z_axis: Vector<S, 3>, x_axis: Vector<S, 3>) -> Self {
        assert_eq!(z_axis.tag(), origin.tag(), "z axis is in a different space");
        assert_eq!(x_axis.tag(), origin.tag(), "x axis is in a different space");
        let z = z_axis
            .coords
            .try_normalize(f64::MIN_POSITIVE)
            .expect("z axis is zero");
        let x = (x_axis.coords - z * z.dot(&x_axis.coords))
            .try_normalize(1e-12 * x_axis.coords.norm())
            .expect("x axis is parallel to z axis");
        let y = z.cross(&x);
        let rotation = Rotation3::from_matrix_unchecked(Matrix3::from_columns(&[x, y, z]));
        Self::new(origin, UnitQuaternion::from_rotation_matrix(&rotation))
    }

    /// The point at local coordinates `local`.
    pub fn point(&self, local: Vector3<f64>) -> Point<S, 3> {
        Point::with_tag(
            self.origin.coords + self.rotation * local,
            self.origin.tag(),
        )
    }

    /// The local vector `local`, in `S`.
    pub fn vector(&self, local: Vector3<f64>) -> Vector<S, 3> {
        Vector::with_tag(self.rotation * local, self.origin.tag())
    }

    /// Local coordinates of `p`. Panics if `p` is in a different instance of `S`.
    pub fn to_local(&self, p: Point<S, 3>) -> Vector3<f64> {
        self.rotation
            .inverse_transform_vector(&(p - self.origin).coords)
    }

    /// Local components of `v`. Panics if `v` is in a different instance of `S`.
    pub fn vector_to_local(&self, v: Vector<S, 3>) -> Vector3<f64> {
        assert_eq!(v.tag(), self.origin.tag(), "vector is in a different space");
        self.rotation.inverse_transform_vector(&v.coords)
    }

    /// Columns are the local x, y and z axes in `S`.
    pub fn rotation_matrix(&self) -> Matrix3<f64> {
        self.rotation.to_rotation_matrix().into_inner()
    }
}

#[cfg(test)]
mod tests {
    use golf_manifold::World;

    use super::*;

    #[test]
    fn from_axes_orthonormalises() {
        let p = Placement::<World>::from_axes(
            Point::new(Vector3::new(1.0, 2.0, 3.0)),
            Vector::new(Vector3::new(0.0, 0.0, 2.0)),
            Vector::new(Vector3::new(1.0, 1.0, 5.0)),
        );
        let m = p.rotation_matrix();
        let x = Vector3::new(1.0, 1.0, 0.0).normalize();
        assert!((m.column(0) - x).norm() < 1e-12);
        assert!((m.column(2) - Vector3::z()).norm() < 1e-12);
        let local = Vector3::new(0.3, -0.7, 1.1);
        assert!((p.to_local(p.point(local)) - local).norm() < 1e-12);
        assert!((p.vector_to_local(p.vector(local)) - local).norm() < 1e-12);
    }

    #[test]
    #[should_panic(expected = "parallel")]
    fn from_axes_rejects_parallel_axes() {
        Placement::<World>::from_axes(
            Point::new(Vector3::zeros()),
            Vector::new(Vector3::z()),
            Vector::new(Vector3::z() * 3.0),
        );
    }
}
