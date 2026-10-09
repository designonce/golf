//! Moving bodies rigidly.

use golf_geom::Transform;
use golf_geom::TransformError;
use golf_manifold::Space;
use nalgebra::Isometry3;

use crate::Body;

impl<S: Space<3>> Body<S> {
    /// The body moved by `motion`, within its own space.
    ///
    /// Vertices, edge curves and face surfaces all move; parametrisations are
    /// kept, so edge ranges and pcurves carry over unchanged.
    pub fn transformed(&self, motion: &Isometry3<f64>) -> Result<Self, TransformError> {
        let mut moved = self.clone();
        moved.transform_in_place(motion)?;
        Ok(moved)
    }

    /// [`Self::transformed`], in place. On error the body is unchanged.
    pub fn transform_in_place(&mut self, motion: &Isometry3<f64>) -> Result<(), TransformError> {
        // Move everything into new storage first, so a failure leaves the body
        // as it was.
        let curves = self
            .edges
            .iter()
            .map(|edge| edge.curve.transformed(motion))
            .collect::<Result<Vec<_>, _>>()?;
        let surfaces = self
            .faces
            .iter()
            .map(|face| face.surface.transformed(motion))
            .collect::<Result<Vec<_>, _>>()?;
        for vertex in &mut self.vertices {
            vertex.point = vertex.point.transformed(motion)?;
        }
        for (edge, curve) in self.edges.iter_mut().zip(curves) {
            edge.curve = curve;
        }
        for (face, surface) in self.faces.iter_mut().zip(surfaces) {
            face.surface = surface;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use golf_geom::AnySurface;
    use golf_geom::Placement;
    use golf_geom::Plane;
    use golf_manifold::Point;
    use golf_manifold::World;
    use nalgebra::Translation3;
    use nalgebra::UnitQuaternion;
    use nalgebra::Vector3;

    use super::*;
    use crate::cuboid;
    use crate::cylinder;

    fn motion() -> Isometry3<f64> {
        Isometry3::from_parts(
            Translation3::new(3.0, -1.0, 2.5),
            UnitQuaternion::from_scaled_axis(Vector3::new(0.4, -0.7, 1.1)),
        )
    }

    fn origin() -> Placement<World> {
        Placement::at(Point::new(Vector3::zeros()))
    }

    #[test]
    fn moved_bodies_stay_valid_with_vertices_moved() {
        for body in [
            cuboid(origin(), Vector3::new(1.0, 2.0, 3.0)),
            cylinder(origin(), 1.5, 4.0),
        ] {
            let moved = body.transformed(&motion()).unwrap();
            assert_eq!(moved.validate(1e-9), Ok(()));
            for ((_, before), (_, after)) in body.vertices().zip(moved.vertices()) {
                let expected = motion().transform_point(&before.point.coords.into()).coords;
                assert!((after.point.coords - expected).norm() < 1e-12);
            }
        }
    }

    #[cfg(feature = "mesh")]
    #[test]
    fn moved_bodies_mesh_to_the_same_volume() {
        use golf_mesh::Tolerance;
        use golf_mesh::mesh;
        let tolerance = Tolerance::new(1e-3, 0.3);
        for body in [
            cuboid(origin(), Vector3::new(1.0, 2.0, 3.0)),
            cylinder(origin(), 1.5, 4.0),
        ] {
            let before = mesh(&body, &tolerance).unwrap();
            let after = mesh(&body.transformed(&motion()).unwrap(), &tolerance).unwrap();
            assert!(after.is_watertight());
            let (v0, v1) = (before.signed_volume(), after.signed_volume());
            assert!((v0 - v1).abs() < 1e-6 * v0, "{v0} vs {v1}");
        }
    }

    #[test]
    fn failure_leaves_the_body_unchanged() {
        let mut body = cuboid(origin(), Vector3::new(1.0, 1.0, 1.0));
        body.faces[2].surface = AnySurface::Custom(Arc::new(Plane::new(origin())));
        let before: Vec<_> = body.vertices().map(|(_, v)| v.point.coords).collect();
        assert!(body.transform_in_place(&motion()).is_err());
        let after: Vec<_> = body.vertices().map(|(_, v)| v.point.coords).collect();
        assert_eq!(before, after);
    }
}
