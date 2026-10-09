//! Meshing bodies with `golf_mesh`.

use golf_geom::AnyCurve;
use golf_geom::AnyCurve2;
use golf_geom::AnySurface;
use golf_manifold::Point;
use golf_manifold::Space;
use golf_mesh::CoedgeData;
use golf_mesh::EdgeData;
use golf_mesh::FaceData;
use golf_mesh::MeshSource;

use crate::Body;
use crate::EdgeId;
use crate::FaceId;
use crate::FaceUv;
use crate::VertexId;

impl<S: Space<3>> MeshSource for Body<S> {
    type Space = S;
    type Curve = AnyCurve<S>;
    type Surface = AnySurface<S>;
    type Pcurve = AnyCurve2<FaceUv<S>>;
    type Vertex = VertexId;
    type Edge = EdgeId;
    type Face = FaceId;

    fn vertices(&self) -> impl Iterator<Item = (VertexId, Point<S, 3>)> + '_ {
        Body::vertices(self).map(|(id, vertex)| (id, vertex.point))
    }

    fn edges(&self) -> impl Iterator<Item = (EdgeId, EdgeData<'_, Self>)> + '_ {
        Body::edges(self).map(|(id, edge)| {
            let data = EdgeData {
                curve: &edge.curve,
                range: edge.range,
                start: edge.start,
                end: edge.end,
            };
            (id, data)
        })
    }

    fn faces(&self) -> impl Iterator<Item = (FaceId, FaceData<'_, Self>)> + '_ {
        Body::faces(self).map(|(id, face)| {
            let loops = face
                .loops
                .iter()
                .map(|l| {
                    l.coedges
                        .iter()
                        .map(|coedge| CoedgeData {
                            edge: coedge.edge,
                            reversed: coedge.reversed,
                            pcurve: coedge.pcurve.as_ref(),
                        })
                        .collect()
                })
                .collect();
            let data = FaceData {
                surface: &face.surface,
                same_sense: face.same_sense,
                loops,
            };
            (id, data)
        })
    }
}

#[cfg(test)]
mod tests {
    use core::f64::consts::FRAC_PI_2;
    use core::f64::consts::PI;

    use golf_geom::Circle;
    use golf_geom::Line;
    use golf_geom::Placement;
    use golf_geom::Plane;
    use golf_geom::Sphere;
    use golf_geom::Torus;
    use golf_manifold::Point;
    use golf_manifold::Vector;
    use golf_manifold::World;
    use golf_mesh::Mesh;
    use golf_mesh::MeshError;
    use golf_mesh::Tolerance;
    use golf_mesh::mesh;
    use nalgebra::Vector2;
    use nalgebra::Vector3;

    use crate::Body;
    use crate::Coedge;
    use crate::FaceId;
    use crate::Loop;
    use crate::cuboid;
    use crate::cylinder;

    fn tolerance() -> Tolerance {
        Tolerance::new(1e-3, 0.3)
    }

    fn tilted() -> Placement<World> {
        Placement::from_axes(
            Point::new(Vector3::new(1.0, -2.0, 0.5)),
            Vector::new(Vector3::new(0.3, -0.2, 1.0)),
            Vector::new(Vector3::new(1.0, 1.0, 0.0)),
        )
    }

    /// Watertight, outward-facing, and enclosing `volume` to within `relative`.
    fn assert_solid(mesh: &Mesh<World, FaceId>, volume: f64, relative: f64) {
        let open = mesh.open_edges();
        assert!(
            open.is_empty(),
            "{} open edges, e.g. {:?}",
            open.len(),
            &open[..open.len().min(5)]
        );
        let v = mesh.signed_volume();
        assert!(
            ((v - volume) / volume).abs() < relative,
            "volume {v}, expected {volume}"
        );
        for t in &mesh.triangles {
            let [a, b, c] = t.vertices.map(|i| mesh.positions[i as usize].coords);
            let facet = (b - a).cross(&(c - a));
            for n in t.normals {
                assert!(
                    n.coords.dot(&facet) > -1e-12,
                    "normal disagrees with winding"
                );
            }
        }
    }

    #[test]
    fn cuboid_meshes_exactly() {
        let mesh = mesh(&cuboid(tilted(), Vector3::new(1.0, 2.0, 3.0)), &tolerance()).unwrap();
        assert_solid(&mesh, 6.0, 1e-12);
        assert!((mesh.area() - 22.0).abs() < 1e-12);
        assert_eq!(mesh.triangles.len(), 12);
    }

    #[test]
    fn cylinder_meshes_watertight() {
        let (radius, height) = (1.5, 4.0);
        let mesh = mesh(&cylinder(tilted(), radius, height), &tolerance()).unwrap();
        assert_solid(&mesh, PI * radius * radius * height, 2e-3);
        // Every vertex is within the chord tolerance of the true cylinder.
        for p in &mesh.positions {
            let local = tilted().to_local(*p);
            let off = (local.xy().norm() - radius)
                .abs()
                .min(local.z.abs())
                .min((local.z - height).abs());
            assert!(off < 1e-9, "{local:?}");
        }
    }

    #[test]
    fn max_length_limits_triangle_edges() {
        let tolerance = tolerance().with_max_length(0.4);
        let mesh = mesh(&cuboid(tilted(), Vector3::new(1.0, 2.0, 3.0)), &tolerance).unwrap();
        assert_solid(&mesh, 6.0, 1e-12);
        for t in &mesh.triangles {
            let [a, b, c] = t.vertices.map(|i| mesh.positions[i as usize].coords);
            for (p, q) in [(a, b), (b, c), (c, a)] {
                assert!((p - q).norm() <= 0.4 + 1e-12);
            }
        }
    }

    #[test]
    fn whole_sphere_meshes_over_its_domain() {
        let mut body = Body::<World>::new();
        let face = body
            .add_face(
                Sphere::new(Point::new(Vector3::new(1.0, 2.0, 3.0)), 2.0),
                true,
                vec![],
            )
            .unwrap();
        body.add_shell(vec![face]).unwrap();
        let mesh = mesh(&body, &tolerance()).unwrap();
        assert_solid(&mesh, 4.0 / 3.0 * PI * 8.0, 3e-3);
    }

    #[test]
    fn torus_meshes_with_both_axes_periodic() {
        let mut body = Body::<World>::new();
        let face = body
            .add_face(Torus::new(tilted(), 3.0, 1.0), true, vec![])
            .unwrap();
        body.add_shell(vec![face]).unwrap();
        let mesh = mesh(&body, &tolerance()).unwrap();
        assert_solid(&mesh, 2.0 * PI * PI * 3.0, 3e-3);
    }

    /// The upper half of a sphere, closed by a disc. The curved face's loop runs
    /// up the seam to the pole and back down, so it crosses the pole's singular
    /// bound with no edge there.
    #[test]
    fn hemisphere_bridges_its_pole() {
        let radius = 2.0;
        let origin = Point::<World, 3>::new(Vector3::zeros());
        let mut body = Body::new();
        let equator_vertex = body
            .add_vertex(Point::new(Vector3::new(-radius, 0.0, 0.0)))
            .unwrap();
        let pole_vertex = body
            .add_vertex(Point::new(Vector3::new(0.0, 0.0, radius)))
            .unwrap();
        let equator = body
            .add_edge(
                Circle::new(Placement::at(origin), radius),
                (-PI, PI),
                equator_vertex,
                equator_vertex,
            )
            .unwrap();
        // The meridian at u = ±π, as a quarter circle from the equator up.
        let meridian = Circle::new(
            Placement::from_axes(
                origin,
                Vector::new(Vector3::y()),
                Vector::new(-Vector3::x()),
            ),
            radius,
        );
        let seam = body
            .add_edge(meridian, (0.0, FRAC_PI_2), equator_vertex, pole_vertex)
            .unwrap();
        let uv_line = |origin: [f64; 2], direction: [f64; 2]| {
            Line::new(
                Point::new(Vector2::from(origin)),
                Vector::new(Vector2::from(direction)),
            )
        };
        let dome = Loop::new(vec![
            Coedge::new(equator, false).with_pcurve(uv_line([0.0, 0.0], [1.0, 0.0])),
            Coedge::new(seam, false).with_pcurve(uv_line([PI, 0.0], [0.0, 1.0])),
            Coedge::new(seam, true).with_pcurve(uv_line([-PI, 0.0], [0.0, 1.0])),
        ]);
        let dome = body
            .add_face(Sphere::new(origin, radius), true, vec![dome])
            .unwrap();
        let disc = Loop::new(vec![Coedge::new(equator, true)]);
        let disc = body
            .add_face(Plane::new(Placement::at(origin)), false, vec![disc])
            .unwrap();
        body.add_shell(vec![dome, disc]).unwrap();
        assert_eq!(body.validate(1e-9), Ok(()));

        let mesh = mesh(&body, &tolerance()).unwrap();
        assert_solid(&mesh, 2.0 / 3.0 * PI * radius.powi(3), 3e-3);
    }

    #[test]
    fn unbounded_face_without_loops_is_rejected() {
        let mut body = Body::<World>::new();
        let face = body.add_face(Plane::new(tilted()), true, vec![]).unwrap();
        body.add_shell(vec![face]).unwrap();
        assert!(matches!(
            mesh(&body, &tolerance()),
            Err(MeshError::UnboundedFace { .. })
        ));
    }

    #[test]
    fn loop_round_a_cylinder_without_a_seam_is_rejected() {
        let mut body = Body::<World>::new();
        let origin = Point::new(Vector3::zeros());
        let start = body
            .add_vertex(Point::new(Vector3::new(-1.0, 0.0, 0.0)))
            .unwrap();
        let circle = body
            .add_edge(
                Circle::new(Placement::at(origin), 1.0),
                (-PI, PI),
                start,
                start,
            )
            .unwrap();
        let side = Loop::new(vec![Coedge::new(circle, false)]);
        body.add_face(
            golf_geom::Cylinder::new(Placement::at(origin), 1.0),
            true,
            vec![side],
        )
        .unwrap();
        assert!(matches!(
            mesh(&body, &tolerance()),
            Err(MeshError::LoopWrapsSurface { .. })
        ));
    }

    #[test]
    fn meshes_in_a_frame_keep_their_tag() {
        use golf_frame::Frame;
        use golf_frame::FrameTree;
        let mut tree = FrameTree::new();
        let frame = tree.add(tree.root(), nalgebra::Isometry3::identity());
        let body = cuboid(
            Placement::<Frame>::at(frame.point(Vector3::zeros())),
            Vector3::new(1.0, 1.0, 1.0),
        );
        let mesh = mesh(&body, &tolerance()).unwrap();
        assert!(mesh.is_watertight());
        assert!(mesh.positions.iter().all(|p| p.tag() == frame));
        assert!(
            mesh.triangles
                .iter()
                .all(|t| t.normals.iter().all(|n| n.tag() == frame))
        );
    }
}
