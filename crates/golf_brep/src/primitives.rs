//! Bodies for simple solids.

use core::f64::consts::PI;
use std::collections::HashMap;

use golf_geom::Circle;
use golf_geom::Cylinder;
use golf_geom::Line;
use golf_geom::Placement;
use golf_geom::Plane;
use golf_manifold::Embedding;
use golf_manifold::Point;
use golf_manifold::Space;
use golf_manifold::Vector;
use nalgebra::Vector2;
use nalgebra::Vector3;

use super::Body;
use super::Coedge;
use super::EdgeId;
use super::Loop;

/// A box with one corner at the placement's origin, extending `size` along its
/// local axes. Faces carry line pcurves.
pub fn cuboid<S: Space<3>>(placement: Placement<S>, size: Vector3<f64>) -> Body<S> {
    let mut body = Body::with_tag(placement.origin.tag());
    // Vertex i is at the corner (x, y, z) = bits (0, 1, 2) of i.
    let corner = |i: usize| {
        let bit = |b: usize| ((i >> b) & 1) as f64;
        Vector3::new(bit(0) * size.x, bit(1) * size.y, bit(2) * size.z)
    };
    let vertices: Vec<_> = (0..8)
        .map(|i| {
            body.add_vertex(placement.point(corner(i)))
                .expect("same space")
        })
        .collect();

    // Each face's corners, anticlockwise seen from outside.
    const FACES: [[usize; 4]; 6] = [
        [0, 2, 3, 1], // z = 0
        [4, 5, 7, 6], // z = 1
        [0, 1, 5, 4], // y = 0
        [2, 6, 7, 3], // y = 1
        [0, 4, 6, 2], // x = 0
        [1, 3, 7, 5], // x = 1
    ];
    // Edges run from the lower-numbered corner to the higher, over [0, 1].
    let mut edges: HashMap<(usize, usize), EdgeId> = HashMap::new();
    let mut faces = Vec::new();
    for corners in FACES {
        let local = corners.map(corner);
        let normal = (local[1] - local[0]).cross(&(local[2] - local[1]));
        let plane = Plane::new(Placement::from_axes(
            placement.point(local[0]),
            placement.vector(normal),
            placement.vector(local[1] - local[0]),
        ));
        let uv = |i: usize| {
            plane
                .project(placement.point(corner(i)), None)
                .expect("planes always project")
                .coords
        };
        let coedges = (0..4)
            .map(|k| {
                let (a, b) = (corners[k], corners[(k + 1) % 4]);
                let (lo, hi) = (a.min(b), a.max(b));
                let edge = *edges.entry((lo, hi)).or_insert_with(|| {
                    let line = Line::new(
                        placement.point(corner(lo)),
                        placement.vector(corner(hi) - corner(lo)),
                    );
                    body.add_edge(line, (0.0, 1.0), vertices[lo], vertices[hi])
                        .expect("vertices exist")
                });
                let pcurve = Line::new(Point::new(uv(lo)), Vector::new(uv(hi) - uv(lo)));
                Coedge::new(edge, a > b).with_pcurve(pcurve)
            })
            .collect();
        faces.push(
            body.add_face(plane, true, vec![Loop::new(coedges)])
                .expect("loop is closed"),
        );
    }
    body.add_shell(faces).expect("faces exist");
    body
}

/// A solid cylinder about the placement's local z axis, from its origin up to
/// `height`.
///
/// Its circular edges start and end on a seam line at local `(-radius, 0)`,
/// which the side face uses twice. The side face carries pcurves; the end
/// discs don't yet, there being no 2D circle to describe them with.
pub fn cylinder<S: Space<3>>(placement: Placement<S>, radius: f64, height: f64) -> Body<S> {
    let mut body = Body::with_tag(placement.origin.tag());
    let top_placement = Placement::new(
        placement.point(Vector3::new(0.0, 0.0, height)),
        placement.rotation,
    );
    let seam_at = |z: f64| placement.point(Vector3::new(-radius, 0.0, z));
    let bottom_vertex = body.add_vertex(seam_at(0.0)).expect("same space");
    let top_vertex = body.add_vertex(seam_at(height)).expect("same space");

    // Circles are parametrised by angle, so ±π is the seam.
    let full_turn = (-PI, PI);
    let bottom = body
        .add_edge(
            Circle::new(placement, radius),
            full_turn,
            bottom_vertex,
            bottom_vertex,
        )
        .expect("vertices exist");
    let top = body
        .add_edge(
            Circle::new(top_placement, radius),
            full_turn,
            top_vertex,
            top_vertex,
        )
        .expect("vertices exist");
    let seam_line = Line::new(seam_at(0.0), placement.vector(Vector3::z()));
    let seam = body
        .add_edge(seam_line, (0.0, height), bottom_vertex, top_vertex)
        .expect("vertices exist");

    // In the side's (angle, height) space the loop runs anticlockwise round the
    // rectangle [-π, π] × [0, height], crossing the seam at u = π going up and at
    // u = -π coming down.
    let uv_line = |origin: [f64; 2], direction: [f64; 2]| {
        Line::new(
            Point::new(Vector2::from(origin)),
            Vector::new(Vector2::from(direction)),
        )
    };
    let side = Loop::new(vec![
        Coedge::new(bottom, false).with_pcurve(uv_line([0.0, 0.0], [1.0, 0.0])),
        Coedge::new(seam, false).with_pcurve(uv_line([PI, 0.0], [0.0, 1.0])),
        Coedge::new(top, true).with_pcurve(uv_line([0.0, height], [1.0, 0.0])),
        Coedge::new(seam, true).with_pcurve(uv_line([-PI, 0.0], [0.0, 1.0])),
    ]);
    let side = body
        .add_face(Cylinder::new(placement, radius), true, vec![side])
        .expect("loop is closed");
    // The bottom plane's normal points up, into the solid, so the face flips it.
    let bottom_face = body
        .add_face(
            Plane::new(placement),
            false,
            vec![Loop::new(vec![Coedge::new(bottom, true)])],
        )
        .expect("loop is closed");
    let top_face = body
        .add_face(
            Plane::new(top_placement),
            true,
            vec![Loop::new(vec![Coedge::new(top, false)])],
        )
        .expect("loop is closed");
    body.add_shell(vec![side, bottom_face, top_face])
        .expect("faces exist");
    body
}

#[cfg(test)]
mod tests {
    use golf_frame::Frame;
    use golf_frame::FrameTree;
    use golf_manifold::Mapping;
    use golf_manifold::Surface;
    use golf_manifold::World;
    use nalgebra::Isometry3;
    use nalgebra::Translation3;
    use nalgebra::UnitQuaternion;

    use super::*;
    use crate::FaceId;
    use crate::TopologyError;

    fn tilted() -> Placement<World> {
        Placement::from_axes(
            Point::new(Vector3::new(1.0, -2.0, 0.5)),
            Vector::new(Vector3::new(0.3, -0.2, 1.0)),
            Vector::new(Vector3::new(1.0, 1.0, 0.0)),
        )
    }

    /// Each face's normal, at the middle of its surface patch, points away from
    /// the body's centroid of vertices.
    fn assert_faces_point_out<S: Space<3>>(body: &Body<S>, inside: Point<S, 3>) {
        for (id, face) in body.faces() {
            let edge = body.edge(face.loops[0].coedges[0].edge);
            let t = Point::new([(edge.range.0 + edge.range.1) / 2.0].into());
            let uv = face.surface.project(edge.curve.apply(t), None).unwrap();
            let normal = face.surface.normal(uv).coords * if face.same_sense { 1.0 } else { -1.0 };
            let outward = (face.surface.apply(uv) - inside).coords;
            assert!(normal.dot(&outward) > 0.0, "{id} points inwards");
        }
    }

    #[test]
    fn cuboid_is_a_valid_closed_box() {
        let size = Vector3::new(1.0, 2.0, 3.0);
        let body = cuboid(tilted(), size);
        assert_eq!(
            (
                body.vertices().len(),
                body.edges().len(),
                body.faces().len()
            ),
            (8, 12, 6)
        );
        assert_eq!(body.validate(1e-9), Ok(()));
        assert_eq!(body.genus(), Some(0));
        assert_faces_point_out(&body, tilted().point(size / 2.0));
        for (_, uses) in body.edge_uses() {
            assert_ne!(uses[0].0, uses[1].0, "box edges join two different faces");
        }
    }

    #[test]
    fn cylinder_is_a_valid_closed_solid() {
        let body = cylinder(tilted(), 1.5, 4.0);
        assert_eq!(
            (
                body.vertices().len(),
                body.edges().len(),
                body.faces().len()
            ),
            (2, 3, 3)
        );
        assert_eq!(body.validate(1e-9), Ok(()));
        assert_eq!(body.genus(), Some(0));
        assert_faces_point_out(&body, tilted().point(Vector3::new(0.0, 0.0, 2.0)));
        // The seam is used twice by the side face.
        let seam_uses = &body.edge_uses()[&EdgeId(2)];
        assert_eq!(seam_uses.len(), 2);
        assert_eq!(seam_uses[0].0, seam_uses[1].0);
    }

    #[test]
    fn bodies_in_a_frame_validate() {
        let mut tree = FrameTree::new();
        let frame = tree.add(
            tree.root(),
            Isometry3::from_parts(
                Translation3::new(5.0, 0.0, 0.0),
                UnitQuaternion::from_scaled_axis(Vector3::y()),
            ),
        );
        let placement: Placement<Frame> = Placement::at(frame.point(Vector3::zeros()));
        assert_eq!(
            cuboid(placement, Vector3::new(1.0, 1.0, 1.0)).validate(1e-9),
            Ok(())
        );
        assert_eq!(cylinder(placement, 1.0, 1.0).validate(1e-9), Ok(()));
        let body = cuboid(placement, Vector3::new(1.0, 1.0, 1.0));
        assert_eq!(body.tag(), frame);
        let mut wrong = Body::<Frame>::with_tag(tree.root());
        assert_eq!(
            wrong.add_vertex(frame.point(Vector3::zeros())),
            Err(TopologyError::WrongSpace)
        );
    }

    #[test]
    fn missing_face_leaves_edges_half_used() {
        let mut body = cuboid(tilted(), Vector3::new(1.0, 1.0, 1.0));
        body.faces.pop();
        body.shells[0].faces.pop();
        let errors = body.validate(1e-9).unwrap_err();
        assert_eq!(errors.len(), 4);
        assert!(errors.iter().all(|e| matches!(
            e,
            TopologyError::EdgeUses {
                forward: 1,
                reversed: 0,
                ..
            } | TopologyError::EdgeUses {
                forward: 0,
                reversed: 1,
                ..
            }
        )));
    }

    #[test]
    fn flipped_face_uses_edges_the_same_way_twice() {
        let mut body = cuboid(tilted(), Vector3::new(1.0, 1.0, 1.0));
        let face = &mut body.faces[2];
        face.same_sense = false;
        face.loops[0].coedges.reverse();
        for coedge in &mut face.loops[0].coedges {
            coedge.reversed = !coedge.reversed;
        }
        let errors = body.validate(1e-9).unwrap_err();
        assert_eq!(errors.len(), 4, "{errors:?}");
    }

    #[test]
    fn moved_vertex_and_stray_edge_are_reported() {
        let mut body = cuboid(tilted(), Vector3::new(1.0, 1.0, 1.0));
        body.vertices[0].point = body.vertices[0].point + Vector::new(Vector3::new(0.0, 0.0, 0.1));
        let errors = body.validate(1e-9).unwrap_err();
        // Three edges meet at the corner.
        let off = errors
            .iter()
            .filter(|e| matches!(e, TopologyError::VertexOffCurve { .. }))
            .count();
        assert_eq!(off, 3, "{errors:?}");

        let mut body = cylinder(tilted(), 1.0, 2.0);
        body.faces[1].surface = Plane::new(Placement::new(
            tilted().point(Vector3::new(0.0, 0.0, 0.5)),
            tilted().rotation,
        ))
        .into();
        let errors = body.validate(1e-9).unwrap_err();
        assert!(
            matches!(
                errors[..],
                [TopologyError::EdgeOffFace {
                    face: FaceId(1),
                    ..
                }]
            ),
            "{errors:?}"
        );
    }

    #[test]
    fn face_outside_any_shell_is_reported() {
        let mut body = cuboid(tilted(), Vector3::new(1.0, 1.0, 1.0));
        body.shells[0].faces.retain(|&f| f != FaceId(3));
        assert_eq!(
            body.validate(1e-9),
            Err(vec![TopologyError::FaceShells {
                face: FaceId(3),
                shells: 0
            }])
        );
    }
}
