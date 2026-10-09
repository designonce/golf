use core::f64::consts::FRAC_PI_2;
use core::f64::consts::TAU;

use golf_brep::Body;
use golf_brep::Coedge;
use golf_brep::EdgeId;
use golf_brep::FaceUv;
use golf_brep::Loop;
use golf_brep::VertexId;
use golf_geom::AnyCurve;
use golf_geom::Circle;
use golf_geom::Cylinder;
use golf_geom::Ellipse;
use golf_geom::Line;
use golf_geom::Placement;
use golf_geom::Plane;
use golf_geom::Torus;
use golf_manifold::Point;
use golf_manifold::Space;
use golf_manifold::Vector;
use golf_sketch::Path;
use golf_sketch::Segment;
use nalgebra::Vector2;
use nalgebra::Vector3;

use crate::error::ModelError;
use crate::extrude::uv_line;

/// Sweeps a circle of `radius` along `path`, drawn in the sketch plane
/// `plane`, into a solid pipe.
///
/// Lines sweep cylinders and arcs sweep tori, so an arc's radius must exceed
/// the pipe's. Corners must be smooth, or sharp between two lines, where the
/// cylinders meet in a mitre (an ellipse); each line must be long enough for
/// the mitres at its ends. An open path is capped with discs at both ends; a
/// closed one makes a ring.
pub fn pipe<S: Space<3>>(
    plane: &Placement<S>,
    path: &Path,
    radius: f64,
) -> Result<Body<S>, ModelError> {
    if !radius.is_finite() || radius <= 0.0 {
        return Err(ModelError::BadTreatment);
    }
    let segments = path.segments();
    let n = segments.len();
    let closed = path.is_closed();
    let joint_count = if closed { n } else { n + 1 };

    // Each joint's point, kind, and the directions the path comes in and goes
    // out (the same at a smooth joint or an open end).
    let joints: Vec<PathJoint> = (0..joint_count)
        .map(|j| {
            let incoming = (closed || j > 0).then(|| tangent(&segments[(j + n - 1) % n], 1.0));
            let outgoing = (closed || j < n).then(|| tangent(&segments[j % n], 0.0));
            let point = if j < n {
                segments[j].start()
            } else {
                segments[n - 1].end()
            };
            let (inward, outward) = (
                incoming.or(outgoing).expect("a joint has a segment"),
                outgoing.or(incoming).expect("a joint has a segment"),
            );
            let smooth = inward.perp(&outward).abs() <= 1e-9 && inward.dot(&outward) > 0.0;
            let lines = |k: usize| matches!(segments[k % n], Segment::Line { .. });
            if smooth {
                Ok(PathJoint {
                    point,
                    incoming: inward,
                    outgoing: outward,
                    mitre: false,
                })
            } else if lines(j + n - 1) && lines(j) && inward.dot(&outward) > -1.0 + 1e-9 {
                Ok(PathJoint {
                    point,
                    incoming: inward,
                    outgoing: outward,
                    mitre: true,
                })
            } else {
                Err(ModelError::UnsupportedCorner { index: j % n })
            }
        })
        .collect::<Result<_, _>>()?;

    // Arcs must be wider than the pipe, lines longer than their mitres reach.
    let reach = |j: &PathJoint| match j.mitre {
        true => radius * (j.incoming.angle(&j.outgoing) / 2.0).tan(),
        false => 0.0,
    };
    for (index, segment) in segments.iter().enumerate() {
        let too_thick = match *segment {
            Segment::Arc { radius: arc, .. } => arc <= radius * (1.0 + 1e-9),
            Segment::Line { .. } => {
                segment.length() < reach(&joints[index]) + reach(&joints[(index + 1) % joint_count])
            }
        };
        if too_thick {
            return Err(ModelError::PipeTooThick { index });
        }
    }

    let up = plane.vector(Vector3::z()).coords;
    let tag = plane.origin.tag();
    let point = |p: Vector2<f64>| plane.point(Vector3::new(p.x, p.y, 0.0));
    let vector = |v: Vector2<f64>| plane.vector(Vector3::new(v.x, v.y, 0.0)).coords;
    let with_tag = |v: Vector3<f64>| Vector::with_tag(v, tag);
    let mut body = Body::with_tag(tag);

    // The seam runs along the top of the pipe: each joint's vertex is there,
    // and its cross-section starts and ends there, its parameter turning from
    // the top towards (tangent × up).
    let vertices: Vec<VertexId> = joints
        .iter()
        .map(|j| body.add_vertex(Point::with_tag(point(j.point).coords + up * radius, tag)))
        .collect::<Result<_, _>>()?;
    let sections: Vec<EdgeId> = joints
        .iter()
        .zip(&vertices)
        .map(|(j, &v)| {
            let (curve, range): (AnyCurve<S>, _) = match j.mitre {
                false => {
                    let frame = Placement::from_axes(
                        point(j.point),
                        with_tag(vector(j.outgoing)),
                        with_tag(up),
                    );
                    (Circle::new(frame, radius).into(), (0.0, TAU))
                }
                // The cylinders meet in the plane bisecting their axes: an
                // ellipse, wider across the bend. Its frame is set so its
                // parameter turns the same way as the cylinders' angle, from
                // the top at 3π/2.
                true => {
                    let bisector = (vector(j.incoming) + vector(j.outgoing)).normalize();
                    let across = vector(j.incoming).cross(&up);
                    let across = (across - bisector * across.dot(&bisector)).normalize();
                    let half_turn = j.incoming.angle(&j.outgoing) / 2.0;
                    let frame = Placement::from_axes(
                        point(j.point),
                        with_tag(across.cross(&-up)),
                        with_tag(across),
                    );
                    let ellipse = Ellipse::new(frame, radius / half_turn.cos(), radius);
                    (ellipse.into(), (3.0 * FRAC_PI_2, 3.0 * FRAC_PI_2 + TAU))
                }
            };
            body.add_edge(curve, range, v, v)
        })
        .collect::<Result<_, _>>()?;

    let mut faces = Vec::new();
    for (i, segment) in segments.iter().enumerate() {
        let (a, b) = (i, (i + 1) % joint_count);
        let top_a = point(joints[a].point).coords + up * radius;
        // The face's parametrisation: where the seam is on its far and near
        // side, and the cross-sections' pcurves (none for a mitre).
        let (surface, seam, seam_range, forward, pcurves): (
            golf_geom::AnySurface<S>,
            AnyCurve<S>,
            _,
            bool,
            Pcurves<S>,
        ) = match *segment {
            Segment::Line { start, end } => {
                let length = (end - start).norm();
                let along = vector((end - start) / length);
                let frame = Placement::from_axes(point(start), with_tag(along), with_tag(up));
                let section = |v: f64, mitre: bool| (!mitre).then(|| uv_line([0.0, v], [1.0, 0.0]));
                (
                    Cylinder::new(frame, radius).into(),
                    Line::new(Point::with_tag(top_a, tag), with_tag(along)).into(),
                    (0.0, length),
                    true,
                    Pcurves {
                        a: section(0.0, joints[a].mitre),
                        b: section(length, joints[b].mitre),
                        far: uv_line([TAU, 0.0], [0.0, 1.0]),
                        near: uv_line([0.0, 0.0], [0.0, 1.0]),
                    },
                )
            }
            Segment::Arc {
                center,
                radius: major,
                start_angle,
                sweep,
            } => {
                // u is the angle round the arc; the tube angle v is π/2 at the
                // top and runs against the section parameter for an
                // anticlockwise arc (tangent × up points out of the arc), with
                // it for a clockwise one.
                let rate = -sweep.signum();
                let frame = Placement::new(point(center), plane.rotation);
                let seam_frame = Placement::new(
                    Point::with_tag(point(center).coords + up * radius, tag),
                    plane.rotation,
                );
                let (from, to) = (start_angle, start_angle + sweep);
                (
                    Torus::new(frame, major, radius).into(),
                    Circle::new(seam_frame, major).into(),
                    (from.min(to), from.max(to)),
                    sweep > 0.0,
                    Pcurves {
                        a: Some(uv_line([from, FRAC_PI_2], [0.0, rate])),
                        b: Some(uv_line([to, FRAC_PI_2], [0.0, rate])),
                        far: uv_line([0.0, FRAC_PI_2 + rate * TAU], [1.0, 0.0]),
                        near: uv_line([0.0, FRAC_PI_2], [1.0, 0.0]),
                    },
                )
            }
        };
        let (start, end) = if forward {
            (vertices[a], vertices[b])
        } else {
            (vertices[b], vertices[a])
        };
        let seam = body.add_edge(seam, seam_range, start, end)?;
        // Round section A, along the seam on its far side, back round section
        // B, and back along the seam on its near side.
        let with = |coedge: Coedge<S>, pcurve: Option<Line<FaceUv<S>, 2>>| match pcurve {
            Some(p) => coedge.with_pcurve(p),
            None => coedge,
        };
        let coedges = vec![
            with(Coedge::new(sections[a], false), pcurves.a),
            Coedge::new(seam, !forward).with_pcurve(pcurves.far),
            with(Coedge::new(sections[b], true), pcurves.b),
            Coedge::new(seam, forward).with_pcurve(pcurves.near),
        ];
        faces.push(body.add_face(surface, true, vec![Loop::new(coedges)])?);
    }

    if !closed {
        // Discs at the ends, facing back along the path at the start and on
        // along it at the end.
        let disc = |j: &PathJoint, direction: Vector2<f64>| {
            Plane::new(Placement::from_axes(
                point(j.point),
                with_tag(vector(direction)),
                with_tag(up),
            ))
        };
        let (first, last) = (&joints[0], &joints[joint_count - 1]);
        faces.push(body.add_face(
            disc(first, first.outgoing),
            false,
            vec![Loop::new(vec![Coedge::new(sections[0], true)])],
        )?);
        faces.push(body.add_face(
            disc(last, last.incoming),
            true,
            vec![Loop::new(vec![Coedge::new(
                sections[joint_count - 1],
                false,
            )])],
        )?);
    }
    body.add_shell(faces)?;
    Ok(body)
}

/// A joint of the path: where segments meet, or an open end.
struct PathJoint {
    point: Vector2<f64>,
    incoming: Vector2<f64>,
    outgoing: Vector2<f64>,
    mitre: bool,
}

/// A pipe face's pcurves: its cross-sections at A and B (none for a mitre),
/// and the seam on its far and near sides.
struct Pcurves<S: Space<3>> {
    a: Option<Line<FaceUv<S>, 2>>,
    b: Option<Line<FaceUv<S>, 2>>,
    far: Line<FaceUv<S>, 2>,
    near: Line<FaceUv<S>, 2>,
}

/// A segment's unit tangent at fraction `s` along it.
fn tangent(segment: &Segment, s: f64) -> Vector2<f64> {
    segment.tangent_at(s)
}
