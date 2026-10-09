use golf_brep::Body;
use golf_brep::Coedge;
use golf_brep::EdgeId;
use golf_brep::FaceUv;
use golf_brep::Loop;
use golf_brep::VertexId;
use golf_geom::AnyCurve;
use golf_geom::AnySurface;
use golf_geom::Circle;
use golf_geom::Cylinder;
use golf_geom::Line;
use golf_geom::Placement;
use golf_geom::Plane;
use golf_manifold::Point;
use golf_manifold::Space;
use golf_manifold::Vector;
use golf_sketch::Region;
use golf_sketch::Segment;
use nalgebra::Vector2;
use nalgebra::Vector3;

use crate::error::ModelError;

/// Extrudes `region`, drawn in the sketch plane `plane`, along the plane's
/// normal by `distance` (backwards if negative) into a solid.
///
/// The solid has a face per segment of each profile (planar for lines,
/// cylindrical for arcs) and two planar caps; holes become through-holes.
pub fn extrude<S: Space<3>>(
    plane: &Placement<S>,
    region: &Region,
    distance: f64,
) -> Result<Body<S>, ModelError> {
    if !distance.is_finite() || distance == 0.0 {
        return Err(ModelError::ZeroDistance);
    }
    // Extrude upwards from whichever end is lower.
    let base = match distance > 0.0 {
        true => *plane,
        false => shifted(plane, distance),
    };
    let height = distance.abs();
    let top = shifted(&base, height);
    let up = base.vector(Vector3::z());

    let mut body = Body::with_tag(base.origin.tag());
    let mut bottom_loops = Vec::new();
    let mut top_loops = Vec::new();
    let mut sides = Vec::new();
    for profile in region.profiles() {
        let segments = profile.segments();
        let n = segments.len();
        let at = |placement: &Placement<S>, p: Vector2<f64>| {
            placement.point(Vector3::new(p.x, p.y, 0.0))
        };
        let bottom_vertices: Vec<VertexId> = segments
            .iter()
            .map(|s| body.add_vertex(at(&base, s.start())))
            .collect::<Result<_, _>>()?;
        let top_vertices: Vec<VertexId> = segments
            .iter()
            .map(|s| body.add_vertex(at(&top, s.start())))
            .collect::<Result<_, _>>()?;
        // Vertical edges up from each vertex, parametrised by height.
        let verticals: Vec<EdgeId> = (0..n)
            .map(|k| {
                let line = Line::new(at(&base, segments[k].start()), up);
                body.add_edge(line, (0.0, height), bottom_vertices[k], top_vertices[k])
            })
            .collect::<Result<_, _>>()?;

        let mut bottom_coedges = Vec::new();
        let mut top_coedges = Vec::new();
        for (k, segment) in segments.iter().enumerate() {
            let (a, b) = (k, (k + 1) % n);
            let swept = SweptSegment::new(segment);
            let edge = |body: &mut Body<S>, placement: &Placement<S>, vertices: &[VertexId]| {
                let (start, end) = match swept.forward {
                    true => (vertices[a], vertices[b]),
                    false => (vertices[b], vertices[a]),
                };
                body.add_edge(swept.curve(placement), swept.range(), start, end)
            };
            let bottom_edge = edge(&mut body, &base, &bottom_vertices)?;
            let top_edge = edge(&mut body, &top, &top_vertices)?;

            // The side face: in its (u, v) the bottom edge runs along v = 0 with
            // u the edge's own parameter, the top along v = height, and the
            // verticals up at u = t_to (B) and down at u = t_from (A).
            let side = Loop::new(vec![
                Coedge::new(bottom_edge, !swept.forward)
                    .with_pcurve(uv_line([0.0, 0.0], [1.0, 0.0])),
                Coedge::new(verticals[b], false)
                    .with_pcurve(uv_line([swept.t_to, 0.0], [0.0, 1.0])),
                Coedge::new(top_edge, swept.forward)
                    .with_pcurve(uv_line([0.0, height], [1.0, 0.0])),
                Coedge::new(verticals[a], true)
                    .with_pcurve(uv_line([swept.t_from, 0.0], [0.0, 1.0])),
            ]);
            let (surface, same_sense) = swept.side(&base);
            sides.push(body.add_face(surface, same_sense, vec![side])?);

            // The caps' uv are sketch coordinates, so a line's pcurve is the
            // line itself.
            let cap_pcurve = swept.cap_pcurve();
            let cap = |edge, reversed| {
                let coedge = Coedge::new(edge, reversed);
                match &cap_pcurve {
                    Some(line) => coedge.with_pcurve(line.clone()),
                    None => coedge,
                }
            };
            top_coedges.push(cap(top_edge, !swept.forward));
            bottom_coedges.push(cap(bottom_edge, swept.forward));
        }
        bottom_coedges.reverse();
        bottom_loops.push(Loop::new(bottom_coedges));
        top_loops.push(Loop::new(top_coedges));
    }

    // The bottom plane's normal points up, into the solid, so the face flips it.
    let bottom = body.add_face(Plane::new(base), false, bottom_loops)?;
    let top_face = body.add_face(Plane::new(top), true, top_loops)?;
    let mut faces = vec![bottom, top_face];
    faces.extend(sides);
    body.add_shell(faces)?;
    Ok(body)
}

/// `placement` moved `distance` along its own z axis.
fn shifted<S: Space<3>>(placement: &Placement<S>, distance: f64) -> Placement<S> {
    Placement::new(
        placement.point(Vector3::new(0.0, 0.0, distance)),
        placement.rotation,
    )
}

pub(crate) fn uv_line<S: Space<3>>(origin: [f64; 2], direction: [f64; 2]) -> Line<FaceUv<S>, 2> {
    Line::new(
        Point::new(Vector2::from(origin)),
        Vector::new(Vector2::from(direction)),
    )
}

/// A profile segment as an edge curve: how it's parametrised, and which way the
/// profile runs along it.
struct SweptSegment<'a> {
    segment: &'a Segment,
    /// The curve parameters where the profile enters and leaves the segment.
    t_from: f64,
    t_to: f64,
    /// Whether the profile runs along increasing parameter.
    forward: bool,
}

impl<'a> SweptSegment<'a> {
    /// Lines are parametrised by length from their start, arcs by angle (about
    /// the sketch's x axis), so arcs run backwards when clockwise.
    fn new(segment: &'a Segment) -> Self {
        match *segment {
            Segment::Line { .. } => Self {
                segment,
                t_from: 0.0,
                t_to: segment.length(),
                forward: true,
            },
            Segment::Arc {
                start_angle, sweep, ..
            } => Self {
                segment,
                t_from: start_angle,
                t_to: start_angle + sweep,
                forward: sweep > 0.0,
            },
        }
    }

    fn range(&self) -> (f64, f64) {
        (self.t_from.min(self.t_to), self.t_from.max(self.t_to))
    }

    /// The segment drawn in `placement`'s xy plane.
    fn curve<S: Space<3>>(&self, placement: &Placement<S>) -> AnyCurve<S> {
        match *self.segment {
            Segment::Line { start, end } => {
                let direction = (end - start).normalize();
                Line::new(
                    placement.point(Vector3::new(start.x, start.y, 0.0)),
                    placement.vector(Vector3::new(direction.x, direction.y, 0.0)),
                )
                .into()
            }
            Segment::Arc { center, radius, .. } => {
                Circle::new(about(placement, center), radius).into()
            }
        }
    }

    /// The surface swept by the segment, and whether its normal already points
    /// out of the solid (to the right of the profile's direction).
    fn side<S: Space<3>>(&self, base: &Placement<S>) -> (AnySurface<S>, bool) {
        match *self.segment {
            Segment::Line { start, end } => {
                // x along the line and z out to its right make y point up, and
                // u the line's own parameter.
                let along = (end - start).normalize();
                let x = base.vector(Vector3::new(along.x, along.y, 0.0));
                let outward = base.vector(Vector3::new(along.y, -along.x, 0.0));
                let origin = base.point(Vector3::new(start.x, start.y, 0.0));
                (
                    Plane::new(Placement::from_axes(origin, outward, x)).into(),
                    true,
                )
            }
            // An anticlockwise arc has the solid towards its centre, so the
            // cylinder's outward normal points out of it.
            Segment::Arc { center, radius, .. } => (
                Cylinder::new(about(base, center), radius).into(),
                self.forward,
            ),
        }
    }

    /// A line segment's pcurve on either cap; arcs have none, there being no 2D
    /// circle pcurve yet.
    fn cap_pcurve<S: Space<3>>(&self) -> Option<Line<FaceUv<S>, 2>> {
        match *self.segment {
            Segment::Line { start, end } => {
                let direction = (end - start).normalize();
                Some(uv_line([start.x, start.y], [direction.x, direction.y]))
            }
            Segment::Arc { .. } => None,
        }
    }
}

/// `placement`'s frame moved to `center` in its xy plane.
fn about<S: Space<3>>(placement: &Placement<S>, center: Vector2<f64>) -> Placement<S> {
    Placement::new(
        placement.point(Vector3::new(center.x, center.y, 0.0)),
        placement.rotation,
    )
}
