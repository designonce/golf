use core::f64::consts::PI;
use core::f64::consts::TAU;
use std::collections::HashMap;

use golf_brep::Body;
use golf_brep::Coedge;
use golf_brep::EdgeId;
use golf_brep::Loop;
use golf_brep::VertexId;
use golf_geom::AnyCurve;
use golf_geom::AnySurface;
use golf_geom::Circle;
use golf_geom::Cone;
use golf_geom::Cylinder;
use golf_geom::Line;
use golf_geom::Placement;
use golf_geom::Plane;
use golf_geom::Sphere;
use golf_geom::Torus;
use golf_manifold::Space;
use golf_sketch::Region;
use golf_sketch::Segment;
use nalgebra::UnitQuaternion;
use nalgebra::Vector2;
use nalgebra::Vector3;

use crate::error::ModelError;
use crate::extrude::uv_line;

/// Revolves `region`, drawn in the sketch plane `plane`, a full turn about the
/// sketch's y axis into a solid.
///
/// The region must lie on the side `x >= 0`, with sketch x the distance from the
/// axis and sketch y the height along it. Segments lying on the axis bound the
/// solid without making a face; each other segment sweeps a plane, cylinder,
/// cone, sphere or torus, and is itself the seam of that face. Holes revolve
/// into internal voids, each its own shell.
pub fn revolve<S: Space<3>>(plane: &Placement<S>, region: &Region) -> Result<Body<S>, ModelError> {
    revolve_by(plane, region, TAU)
}

/// [`revolve`] through `angle` radians (anticlockwise about the sketch's y
/// axis, seen from its tip), in `(0, 2π]`.
///
/// Short of a full turn the solid is closed by two planar caps: the region at
/// angle 0 and its copy at `angle`. Each revolved face is then bounded by the
/// segment at both angles rather than a seam, and holes become tunnels rather
/// than voids, so everything is one shell.
pub fn revolve_by<S: Space<3>>(
    plane: &Placement<S>,
    region: &Region,
    angle: f64,
) -> Result<Body<S>, ModelError> {
    if !(angle > 0.0 && angle <= TAU * (1.0 + 1e-12)) {
        return Err(ModelError::BadAngle(angle));
    }
    let full = angle >= TAU * (1.0 - 1e-12);
    let angle = if full { TAU } else { angle };
    // The frame the surfaces of revolution are placed in: local x is the
    // sketch's x (the radial direction at angle 0), local z the axis. `end` is
    // the same frame turned through the angle, for the far side of a partial
    // revolution.
    let frame = Placement::from_axes(
        plane.origin,
        plane.vector(Vector3::y()),
        plane.vector(Vector3::x()),
    );
    let axis = nalgebra::Unit::new_normalize(frame.vector(Vector3::z()).coords);
    let turn = UnitQuaternion::from_axis_angle(&axis, angle);
    let end = Placement::new(frame.origin, turn * frame.rotation);
    let end_plane = Placement::new(plane.origin, turn * plane.rotation);
    let on_axis = |x: f64| x.abs() <= 1e-12 * (1.0 + region_size(region));

    let mut body = Body::with_tag(plane.origin.tag());
    let mut all_faces = Vec::new();
    let mut start_caps = Vec::new();
    let mut end_caps = Vec::new();
    for profile in region.profiles() {
        let segments: Vec<Segment> = profile.segments().iter().map(normalise).collect();
        for (index, segment) in segments.iter().enumerate() {
            if closest_to_axis(segment) < -1e-12 * (1.0 + region_size(region)) {
                return Err(ModelError::CrossesAxis { index });
            }
        }
        let n = segments.len();
        let mut corners = Corners {
            frame: &frame,
            end: &end,
            angle,
            full,
            on_axis: &on_axis,
            vertices: HashMap::new(),
            end_vertices: HashMap::new(),
            arcs: HashMap::new(),
        };
        let mut faces = Vec::new();
        let mut start_cap = Vec::new();
        let mut end_cap = Vec::new();

        for (index, segment) in segments.iter().enumerate() {
            let (a, b) = (index, (index + 1) % n);
            let (pa, pb) = (segment.start(), segment.end());
            let along_axis = matches!(segment, Segment::Line { .. }) && on_axis(pa.x) && on_axis(pb.x);
            let perpendicular = matches!(segment, &Segment::Line { start, end }
                if (end.y - start.y).abs() <= 1e-12 * (end - start).norm());
            if full && along_axis {
                continue;
            }
            let arc_a = corners.arc(&mut body, a, pa)?;
            let arc_b = corners.arc(&mut body, b, pb)?;
            // A partial revolution needs the segment as an edge at both angles
            // (one edge if it lies on the axis), whatever its face.
            let partial_edges = match full {
                true => None,
                false => {
                    let (va, vb) = (corners.vertex(&mut body, a, pa)?, corners.vertex(&mut body, b, pb)?);
                    let start_edge = profile_edge(&mut body, segment, &frame, |forward| {
                        if forward { (va, vb) } else { (vb, va) }
                    })?;
                    let end_edge = match along_axis {
                        true => start_edge,
                        false => {
                            let (va_end, vb_end) =
                                (corners.end_vertex(&mut body, a, pa)?, corners.end_vertex(&mut body, b, pb)?);
                            profile_edge(&mut body, segment, &end, |forward| {
                                if forward { (va_end, vb_end) } else { (vb_end, va_end) }
                            })?
                        }
                    };
                    let forward = segment_forward(segment);
                    let cap_pcurve = cap_pcurve::<S>(segment);
                    let cap = |edge, reversed| match &cap_pcurve {
                        Some(line) => Coedge::new(edge, reversed).with_pcurve(line.clone()),
                        None => Coedge::new(edge, reversed),
                    };
                    start_cap.push(cap(start_edge, !forward));
                    end_cap.push(cap(end_edge, forward));
                    Some((start_edge, end_edge, forward))
                }
            };
            if along_axis {
                continue;
            }

            let face = if perpendicular {
                // A disc, annulus or sector perpendicular to the axis. Its
                // plane's normal is the axis; the solid is to the profile's
                // left, so the face points along the axis when the profile runs
                // inwards. Arc A runs forwards and arc B back, as on the curved
                // faces; a partial one is closed by the segment at each angle.
                let surface = Plane::new(Placement::new(
                    frame.point(Vector3::new(0.0, 0.0, pa.y)),
                    frame.rotation,
                ));
                let same_sense = pb.x < pa.x;
                let loops = match partial_edges {
                    None => [arc_a.map(|e| (e, false)), arc_b.map(|e| (e, true))]
                        .into_iter()
                        .flatten()
                        .map(|(e, reversed)| Loop::new(vec![Coedge::new(e, reversed)]))
                        .collect(),
                    Some((start_edge, end_edge, forward)) => {
                        let mut coedges = Vec::new();
                        coedges.extend(arc_a.map(|e| Coedge::new(e, false)));
                        coedges.push(Coedge::new(end_edge, !forward));
                        coedges.extend(arc_b.map(|e| Coedge::new(e, true)));
                        coedges.push(Coedge::new(start_edge, forward));
                        vec![Loop::new(coedges)]
                    }
                };
                body.add_face(surface, same_sense, loops)?
            } else {
                let swept = Revolved::new(index, segment, &frame)?;
                // In the face's (u, v), u the angle of revolution: arc A along
                // v_a, up the segment's edge at u = angle, back along arc B at
                // v_b, down its edge at u = 0. That runs anticlockwise about the
                // surface normal exactly when v increases from A to B, which is
                // when the normal points out of the solid. A full revolution
                // uses one seam edge for both.
                let (start_edge, end_edge) = match partial_edges {
                    Some((start_edge, end_edge, _)) => (start_edge, end_edge),
                    None => {
                        let (va, vb) = (corners.vertex(&mut body, a, pa)?, corners.vertex(&mut body, b, pb)?);
                        let (start, finish) = if swept.forward { (va, vb) } else { (vb, va) };
                        let seam = body.add_edge(swept.curve.clone(), swept.range(), start, finish)?;
                        (seam, seam)
                    }
                };
                let (va_param, vb_param) = (swept.v(swept.t_from), swept.v(swept.t_to));
                let mut coedges = Vec::new();
                if let Some(e) = arc_a {
                    coedges.push(Coedge::new(e, false).with_pcurve(uv_line([0.0, va_param], [1.0, 0.0])));
                }
                coedges.push(Coedge::new(end_edge, !swept.forward).with_pcurve(swept.seam_pcurve(angle)));
                if let Some(e) = arc_b {
                    coedges.push(Coedge::new(e, true).with_pcurve(uv_line([0.0, vb_param], [1.0, 0.0])));
                }
                coedges.push(Coedge::new(start_edge, swept.forward).with_pcurve(swept.seam_pcurve(0.0)));
                body.add_face(swept.surface, vb_param > va_param, vec![Loop::new(coedges)])?
            };
            faces.push(face);
        }
        match full {
            true => {
                body.add_shell(faces)?;
            }
            false => {
                all_faces.extend(faces);
                start_caps.push(Loop::new(start_cap));
                end_cap.reverse();
                end_caps.push(Loop::new(end_cap));
            }
        }
    }
    if !full {
        // The region itself closes the solid at angle 0, its normal (the
        // sketch's) pointing away from the turn; its copy closes it at the end,
        // facing the other way.
        all_faces.push(body.add_face(Plane::new(*plane), true, start_caps)?);
        all_faces.push(body.add_face(Plane::new(end_plane), false, end_caps)?);
        body.add_shell(all_faces)?;
    }
    Ok(body)
}

/// A profile segment as an edge in the profile plane of `frame`, from the
/// vertices `along` gives for its direction along the curve.
fn profile_edge<S: Space<3>>(
    body: &mut Body<S>,
    segment: &Segment,
    frame: &Placement<S>,
    along: impl Fn(bool) -> (VertexId, VertexId),
) -> Result<EdgeId, ModelError> {
    let in_plane = |p: Vector2<f64>| frame.point(Vector3::new(p.x, 0.0, p.y));
    let (curve, range): (AnyCurve<S>, _) = match *segment {
        Segment::Line { start, end } => {
            let length = (end - start).norm();
            let direction = (end - start) / length;
            let line = Line::new(in_plane(start), frame.vector(Vector3::new(direction.x, 0.0, direction.y)));
            (line.into(), (0.0, length))
        }
        Segment::Arc {
            center,
            radius,
            start_angle,
            sweep,
        } => {
            let circle = Circle::new(meridional(frame, center), radius);
            let (t0, t1) = (start_angle, start_angle + sweep);
            (circle.into(), (t0.min(t1), t0.max(t1)))
        }
    };
    let (start, end) = along(segment_forward(segment));
    Ok(body.add_edge(curve, range, start, end)?)
}

/// The frame of an arc's circle in `frame`'s profile plane: x radial, y along
/// the axis, so the circle's angle is the sketch angle.
fn meridional<S: Space<3>>(frame: &Placement<S>, center: Vector2<f64>) -> Placement<S> {
    Placement::from_axes(
        frame.point(Vector3::new(center.x, 0.0, center.y)),
        frame.vector(-Vector3::y()),
        frame.vector(Vector3::x()),
    )
}

/// Whether a segment's curve runs the same way as the profile: lines always,
/// arcs when anticlockwise.
fn segment_forward(segment: &Segment) -> bool {
    match *segment {
        Segment::Line { .. } => true,
        Segment::Arc { sweep, .. } => sweep > 0.0,
    }
}

/// A line segment's pcurve on a cap, whose uv are sketch coordinates; arcs have
/// none, there being no 2D circle pcurve yet.
fn cap_pcurve<S: Space<3>>(segment: &Segment) -> Option<Line<golf_brep::FaceUv<S>, 2>> {
    match *segment {
        Segment::Line { start, end } => {
            let direction = (end - start).normalize();
            Some(uv_line([start.x, start.y], [direction.x, direction.y]))
        }
        Segment::Arc { .. } => None,
    }
}

/// A profile's corners as they're needed: their vertices (at both angles of a
/// partial revolution), and the arcs off-axis ones sweep. Made on demand, so a
/// corner on the axis that only bounds a disc leaves no stray vertex.
struct Corners<'a, S: Space<3>, F: Fn(f64) -> bool> {
    frame: &'a Placement<S>,
    end: &'a Placement<S>,
    angle: f64,
    full: bool,
    on_axis: &'a F,
    vertices: HashMap<usize, VertexId>,
    end_vertices: HashMap<usize, VertexId>,
    arcs: HashMap<usize, EdgeId>,
}

impl<S: Space<3>, F: Fn(f64) -> bool> Corners<'_, S, F> {
    /// The vertex at corner `k`, which is at `p` in the sketch, at angle 0.
    fn vertex(&mut self, body: &mut Body<S>, k: usize, p: Vector2<f64>) -> Result<VertexId, ModelError> {
        if let Some(&v) = self.vertices.get(&k) {
            return Ok(v);
        }
        let v = body.add_vertex(self.frame.point(Vector3::new(p.x, 0.0, p.y)))?;
        self.vertices.insert(k, v);
        Ok(v)
    }

    /// The vertex at corner `k` at the end of the revolution: the same vertex
    /// for a full turn, or on the axis.
    fn end_vertex(&mut self, body: &mut Body<S>, k: usize, p: Vector2<f64>) -> Result<VertexId, ModelError> {
        if self.full || (self.on_axis)(p.x) {
            return self.vertex(body, k, p);
        }
        if let Some(&v) = self.end_vertices.get(&k) {
            return Ok(v);
        }
        let v = body.add_vertex(self.end.point(Vector3::new(p.x, 0.0, p.y)))?;
        self.end_vertices.insert(k, v);
        Ok(v)
    }

    /// The arc corner `k` sweeps, parametrised by the angle of revolution;
    /// `None` on the axis.
    fn arc(&mut self, body: &mut Body<S>, k: usize, p: Vector2<f64>) -> Result<Option<EdgeId>, ModelError> {
        if (self.on_axis)(p.x) {
            return Ok(None);
        }
        if let Some(&e) = self.arcs.get(&k) {
            return Ok(Some(e));
        }
        let start = self.vertex(body, k, p)?;
        let finish = self.end_vertex(body, k, p)?;
        let centre = Placement::new(self.frame.point(Vector3::new(0.0, 0.0, p.y)), self.frame.rotation);
        let e = body.add_edge(Circle::new(centre, p.x), (0.0, self.angle), start, finish)?;
        self.arcs.insert(k, e);
        Ok(Some(e))
    }
}

/// A profile segment revolved: its surface, and the segment as that surface's
/// seam edge at angle 0.
struct Revolved<S: Space<3>> {
    surface: AnySurface<S>,
    curve: AnyCurve<S>,
    t_from: f64,
    t_to: f64,
    forward: bool,
    /// The surface's v along the seam is `v_origin + v_rate · t`.
    v_origin: f64,
    v_rate: f64,
}

impl<S: Space<3>> Revolved<S> {
    fn new(index: usize, segment: &Segment, frame: &Placement<S>) -> Result<Self, ModelError> {
        let in_profile_plane = |p: Vector2<f64>| frame.point(Vector3::new(p.x, 0.0, p.y));
        let on_axis_at = |height: f64| {
            Placement::new(frame.point(Vector3::new(0.0, 0.0, height)), frame.rotation)
        };
        Ok(match *segment {
            // Lines are parametrised by length from their start; the surfaces'
            // v is height, measured from the frame's origin.
            Segment::Line { start, end } => {
                let length = (end - start).norm();
                let direction = (end - start) / length;
                let curve = Line::new(
                    in_profile_plane(start),
                    frame.vector(Vector3::new(direction.x, 0.0, direction.y)),
                );
                let surface: AnySurface<S> = if direction.x.abs() <= 1e-12 {
                    Cylinder::new(*frame, start.x).into()
                } else {
                    // Radius at height 0, and the angle from the axis.
                    let slope = direction.x / direction.y;
                    Cone::new(*frame, start.x - start.y * slope, slope.atan()).into()
                };
                Self {
                    surface,
                    curve: curve.into(),
                    t_from: 0.0,
                    t_to: length,
                    forward: true,
                    v_origin: start.y,
                    v_rate: direction.y,
                }
            }
            // Arcs are parametrised by their angle in the profile plane, which
            // is the sphere's latitude and the torus's minor angle.
            Segment::Arc {
                center,
                radius,
                start_angle,
                sweep,
            } => {
                let circle_frame = meridional(frame, center);
                let surface: AnySurface<S> = if center.x.abs() <= 1e-12 * (1.0 + radius) {
                    Sphere::from_placement(on_axis_at(center.y), radius).into()
                } else if center.x > radius {
                    Torus::new(on_axis_at(center.y), center.x, radius).into()
                } else {
                    return Err(ModelError::SpindleTorus { index });
                };
                Self {
                    surface,
                    curve: Circle::new(circle_frame, radius).into(),
                    t_from: start_angle,
                    t_to: start_angle + sweep,
                    forward: sweep > 0.0,
                    v_origin: 0.0,
                    v_rate: 1.0,
                }
            }
        })
    }

    fn range(&self) -> (f64, f64) {
        (self.t_from.min(self.t_to), self.t_from.max(self.t_to))
    }

    fn v(&self, t: f64) -> f64 {
        self.v_origin + self.v_rate * t
    }

    /// The seam's pcurve on the side of the face at angle `u`.
    fn seam_pcurve(&self, u: f64) -> Line<golf_brep::FaceUv<S>, 2> {
        uv_line([u, self.v_origin], [0.0, self.v_rate])
    }
}

/// Arcs with their start angle in (-π, π], so a right-half arc's angles are
/// latitudes in [-π/2, π/2].
fn normalise(segment: &Segment) -> Segment {
    match *segment {
        Segment::Arc {
            center,
            radius,
            start_angle,
            sweep,
        } => Segment::Arc {
            center,
            radius,
            start_angle: PI - (PI - start_angle).rem_euclid(TAU),
            sweep,
        },
        line => line,
    }
}

/// The smallest sketch x a segment reaches.
fn closest_to_axis(segment: &Segment) -> f64 {
    match *segment {
        Segment::Line { start, end } => start.x.min(end.x),
        Segment::Arc {
            center,
            radius,
            start_angle,
            sweep,
        } => {
            // The circle's leftmost point is at angle π; is it on the arc?
            let (lo, hi) = (
                start_angle.min(start_angle + sweep),
                start_angle.max(start_angle + sweep),
            );
            let reaches_left = (lo..=hi).contains(&(PI + TAU * ((lo - PI) / TAU).ceil()));
            match reaches_left {
                true => center.x - radius,
                false => segment.start().x.min(segment.end().x),
            }
        }
    }
}

fn region_size(region: &Region) -> f64 {
    region
        .profiles()
        .flat_map(|p| p.segments())
        .map(|s| s.start().norm() + s.length())
        .fold(0.0, f64::max)
}
