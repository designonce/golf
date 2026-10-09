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
use nalgebra::Vector2;
use nalgebra::Vector3;

use crate::error::ModelError;
use crate::extrude::uv_line;
use crate::region::Region;
use crate::segment::Segment;

/// Revolves `region`, drawn in the sketch plane `plane`, a full turn about the
/// sketch's y axis into a solid.
///
/// The region must lie on the side `x >= 0`, with sketch x the distance from the
/// axis and sketch y the height along it. Segments lying on the axis bound the
/// solid without making a face; each other segment sweeps a plane, cylinder,
/// cone, sphere or torus, and is itself the seam of that face. Holes revolve
/// into internal voids, each its own shell.
pub fn revolve<S: Space<3>>(plane: &Placement<S>, region: &Region) -> Result<Body<S>, ModelError> {
    // The frame the surfaces of revolution are placed in: local x is the
    // sketch's x (the radial direction at angle 0), local z the axis.
    let frame = Placement::from_axes(
        plane.origin,
        plane.vector(Vector3::y()),
        plane.vector(Vector3::x()),
    );
    let on_axis = |x: f64| x.abs() <= 1e-12 * (1.0 + region_size(region));

    let mut body = Body::with_tag(plane.origin.tag());
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
            on_axis: &on_axis,
            vertices: HashMap::new(),
            circles: HashMap::new(),
        };
        let mut faces = Vec::new();

        for (index, segment) in segments.iter().enumerate() {
            let (a, b) = (index, (index + 1) % n);
            let (pa, pb) = (segment.start(), segment.end());
            if matches!(segment, Segment::Line { .. }) && on_axis(pa.x) && on_axis(pb.x) {
                continue;
            }
            let circle_a = corners.circle(&mut body, a, pa)?;
            let circle_b = corners.circle(&mut body, b, pb)?;

            let face = match segment {
                // Perpendicular to the axis: a disc or annulus, bounded by the
                // circles alone. Its plane's normal is the axis; the solid is to
                // the profile's left, so the face points along the axis when the
                // profile runs inwards.
                &Segment::Line { start, end }
                    if (end.y - start.y).abs() <= 1e-12 * (end - start).norm() =>
                {
                    let surface = Plane::new(Placement::new(
                        frame.point(Vector3::new(0.0, 0.0, start.y)),
                        frame.rotation,
                    ));
                    let loops = [circle_a.map(|e| (e, false)), circle_b.map(|e| (e, true))]
                        .into_iter()
                        .flatten()
                        .map(|(e, reversed)| Loop::new(vec![Coedge::new(e, reversed)]))
                        .collect();
                    body.add_face(surface, end.x < start.x, loops)?
                }
                _ => {
                    let swept = Revolved::new(index, segment, &frame)?;
                    let (va, vb) = (
                        corners.vertex(&mut body, a, pa)?,
                        corners.vertex(&mut body, b, pb)?,
                    );
                    let (start, end) = match swept.forward {
                        true => (va, vb),
                        false => (vb, va),
                    };
                    let seam = body.add_edge(swept.curve.clone(), swept.range(), start, end)?;
                    // In the face's (u, v), u the angle of revolution: circle A
                    // along v_a, up the seam at u = 2π, back along circle B at v_b,
                    // down the seam at u = 0. That runs anticlockwise about the
                    // surface normal exactly when v increases from A to B, which
                    // is when the normal points out of the solid.
                    let (va_param, vb_param) = (swept.v(swept.t_from), swept.v(swept.t_to));
                    let mut coedges = Vec::new();
                    if let Some(e) = circle_a {
                        coedges.push(
                            Coedge::new(e, false).with_pcurve(uv_line([0.0, va_param], [1.0, 0.0])),
                        );
                    }
                    coedges.push(
                        Coedge::new(seam, !swept.forward).with_pcurve(swept.seam_pcurve(TAU)),
                    );
                    if let Some(e) = circle_b {
                        coedges.push(
                            Coedge::new(e, true).with_pcurve(uv_line([0.0, vb_param], [1.0, 0.0])),
                        );
                    }
                    coedges
                        .push(Coedge::new(seam, swept.forward).with_pcurve(swept.seam_pcurve(0.0)));
                    body.add_face(swept.surface, vb_param > va_param, vec![Loop::new(coedges)])?
                }
            };
            faces.push(face);
        }
        body.add_shell(faces)?;
    }
    Ok(body)
}

/// A profile's corners as they're needed: their vertices, and the circles
/// off-axis ones sweep. Made on demand, so a corner on the axis that only bounds
/// a disc leaves no stray vertex.
struct Corners<'a, S: Space<3>, F: Fn(f64) -> bool> {
    frame: &'a Placement<S>,
    on_axis: &'a F,
    vertices: HashMap<usize, VertexId>,
    circles: HashMap<usize, EdgeId>,
}

impl<S: Space<3>, F: Fn(f64) -> bool> Corners<'_, S, F> {
    /// The vertex at corner `k`, which is at `p` in the sketch.
    fn vertex(
        &mut self,
        body: &mut Body<S>,
        k: usize,
        p: Vector2<f64>,
    ) -> Result<VertexId, ModelError> {
        if let Some(&v) = self.vertices.get(&k) {
            return Ok(v);
        }
        let v = body.add_vertex(self.frame.point(Vector3::new(p.x, 0.0, p.y)))?;
        self.vertices.insert(k, v);
        Ok(v)
    }

    /// The circle corner `k` sweeps, parametrised by the angle of revolution;
    /// `None` on the axis.
    fn circle(
        &mut self,
        body: &mut Body<S>,
        k: usize,
        p: Vector2<f64>,
    ) -> Result<Option<EdgeId>, ModelError> {
        if (self.on_axis)(p.x) {
            return Ok(None);
        }
        if let Some(&e) = self.circles.get(&k) {
            return Ok(Some(e));
        }
        let v = self.vertex(body, k, p)?;
        let centre = Placement::new(
            self.frame.point(Vector3::new(0.0, 0.0, p.y)),
            self.frame.rotation,
        );
        let e = body.add_edge(Circle::new(centre, p.x), (0.0, TAU), v, v)?;
        self.circles.insert(k, e);
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
                let circle_frame = Placement::from_axes(
                    in_profile_plane(center),
                    frame.vector(-Vector3::y()),
                    frame.vector(Vector3::x()),
                );
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
