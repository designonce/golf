//! Extrusions built as a stack of rings: the profile inset by some amount at a
//! series of heights, joined by bands that are ruled (straight) or round
//! (quarter-circle fillets).
//!
//! A plain extrusion is two rings joined by a vertical ruled band; a draft
//! insets the top ring; a chamfered end adds a ruled band out to an inset cap
//! ring, and a filleted end a round one.

use core::f64::consts::FRAC_PI_2;
use core::f64::consts::PI;

use golf_brep::Body;
use golf_brep::Coedge;
use golf_brep::EdgeId;
use golf_brep::FaceUv;
use golf_brep::Loop;
use golf_brep::VertexId;
use golf_geom::AnyCurve;
use golf_geom::AnySurface;
use golf_geom::Circle;
use golf_geom::Cone;
use golf_geom::Cylinder;
use golf_geom::Ellipse;
use golf_geom::Line;
use golf_geom::Placement;
use golf_geom::Plane;
use golf_geom::Sphere;
use golf_geom::Torus;
use golf_manifold::Point;
use golf_manifold::Space;
use golf_manifold::Vector;
use golf_sketch::Profile;
use golf_sketch::Region;
use golf_sketch::Segment;
use nalgebra::Vector2;
use nalgebra::Vector3;

use crate::error::ModelError;
use crate::extrude::uv_line;

/// The profile inset by `inset` (into the material) at height `z`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Ring {
    pub z: f64,
    pub inset: f64,
}

/// How two consecutive rings are joined.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum Band {
    /// Straight from one ring to the next.
    Ruled,
    /// A quarter-round fillet between a side ring (outermost) and a cap ring
    /// inset and raised (or lowered) from it by the fillet radius.
    Round,
}

/// Builds the solid whose rings (in increasing height) are joined by `bands`
/// (one fewer), capped at both ends, from `region` in `plane`.
pub(crate) fn build<S: Space<3>>(
    plane: &Placement<S>,
    region: &Region,
    rings: &[Ring],
    bands: &[Band],
) -> Result<Body<S>, ModelError> {
    debug_assert_eq!(rings.len(), bands.len() + 1);
    let mut body = Body::with_tag(plane.origin.tag());
    let mut faces = Vec::new();
    let mut bottom_loops = Vec::new();
    let mut top_loops = Vec::new();
    for profile in region.profiles() {
        let shape = Shape::new(profile)?;
        let mut stack = Stack {
            plane,
            shape: &shape,
            body: &mut body,
        };
        // Each ring's inset segments, vertices and edges.
        let insets: Vec<Vec<Inset>> = rings
            .iter()
            .enumerate()
            .map(|(i, ring)| {
                // Only the cap ring of a round band may collapse an arc to a
                // point (a sphere corner).
                let may_collapse =
                    (i > 0 && bands[i - 1] == Band::Round && rings[i].inset > rings[i - 1].inset)
                        || (i + 1 < rings.len()
                            && bands[i] == Band::Round
                            && rings[i].inset > rings[i + 1].inset);
                shape.inset(ring.inset, may_collapse)
            })
            .collect::<Result<_, _>>()?;
        let ring_edges: Vec<RingEdges> = rings
            .iter()
            .zip(&insets)
            .map(|(ring, inset)| stack.ring(ring.z, ring.inset, inset))
            .collect::<Result<_, _>>()?;
        for (i, band) in bands.iter().enumerate() {
            let (lo, hi) = (i, i + 1);
            faces.extend(stack.band(
                *band,
                (&rings[lo], &insets[lo], &ring_edges[lo]),
                (&rings[hi], &insets[hi], &ring_edges[hi]),
            )?);
        }
        let (first, last) = (&ring_edges[0], &ring_edges[rings.len() - 1]);
        if let Some(l) = stack.cap(&insets[0], first, false) {
            bottom_loops.push(l);
        }
        if let Some(l) = stack.cap(&insets[rings.len() - 1], last, true) {
            top_loops.push(l);
        }
    }
    let at = |z: f64| Placement::new(plane.point(Vector3::new(0.0, 0.0, z)), plane.rotation);
    // The bottom cap's plane points up, into the solid, so the face flips it.
    if !bottom_loops.is_empty() {
        faces.push(body.add_face(Plane::new(at(rings[0].z)), false, bottom_loops)?);
    }
    if !top_loops.is_empty() {
        faces.push(body.add_face(Plane::new(at(rings[rings.len() - 1].z)), true, top_loops)?);
    }
    body.add_shell(faces)?;
    Ok(body)
}

/// How two consecutive segments meet.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Joint {
    /// Tangent: offsets meet along the shared normal.
    Smooth,
    /// Two lines at an angle: offsets meet where the offset lines cross.
    Mitre,
}

/// A profile with what its insets need: each joint's kind and the segments'
/// directions.
struct Shape<'a> {
    segments: &'a [Segment],
    /// Joint `k` is where segment `k` starts.
    joints: Vec<Joint>,
}

impl<'a> Shape<'a> {
    fn new(profile: &'a Profile) -> Result<Self, ModelError> {
        let segments = profile.segments();
        let n = segments.len();
        let joints = (0..n)
            .map(|k| {
                let (before, after) = (&segments[(k + n - 1) % n], &segments[k]);
                let (incoming, outgoing) = (tangent(before, 1.0), tangent(after, 0.0));
                if incoming.perp(&outgoing).abs() <= 1e-9 && incoming.dot(&outgoing) > 0.0 {
                    Ok(Joint::Smooth)
                } else if matches!(before, Segment::Line { .. })
                    && matches!(after, Segment::Line { .. })
                {
                    Ok(Joint::Mitre)
                } else {
                    Err(ModelError::UnsupportedCorner { index: k })
                }
            })
            .collect::<Result<_, _>>()?;
        Ok(Self { segments, joints })
    }

    /// Joint `k` moved `inset` into the material.
    fn joint(&self, k: usize, inset: f64) -> Vector2<f64> {
        let n = self.segments.len();
        let point = self.segments[k].start();
        let normal_in = |s: &Segment, at: f64| left(tangent(s, at));
        match self.joints[k] {
            Joint::Smooth => point + normal_in(&self.segments[k], 0.0) * inset,
            Joint::Mitre => {
                let (a, b) = (
                    normal_in(&self.segments[(k + n - 1) % n], 1.0),
                    normal_in(&self.segments[k], 0.0),
                );
                point + (a + b) * (inset / (1.0 + a.dot(&b)))
            }
        }
    }

    /// The profile's segments inset by `inset`. An arc shrunk to nothing is
    /// `Collapsed` if `may_collapse`, an error otherwise; a segment turned
    /// inside out is always an error.
    fn inset(&self, inset: f64, may_collapse: bool) -> Result<Vec<Inset>, ModelError> {
        let n = self.segments.len();
        let size = self
            .segments
            .iter()
            .map(|s| s.start().norm() + s.length())
            .fold(1.0, f64::max);
        self.segments
            .iter()
            .enumerate()
            .map(|(k, segment)| {
                let (start, end) = (self.joint(k, inset), self.joint((k + 1) % n, inset));
                match *segment {
                    Segment::Line { start: a, end: b } => {
                        let direction = (b - a).normalize();
                        let length = (end - start).dot(&direction);
                        if length <= 1e-9 * size {
                            return Err(ModelError::InsetCollapses { index: k });
                        }
                        Ok(Inset::Line {
                            start,
                            direction,
                            length,
                        })
                    }
                    Segment::Arc {
                        center,
                        radius,
                        start_angle,
                        sweep,
                    } => {
                        // The material is towards the centre of an anticlockwise
                        // arc, away from it for a clockwise one.
                        let shrunk = radius - sweep.signum() * inset;
                        if shrunk.abs() <= 1e-9 * size && may_collapse {
                            Ok(Inset::Collapsed)
                        } else if shrunk <= 1e-9 * size {
                            Err(ModelError::InsetCollapses { index: k })
                        } else {
                            Ok(Inset::Arc {
                                center,
                                radius: shrunk,
                                start_angle,
                                sweep,
                            })
                        }
                    }
                }
            })
            .collect()
    }
}

/// A segment inset.
#[derive(Clone, Copy, Debug)]
enum Inset {
    Line {
        start: Vector2<f64>,
        direction: Vector2<f64>,
        length: f64,
    },
    Arc {
        center: Vector2<f64>,
        radius: f64,
        start_angle: f64,
        sweep: f64,
    },
    /// An arc shrunk to its centre.
    Collapsed,
}

impl Inset {
    /// The edge curve's parameters where the profile enters and leaves it.
    fn params(&self) -> (f64, f64) {
        match *self {
            Self::Line { length, .. } => (0.0, length),
            Self::Arc {
                start_angle, sweep, ..
            } => (start_angle, start_angle + sweep),
            Self::Collapsed { .. } => (0.0, 0.0),
        }
    }

    fn forward(&self) -> bool {
        let (from, to) = self.params();
        to > from
    }
}

/// A ring's vertices (one per joint, shared across a collapsed arc) and edges
/// (none for a collapsed arc).
struct RingEdges {
    vertices: Vec<VertexId>,
    edges: Vec<Option<EdgeId>>,
}

struct Stack<'a, 'b, S: Space<3>> {
    plane: &'a Placement<S>,
    shape: &'a Shape<'b>,
    body: &'a mut Body<S>,
}

impl<S: Space<3>> Stack<'_, '_, S> {
    fn point(&self, p: Vector2<f64>, z: f64) -> Point<S, 3> {
        self.plane.point(Vector3::new(p.x, p.y, z))
    }

    fn vector(&self, v: Vector3<f64>) -> Vector<S, 3> {
        self.plane.vector(v)
    }

    /// `plane`'s frame moved to `p` at height `z`.
    fn frame_at(&self, p: Vector2<f64>, z: f64) -> Placement<S> {
        Placement::new(self.point(p, z), self.plane.rotation)
    }

    fn ring(&mut self, z: f64, inset: f64, insets: &[Inset]) -> Result<RingEdges, ModelError> {
        let n = insets.len();
        // Joints either side of a collapsed arc are one point, so one vertex:
        // walk round from a joint that starts a group, starting a new vertex
        // only after a segment that didn't collapse.
        let collapsed = |k: usize| matches!(insets[k], Inset::Collapsed);
        let first = (0..n).find(|&k| !collapsed((k + n - 1) % n)).unwrap_or(0);
        let mut vertices: Vec<Option<VertexId>> = vec![None; n];
        let mut current = None;
        for step in 0..n {
            let k = (first + step) % n;
            let joins_previous = step > 0 && collapsed((k + n - 1) % n);
            if !joins_previous || current.is_none() {
                current = Some(
                    self.body
                        .add_vertex(self.point(self.shape.joint(k, inset), z))?,
                );
            }
            vertices[k] = current;
        }
        let vertices: Vec<VertexId> = vertices
            .into_iter()
            .map(|v| v.expect("every joint visited"))
            .collect();
        let edges = (0..n)
            .map(|k| {
                let (a, b) = (vertices[k], vertices[(k + 1) % n]);
                let inset = &insets[k];
                let (curve, range): (AnyCurve<S>, (f64, f64)) = match *inset {
                    Inset::Collapsed => return Ok(None),
                    Inset::Line {
                        start,
                        direction,
                        length,
                    } => (
                        Line::new(
                            self.point(start, z),
                            self.vector(Vector3::new(direction.x, direction.y, 0.0)),
                        )
                        .into(),
                        (0.0, length),
                    ),
                    Inset::Arc { center, radius, .. } => {
                        let (from, to) = inset.params();
                        (
                            Circle::new(self.frame_at(center, z), radius).into(),
                            (from.min(to), from.max(to)),
                        )
                    }
                };
                let (start, end) = if inset.forward() { (a, b) } else { (b, a) };
                Ok(Some(self.body.add_edge(curve, range, start, end)?))
            })
            .collect::<Result<_, ModelError>>()?;
        Ok(RingEdges { vertices, edges })
    }

    /// The cap closing a ring, or `None` if every segment collapsed. The top
    /// runs the profile's way; the bottom, facing down, runs back.
    fn cap(&self, insets: &[Inset], ring: &RingEdges, top: bool) -> Option<Loop<S>> {
        let mut coedges: Vec<Coedge<S>> = insets
            .iter()
            .zip(&ring.edges)
            .filter_map(|(inset, edge)| {
                let edge = (*edge)?;
                let coedge = Coedge::new(edge, top != inset.forward());
                Some(match *inset {
                    Inset::Line {
                        start, direction, ..
                    } => {
                        coedge.with_pcurve(uv_line([start.x, start.y], [direction.x, direction.y]))
                    }
                    _ => coedge,
                })
            })
            .collect();
        if coedges.is_empty() {
            return None;
        }
        if !top {
            coedges.reverse();
        }
        Some(Loop::new(coedges))
    }

    /// The faces of the band between rings `lo` and `hi`, one per segment.
    fn band(
        &mut self,
        band: Band,
        (lo_ring, lo, lo_edges): (&Ring, &[Inset], &RingEdges),
        (hi_ring, hi, hi_edges): (&Ring, &[Inset], &RingEdges),
    ) -> Result<Vec<golf_brep::FaceId>, ModelError> {
        let n = lo.len();
        // Joint edges up the band, stored from `lo` to `hi` for a ruled band and
        // from the side ring to the cap ring for a round one.
        let cap_is_hi = hi_ring.inset > lo_ring.inset;
        let joints: Vec<(EdgeId, Option<JointPcurves>)> = (0..n)
            .map(|k| self.joint_edge(band, k, (lo_ring, lo_edges), (hi_ring, hi_edges), cap_is_hi))
            .collect::<Result<_, _>>()?;
        let joint_runs_up = band == Band::Ruled || cap_is_hi;

        (0..n)
            .map(|k| {
                let next = (k + 1) % n;
                let (surface, same_sense, pcurves) =
                    self.band_surface(band, k, (lo_ring, &lo[k]), (hi_ring, &hi[k]), cap_is_hi)?;
                // Along lo from A to B, up joint B, back along hi, down joint A.
                let mut coedges = Vec::new();
                if let Some(e) = lo_edges.edges[k] {
                    coedges.push(with(Coedge::new(e, !lo[k].forward()), pcurves.lo.clone()));
                }
                coedges.push(with(
                    Coedge::new(joints[next].0, !joint_runs_up),
                    pcurves.joint(&joints[next].1, true),
                ));
                if let Some(e) = hi_edges.edges[k] {
                    coedges.push(with(Coedge::new(e, hi[k].forward()), pcurves.hi.clone()));
                }
                coedges.push(with(
                    Coedge::new(joints[k].0, joint_runs_up),
                    pcurves.joint(&joints[k].1, false),
                ));
                Ok(self
                    .body
                    .add_face(surface, same_sense, vec![Loop::new(coedges)])?)
            })
            .collect()
    }

    /// The edge up joint `k` of a band.
    fn joint_edge(
        &mut self,
        band: Band,
        k: usize,
        (lo_ring, lo_edges): (&Ring, &RingEdges),
        (hi_ring, hi_edges): (&Ring, &RingEdges),
        cap_is_hi: bool,
    ) -> Result<(EdgeId, Option<JointPcurves>), ModelError> {
        let (lo_point, hi_point) = (
            self.point(self.shape.joint(k, lo_ring.inset), lo_ring.z),
            self.point(self.shape.joint(k, hi_ring.inset), hi_ring.z),
        );
        let (lo_v, hi_v) = (lo_edges.vertices[k], hi_edges.vertices[k]);
        match band {
            Band::Ruled => {
                let line = Line::new(
                    lo_point,
                    Vector::with_tag(hi_point.coords - lo_point.coords, lo_point.tag()),
                );
                Ok((self.body.add_edge(line, (0.0, 1.0), lo_v, hi_v)?, None))
            }
            Band::Round => {
                // From the side ring (outermost) to the cap ring, about the
                // axis ring: the side ring's height and the cap ring's inset.
                let (side, cap, side_v, cap_v) = match cap_is_hi {
                    true => (lo_ring, hi_ring, lo_v, hi_v),
                    false => (hi_ring, lo_ring, hi_v, lo_v),
                };
                let radius = cap.inset - side.inset;
                let axis = self.point(self.shape.joint(k, cap.inset), side.z);
                let side_point = self.point(self.shape.joint(k, side.inset), side.z);
                let outward = side_point.coords - axis.coords;
                let vertical = self.vector(Vector3::new(0.0, 0.0, (cap.z - side.z).signum()));
                let x = Vector::with_tag(outward.normalize(), axis.tag());
                let frame = Placement::from_axes(
                    axis,
                    Vector::with_tag(x.coords.cross(&vertical.coords), axis.tag()),
                    x,
                );
                let curve: AnyCurve<S> = match self.shape.joints[k] {
                    Joint::Smooth => Circle::new(frame, radius).into(),
                    // Two equal cylinders whose axes cross meet in an ellipse in
                    // the mitre plane: across the mitre to the side corner, and
                    // up the fillet radius.
                    Joint::Mitre => Ellipse::new(frame, outward.norm(), radius).into(),
                };
                let pcurves = (self.shape.joints[k] == Joint::Smooth).then_some(JointPcurves);
                Ok((
                    self.body.add_edge(curve, (0.0, FRAC_PI_2), side_v, cap_v)?,
                    pcurves,
                ))
            }
        }
    }

    /// A band face's surface, sense, and its edges' pcurves.
    fn band_surface(
        &self,
        band: Band,
        k: usize,
        (lo_ring, lo): (&Ring, &Inset),
        (hi_ring, hi): (&Ring, &Inset),
        cap_is_hi: bool,
    ) -> Result<(AnySurface<S>, bool, BandPcurves<S>), ModelError> {
        let segment = &self.shape.segments[k];
        let rise = hi_ring.z - lo_ring.z;
        match (band, segment) {
            (Band::Ruled, Segment::Line { start, end }) => {
                // The plane through both inset lines, x along them and its normal
                // out of the material (right of the profile, tilted by the inset).
                let along = (end - start).normalize();
                let across = left(along) * (hi_ring.inset - lo_ring.inset);
                let x = Vector3::new(along.x, along.y, 0.0);
                let w = Vector3::new(across.x, across.y, rise);
                let origin = self.point(self.shape.joint(k, lo_ring.inset), lo_ring.z);
                let frame = Placement::from_axes(origin, self.vector(x.cross(&w)), self.vector(x));
                let uv = |p: Vector2<f64>, z: f64| {
                    let local = frame.to_local(self.point(p, z));
                    [local.x, local.y]
                };
                let line_pcurve = |inset: &Inset, z: f64| match *inset {
                    Inset::Line { start, .. } => Some(uv_line(uv(start, z), [1.0, 0.0])),
                    _ => None,
                };
                let joint_at = |j: usize| {
                    let (a, b) = (
                        uv(self.shape.joint(j, lo_ring.inset), lo_ring.z),
                        uv(self.shape.joint(j, hi_ring.inset), hi_ring.z),
                    );
                    uv_line(a, [b[0] - a[0], b[1] - a[1]])
                };
                let n = self.shape.segments.len();
                Ok((
                    Plane::new(frame).into(),
                    true,
                    BandPcurves {
                        lo: line_pcurve(lo, lo_ring.z),
                        hi: line_pcurve(hi, hi_ring.z),
                        joints: [Some(joint_at(k)), Some(joint_at((k + 1) % n))],
                        round: None,
                    },
                ))
            }
            (Band::Ruled, &Segment::Arc { center, .. }) => {
                let (Inset::Arc { radius: r_lo, .. }, Inset::Arc { radius: r_hi, .. }) = (lo, hi)
                else {
                    return Err(ModelError::InsetCollapses { index: k });
                };
                let frame = self.frame_at(center, lo_ring.z);
                let surface: AnySurface<S> = match (r_hi - r_lo).abs() <= 1e-12 * r_lo.max(1.0) {
                    true => Cylinder::new(frame, *r_lo).into(),
                    false => Cone::new(frame, *r_lo, ((r_hi - r_lo) / rise).atan()).into(),
                };
                let (from, to) = lo.params();
                Ok((
                    surface,
                    lo.forward(),
                    BandPcurves {
                        lo: Some(uv_line([0.0, 0.0], [1.0, 0.0])),
                        hi: Some(uv_line([0.0, rise], [1.0, 0.0])),
                        joints: [
                            Some(uv_line([from, 0.0], [0.0, rise])),
                            Some(uv_line([to, 0.0], [0.0, rise])),
                        ],
                        round: None,
                    },
                ))
            }
            (Band::Round, Segment::Line { start, end }) => {
                let (side, cap) = if cap_is_hi {
                    (lo_ring, hi_ring)
                } else {
                    (hi_ring, lo_ring)
                };
                let radius = cap.inset - side.inset;
                let along = (end - start).normalize();
                // u from outward (0) to vertical (π/2), v along x × y.
                let axis = self.point(self.shape.joint(k, cap.inset), side.z);
                let x = -left(along);
                let x = Vector3::new(x.x, x.y, 0.0);
                let y = Vector3::new(0.0, 0.0, (cap.z - side.z).signum());
                let z = x.cross(&y);
                let frame = Placement::from_axes(axis, self.vector(z), self.vector(x));
                let v_of = |p: Vector2<f64>, height: f64| frame.to_local(self.point(p, height)).z;
                let rate = Vector3::new(along.x, along.y, 0.0).dot(&z);
                let line_pcurve = |inset: &Inset, ring: &Ring| match *inset {
                    Inset::Line { start, .. } => {
                        let u = if ring.inset == side.inset {
                            0.0
                        } else {
                            FRAC_PI_2
                        };
                        Some(uv_line([u, v_of(start, ring.z)], [0.0, rate]))
                    }
                    _ => None,
                };
                let n = self.shape.segments.len();
                let joint_v = |j: usize| v_of(self.shape.joint(j, cap.inset), side.z);
                Ok((
                    Cylinder::new(frame, radius).into(),
                    true,
                    BandPcurves {
                        lo: line_pcurve(lo, lo_ring),
                        hi: line_pcurve(hi, hi_ring),
                        joints: [None, None],
                        round: Some(RoundPcurves {
                            v: [joint_v(k), joint_v((k + 1) % n)],
                            u_of_s: (0.0, 1.0),
                            vertical_is_u: false,
                        }),
                    },
                ))
            }
            (Band::Round, &Segment::Arc { center, sweep, .. }) => {
                let (side, cap) = if cap_is_hi {
                    (lo_ring, hi_ring)
                } else {
                    (hi_ring, lo_ring)
                };
                let (side_inset, cap_inset) = if cap_is_hi { (lo, hi) } else { (hi, lo) };
                let radius = cap.inset - side.inset;
                let frame = self.frame_at(center, side.z);
                let Inset::Arc {
                    radius: side_radius,
                    ..
                } = side_inset
                else {
                    return Err(ModelError::InsetCollapses { index: k });
                };
                // The tube's centre circle: the side arc moved in by the radius.
                let major = side_radius - sweep.signum() * radius;
                let up = (cap.z - side.z).signum();
                let surface: AnySurface<S> = match cap_inset {
                    Inset::Collapsed => Sphere::from_placement(frame, radius).into(),
                    _ if major > radius * (1.0 + 1e-9) => Torus::new(frame, major, radius).into(),
                    _ => return Err(ModelError::FilletTooLarge { index: k }),
                };
                // The tube angle v along a joint's quarter circle s ∈ [0, π/2]:
                // the side point is outward of the material, radially out for an
                // anticlockwise arc and in for a clockwise one.
                let (v0, dv) = match (sweep > 0.0, up > 0.0) {
                    (true, true) => (0.0, 1.0),
                    (true, false) => (0.0, -1.0),
                    (false, true) => (PI, -1.0),
                    (false, false) => (-PI, 1.0),
                };
                let ring_pcurve = |inset: &Inset, is_side: bool| match inset {
                    Inset::Arc { .. } => Some(uv_line(
                        [0.0, if is_side { v0 } else { v0 + dv * FRAC_PI_2 }],
                        [1.0, 0.0],
                    )),
                    _ => None,
                };
                let (from, to) = side_inset.params();
                Ok((
                    surface,
                    true,
                    BandPcurves {
                        lo: ring_pcurve(lo, cap_is_hi),
                        hi: ring_pcurve(hi, !cap_is_hi),
                        joints: [None, None],
                        round: Some(RoundPcurves {
                            v: [from, to],
                            u_of_s: (v0, dv),
                            vertical_is_u: true,
                        }),
                    },
                ))
            }
        }
    }
}

/// Marks a round band's joint as a circle, whose parameter maps linearly onto
/// the face (mitre ellipses don't).
#[derive(Clone, Copy)]
struct JointPcurves;

/// A round band face's joint pcurves: the joint circle's parameter s maps to
/// `offset + rate · s` on one surface axis, the other fixed per joint.
struct RoundPcurves {
    /// The fixed coordinate at joints A and B.
    v: [f64; 2],
    u_of_s: (f64, f64),
    /// Whether s maps onto the surface's second coordinate (torus/sphere tube
    /// angle) rather than its first (cylinder angle).
    vertical_is_u: bool,
}

/// Pcurves for a band face's edges.
struct BandPcurves<S: Space<3>> {
    lo: Option<Line<FaceUv<S>, 2>>,
    hi: Option<Line<FaceUv<S>, 2>>,
    /// Ruled bands' joint pcurves at A and B.
    joints: [Option<Line<FaceUv<S>, 2>>; 2],
    round: Option<RoundPcurves>,
}

impl<S: Space<3>> BandPcurves<S> {
    /// The pcurve of the joint at B (`at_b`) or A. Round joints only have one
    /// if they're circles.
    fn joint(&self, circle: &Option<JointPcurves>, at_b: bool) -> Option<Line<FaceUv<S>, 2>> {
        let index = usize::from(at_b);
        match &self.round {
            None => self.joints[index].clone(),
            Some(round) => {
                circle.as_ref()?;
                let fixed = round.v[index];
                let (offset, rate) = round.u_of_s;
                Some(match round.vertical_is_u {
                    // Torus/sphere: u fixed at the joint's angle, v = offset + rate s.
                    true => uv_line([fixed, offset], [0.0, rate]),
                    // Cylinder: u = s, v fixed.
                    false => uv_line([offset, fixed], [rate, 0.0]),
                })
            }
        }
    }
}

fn with<S: Space<3>>(coedge: Coedge<S>, pcurve: Option<Line<FaceUv<S>, 2>>) -> Coedge<S> {
    match pcurve {
        Some(p) => coedge.with_pcurve(p),
        None => coedge,
    }
}

/// A segment's unit tangent at fraction `s` along it.
fn tangent(segment: &Segment, s: f64) -> Vector2<f64> {
    match *segment {
        Segment::Line { start, end } => (end - start).normalize(),
        Segment::Arc {
            start_angle, sweep, ..
        } => {
            let angle = start_angle + sweep * s;
            Vector2::new(-angle.sin(), angle.cos()) * sweep.signum()
        }
    }
}

/// `v` turned a quarter anticlockwise: the material side of a profile.
fn left(v: Vector2<f64>) -> Vector2<f64> {
    Vector2::new(-v.y, v.x)
}
