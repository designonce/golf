use golf_brep::Body;
use golf_brep::Coedge;
use golf_brep::EdgeId;
use golf_brep::Loop;
use golf_brep::VertexId;
use golf_geom::Line;
use golf_geom::NurbsCurve;
use golf_geom::NurbsSurface;
use golf_geom::Placement;
use golf_geom::Plane;
use golf_manifold::Space;
use golf_manifold::Vector;
use golf_sketch::Profile;
use golf_sketch::Region;
use golf_sketch::Segment;
use nalgebra::Vector2;
use nalgebra::Vector3;

use crate::error::ModelError;
use crate::extrude::uv_line;

/// Lofts a ruled solid between two regions, each drawn in its own sketch
/// plane.
///
/// The regions must match: the same number of holes, and each profile the same
/// number of segments as its counterpart. Segment `i` of one profile is joined
/// to segment `i` of the other by a ruled surface (an exact rational NURBS, so
/// arcs stay arcs), the matching vertices by straight edges, and each region
/// closes its end as a planar cap. Both sketch planes must face the same way
/// along the loft (from one region towards the other, or both against it).
pub fn loft<S: Space<3>>(
    (bottom_plane, bottom): (&Placement<S>, &Region),
    (top_plane, top): (&Placement<S>, &Region),
) -> Result<Body<S>, ModelError> {
    let bottom_profiles: Vec<&Profile> = bottom.profiles().collect();
    let top_profiles: Vec<&Profile> = top.profiles().collect();
    let mismatched = bottom_profiles.len() != top_profiles.len()
        || bottom_profiles
            .iter()
            .zip(&top_profiles)
            .any(|(b, t)| b.segments().len() != t.segments().len());
    if mismatched {
        return Err(ModelError::LoftMismatch);
    }

    // Which way the loft runs, against each plane's normal.
    let centre = |plane: &Placement<S>, region: &Region| {
        let points: Vec<Vector2<f64>> = region
            .outer()
            .segments()
            .iter()
            .map(Segment::start)
            .collect();
        let mean = points.iter().sum::<Vector2<f64>>() / points.len() as f64;
        plane.point(Vector3::new(mean.x, mean.y, 0.0)).coords
    };
    let direction = centre(top_plane, top) - centre(bottom_plane, bottom);
    let normal = |plane: &Placement<S>| plane.vector(Vector3::z()).coords;
    let (along_bottom, along_top) = (
        direction.dot(&normal(bottom_plane)),
        direction.dot(&normal(top_plane)),
    );
    if along_bottom == 0.0 || along_bottom.signum() != along_top.signum() {
        return Err(ModelError::LoftFacing);
    }
    // When the profiles (anticlockwise about their planes' normals) run
    // anticlockwise about the loft direction, the ruled surfaces' normals point
    // out of the solid.
    let outward = along_bottom > 0.0;

    let tag = bottom_plane.origin.tag();
    let mut body = Body::with_tag(tag);
    let mut sides = Vec::new();
    let mut bottom_loops = Vec::new();
    let mut top_loops = Vec::new();
    for (b_profile, t_profile) in bottom_profiles.iter().zip(&top_profiles) {
        let (b_segments, t_segments) = (b_profile.segments(), t_profile.segments());
        let n = b_segments.len();
        let at = |plane: &Placement<S>, p: Vector2<f64>| plane.point(Vector3::new(p.x, p.y, 0.0));
        let b_vertices: Vec<VertexId> = b_segments
            .iter()
            .map(|s| body.add_vertex(at(bottom_plane, s.start())))
            .collect::<Result<_, _>>()?;
        let t_vertices: Vec<VertexId> = t_segments
            .iter()
            .map(|s| body.add_vertex(at(top_plane, s.start())))
            .collect::<Result<_, _>>()?;
        let rails: Vec<EdgeId> = (0..n)
            .map(|k| {
                let from = at(bottom_plane, b_segments[k].start());
                let to = at(top_plane, t_segments[k].start());
                body.add_edge(
                    Line::new(from, Vector::with_tag(to.coords - from.coords, tag)),
                    (0.0, 1.0),
                    b_vertices[k],
                    t_vertices[k],
                )
            })
            .collect::<Result<_, _>>()?;

        let mut b_coedges = Vec::new();
        let mut t_coedges = Vec::new();
        for k in 0..n {
            let next = (k + 1) % n;
            let (b_curve, t_curve) = compatible(
                segment_nurbs(bottom_plane, &b_segments[k]),
                segment_nurbs(top_plane, &t_segments[k]),
            )?;
            let b_edge = body.add_edge(
                NurbsCurve::with_tag(b_curve.clone(), tag),
                (0.0, 1.0),
                b_vertices[k],
                b_vertices[next],
            )?;
            let t_edge = body.add_edge(
                NurbsCurve::with_tag(t_curve.clone(), tag),
                (0.0, 1.0),
                t_vertices[k],
                t_vertices[next],
            )?;

            // The ruled face, u along the profiles and v from bottom to top: the
            // bottom edge along v = 0, up the next rail at u = 1, back along the
            // top at v = 1, down this rail at u = 0.
            let surface = ruled(&b_curve, &t_curve)?;
            let mut side = vec![
                Coedge::new(b_edge, false).with_pcurve(uv_line([0.0, 0.0], [1.0, 0.0])),
                Coedge::new(rails[next], false).with_pcurve(uv_line([1.0, 0.0], [0.0, 1.0])),
                Coedge::new(t_edge, true).with_pcurve(uv_line([0.0, 1.0], [1.0, 0.0])),
                Coedge::new(rails[k], true).with_pcurve(uv_line([0.0, 0.0], [0.0, 1.0])),
            ];
            if !outward {
                side = reversed(side);
            }
            sides.push(body.add_face(
                NurbsSurface::with_tag(surface, tag),
                outward,
                vec![Loop::new(side)],
            )?);
            b_coedges.push(Coedge::new(b_edge, false));
            t_coedges.push(Coedge::new(t_edge, false));
        }
        // Each cap runs its edges the other way from the sides.
        bottom_loops.push(Loop::new(match outward {
            true => reversed(b_coedges),
            false => b_coedges,
        }));
        top_loops.push(Loop::new(match outward {
            true => t_coedges,
            false => reversed(t_coedges),
        }));
    }
    // The bottom cap faces away from the loft, the top along it.
    let mut faces = vec![
        body.add_face(Plane::new(*bottom_plane), !outward, bottom_loops)?,
        body.add_face(Plane::new(*top_plane), outward, top_loops)?,
    ];
    faces.extend(sides);
    body.add_shell(faces)?;
    Ok(body)
}

/// A loop of coedges traversed the other way.
fn reversed<S: Space<3>>(coedges: Vec<Coedge<S>>) -> Vec<Coedge<S>> {
    coedges
        .into_iter()
        .rev()
        .map(|mut c| {
            c.reversed = !c.reversed;
            c
        })
        .collect()
}

/// A segment drawn in `plane`, as a NURBS on [0, 1] running the segment's way:
/// a line as a straight Bézier, an arc as an exact rational quadratic.
fn segment_nurbs<S: Space<3>>(
    plane: &Placement<S>,
    segment: &Segment,
) -> golf_nurbs::NurbsCurve<3> {
    let at = |p: Vector2<f64>| plane.point(Vector3::new(p.x, p.y, 0.0)).coords;
    match *segment {
        Segment::Line { start, end } => {
            golf_nurbs::NurbsCurve::bezier(vec![at(start), at(end)]).expect("two points")
        }
        Segment::Arc {
            center,
            radius,
            start_angle,
            sweep,
        } => {
            let (x, y) = (
                plane.vector(Vector3::x()).coords,
                plane.vector(Vector3::y()).coords,
            );
            let anticlockwise = |from: f64, to: f64| {
                golf_nurbs::NurbsCurve::circular_arc(at(center), x, y, radius, from, to)
            };
            match sweep > 0.0 {
                true => anticlockwise(start_angle, start_angle + sweep),
                false => anticlockwise(start_angle + sweep, start_angle).reversed(),
            }
        }
    }
}

/// The two curves at a common degree and on common knots, ready to be rows
/// of one control grid.
fn compatible(
    a: golf_nurbs::NurbsCurve<3>,
    b: golf_nurbs::NurbsCurve<3>,
) -> Result<(golf_nurbs::NurbsCurve<3>, golf_nurbs::NurbsCurve<3>), ModelError> {
    let degree = a.degree().max(b.degree());
    let raise = |c: golf_nurbs::NurbsCurve<3>| {
        let by = degree - c.degree();
        c.elevate_degree(by)
    };
    let (mut a, mut b) = (raise(a).map_err(nurbs)?, raise(b).map_err(nurbs)?);
    // Each interior knot at the larger of its two multiplicities.
    let interior = |c: &golf_nurbs::NurbsCurve<3>| -> Vec<(f64, usize)> {
        let (lo, hi) = c.domain();
        c.knots()
            .distinct()
            .filter(|&(k, _)| lo < k && k < hi)
            .collect()
    };
    for (knot, multiplicity) in interior(&a).into_iter().chain(interior(&b)) {
        for curve in [&mut a, &mut b] {
            let missing = multiplicity.saturating_sub(curve.knots().multiplicity(knot));
            if missing > 0 {
                *curve = curve.insert_knot(knot, missing).map_err(nurbs)?;
            }
        }
    }
    Ok((a, b))
}

/// The ruled surface joining two compatible curves: u along them, v from `a`
/// (v = 0) to `b` (v = 1).
fn ruled(
    a: &golf_nurbs::NurbsCurve<3>,
    b: &golf_nurbs::NurbsCurve<3>,
) -> Result<golf_nurbs::NurbsSurface<3>, ModelError> {
    let count = a.control_point_count();
    golf_nurbs::NurbsSurface::new(
        (a.degree(), a.knots().knots().to_vec()),
        (1, vec![0.0, 0.0, 1.0, 1.0]),
        (0..count)
            .map(|i| vec![a.control_point(i), b.control_point(i)])
            .collect(),
        Some((0..count).map(|i| vec![a.weight(i), b.weight(i)]).collect()),
    )
    .map_err(nurbs)
}

fn nurbs(error: golf_nurbs::NurbsError) -> ModelError {
    ModelError::Nurbs(error.to_string())
}
