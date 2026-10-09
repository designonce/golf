//! Cutting faces open along seams, where their loops wrap round a periodic
//! surface instead of enclosing part of it.

use golf_brep::Body;
use golf_brep::Coedge;
use golf_brep::FaceId;
use golf_brep::FaceUv;
use golf_brep::Loop;
use golf_brep::VertexId;
use golf_geom::AnyCurve;
use golf_geom::AnyCurve2;
use golf_geom::AnySurface;
use golf_geom::Circle;
use golf_geom::Line;
use golf_geom::NurbsCurve;
use golf_geom::Placement;
use golf_manifold::Domain;
use golf_manifold::Embedding;
use golf_manifold::Mapping;
use golf_manifold::Point;
use golf_manifold::Space;
use golf_manifold::Vector;
use nalgebra::Vector2;
use nalgebra::Vector3;

/// Where a loop goes in parameter space, from its pcurves: each coedge's
/// start, and how many periods it winds round each axis.
struct Winding {
    starts: Vec<(VertexId, Vector2<f64>)>,
    wraps: [i64; 2],
}

/// Why a face couldn't be cut open.
pub(crate) type SeamError = &'static str;

/// Whether any loop of `face` winds round its surface.
pub(crate) fn needs_seam<S: Space<3>>(body: &Body<S>, face: FaceId) -> bool {
    let domain = body.face(face).surface.domain();
    (0..body.face(face).loops.len())
        .filter_map(|l| winding(body, face, l, &domain))
        .any(|w| w.wraps != [0, 0])
}

/// Cuts `face` open with a seam edge (used once each way) joining the loops
/// that wrap round its surface: two loops bounding a band, or one loop and the
/// pole it surrounds. Splits an edge where the seam meets a loop between
/// vertices. Leaves the new coedges without pcurves, for healing.
pub(crate) fn insert_seam<S: Space<3>>(body: &mut Body<S>, face: FaceId) -> Result<(), SeamError> {
    let surface = body.face(face).surface.clone();
    let domain = surface.domain();
    let same_sense = body.face(face).same_sense;
    let loop_count = body.face(face).loops.len();
    let windings: Vec<Winding> = (0..loop_count)
        .map(|l| winding(body, face, l, &domain).ok_or("some of its coedges have no pcurve"))
        .collect::<Result<_, _>>()?;
    let wrapping: Vec<usize> = (0..loop_count)
        .filter(|&l| windings[l].wraps != [0, 0])
        .collect();
    // The axis wound round, which the seam crosses; it runs along the other.
    let k = match windings[wrapping[0]].wraps {
        [w, 0] if w.abs() == 1 => 0,
        [0, w] if w.abs() == 1 => 1,
        [_, 0] | [0, _] => return Err("a loop winds round more than once"),
        _ => return Err("a loop winds round both ways"),
    };
    let j = 1 - k;
    let period = domain.axes[k].max - domain.axes[k].min;
    let first = wrapping[0];
    // Walking the first loop the way it winds, the face is on its left in
    // parameter space (right, if the face is against its surface): that way
    // along the seam's axis is into the face.
    let winds = windings[first].wraps[k].signum() as f64;
    let sense = if same_sense { 1.0 } else { -1.0 };
    let inward = sense * if k == 0 { winds } else { -winds };

    let (from, c) = windings[first].starts[0];
    let c = c[k];
    let from_uv = windings[first].starts[0].1;
    let (to, to_j, second) = match wrapping.as_slice() {
        [_] => {
            // A cap: the seam runs to the pole the face surrounds.
            let axis = &domain.axes[j];
            let pole = match inward > 0.0 {
                true if axis.singular_at_max => axis.max,
                false if axis.singular_at_min => axis.min,
                _ => return Err("one loop winds round, with no pole for it to surround"),
            };
            let mut uv = Vector2::zeros();
            uv[k] = c;
            uv[j] = pole;
            let vertex = body
                .add_vertex(surface.apply(Point::new(uv)))
                .map_err(|_| "its pole couldn't be added")?;
            (vertex, pole, None)
        }
        [_, other] if windings[*other].wraps[k] == -windings[first].wraps[k] => {
            let (vertex, uv) = meet(body, face, *other, &windings[*other], k, c, period)?;
            (vertex, uv[j], Some(*other))
        }
        [_, other] if windings[*other].wraps == windings[first].wraps => {
            return Err("its two winding loops wind the same way");
        }
        [_, _] => return Err("its two winding loops wind round different axes"),
        _ => return Err("more than two of its loops wind round"),
    };

    // Along the seam's axis from the first loop into the face, a period at
    // most round a periodic one.
    let mut to_j = to_j;
    if let Some(p) = domain.axes[j].period() {
        to_j = from_uv[j] + inward * (inward * (to_j - from_uv[j])).rem_euclid(p);
    }
    if inward * (to_j - from_uv[j]) <= 0.0 {
        return Err("the seam would leave the face");
    }
    let curve = iso_curve(&surface, k, c).ok_or("its surface has no exact seam curve")?;
    let (lo, hi) = (from_uv[j].min(to_j), from_uv[j].max(to_j));
    let (start, end) = if from_uv[j] < to_j {
        (from, to)
    } else {
        (to, from)
    };
    let seam = body
        .add_edge(curve, (lo, hi), start, end)
        .map_err(|_| "the seam edge couldn't be added")?;

    // The first loop from where the seam leaves it, along the seam, the other
    // loop from where the seam meets it, and back.
    let loops = &body.face(face).loops;
    let rotated = |l: usize, at: VertexId| -> Vec<Coedge<S>> {
        let coedges = &loops[l].coedges;
        let i = coedges
            .iter()
            .position(|c| body.coedge_start(c) == at)
            .expect("a coedge starts at the seam's vertex");
        coedges[i..].iter().chain(&coedges[..i]).cloned().collect()
    };
    let mut merged = rotated(first, from);
    merged.push(Coedge::new(seam, start != from));
    if let Some(second) = second {
        merged.extend(rotated(second, to));
    }
    merged.push(Coedge::new(seam, start == from));
    let mut new_loops = vec![Loop::new(merged)];
    new_loops.extend(
        (0..loop_count)
            .filter(|l| !wrapping.contains(l))
            .map(|l| loops[l].clone()),
    );
    body.set_loops(face, new_loops)
        .map_err(|_| "the cut loop doesn't close")
}

/// The vertex of loop `l` at `c` on axis `k` (modulo `period`), splitting the
/// edge that crosses there if no vertex is; and its uv.
fn meet<S: Space<3>>(
    body: &mut Body<S>,
    face: FaceId,
    l: usize,
    winding: &Winding,
    k: usize,
    c: f64,
    period: f64,
) -> Result<(VertexId, Vector2<f64>), SeamError> {
    let near = |x: f64| {
        ((x - c) / period - ((x - c) / period).round()).abs() * period <= 1e-9 * (1.0 + c.abs())
    };
    if let Some(&(vertex, uv)) = winding.starts.iter().find(|(_, uv)| near(uv[k])) {
        return Ok((vertex, uv));
    }
    let coedges = body.face(face).loops[l].coedges.clone();
    for coedge in coedges {
        let Some(pcurve) = coedge.pcurve.clone() else {
            continue;
        };
        let range = body.edge(coedge.edge).range;
        if let Some(t) = crossing(&pcurve, range, k, c, period) {
            let uv = pcurve.apply(Point::new([t].into())).coords;
            let (vertex, _) = body
                .split_edge(coedge.edge, t)
                .map_err(|_| "an edge couldn't be split for the seam")?;
            return Ok((vertex, uv));
        }
    }
    Err("the seam meets the other loop nowhere")
}

/// A parameter strictly inside `range` where `pcurve` crosses `c` (modulo
/// `period`) on axis `k`.
fn crossing<S: Space<3>>(
    pcurve: &AnyCurve2<FaceUv<S>>,
    (a, b): (f64, f64),
    k: usize,
    c: f64,
    period: f64,
) -> Option<f64> {
    const SAMPLES: usize = 64;
    let at = |t: f64| pcurve.apply(Point::new([t].into())).coords[k];
    let t = |i: usize| a + (b - a) * i as f64 / SAMPLES as f64;
    let margin = 1e-9 * (b - a);
    let inside = |t: f64| (t > a + margin && t < b - margin).then_some(t);
    for i in 0..SAMPLES {
        let (t0, t1) = (t(i), t(i + 1));
        let (x0, x1) = (at(t0), at(t1));
        let (lo, hi) = (x0.min(x1), x0.max(x1));
        let target = c + ((lo - c) / period).ceil() * period;
        // A sample can land on it exactly.
        if x0 == target {
            return inside(t0);
        }
        if !(lo <= target && target <= hi) {
            continue;
        }
        if x1 == target {
            return inside(t1);
        }
        // Bisect for it.
        let (mut l, mut h) = (t0, t1);
        let rising = x1 > x0;
        for _ in 0..80 {
            let m = (l + h) / 2.0;
            if (at(m) < target) == rising {
                l = m;
            } else {
                h = m;
            }
        }
        return inside((l + h) / 2.0);
    }
    None
}

/// Where loop `l` goes in parameter space, if all its coedges have pcurves.
fn winding<S: Space<3>>(
    body: &Body<S>,
    face: FaceId,
    l: usize,
    domain: &Domain<2>,
) -> Option<Winding> {
    let coedges = &body.face(face).loops[l].coedges;
    let mut starts = Vec::with_capacity(coedges.len());
    let mut ends = Vec::with_capacity(coedges.len());
    for c in coedges {
        let pcurve = c.pcurve.as_ref()?;
        let range = body.edge(c.edge).range;
        let (t0, t1) = if c.reversed {
            (range.1, range.0)
        } else {
            range
        };
        let at = |t: f64| pcurve.apply(Point::new([t].into())).coords;
        starts.push((body.coedge_start(c), at(t0)));
        ends.push(at(t1));
    }
    let n = starts.len();
    let mut net = Vector2::zeros();
    for i in 0..n {
        net += ends[i] - starts[i].1;
        // The step to the next coedge, the short way round.
        let (end, next) = (ends[i], starts[(i + 1) % n].1);
        let at_pole = domain.is_singular(Point::<FaceUv<S>, 2>::new(end), 1e-9);
        for (k, axis) in domain.axes.iter().enumerate() {
            let step = next[k] - end[k];
            net[k] += match axis.period() {
                // Along a pole the loop goes as far as it goes: a whole
                // period, round a cone's apex.
                Some(_) if at_pole => step,
                Some(p) => step - (step / p).round() * p,
                None => step,
            };
        }
    }
    let wraps = [0, 1].map(|k| match domain.axes[k].period() {
        Some(p) => (net[k] / p).round() as i64,
        None => 0,
    });
    Some(Winding { starts, wraps })
}

/// The curve of `surface` along axis `1 - k` at `c` on axis `k`, parametrised
/// by that axis's parameter, where the surface has one exactly.
fn iso_curve<S: Space<3>>(surface: &AnySurface<S>, k: usize, c: f64) -> Option<AnyCurve<S>> {
    let tag = surface.apply(Point::new(Vector2::zeros())).tag();
    let point = |v: Vector3<f64>| Point::with_tag(v, tag);
    let vector = |v: Vector3<f64>| Vector::with_tag(v, tag);
    let radial = |p: &Placement<S>| p.vector(Vector3::new(c.cos(), c.sin(), 0.0)).coords;
    // A circle through its placement's x at parameter 0 and towards `y` at
    // a quarter turn.
    let circle =
        |centre: Vector3<f64>, x: Vector3<f64>, y: Vector3<f64>, radius: f64| -> AnyCurve<S> {
            Circle::new(
                Placement::from_axes(point(centre), vector(x.cross(&y)), vector(x)),
                radius,
            )
            .into()
        };
    let curve: AnyCurve<S> = match (surface, k) {
        (AnySurface::Cylinder(_) | AnySurface::Cone(_), 0) => {
            let uv = |v: f64| Point::new(Vector2::new(c, v));
            let origin = surface.apply(uv(0.0));
            let along = surface.jacobian(uv(0.0)).column(1).into_owned();
            Line::new(origin, vector(along)).into()
        }
        (AnySurface::Sphere(s), 0) => {
            let z = s.placement.vector(Vector3::z()).coords;
            circle(s.placement.origin.coords, radial(&s.placement), z, s.radius)
        }
        (AnySurface::Torus(t), 0) => {
            let r = radial(&t.placement);
            let z = t.placement.vector(Vector3::z()).coords;
            circle(
                t.placement.origin.coords + r * t.major_radius,
                r,
                z,
                t.minor_radius,
            )
        }
        (AnySurface::Torus(t), 1) => {
            let p = &t.placement;
            let centre =
                p.origin.coords + p.vector(Vector3::z()).coords * (t.minor_radius * c.sin());
            let radius = t.major_radius + t.minor_radius * c.cos();
            circle(
                centre,
                p.vector(Vector3::x()).coords,
                p.vector(Vector3::y()).coords,
                radius,
            )
        }
        (AnySurface::Nurbs(n), _) => {
            let g = n.geometry();
            let iso = if k == 0 {
                g.isocurve_u(c)
            } else {
                g.isocurve_v(c)
            };
            NurbsCurve::with_tag(iso, tag).into()
        }
        _ => return None,
    };
    // It must be the surface's own curve, at the surface's own parameter.
    let domain = surface.domain();
    let axis = domain.axes[1 - k];
    let (lo, hi) = (axis.min.max(-1e3), axis.max.min(1e3));
    let agrees = (0..=8).all(|i| {
        let v = lo + (hi - lo) * i as f64 / 8.0;
        let mut uv = Vector2::zeros();
        uv[k] = c;
        uv[1 - k] = v;
        let on_surface = surface.apply(Point::new(uv)).coords;
        let on_curve = curve.apply(Point::new([v].into())).coords;
        (on_surface - on_curve).norm() <= 1e-9 * (1.0 + on_surface.norm())
    });
    agrees.then_some(curve)
}
