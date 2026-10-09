//! Triangulating a face inside its boundary, in parameter space.

use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::VecDeque;

use golf_manifold::Domain;
use golf_manifold::Embedding;
use golf_manifold::Mapping;
use golf_manifold::Point;
use golf_manifold::Surface;
use golf_manifold::Vector;
use nalgebra::Vector2;
use nalgebra::Vector3;
use spade::ConstrainedDelaunayTriangulation;
use spade::Point2;
use spade::Triangulation;
use spade::handles::FixedFaceHandle;
use spade::handles::FixedVertexHandle;
use spade::handles::InnerTag;

use crate::boundary::BoundaryPoint;
use crate::error::MeshError;
use crate::mesh::Mesh;
use crate::mesh::Triangle;
use crate::sample::distance_to_segment;
use crate::source::MeshErrorOf;
use crate::source::MeshSource;
use crate::tolerance::Tolerance;

type Cdt = ConstrainedDelaunayTriangulation<Point2<f64>>;

/// Rounds of splitting triangles that miss the tolerance.
const MAX_ROUNDS: usize = 32;
/// A cap on interior vertices per face, against tolerances too tight to meet.
const MAX_INTERIOR_VERTICES: usize = 200_000;
/// A cap on the seeded grid's points per face.
const MAX_SEEDS: usize = 100_000;
/// The most one parameter axis is stretched against the other, so triangles
/// can be this many times longer along the surface's straighter direction.
const MAX_STRETCH: f64 = 64.0;

/// Triangulates `face` inside `loops` and appends its triangles (and any new
/// interior vertices) to `mesh`.
pub(crate) fn triangulate<M: MeshSource>(
    face_id: M::Face,
    surface: &M::Surface,
    same_sense: bool,
    loops: &[Vec<BoundaryPoint>],
    mesh: &mut Mesh<M::Space, M::Face>,
    tolerance: &Tolerance,
) -> Result<(), MeshErrorOf<M>> {
    let fail = |error: spade::InsertionError| MeshError::Triangulation {
        face: face_id,
        message: error.to_string(),
    };

    // Delaunay triangles are only well shaped if parameter distances are
    // roughly surface distances, so scale each axis by its average speed. Then
    // stretch the axis needing the shorter steps (by the curvature), so
    // triangles come out long where the surface is straight, as along a
    // cylinder: well shaped in steps rather than in millimetres.
    let speed = axis_scale::<M>(surface, loops);
    let scale = match Bounds::of(loops) {
        Some(bounds) => {
            let step = |k: usize| {
                bounds.range[k]
                    / steps_along::<M>(surface, &bounds, k, &seed_tolerance(tolerance)) as f64
                    * speed[k]
            };
            let (su, sv) = (step(0), step(1));
            let finest = su.min(sv);
            Vector2::new(
                speed.x * (finest / su).max(1.0 / MAX_STRETCH),
                speed.y * (finest / sv).max(1.0 / MAX_STRETCH),
            )
        }
        None => speed,
    };
    let to_cdt = |uv: Vector2<f64>| cdt_point(uv.component_mul(&scale));
    let to_uv = |p: Point2<f64>| Vector2::new(p.x / scale.x, p.y / scale.y);

    let mut cdt = Cdt::new();
    let mut boundary_vertex: HashMap<FixedVertexHandle, u32> = HashMap::new();
    // A mesh vertex reached more than once (where a loop touches itself, or
    // touches another) is one point, though its parameters may differ by a
    // rounding error each time; two points that close would make constraints
    // cross. (At a pole one vertex has parameters far apart: those stay.)
    let extent = {
        let uvs = loops.iter().flatten().map(|p| p.uv);
        let lo = uvs
            .clone()
            .fold(Vector2::repeat(f64::INFINITY), |a, b| a.inf(&b));
        let hi = uvs.fold(Vector2::repeat(f64::NEG_INFINITY), |a, b| a.sup(&b));
        (hi - lo).norm()
    };
    let mut inserted: HashMap<u32, Vec<(Vector2<f64>, FixedVertexHandle)>> = HashMap::new();
    for l in loops {
        let handles = l
            .iter()
            .map(|p| {
                let close = |uv: &Vector2<f64>| {
                    (uv - p.uv).norm() <= 1e-9 * (1.0 + p.uv.norm()) + 1e-6 * extent
                };
                let earlier = inserted
                    .get(&p.vertex)
                    .and_then(|seen| seen.iter().find(|(uv, _)| close(uv)).map(|&(_, h)| h));
                let handle = match earlier {
                    Some(handle) => handle,
                    None => cdt.insert(to_cdt(p.uv)).map_err(fail)?,
                };
                inserted.entry(p.vertex).or_default().push((p.uv, handle));
                boundary_vertex.entry(handle).or_insert(p.vertex);
                Ok(handle)
            })
            .collect::<Result<Vec<_>, MeshErrorOf<M>>>()?;
        for (i, &a) in handles.iter().enumerate() {
            let b = handles[(i + 1) % handles.len()];
            if a == b
                || cdt
                    .get_edge_from_neighbors(a, b)
                    .is_some_and(|e| e.is_constraint_edge())
            {
                continue;
            }
            if cdt.try_add_constraint(a, b).is_empty() {
                return Err(MeshError::SelfIntersectingBoundary { face: face_id });
            }
        }
    }

    seed_interior::<M>(&mut cdt, surface, loops, scale, tolerance).map_err(fail)?;

    let position = |cdt: &Cdt, handle: FixedVertexHandle| -> Vector3<f64> {
        match boundary_vertex.get(&handle) {
            Some(&v) => mesh.positions[v as usize].coords,
            None => {
                surface
                    .apply(Point::new(to_uv(cdt.vertex(handle).position())))
                    .coords
            }
        }
    };

    let domain = surface.domain();
    // Triangles found to meet the tolerance, by their corners, which stay
    // the same while nothing is inserted inside them; so each is checked once.
    let mut settled: HashSet<[FixedVertexHandle; 3]> = HashSet::new();
    for _ in 0..MAX_ROUNDS {
        let inside = inside_faces(&cdt);
        let mut splits = Vec::new();
        for triangle in cdt.inner_faces() {
            if !inside.contains(&triangle.fix()) {
                continue;
            }
            let handles = triangle.vertices().map(|v| v.fix());
            let mut key = handles;
            key.sort();
            if settled.contains(&key) {
                continue;
            }
            let Some(uvs) =
                lift_singular(triangle.vertices().map(|v| to_uv(v.position())), &domain)
            else {
                // Every corner on one singular bound: collapsed, and dropped.
                continue;
            };
            let corners = handles.map(|h| position(&cdt, h));
            match meets_tolerance::<M>(surface, &triangle, uvs, corners, tolerance) {
                true => {
                    settled.insert(key);
                }
                false => splits.push(to_cdt(split_point(&triangle, handles, uvs, corners))),
            }
        }
        if splits.is_empty() || cdt.num_vertices() > boundary_vertex.len() + MAX_INTERIOR_VERTICES {
            break;
        }
        for split in splits {
            cdt.insert(split).map_err(fail)?;
        }
    }

    let sense = if same_sense { 1.0 } else { -1.0 };
    let mut interior_vertex: HashMap<FixedVertexHandle, u32> = HashMap::new();
    let inside = inside_faces(&cdt);
    for triangle in cdt.inner_faces() {
        if !inside.contains(&triangle.fix()) {
            continue;
        }
        let uvs = triangle.vertices().map(|v| to_uv(v.position()));
        let mut corners: [(u32, Vector2<f64>); 3] = [(0, Vector2::zeros()); 3];
        for (corner, (v, uv)) in corners
            .iter_mut()
            .zip(triangle.vertices().into_iter().zip(uvs))
        {
            let handle = v.fix();
            let index = match boundary_vertex.get(&handle) {
                Some(&index) => index,
                None => *interior_vertex.entry(handle).or_insert_with(|| {
                    mesh.positions.push(surface.apply(Point::new(uv)));
                    (mesh.positions.len() - 1) as u32
                }),
            };
            *corner = (index, uv);
        }
        let [a, b, c] = corners.map(|(i, _)| i);
        if a == b || b == c || c == a {
            // Collapsed onto a pole.
            continue;
        }
        // Anticlockwise in uv is anticlockwise about the surface normal.
        if !same_sense {
            corners.swap(1, 2);
        }
        let [p0, p1, p2] = corners.map(|(i, _)| mesh.positions[i as usize].coords);
        let flat = (p1 - p0).cross(&(p2 - p0)).normalize();
        let normals = corners.map(|(_, uv)| {
            let n = surface.normal(Point::new(uv)).coords * sense;
            let n = if n.iter().all(|x| x.is_finite()) {
                n
            } else {
                flat
            };
            Vector::with_tag(n, mesh.positions[0].tag())
        });
        mesh.triangles.push(Triangle {
            vertices: corners.map(|(i, _)| i),
            face: face_id,
            normals,
        });
    }
    Ok(())
}

/// Where, in parameter space, to split a triangle that misses the tolerance:
/// the midpoint of its longest edge (in 3D), which keeps Delaunay triangles well
/// shaped, unless that edge is on the boundary, which neighbouring faces share;
/// then its centroid. `uvs` are the corners after [`lift_singular`].
fn split_point(
    triangle: &spade::handles::FaceHandle<'_, InnerTag, Point2<f64>, (), spade::CdtEdge<()>, ()>,
    handles: [FixedVertexHandle; 3],
    uvs: [Vector2<f64>; 3],
    corners: [Vector3<f64>; 3],
) -> Vector2<f64> {
    let corner =
        |handle: FixedVertexHandle| handles.iter().position(|&h| h == handle).expect("a corner");
    let longest = triangle
        .adjacent_edges()
        .into_iter()
        .map(|edge| (edge, corner(edge.from().fix()), corner(edge.to().fix())))
        .max_by(|a, b| {
            let length = |&(_, i, j): &(_, usize, usize)| (corners[j] - corners[i]).norm_squared();
            length(a).total_cmp(&length(b))
        })
        .expect("three edges");
    match longest {
        (edge, _, _) if edge.is_constraint_edge() => (uvs[0] + uvs[1] + uvs[2]) / 3.0,
        (_, i, j) => (uvs[i] + uvs[j]) / 2.0,
    }
}

/// The corners of a triangle with those on a singular bound moved along it to
/// the mean of the others' coordinate there, or `None` if all three are on it.
///
/// At a pole every longitude is the same point, so a corner's longitude is
/// arbitrary; the triangle's edges to it really run along the other corners'
/// meridians. Lifting makes centroids and midpoints land on those meridians,
/// rather than on whatever meridian the arbitrary longitude gives.
fn lift_singular(uvs: [Vector2<f64>; 3], domain: &Domain<2>) -> Option<[Vector2<f64>; 3]> {
    let mut lifted = uvs;
    for (k, axis) in domain.axes.iter().enumerate() {
        let along = 1 - k;
        let singular = uvs.map(|uv| axis.is_singular(uv[k], 1e-12 * (1.0 + uv[k].abs())));
        if !singular.contains(&true) {
            continue;
        }
        let free: Vec<f64> = (0..3)
            .filter(|&i| !singular[i])
            .map(|i| uvs[i][along])
            .collect();
        if free.is_empty() {
            return None;
        }
        let mean = free.iter().sum::<f64>() / free.len() as f64;
        for i in (0..3).filter(|&i| singular[i]) {
            lifted[i][along] = mean;
        }
    }
    Some(lifted)
}

/// Whether a triangle with corners at `uvs` (raised to `corners`) follows the
/// surface closely enough: the surface at its centroid is within the chord of
/// its plane, at the middle of each interior edge within the chord of that
/// edge, and its normal within the angle of the surface's at each corner.
/// Boundary edges are the edge sampler's business.
///
/// A triangle no longer than the chord is close enough whatever its normal,
/// which near a pole or a sharp boundary corner may never settle.
fn meets_tolerance<M: MeshSource>(
    surface: &M::Surface,
    triangle: &spade::handles::FaceHandle<'_, InnerTag, Point2<f64>, (), spade::CdtEdge<()>, ()>,
    uvs: [Vector2<f64>; 3],
    corners: [Vector3<f64>; 3],
    tolerance: &Tolerance,
) -> bool {
    let raise = |uv: Vector2<f64>| surface.apply(Point::new(uv)).coords;
    let normal = (corners[1] - corners[0]).cross(&(corners[2] - corners[0]));
    let longest = (0..3)
        .map(|i| (corners[(i + 1) % 3] - corners[i]).norm())
        .fold(0.0, f64::max);
    if longest <= tolerance.chord || normal.norm() <= 1e-12 * longest * longest {
        // Too small to miss, or collapsed onto a line (a pole, or a straight
        // seam meeting an apex), where its normal means nothing.
        return true;
    }
    let centroid = (uvs[0] + uvs[1] + uvs[2]) / 3.0;
    if (raise(centroid) - corners[0])
        .dot(&normal.normalize())
        .abs()
        > tolerance.chord
    {
        return false;
    }
    let edges = triangle.adjacent_edges();
    for i in 0..3 {
        let j = (i + 1) % 3;
        if edges[i].is_constraint_edge() {
            continue;
        }
        let middle = raise((uvs[i] + uvs[j]) / 2.0);
        if distance_to_segment(middle, corners[i], corners[j]) > tolerance.chord {
            return false;
        }
    }
    if let Some(max) = tolerance.max_length {
        if (0..3).any(|i| (corners[i] - corners[(i + 1) % 3]).norm() > max) {
            return false;
        }
    }
    // The facet must face the way the surface does at each corner; this also
    // catches facets folded over, which face the opposite way.
    let facet = normal.normalize();
    uvs.iter()
        .map(|&uv| surface.normal(Point::new(uv)).coords)
        .filter(|n| n.iter().all(|x| x.is_finite()))
        .all(|n| facet.angle(&n) <= tolerance.angle)
}

/// Per-axis scale making parameter distances roughly surface distances: the
/// mean length of each partial derivative over the boundary.
fn axis_scale<M: MeshSource>(surface: &M::Surface, loops: &[Vec<BoundaryPoint>]) -> Vector2<f64> {
    let mut sum = Vector2::<f64>::zeros();
    let mut count = Vector2::<f64>::zeros();
    for p in loops.iter().flatten() {
        let j = surface.jacobian(Point::new(p.uv));
        for k in 0..2 {
            let speed = j.column(k).norm();
            if speed.is_finite() && speed > 0.0 {
                sum[k] += speed;
                count[k] += 1.0;
            }
        }
    }
    sum.zip_map(&count, |s, n| if n > 0.0 { s / n } else { 1.0 })
}

/// The triangles inside the boundary, found by walking out from the convex
/// hull and flipping between outside and inside at each constraint edge.
///
/// Unlike testing a point of each triangle against the loops, this classifies
/// slivers between nearly collinear boundary points by how they connect, not by
/// which way rounding tipped them.
fn inside_faces(cdt: &Cdt) -> HashSet<FixedFaceHandle<InnerTag>> {
    let mut inside: HashMap<FixedFaceHandle<InnerTag>, bool> = HashMap::new();
    let mut queue = VecDeque::new();
    for hull_edge in cdt.convex_hull() {
        for edge in [hull_edge, hull_edge.rev()] {
            if let Some(face) = edge.face().as_inner() {
                if inside
                    .insert(face.fix(), hull_edge.is_constraint_edge())
                    .is_none()
                {
                    queue.push_back(face.fix());
                }
            }
        }
    }
    while let Some(face) = queue.pop_front() {
        let parity = inside[&face];
        for edge in cdt.face(face).adjacent_edges() {
            if let Some(next) = edge.rev().face().as_inner() {
                inside.entry(next.fix()).or_insert_with(|| {
                    queue.push_back(next.fix());
                    parity ^ edge.is_constraint_edge()
                });
            }
        }
    }
    inside
        .into_iter()
        .filter(|&(_, inside)| inside)
        .map(|(face, _)| face)
        .collect()
}

/// Inserts a grid of interior points, spaced as the surface's curvature needs,
/// so refinement starts close to the tolerance instead of working inwards from
/// the boundary. Points outside the loops, or crowding the boundary, are left
/// out.
fn seed_interior<M: MeshSource>(
    cdt: &mut Cdt,
    surface: &M::Surface,
    loops: &[Vec<BoundaryPoint>],
    scale: Vector2<f64>,
    tolerance: &Tolerance,
) -> Result<(), spade::InsertionError> {
    let Some(bounds) = Bounds::of(loops) else {
        return Ok(());
    };
    let (lo, range) = (bounds.lo, bounds.range);
    let uvs: Vec<Vector2<f64>> = loops.iter().flatten().map(|p| p.uv).collect();
    let tolerance = &seed_tolerance(tolerance);
    let steps_along = |k: usize, fixed: f64| steps_at::<M>(surface, &bounds, k, fixed, tolerance);

    // Rows of constant v, as many as the curves along v need; then each row
    // with as many points as its own curve along u needs, so rows that are
    // short on the surface (near a pole, or inside a torus) get fewer.
    let mut rows = [0.1, 0.3, 0.5, 0.7, 0.9]
        .map(|fraction| steps_along(1, lo.x + range.x * fraction))
        .into_iter()
        .max()
        .unwrap_or(1);
    let mut columns: Vec<usize> = (1..rows)
        .map(|j| steps_along(0, lo.y + range.y * j as f64 / rows as f64))
        .collect();
    let total: usize = columns.iter().sum();
    if total > MAX_SEEDS {
        let shrink = (MAX_SEEDS as f64 / total as f64).sqrt();
        rows = ((rows as f64 * shrink) as usize).max(1);
        columns = (1..rows)
            .map(|j| {
                ((steps_along(0, lo.y + range.y * j as f64 / rows as f64) as f64 * shrink) as usize)
                    .max(1)
            })
            .collect();
    }

    // Boundary points bucketed by row-sized cells (in scaled parameters), to
    // keep seeds at least half a spacing off the boundary.
    let scaled = |uv: Vector2<f64>| uv.component_mul(&scale);
    let row_spacing = range.y / rows as f64 * scale.y;
    let cell = |p: Vector2<f64>| {
        (
            (p.x / row_spacing).floor() as i64,
            (p.y / row_spacing).floor() as i64,
        )
    };
    let mut near: HashMap<(i64, i64), Vec<Vector2<f64>>> = HashMap::new();
    for &uv in &uvs {
        near.entry(cell(scaled(uv))).or_default().push(scaled(uv));
    }
    let crowded = |p: Vector2<f64>, spacing: f64| {
        let (i, j) = cell(p);
        (i - 1..=i + 1)
            .flat_map(|a| (j - 1..=j + 1).map(move |b| (a, b)))
            .any(|key| {
                near.get(&key)
                    .is_some_and(|points| points.iter().any(|q| (q - p).norm() < spacing / 2.0))
            })
    };
    let polygons: Vec<Vec<Vector2<f64>>> = loops
        .iter()
        .map(|l| l.iter().map(|p| p.uv).collect())
        .collect();
    for (j, &count) in (1..rows).zip(&columns) {
        let v = lo.y + range.y * j as f64 / rows as f64;
        let spacing = (range.x / count as f64 * scale.x).min(row_spacing);
        for i in 1..count {
            let uv = Vector2::new(lo.x + range.x * i as f64 / count as f64, v);
            if inside_even_odd(&polygons, uv) && !crowded(scaled(uv), spacing) {
                cdt.insert(cdt_point(uv.component_mul(&scale)))?;
            }
        }
    }
    Ok(())
}

/// The parameter box around a face's loops, if it has area.
struct Bounds {
    lo: Vector2<f64>,
    range: Vector2<f64>,
}

impl Bounds {
    fn of(loops: &[Vec<BoundaryPoint>]) -> Option<Self> {
        let uvs = loops.iter().flatten().map(|p| p.uv);
        let lo = uvs
            .clone()
            .fold(Vector2::repeat(f64::INFINITY), |a, b| a.inf(&b));
        let hi = uvs.fold(Vector2::repeat(f64::NEG_INFINITY), |a, b| a.sup(&b));
        let range = hi - lo;
        (range.x > 0.0 && range.y > 0.0).then_some(Self { lo, range })
    }
}

/// The tolerance seeds are spaced for: half the chord, as a cell's diagonal is
/// longer than either side, and its triangles sag about twice as much.
fn seed_tolerance(tolerance: &Tolerance) -> Tolerance {
    Tolerance {
        chord: tolerance.chord / 2.0,
        ..*tolerance
    }
}

/// The most steps [`steps_at`] takes along axis `k`, over a few lines across
/// the bounds.
fn steps_along<M: MeshSource>(
    surface: &M::Surface,
    bounds: &Bounds,
    k: usize,
    tolerance: &Tolerance,
) -> usize {
    [0.1, 0.3, 0.5, 0.7, 0.9]
        .map(|f| {
            steps_at::<M>(
                surface,
                bounds,
                k,
                bounds.lo[1 - k] + bounds.range[1 - k] * f,
                tolerance,
            )
        })
        .into_iter()
        .max()
        .unwrap_or(1)
}

/// How many steps the curve sampler takes along the isoparametric curve in
/// axis `k` at `fixed` in the other axis, across the bounds, at its finest.
fn steps_at<M: MeshSource>(
    surface: &M::Surface,
    bounds: &Bounds,
    k: usize,
    fixed: f64,
    tolerance: &Tolerance,
) -> usize {
    let (lo, range) = (bounds.lo, bounds.range);
    let other = 1 - k;
    let uv = |t: f64| {
        let mut uv = Vector2::zeros();
        uv[k] = t;
        uv[other] = fixed;
        uv
    };
    let point = |t: f64| surface.apply(Point::new(uv(t))).coords;
    let tangent = |t: f64| -> Vector3<f64> { surface.jacobian(Point::new(uv(t))).column(k).into() };
    let hi = lo[k] + range[k];
    let closed = (point(lo[k]) - point(hi)).norm() <= tolerance.chord;
    let ts = crate::sample::sample(
        point,
        tangent,
        (lo[k], hi),
        if closed { 3 } else { 1 },
        tolerance,
    );
    // The sampler halves segments, so its steps can be up to twice as fine as
    // needed. Sag grows with the square of the step, so each segment's own sag
    // says how long it could have been.
    let step = ts
        .windows(2)
        .map(|w| {
            let (a, b) = (w[0], w[1]);
            let sag = distance_to_segment(point((a + b) / 2.0), point(a), point(b));
            match sag > 0.0 {
                true => (b - a) * (tolerance.chord / sag).sqrt(),
                false => range[k],
            }
        })
        .fold(range[k], f64::min);
    (range[k] / step).ceil().max(1.0) as usize
}

/// Whether `p` is inside the polygons by the even-odd rule. Only for points
/// well clear of the boundary; triangles are classified by [`inside_faces`].
fn inside_even_odd(polygons: &[Vec<Vector2<f64>>], p: Vector2<f64>) -> bool {
    let mut inside = false;
    for polygon in polygons {
        for i in 0..polygon.len() {
            let (a, b) = (polygon[i], polygon[(i + 1) % polygon.len()]);
            if (a.y > p.y) != (b.y > p.y) {
                let x = a.x + (p.y - a.y) / (b.y - a.y) * (b.x - a.x);
                if p.x < x {
                    inside = !inside;
                }
            }
        }
    }
    inside
}

/// A point for the triangulation, with negligible coordinates made exactly
/// zero: rotations leave values like 1e-50 where zero belongs, and spade
/// rejects anything that small that isn't zero.
fn cdt_point(p: Vector2<f64>) -> Point2<f64> {
    let snap = |x: f64| if x.abs() < 1e-30 { 0.0 } else { x };
    Point2::new(snap(p.x), snap(p.y))
}
