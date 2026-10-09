//! A face's boundary in its surface's parameter space.

use std::collections::HashMap;

use golf_manifold::Domain;
use golf_manifold::Embedding;
use golf_manifold::Mapping;
use golf_manifold::Point;
use nalgebra::Vector2;
use nalgebra::Vector3;

use crate::error::MeshError;
use crate::mesh::Mesh;
use crate::sample::sample;
use crate::source::FaceData;
use crate::source::MeshErrorOf;
use crate::source::MeshSource;
use crate::tolerance::Tolerance;

/// A boundary point: where it is in the face's parameter space, and which mesh
/// vertex it raises to.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct BoundaryPoint {
    pub uv: Vector2<f64>,
    pub vertex: u32,
}

/// One sample of an edge: its curve parameter and mesh vertex.
#[derive(Clone, Copy, Debug)]
pub(crate) struct EdgeSample {
    pub t: f64,
    pub vertex: u32,
}

/// The face's loops as closed polygons in parameter space (the last point
/// joins the first), each point tied to a mesh vertex.
pub(crate) fn face_loops<M: MeshSource>(
    face_id: M::Face,
    face: &FaceData<'_, M>,
    edge_samples: &HashMap<M::Edge, Vec<EdgeSample>>,
    mesh: &mut Mesh<M::Space, M::Face>,
    tolerance: &Tolerance,
) -> Result<Vec<Vec<BoundaryPoint>>, MeshErrorOf<M>> {
    let domain = face.surface.domain();
    if face.loops.is_empty() {
        return domain_rectangle::<M>(face_id, face.surface, &domain, mesh, tolerance)
            .map(|l| vec![l]);
    }
    let mut loops = Vec::with_capacity(face.loops.len());
    for (loop_index, l) in face.loops.iter().enumerate() {
        let mut points: Vec<BoundaryPoint> = Vec::new();
        for coedge in l {
            let samples = &edge_samples[&coedge.edge];
            let ordered: Vec<EdgeSample> = match coedge.reversed {
                true => samples.iter().rev().copied().collect(),
                false => samples.clone(),
            };
            // Each coedge's last sample is the next one's first vertex, but not
            // necessarily the same parameters: at a pole the two differ. Exact
            // repeats are dropped below.
            for s in &ordered {
                let uv = match &coedge.pcurve {
                    Some(pcurve) => pcurve.apply(Point::new([s.t].into())).coords,
                    None => {
                        let hint = points.last().map(|p| Point::new(p.uv));
                        let p = mesh.positions[s.vertex as usize];
                        let q = face.surface.project(p, hint).map_err(|source| {
                            MeshError::Projection {
                                face: face_id,
                                edge: coedge.edge,
                                source,
                            }
                        })?;
                        match hint {
                            Some(hint) => domain.unwrap_near(q, hint).coords,
                            None => q.coords,
                        }
                    }
                };
                points.push(BoundaryPoint {
                    uv,
                    vertex: s.vertex,
                });
            }
        }
        // The same vertex reached through a pcurve and through projection can
        // land a rounding error apart; those are one point. (At a pole one
        // vertex has genuinely different parameters, far apart.)
        let same = |a: &BoundaryPoint, b: &BoundaryPoint| {
            a.uv == b.uv
                || (a.vertex == b.vertex && (a.uv - b.uv).norm() <= 1e-9 * (1.0 + a.uv.norm()))
        };
        points.dedup_by(|b, a| same(a, b));
        if points.len() > 1 && same(&points[0], &points[points.len() - 1]) {
            points.pop();
        }
        let points = bridge_singular_bounds(&points, &domain);
        if wraps(&points, &domain) {
            return Err(MeshError::LoopWrapsSurface {
                face: face_id,
                loop_index,
            });
        }
        loops.push(points);
    }
    Ok(loops)
}

/// Whether a closed polygon winds round a periodic axis rather than enclosing
/// a region: its steps, each taken the short way round, don't sum to zero.
fn wraps(points: &[BoundaryPoint], domain: &Domain<2>) -> bool {
    domain.axes.iter().enumerate().any(|(k, axis)| {
        let Some(period) = axis.period() else {
            return false;
        };
        let total: f64 = (0..points.len())
            .map(|i| {
                let (a, b) = (points[i].uv[k], points[(i + 1) % points.len()].uv[k]);
                axis.unwrap_near(b, a) - a
            })
            .sum();
        total.abs() > period / 2.0
    })
}

/// Fills in the stretches where a loop runs along a singular bound between two
/// points (crossing a pole, say), which no edge covers. The added points all
/// raise to the same vertex as the pole itself.
fn bridge_singular_bounds(points: &[BoundaryPoint], domain: &Domain<2>) -> Vec<BoundaryPoint> {
    let mut out = Vec::with_capacity(points.len());
    for i in 0..points.len() {
        let (a, b) = (points[i], points[(i + 1) % points.len()]);
        out.push(a);
        for (k, axis) in domain.axes.iter().enumerate() {
            let along = 1 - k;
            let tolerance = 1e-9 * (1.0 + a.uv[k].abs());
            let on_same_bound =
                axis.is_singular(a.uv[k], tolerance) && (a.uv[k] - b.uv[k]).abs() <= tolerance;
            let gap = b.uv[along] - a.uv[along];
            if !on_same_bound || gap.abs() <= tolerance {
                continue;
            }
            let span = domain.axes[along]
                .period()
                .unwrap_or(domain.axes[along].max - domain.axes[along].min);
            let steps = ((gap.abs() / span * 16.0).ceil() as usize).max(1);
            for s in 1..steps {
                let mut uv = a.uv;
                uv[along] += gap * s as f64 / steps as f64;
                out.push(BoundaryPoint {
                    uv,
                    vertex: a.vertex,
                });
            }
        }
    }
    out
}

/// The boundary of a face with no loops: its whole domain, with the sides a
/// periodic axis joins sharing vertices, and each singular side one vertex.
fn domain_rectangle<M: MeshSource>(
    face_id: M::Face,
    surface: &M::Surface,
    domain: &Domain<2>,
    mesh: &mut Mesh<M::Space, M::Face>,
    tolerance: &Tolerance,
) -> Result<Vec<BoundaryPoint>, MeshErrorOf<M>> {
    let [u_axis, v_axis] = domain.axes;
    let bounds = [u_axis.min, u_axis.max, v_axis.min, v_axis.max];
    if bounds.iter().any(|b| !b.is_finite()) {
        return Err(MeshError::UnboundedFace { face: face_id });
    }
    let (u0, u1, v0, v1) = (u_axis.min, u_axis.max, v_axis.min, v_axis.max);

    // Samples along an isoparametric side: u varying at fixed v, or the reverse.
    let side = |fixed: f64, along_u: bool| {
        let uv = move |t: f64| match along_u {
            true => Vector2::new(t, fixed),
            false => Vector2::new(fixed, t),
        };
        let direction = match along_u {
            true => Vector2::x(),
            false => Vector2::y(),
        };
        let point = |t: f64| surface.apply(Point::new(uv(t))).coords;
        let tangent = |t: f64| -> Vector3<f64> { surface.jacobian(Point::new(uv(t))) * direction };
        let range = if along_u { (u0, u1) } else { (v0, v1) };
        sample(point, tangent, range, 4, tolerance)
    };
    let bottom = side(v0, true);
    let top = if v_axis.periodic {
        bottom.clone()
    } else {
        side(v1, true)
    };
    let left = side(u0, false);
    let right = if u_axis.periodic {
        left.clone()
    } else {
        side(u1, false)
    };

    // Points that are the same place on the surface share a vertex: identify
    // the ends of periodic axes, and collapse singular sides.
    let mut vertices: HashMap<(u64, u64), u32> = HashMap::new();
    let mut point = |uv: Vector2<f64>| {
        let mut key = uv;
        for k in 0..2 {
            key[k] = domain.axes[k].wrap(key[k]);
            if domain.axes[k].is_singular(uv[k], 0.0) {
                key[1 - k] = domain.axes[1 - k].min;
            }
        }
        let vertex = *vertices
            .entry((key.x.to_bits(), key.y.to_bits()))
            .or_insert_with(|| {
                mesh.positions.push(surface.apply(Point::new(uv)));
                (mesh.positions.len() - 1) as u32
            });
        BoundaryPoint { uv, vertex }
    };

    let mut loop_points = Vec::new();
    for &u in &bottom[..bottom.len() - 1] {
        loop_points.push(point(Vector2::new(u, v0)));
    }
    for &v in &right[..right.len() - 1] {
        loop_points.push(point(Vector2::new(u1, v)));
    }
    for &u in top[1..].iter().rev() {
        loop_points.push(point(Vector2::new(u, v1)));
    }
    for &v in left[1..].iter().rev() {
        loop_points.push(point(Vector2::new(u0, v)));
    }
    Ok(loop_points)
}
