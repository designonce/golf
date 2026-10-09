use std::collections::HashMap;

use golf_manifold::Embedding;
use golf_manifold::Mapping;
use golf_manifold::Point;
use golf_manifold::Surface;
use nalgebra::Vector3;

use crate::boundary::EdgeSample;
use crate::boundary::face_loops;
use crate::error::MeshError;
use crate::mesh::Mesh;
use crate::sample::sample_with;
use crate::source::MeshErrorOf;
use crate::source::MeshOf;
use crate::source::MeshSource;
use crate::tolerance::Tolerance;
use crate::triangulate::triangulate;

/// Tessellates every face of `source` to `tolerance`.
///
/// Edges are sampled once and shared by the faces either side, so a closed
/// source gives a watertight mesh. Where a face's boundary, so sampled,
/// crosses itself (a hole close to the outline, whose chords cut across it),
/// its edges are sampled more finely and the source meshed again.
pub fn mesh<M: MeshSource>(source: &M, tolerance: &Tolerance) -> Result<MeshOf<M>, MeshErrorOf<M>> {
    const RETRIES: usize = 4;
    let mut finer: HashMap<M::Edge, f64> = HashMap::new();
    let mut attempt = 0;
    loop {
        match mesh_with(source, tolerance, &finer) {
            Err(MeshError::SelfIntersectingBoundary { face }) if attempt < RETRIES => {
                attempt += 1;
                let edges = source
                    .faces()
                    .find(|(id, _)| *id == face)
                    .map(|(_, f)| f.loops.iter().flatten().map(|c| c.edge).collect::<Vec<_>>())
                    .unwrap_or_default();
                for edge in edges {
                    *finer.entry(edge).or_insert(1.0) /= 4.0;
                }
            }
            result => return result,
        }
    }
}

/// [`mesh`], with some edges' chord tolerance scaled down.
fn mesh_with<M: MeshSource>(
    source: &M,
    tolerance: &Tolerance,
    finer: &HashMap<M::Edge, f64>,
) -> Result<MeshOf<M>, MeshErrorOf<M>> {
    let mut mesh = Mesh::default();
    let vertex_index: HashMap<M::Vertex, u32> = source
        .vertices()
        .map(|(id, point)| {
            mesh.positions.push(point);
            (id, (mesh.positions.len() - 1) as u32)
        })
        .collect();

    // Half the chord, so the faces either side can meet the tolerance along
    // their boundary: an edge's sag shows up in their triangles too.
    let edge_tolerance = Tolerance {
        chord: tolerance.chord / 2.0,
        ..*tolerance
    };
    // The faces either side of each edge, through its pcurves: their triangles
    // along the edge are flat, so the edge must be sampled finely enough for
    // the surface normals to turn no more than the angle between samples.
    let faces: Vec<_> = source.faces().collect();
    let mut sides: HashMap<M::Edge, Vec<Side<'_, M>>> = HashMap::new();
    for (_, face) in &faces {
        let sense = if face.same_sense { 1.0 } else { -1.0 };
        for coedge in face.loops.iter().flatten() {
            if let Some(pcurve) = coedge.pcurve {
                sides.entry(coedge.edge).or_default().push(Side {
                    surface: face.surface,
                    pcurve,
                    sense,
                });
            }
        }
    }
    let mut edge_samples: HashMap<M::Edge, Vec<EdgeSample>> = HashMap::new();
    for (id, edge) in source.edges() {
        let sides = sides.get(&id).map_or(&[][..], Vec::as_slice);
        // A normal at (or a rounding error past) a pole means nothing.
        let normal = |side: &Side<'_, M>, t: f64| {
            let uv = side.pcurve.apply(Point::new([t].into()));
            let at_pole = side.surface.domain().is_singular(uv, 1e-9);
            match at_pole {
                true => Vector3::repeat(f64::NAN),
                false => side.surface.normal(uv).coords * side.sense,
            }
        };
        let flat_enough = |a: f64, b: f64| {
            sides.iter().all(|side| {
                let (na, nb) = (normal(side, a), normal(side, b));
                // A pole's normal is undefined; the edge sampler can't help there.
                !(na.iter().chain(&nb).all(|x| x.is_finite())) || na.angle(&nb) <= tolerance.angle
            })
        };
        let at = |t: f64| Point::new([t].into());
        let point = |t: f64| edge.curve.apply(at(t)).coords;
        let tangent = |t: f64| -> Vector3<f64> { edge.curve.jacobian(at(t)).column(0).into() };
        // A closed edge needs a few segments to begin with, or its two ends,
        // being one point, would look like a settled chord.
        let initial = if edge.start == edge.end { 3 } else { 1 };
        let ts = sample_with(
            point,
            tangent,
            edge.range,
            initial,
            &Tolerance {
                chord: edge_tolerance.chord * finer.get(&id).copied().unwrap_or(1.0),
                ..edge_tolerance
            },
            flat_enough,
        );
        let last = ts.len() - 1;
        let samples = ts
            .iter()
            .enumerate()
            .map(|(i, &t)| {
                let vertex = match i {
                    0 => vertex_index[&edge.start],
                    i if i == last => vertex_index[&edge.end],
                    _ => {
                        mesh.positions.push(edge.curve.apply(at(t)));
                        (mesh.positions.len() - 1) as u32
                    }
                };
                EdgeSample { t, vertex }
            })
            .collect();
        edge_samples.insert(id, samples);
    }

    for (face_id, face) in faces {
        let loops = face_loops::<M>(face_id, &face, &edge_samples, &mut mesh, tolerance)?;
        // spade can panic on nearly degenerate input; that face then fails
        // alone, as a triangulation error, rather than taking the caller down.
        // The mesh is only appended to, so a failed face leaves no partial
        // state that matters: the error is returned.
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            triangulate::<M>(
                face_id,
                face.surface,
                face.same_sense,
                &loops,
                &mut mesh,
                tolerance,
            )
        }));
        match result {
            Ok(done) => done?,
            Err(panic) => {
                let message = panic
                    .downcast_ref::<&str>()
                    .map(|s| s.to_string())
                    .or_else(|| panic.downcast_ref::<String>().cloned())
                    .unwrap_or_else(|| "the triangulator panicked".to_string());
                return Err(MeshError::Triangulation {
                    face: face_id,
                    message,
                });
            }
        }
    }
    Ok(mesh)
}

/// A face using an edge, as the edge sampler sees it.
struct Side<'a, M: MeshSource + ?Sized> {
    surface: &'a M::Surface,
    pcurve: &'a M::Pcurve,
    /// The face normal's sign against the surface's.
    sense: f64,
}

#[cfg(test)]
mod tests {
    use core::f64::consts::PI;

    use golf_geom::Line;
    use golf_geom::Placement;
    use golf_geom::Sphere;
    use golf_geom::Torus;
    use golf_manifold::Embedding;
    use golf_manifold::Space;
    use golf_manifold::Uv;
    use golf_manifold::World;

    use super::*;
    use crate::source::EdgeData;
    use crate::source::FaceData;

    /// The smallest source: one face covering a whole closed surface. Not a
    /// B-rep at all, which is the point.
    struct Whole<G>(G);

    impl<G> MeshSource for Whole<G>
    where
        G: Embedding<2, 3, To = World, From = Uv<G>>,
        Uv<G>: Space<2, Tag = ()>,
    {
        type Space = World;
        type Curve = Line<World, 3>;
        type Surface = G;
        type Pcurve = Line<Uv<G>, 2>;
        type Vertex = ();
        type Edge = ();
        type Face = ();

        fn vertices(&self) -> impl Iterator<Item = ((), Point<World, 3>)> + '_ {
            core::iter::empty()
        }

        fn edges(&self) -> impl Iterator<Item = ((), EdgeData<'_, Self>)> + '_ {
            core::iter::empty()
        }

        fn faces(&self) -> impl Iterator<Item = ((), FaceData<'_, Self>)> + '_ {
            core::iter::once((
                (),
                FaceData {
                    surface: &self.0,
                    same_sense: true,
                    loops: Vec::new(),
                },
            ))
        }
    }

    fn tolerance() -> Tolerance {
        Tolerance::new(1e-3, 0.3)
    }

    #[test]
    fn whole_sphere_is_watertight() {
        let sphere = Sphere::new(Point::new(Vector3::new(1.0, 2.0, 3.0)), 2.0);
        let mesh = mesh(&Whole(sphere), &tolerance()).unwrap();
        assert!(mesh.is_watertight());
        let exact = 4.0 / 3.0 * PI * 8.0;
        assert!((mesh.signed_volume() - exact).abs() / exact < 3e-3);
    }

    #[test]
    fn whole_torus_is_watertight() {
        let torus = Torus::new(Placement::at(Point::new(Vector3::zeros())), 3.0, 1.0);
        let mesh = mesh(&Whole(torus), &tolerance()).unwrap();
        assert!(mesh.is_watertight());
        let exact = 2.0 * PI * PI * 3.0;
        assert!((mesh.signed_volume() - exact).abs() / exact < 3e-3);
    }
}
