use std::collections::HashMap;

use golf_manifold::Mapping;
use golf_manifold::Point;
use nalgebra::Vector3;

use crate::boundary::EdgeSample;
use crate::boundary::face_loops;
use crate::mesh::Mesh;
use crate::sample::sample;
use crate::source::MeshErrorOf;
use crate::source::MeshOf;
use crate::source::MeshSource;
use crate::tolerance::Tolerance;
use crate::triangulate::triangulate;

/// Tessellates every face of `source` to `tolerance`.
///
/// Edges are sampled once and shared by the faces either side, so a closed
/// source gives a watertight mesh.
pub fn mesh<M: MeshSource>(source: &M, tolerance: &Tolerance) -> Result<MeshOf<M>, MeshErrorOf<M>> {
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
    let mut edge_samples: HashMap<M::Edge, Vec<EdgeSample>> = HashMap::new();
    for (id, edge) in source.edges() {
        let at = |t: f64| Point::new([t].into());
        let point = |t: f64| edge.curve.apply(at(t)).coords;
        let tangent = |t: f64| -> Vector3<f64> { edge.curve.jacobian(at(t)).column(0).into() };
        // A closed edge needs a few segments to begin with, or its two ends,
        // being one point, would look like a settled chord.
        let initial = if edge.start == edge.end { 3 } else { 1 };
        let ts = sample(point, tangent, edge.range, initial, &edge_tolerance);
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

    for (face_id, face) in source.faces() {
        let loops = face_loops::<M>(face_id, &face, &edge_samples, &mut mesh, tolerance)?;
        triangulate::<M>(
            face_id,
            face.surface,
            face.same_sense,
            &loops,
            &mut mesh,
            tolerance,
        )?;
    }
    Ok(mesh)
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
