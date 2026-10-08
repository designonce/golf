use std::collections::HashMap;

use golf_brep::FaceId;
use golf_manifold::Point;
use golf_manifold::Space;
use golf_manifold::Vector;

/// A triangle mesh in space `S`.
#[derive(Clone, Debug)]
pub struct Mesh<S: Space<3>> {
    pub positions: Vec<Point<S, 3>>,
    pub triangles: Vec<Triangle<S>>,
}

/// One facet: three indices into [`Mesh::positions`], anticlockwise seen from
/// the side its face's normal points to.
#[derive(Clone, Copy, Debug)]
pub struct Triangle<S: Space<3>> {
    pub vertices: [u32; 3],
    /// The face it tessellates.
    pub face: FaceId,
    /// The face's unit normal at each corner.
    pub normals: [Vector<S, 3>; 3],
}

impl<S: Space<3>> Default for Mesh<S> {
    fn default() -> Self {
        Self {
            positions: Vec::new(),
            triangles: Vec::new(),
        }
    }
}

impl<S: Space<3>> Mesh<S> {
    /// The directed edges `[from, to]` not matched by exactly one triangle using
    /// them `[to, from]`: the mesh's holes and non-manifold seams.
    pub fn open_edges(&self) -> Vec<[u32; 2]> {
        let mut uses: HashMap<[u32; 2], usize> = HashMap::new();
        for triangle in &self.triangles {
            let [a, b, c] = triangle.vertices;
            for edge in [[a, b], [b, c], [c, a]] {
                *uses.entry(edge).or_default() += 1;
            }
        }
        let mut open: Vec<[u32; 2]> = uses
            .iter()
            .filter(|&(&[a, b], &n)| n != 1 || uses.get(&[b, a]) != Some(&1))
            .map(|(&edge, _)| edge)
            .collect();
        open.sort_unstable();
        open
    }

    /// Whether every edge is shared by exactly two triangles, running it in
    /// opposite directions.
    pub fn is_watertight(&self) -> bool {
        self.open_edges().is_empty()
    }

    /// The volume enclosed, positive when the triangles face outwards. Only
    /// meaningful for a watertight mesh.
    pub fn signed_volume(&self) -> f64 {
        self.triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.vertices.map(|i| self.positions[i as usize].coords);
                a.dot(&b.cross(&c)) / 6.0
            })
            .sum()
    }

    pub fn area(&self) -> f64 {
        self.triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.vertices.map(|i| self.positions[i as usize].coords);
                (b - a).cross(&(c - a)).norm() / 2.0
            })
            .sum()
    }
}
