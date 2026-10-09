use bevy::asset::RenderAssetUsages;
use bevy::mesh::PrimitiveTopology;
use bevy::prelude::*;
use golf_color::Color as GolfColor;
use golf_manifold::Space;

/// Meshes to show, each face coloured.
#[derive(Clone, Debug, Default)]
pub struct Scene {
    meshes: Vec<Shown>,
}

/// A mesh ready for Bevy: one vertex per triangle corner, so each keeps its
/// face's normal and colour.
#[derive(Clone, Debug, Default)]
struct Shown {
    positions: Vec<[f32; 3]>,
    normals: Vec<[f32; 3]>,
    colors: Vec<[f32; 4]>,
    transparent: bool,
}

impl Scene {
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a mesh, each triangle coloured by `color` of its face.
    pub fn add_mesh<S: Space<3>, F>(
        &mut self,
        mesh: &golf_mesh::Mesh<S, F>,
        color: impl Fn(&F) -> GolfColor,
    ) -> &mut Self {
        let mut shown = Shown::default();
        for triangle in &mesh.triangles {
            let face_color = color(&triangle.face);
            let [r, g, b] = face_color.unit_rgb().map(|c| c as f32);
            let rgba = Color::srgba(r, g, b, face_color.opacity() as f32);
            shown.transparent |= rgba.alpha() < 1.0;
            let linear = rgba.to_linear().to_f32_array();
            for (&vertex, normal) in triangle.vertices.iter().zip(&triangle.normals) {
                let p = mesh.positions[vertex as usize].coords;
                shown.positions.push([p.x as f32, p.y as f32, p.z as f32]);
                let n = normal.coords;
                shown.normals.push([n.x as f32, n.y as f32, n.z as f32]);
                shown.colors.push(linear);
            }
        }
        self.meshes.push(shown);
        self
    }

    pub fn is_empty(&self) -> bool {
        self.meshes.iter().all(|m| m.positions.is_empty())
    }

    /// The smallest box holding every mesh, as its least and greatest corners.
    pub fn bounds(&self) -> Option<(Vec3, Vec3)> {
        let mut points = self
            .meshes
            .iter()
            .flat_map(|m| &m.positions)
            .map(|&p| Vec3::from(p));
        let first = points.next()?;
        Some(points.fold((first, first), |(lo, hi), p| (lo.min(p), hi.max(p))))
    }

    /// Each mesh as a Bevy mesh and the material to show it with.
    pub(crate) fn bevy_meshes(&self) -> impl Iterator<Item = (Mesh, StandardMaterial)> + '_ {
        self.meshes.iter().map(|shown| {
            let mesh = Mesh::new(
                PrimitiveTopology::TriangleList,
                RenderAssetUsages::RENDER_WORLD,
            )
            .with_inserted_attribute(Mesh::ATTRIBUTE_POSITION, shown.positions.clone())
            .with_inserted_attribute(Mesh::ATTRIBUTE_NORMAL, shown.normals.clone())
            .with_inserted_attribute(Mesh::ATTRIBUTE_COLOR, shown.colors.clone());
            // White, so the vertex colours show as they are.
            let material = StandardMaterial {
                base_color: Color::WHITE,
                perceptual_roughness: 0.6,
                alpha_mode: match shown.transparent {
                    true => AlphaMode::Blend,
                    false => AlphaMode::Opaque,
                },
                double_sided: true,
                cull_mode: None,
                ..default()
            };
            (mesh, material)
        })
    }
}

#[cfg(test)]
mod tests {
    use golf_geom::Placement;
    use golf_manifold::Point;
    use golf_manifold::World;
    use golf_mesh::Tolerance;
    use golf_model::primitives;
    use nalgebra::Vector3;

    use super::*;

    #[test]
    fn meshes_become_coloured_corners() {
        let origin = Placement::<World>::at(Point::new(Vector3::zeros()));
        let cuboid = primitives::cuboid(&origin, Vector3::new(1.0, 2.0, 3.0)).unwrap();
        let mesh = golf_mesh::mesh(&cuboid, &Tolerance::new(1e-3, 0.5)).unwrap();
        let mut scene = Scene::new();
        scene.add_mesh(&mesh, |_| GolfColor::rgba(255, 0, 0, 128));
        let shown = &scene.meshes[0];
        assert_eq!(shown.positions.len(), 3 * mesh.triangles.len());
        assert!(shown.transparent);
        assert_eq!(shown.colors[0][..3], [1.0, 0.0, 0.0]);
        let (lo, hi) = scene.bounds().unwrap();
        assert!((hi - lo - Vec3::new(1.0, 2.0, 3.0)).length() < 1e-6);
    }
}
