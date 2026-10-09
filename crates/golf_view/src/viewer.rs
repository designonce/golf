use bevy::asset::RenderAssetUsages;
use bevy::mesh::PrimitiveTopology;
use bevy::prelude::*;
use golf_color::Color as GolfColor;
use golf_manifold::Space;

use crate::orbit::Orbit;
use crate::orbit::orbit;

/// Meshes to show, gathered before opening the window.
#[derive(Clone, Debug, Default)]
pub struct Viewer {
    title: String,
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

#[derive(Resource)]
struct Scene(Vec<Shown>);

impl Viewer {
    pub fn new(title: impl Into<String>) -> Self {
        Self {
            title: title.into(),
            meshes: Vec::new(),
        }
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

    /// Opens a window showing the meshes, until it's closed.
    pub fn run(self) -> AppExit {
        let window = Window {
            title: self.title.clone(),
            ..default()
        };
        let (lo, hi) = self
            .bounds()
            .unwrap_or((Vec3::splat(-1.0), Vec3::splat(1.0)));
        let camera = Orbit {
            target: (lo + hi) / 2.0,
            yaw: -1.0,
            pitch: 0.5,
            size: (hi - lo).length().max(1e-3) * 1.1,
            // Far enough to keep everything in view while panning.
            reach: (hi - lo).length().max(1e-3) * 10.0,
        };
        App::new()
            .add_plugins(DefaultPlugins.set(WindowPlugin {
                primary_window: Some(window),
                ..default()
            }))
            .insert_resource(ClearColor(Color::srgb(0.92, 0.93, 0.95)))
            .insert_resource(GlobalAmbientLight {
                brightness: 200.0,
                ..default()
            })
            .insert_resource(Scene(self.meshes))
            .insert_resource(StartAt(camera))
            .add_systems(Startup, setup)
            .add_systems(Update, orbit)
            .run()
    }
}

#[derive(Resource)]
struct StartAt(Orbit);

fn setup(
    mut commands: Commands,
    scene: Res<Scene>,
    start: Res<StartAt>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    for shown in &scene.0 {
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
        commands.spawn((
            Mesh3d(meshes.add(mesh)),
            MeshMaterial3d(materials.add(material)),
        ));
    }
    // A light at the camera, shining where it looks.
    commands.spawn((
        Camera3d::default(),
        start.0.transform(),
        start.0.projection(),
        start.0,
        children![(
            DirectionalLight {
                illuminance: 8000.0,
                ..default()
            },
            Transform::default(),
        )],
    ));
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
        let mut viewer = Viewer::new("test");
        viewer.add_mesh(&mesh, |_| GolfColor::rgba(255, 0, 0, 128));
        let shown = &viewer.meshes[0];
        assert_eq!(shown.positions.len(), 3 * mesh.triangles.len());
        assert!(shown.transparent);
        assert_eq!(shown.colors[0][..3], [1.0, 0.0, 0.0]);
        let (lo, hi) = viewer.bounds().unwrap();
        assert!((hi - lo - Vec3::new(1.0, 2.0, 3.0)).length() < 1e-6);
    }
}
