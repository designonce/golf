use bevy::prelude::*;

use crate::orbit::Orbit;
use crate::orbit::orbit;
use crate::scene::Scene;

/// Shows a [`Scene`] with an orbiting camera lit from where it looks: put a
/// scene in [`ShowScene`] to show it in place of the last, framed. F frames
/// the scene again.
pub struct ViewerPlugin;

/// The scene to show next; taken once it's shown.
#[derive(Resource, Default)]
pub struct ShowScene(pub Option<Scene>);

/// The scene's meshes, to remove when another replaces it.
#[derive(Component)]
struct SceneMesh;

/// The bounds of the scene shown, for framing it again.
#[derive(Resource, Default)]
struct Shown(Option<(Vec3, Vec3)>);

impl Plugin for ViewerPlugin {
    fn build(&self, app: &mut App) {
        app.insert_resource(ClearColor(Color::srgb(0.92, 0.93, 0.95)))
            .insert_resource(GlobalAmbientLight {
                brightness: 200.0,
                ..default()
            })
            .init_resource::<ShowScene>()
            .init_resource::<Shown>()
            .add_systems(Startup, camera)
            .add_systems(Update, (show, frame_on_key, orbit).chain());
    }
}

fn camera(mut commands: Commands) {
    let start = framed(
        Orbit {
            target: Vec3::ZERO,
            yaw: -1.0,
            pitch: 0.5,
            size: 1.0,
            reach: 1.0,
        },
        (Vec3::splat(-1.0), Vec3::splat(1.0)),
    );
    // A light at the camera, shining where it looks.
    commands.spawn((
        Camera3d::default(),
        start.transform(),
        start.projection(),
        start,
        children![(
            DirectionalLight {
                illuminance: 8000.0,
                ..default()
            },
            Transform::default(),
        )],
    ));
}

/// `orbit` looking at the middle of `(lo, hi)`, with all of it in view.
fn framed(orbit: Orbit, (lo, hi): (Vec3, Vec3)) -> Orbit {
    let diagonal = (hi - lo).length().max(1e-3);
    Orbit {
        target: (lo + hi) / 2.0,
        size: diagonal * 1.1,
        // Far enough to keep everything in view while panning.
        reach: diagonal * 10.0,
        ..orbit
    }
}

fn show(
    mut commands: Commands,
    mut next: ResMut<ShowScene>,
    mut shown: ResMut<Shown>,
    old: Query<Entity, With<SceneMesh>>,
    mut cameras: Query<&mut Orbit>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    let Some(scene) = next.0.take() else { return };
    for entity in &old {
        commands.entity(entity).despawn();
    }
    for (mesh, material) in scene.bevy_meshes() {
        commands.spawn((
            Mesh3d(meshes.add(mesh)),
            MeshMaterial3d(materials.add(material)),
            SceneMesh,
        ));
    }
    shown.0 = scene.bounds();
    if let Some(bounds) = shown.0 {
        for mut orbit in &mut cameras {
            *orbit = framed(*orbit, bounds);
        }
    }
}

fn frame_on_key(
    keys: Res<ButtonInput<KeyCode>>,
    shown: Res<Shown>,
    mut cameras: Query<&mut Orbit>,
) {
    if let (true, Some(bounds)) = (keys.just_pressed(KeyCode::KeyF), shown.0) {
        for mut orbit in &mut cameras {
            *orbit = framed(*orbit, bounds);
        }
    }
}
