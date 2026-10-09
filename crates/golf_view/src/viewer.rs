use bevy::prelude::*;
use golf_color::Color as GolfColor;
use golf_manifold::Space;

use crate::plugin::ShowScene;
use crate::plugin::ViewerPlugin;
use crate::scene::Scene;

/// Meshes to show, gathered before opening the window.
#[derive(Clone, Debug, Default)]
pub struct Viewer {
    title: String,
    scene: Scene,
}

impl Viewer {
    pub fn new(title: impl Into<String>) -> Self {
        Self {
            title: title.into(),
            scene: Scene::new(),
        }
    }

    /// Adds a mesh, each triangle coloured by `color` of its face.
    pub fn add_mesh<S: Space<3>, F>(
        &mut self,
        mesh: &golf_mesh::Mesh<S, F>,
        color: impl Fn(&F) -> GolfColor,
    ) -> &mut Self {
        self.scene.add_mesh(mesh, color);
        self
    }

    /// The smallest box holding every mesh, as its least and greatest corners.
    pub fn bounds(&self) -> Option<(Vec3, Vec3)> {
        self.scene.bounds()
    }

    /// Opens a window showing the meshes, until it's closed.
    pub fn run(self) -> AppExit {
        App::new()
            .add_plugins(DefaultPlugins.set(WindowPlugin {
                primary_window: Some(Window {
                    title: self.title,
                    ..default()
                }),
                ..default()
            }))
            .add_plugins(ViewerPlugin)
            .insert_resource(ShowScene(Some(self.scene)))
            .run()
    }
}
