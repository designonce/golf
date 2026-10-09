use bevy::camera::ScalingMode;
use bevy::input::mouse::AccumulatedMouseMotion;
use bevy::input::mouse::AccumulatedMouseScroll;
use bevy::input::mouse::MouseScrollUnit;
use bevy::prelude::*;

/// An orthographic camera orbiting `target` with z up: turned `yaw` about z
/// from the +x side, raised `pitch` above the xy plane, and showing at least
/// `size` of the scene across the window, both ways.
#[derive(Component, Clone, Copy, Debug, PartialEq)]
pub struct Orbit {
    pub target: Vec3,
    pub yaw: f32,
    pub pitch: f32,
    pub size: f32,
    /// How far the camera stands back from the target, so everything within
    /// this of the target is in view.
    pub reach: f32,
}

impl Orbit {
    /// The camera's placement: `reach` back from the target, looking at it.
    pub fn transform(&self) -> Transform {
        let (sin_yaw, cos_yaw) = self.yaw.sin_cos();
        let (sin_pitch, cos_pitch) = self.pitch.sin_cos();
        let offset = Vec3::new(cos_pitch * cos_yaw, cos_pitch * sin_yaw, sin_pitch);
        Transform::from_translation(self.target + offset * self.reach)
            .looking_at(self.target, Vec3::Z)
    }

    /// The camera's projection, showing `size` and all of `reach` either side
    /// of the target.
    pub fn projection(&self) -> Projection {
        Projection::Orthographic(OrthographicProjection {
            near: 0.0,
            far: 2.0 * self.reach,
            scaling_mode: ScalingMode::AutoMin {
                min_width: self.size,
                min_height: self.size,
            },
            ..OrthographicProjection::default_3d()
        })
    }
}

/// Turns, pans and zooms orbiting cameras with the mouse.
pub(crate) fn orbit(
    buttons: Res<ButtonInput<MouseButton>>,
    motion: Res<AccumulatedMouseMotion>,
    scroll: Res<AccumulatedMouseScroll>,
    windows: Query<&Window>,
    mut cameras: Query<(&mut Orbit, &mut Transform, &mut Projection)>,
) {
    // Pan so the scene follows the cursor: `size` spans the window's shorter
    // side.
    let pixels = windows
        .iter()
        .next()
        .map_or(720.0, |w| w.width().min(w.height()).max(1.0));
    for (mut orbit, mut transform, mut projection) in &mut cameras {
        if buttons.pressed(MouseButton::Left) {
            orbit.yaw -= motion.delta.x * 0.01;
            orbit.pitch = (orbit.pitch + motion.delta.y * 0.01).clamp(-1.55, 1.55);
        }
        if buttons.pressed(MouseButton::Right) {
            let scale = orbit.size / pixels;
            orbit.target +=
                (transform.right() * -motion.delta.x + transform.up() * motion.delta.y) * scale;
        }
        let lines = match scroll.unit {
            MouseScrollUnit::Line => scroll.delta.y,
            MouseScrollUnit::Pixel => scroll.delta.y / 40.0,
        };
        orbit.size *= (-lines * 0.1).exp();
        *transform = orbit.transform();
        *projection = orbit.projection();
    }
}
