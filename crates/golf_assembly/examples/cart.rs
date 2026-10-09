//! A cart, as an assembly with shared parts and sub-assemblies, written as a
//! STEP assembly.
//!
//! ```text
//! cargo run -p golf_assembly --example cart -- cart.step
//! cargo run -p golf_assembly --example cart -- cart_flat.step --flat
//! cargo run -p golf_assembly --example cart -- --view
//! ```
//!
//! The cart is a tray, a handle and two axle groups (an axle and two wheels),
//! so the wheel is defined once and placed four times; the back axle is a copy
//! of the front, coloured differently. With `--flat`, every
//! placed part is written separately, moved into place, instead. With
//! `--view`, the cart is shown in a window rather than written.

use core::f64::consts::FRAC_PI_2;
use std::error::Error;

use golf_assembly::Assembly;
use golf_color::Color;
use golf_export_step::StepFile;
use golf_export_step::StepOptions;
use golf_geom::Placement;
use golf_manifold::Point;
use golf_manifold::World;
use golf_mesh::Tolerance;
use golf_mesh::mesh;
use golf_model::EdgeTreatment;
use golf_model::ExtrudeEnds;
use golf_model::Profile;
use golf_model::extrude_hollow;
use golf_model::extrude_with;
use golf_model::pipe;
use golf_model::primitives;
use golf_view::Viewer;
use nalgebra::Isometry3;
use nalgebra::Translation3;
use nalgebra::UnitQuaternion;
use nalgebra::Vector2;
use nalgebra::Vector3;

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let flag = |name: &str| args.iter().any(|a| a == name);
    let path = args
        .iter()
        .find(|a| !a.starts_with("--"))
        .cloned()
        .unwrap_or_else(|| "cart.step".to_string());

    let cart = cart()?;
    if flag("--view") {
        let mut viewer = Viewer::new("cart");
        for (_, body) in cart.posed_bodies()? {
            let mesh = mesh(&body, &Tolerance::new(0.05, 0.3))?;
            viewer.add_mesh(&mesh, |&face| {
                body.face_color(face).unwrap_or(Color::rgb(170, 175, 180))
            });
        }
        viewer.run();
        return Ok(());
    }
    let mut file = StepFile::new(StepOptions {
        name: "cart".to_string(),
        ..StepOptions::default()
    })?;
    match flag("--flat") {
        false => {
            file.add_assembly(&cart)?;
        }
        true => {
            for (name, body) in cart.posed_bodies()? {
                file.add_body(&name, &body)?;
            }
        }
    }
    std::fs::write(&path, file.finish())?;
    for occurrence in cart.occurrences() {
        let at = occurrence.placement * nalgebra::Point3::origin();
        println!(
            "{:<20} at ({:6.1}, {:6.1}, {:6.1})",
            cart.path_name(occurrence.frame),
            at.x,
            at.y,
            at.z
        );
    }
    println!("wrote {path}");
    Ok(())
}

fn cart() -> Result<Assembly<World>, Box<dyn Error>> {
    let origin = Placement::<World>::at(Point::new(Vector3::zeros()));
    let mut cart = Assembly::new("cart");

    // Parts, each in its own coordinates.
    let tray_outline =
        Profile::rounded_rectangle(Vector2::zeros(), Vector2::new(120.0, 70.0), 8.0)?;
    let tray = cart.add_part(
        "tray",
        extrude_hollow(&origin, &tray_outline, 30.0, 3.0, Some(4.0))?,
    );
    let tyre = Profile::circle(Vector2::zeros(), 15.0)?;
    let ends = ExtrudeEnds::both(EdgeTreatment::Fillet(2.5));
    let wheel = cart.add_part(
        "wheel",
        extrude_with(&origin, &tyre.into(), 8.0, ends)?.with_color(Color::rgb(30, 30, 30)),
    );
    let axle = cart.add_part("axle", primitives::cylinder(&origin, 2.5, 69.0)?);
    let handle_path = Profile::builder(Vector2::zeros())
        .line_to(Vector2::new(0.0, 45.0))
        .tangent_arc_to(Vector2::new(-20.0, 65.0))
        .line_to(Vector2::new(-60.0, 65.0))
        .path()?;
    let handle = cart.add_part(
        "handle",
        pipe(&origin, &handle_path, 3.0)?.with_color(Color::rgb(200, 30, 30)),
    );

    // The front axle group, under the tray: the axle along y, between a wheel
    // at each end turned to face along it (clear of them by half a
    // millimetre). The back axle is a copy.
    let root = cart.root();
    let along_y = |y: f64| {
        Isometry3::from_parts(
            Translation3::new(0.0, y, 0.0),
            UnitQuaternion::from_scaled_axis(Vector3::x() * -FRAC_PI_2),
        )
    };
    let front = cart.add_group(
        root,
        "axle assembly",
        Isometry3::translation(20.0, 0.0, 12.0),
    )?;
    cart.add_instance(front, axle, along_y(0.5))?;
    cart.add_instance(front, wheel, along_y(-8.5))?;
    cart.add_instance(front, wheel, along_y(70.5))?;
    let back = cart.copy_group(front, root, Isometry3::translation(100.0, 0.0, 12.0))?;
    // The back axle in steel blue, wheels and all.
    cart.set_color(back, Some(Color::rgb(70, 110, 160)))?;

    // The tray resting just above the axles, the handle standing behind it:
    // its sketch stood up, x along the cart and y upwards.
    let tray_frame = cart.add_instance(root, tray, Isometry3::translation(0.0, 0.0, 15.0))?;
    cart.set_color(tray_frame, Some(Color::rgb(240, 200, 40)))?;
    let upright = UnitQuaternion::from_scaled_axis(Vector3::x() * FRAC_PI_2);
    cart.add_instance(
        root,
        handle,
        Isometry3::from_parts(Translation3::new(124.0, 35.0, 46.0), upright),
    )?;
    Ok(cart)
}
