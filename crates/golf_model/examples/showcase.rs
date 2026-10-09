//! Builds one of everything `golf_model` can make, and writes them all to a
//! STEP file.
//!
//! ```text
//! cargo run -p golf_model --example showcase -- showcase.step
//! ```
//!
//! Each part is checked (topology valid, mesh watertight) and summarised.
//! Lengths are millimetres.

use core::f64::consts::PI;
use core::f64::consts::TAU;
use std::error::Error;

use golf_brep::Body;
use golf_export_step::StepFile;
use golf_export_step::StepOptions;
use golf_geom::Placement;
use golf_manifold::Point;
use golf_manifold::Vector;
use golf_manifold::World;
use golf_mesh::Tolerance;
use golf_mesh::mesh;
use golf_model::Profile;
use golf_model::Region;
use golf_model::Segment;
use golf_model::extrude;
use golf_model::primitives;
use golf_model::primitives::axial_sketch;
use golf_model::revolve;
use nalgebra::Vector2;
use nalgebra::Vector3;

fn main() -> Result<(), Box<dyn Error>> {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "showcase.step".to_string());

    let parts: Vec<(&str, Body<World>)> = vec![
        ("bracket", bracket(at(0.0, 0.0, 0.0))?),
        ("star prism", star(at(140.0, 30.0, 0.0))?),
        ("pulley", pulley(at(0.0, 120.0, 0.0))?),
        ("vase", vase(at(110.0, 120.0, 0.0))?),
        ("hollow ball", hollow_ball(at(210.0, 120.0, 25.0))?),
        (
            "box",
            primitives::cuboid(&at(0.0, 220.0, 0.0), Vector3::new(30.0, 20.0, 10.0))?,
        ),
        (
            "cylinder",
            primitives::cylinder(&at(60.0, 230.0, 0.0), 10.0, 30.0)?,
        ),
        (
            "cone",
            primitives::cone(&at(100.0, 230.0, 0.0), 12.0, 0.0, 30.0)?,
        ),
        (
            "frustum",
            primitives::cone(&at(140.0, 230.0, 0.0), 14.0, 6.0, 20.0)?,
        ),
        ("sphere", primitives::sphere(&at(180.0, 230.0, 15.0), 15.0)?),
        (
            "torus",
            primitives::torus(&tilted(230.0, 230.0, 15.0), 15.0, 5.0)?,
        ),
    ];

    println!(
        "{:<12} {:>6} {:>6} {:>9} {:>6} {:>12}",
        "part", "faces", "edges", "vertices", "genus", "volume mm³"
    );
    let mut file = StepFile::new(StepOptions {
        name: "golf showcase".to_string(),
        ..StepOptions::default()
    })?;
    for (name, body) in &parts {
        if let Err(errors) = body.validate(1e-9) {
            return Err(format!("{name} is invalid: {errors:?}").into());
        }
        let tessellation = mesh(body, &Tolerance::new(0.01, 0.3))?;
        if !tessellation.is_watertight() {
            return Err(format!("{name} doesn't mesh watertight").into());
        }
        let genus = body.genus().map_or("-".to_string(), |g| g.to_string());
        println!(
            "{:<12} {:>6} {:>6} {:>9} {:>6} {:>12.1}",
            name,
            body.faces().len(),
            body.edges().len(),
            body.vertices().len(),
            genus,
            tessellation.signed_volume()
        );
        file.add_body(name, body)?;
    }
    std::fs::write(&path, file.finish())?;
    println!("wrote {} parts to {path}", parts.len());
    Ok(())
}

/// A placement at `(x, y, z)` with the world's axes.
fn at(x: f64, y: f64, z: f64) -> Placement<World> {
    Placement::at(Point::new(Vector3::new(x, y, z)))
}

/// A placement at `(x, y, z)` tilted off the world's axes.
fn tilted(x: f64, y: f64, z: f64) -> Placement<World> {
    Placement::from_axes(
        Point::new(Vector3::new(x, y, z)),
        Vector::new(Vector3::new(0.4, -0.3, 1.0)),
        Vector::new(Vector3::x()),
    )
}

/// A plate with rounded corners, two bolt holes and a slot, extruded.
///
/// Uses the profile builder (lines and arcs), circular and stadium-shaped
/// holes, and a region with several holes.
fn bracket(placement: Placement<World>) -> Result<Body<World>, Box<dyn Error>> {
    let (w, h, r) = (100.0, 60.0, 8.0);
    let outline = Profile::builder(Vector2::new(r, 0.0))
        .line_to(Vector2::new(w - r, 0.0))
        .arc_to(Vector2::new(w, r), Vector2::new(w - r, r), true)
        .line_to(Vector2::new(w, h - r))
        .arc_to(Vector2::new(w - r, h), Vector2::new(w - r, h - r), true)
        .line_to(Vector2::new(r, h))
        .arc_to(Vector2::new(0.0, h - r), Vector2::new(r, h - r), true)
        .line_to(Vector2::new(0.0, r))
        .arc_to(Vector2::new(r, 0.0), Vector2::new(r, r), true)
        .close()?;
    let bolt = |x: f64, y: f64| Profile::circle(Vector2::new(x, y), 5.0);
    // A slot: two semicircles joined by straight sides.
    let slot = Profile::builder(Vector2::new(40.0, 25.0))
        .line_to(Vector2::new(70.0, 25.0))
        .arc_to(Vector2::new(70.0, 35.0), Vector2::new(70.0, 30.0), true)
        .line_to(Vector2::new(40.0, 35.0))
        .arc_to(Vector2::new(40.0, 25.0), Vector2::new(40.0, 30.0), true)
        .close()?;
    let region = Region::new(outline, vec![bolt(15.0, 15.0)?, bolt(15.0, 45.0)?, slot]);
    Ok(extrude(&placement, &region, 8.0)?)
}

/// A five-pointed star, extruded backwards (below its sketch plane).
fn star(placement: Placement<World>) -> Result<Body<World>, Box<dyn Error>> {
    let points: Vec<Vector2<f64>> = (0..10)
        .map(|i| {
            let angle = PI / 2.0 + TAU * i as f64 / 10.0;
            let radius = if i % 2 == 0 { 30.0 } else { 12.0 };
            Vector2::new(angle.cos(), angle.sin()) * radius
        })
        .collect();
    Ok(extrude(
        &placement,
        &Profile::polygon(&points)?.into(),
        -15.0,
    )?)
}

/// A pulley: a wheel with a bore through it and a semicircular groove round
/// its rim, revolved about the placement's z axis.
///
/// The profile is off the axis (so the solid has a bore) and has a concave arc
/// (so the groove is a torus with its normal flipped).
fn pulley(placement: Placement<World>) -> Result<Body<World>, Box<dyn Error>> {
    let (bore, rim, width, groove) = (8.0, 35.0, 20.0, 6.0);
    let profile = Profile::builder(Vector2::new(bore, 0.0))
        .line_to(Vector2::new(rim, 0.0))
        .line_to(Vector2::new(rim, width / 2.0 - groove))
        .arc_to(
            Vector2::new(rim, width / 2.0 + groove),
            Vector2::new(rim, width / 2.0),
            false,
        )
        .line_to(Vector2::new(rim, width))
        .line_to(Vector2::new(bore, width))
        .close()?;
    Ok(revolve(&axial_sketch(&placement), &profile.into())?)
}

/// A vase: a flat base, a belly bulging out (an arc centred on the axis, so a
/// sphere), a neck pinching in (a concave arc, so a torus), and a flared lip,
/// closed down the axis.
fn vase(placement: Placement<World>) -> Result<Body<World>, Box<dyn Error>> {
    let profile = Profile::builder(Vector2::new(0.0, 0.0))
        .line_to(Vector2::new(20.0, 0.0))
        .arc_to(Vector2::new(20.0, 50.0), Vector2::new(0.0, 25.0), true)
        .arc_to(Vector2::new(20.0, 70.0), Vector2::new(30.0, 60.0), false)
        .line_to(Vector2::new(24.0, 74.0))
        .line_to(Vector2::new(0.0, 74.0))
        .close()?;
    Ok(revolve(&axial_sketch(&placement), &profile.into())?)
}

/// A ball with a ring-shaped cavity inside it: a half disc with a circular hole,
/// revolved, so the hole becomes a toroidal void (a second shell).
fn hollow_ball(placement: Placement<World>) -> Result<Body<World>, Box<dyn Error>> {
    let radius = 25.0;
    let half_disc = Profile::new(vec![
        Segment::Arc {
            center: Vector2::zeros(),
            radius,
            start_angle: -PI / 2.0,
            sweep: PI,
        },
        Segment::Line {
            start: Vector2::new(0.0, radius),
            end: Vector2::new(0.0, -radius),
        },
    ])?;
    let cavity = Profile::circle(Vector2::new(12.0, 0.0), 6.0)?;
    Ok(revolve(
        &axial_sketch(&placement),
        &Region::new(half_disc, vec![cavity]),
    )?)
}
