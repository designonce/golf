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
use golf_model::EdgeTreatment;
use golf_model::ExtrudeEnds;
use golf_model::Path;
use golf_model::Profile;
use golf_model::Region;
use golf_model::Segment;
use golf_model::circular_pattern;
use golf_model::extrude;
use golf_model::extrude_drafted;
use golf_model::extrude_hollow;
use golf_model::extrude_with;
use golf_model::loft;
use golf_model::pipe;
use golf_model::primitives;
use golf_model::primitives::axial_sketch;
use golf_model::revolve;
use golf_model::revolve_about;
use golf_model::revolve_by;
use nalgebra::Vector2;
use nalgebra::Vector3;

type Part = (String, Body<World>);

fn main() -> Result<(), Box<dyn Error>> {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "showcase.step".to_string());

    // Rows of parts, 120 mm apart.
    let mut parts: Vec<Part> = Vec::new();
    let mut add = |name: &str, body: Body<World>| parts.push((name.to_string(), body));

    // Sketching and extruding.
    add("bracket", bracket(at(0.0, 0.0, 0.0))?);
    add("star prism", star(at(160.0, 30.0, 0.0))?);
    add("drafted boss", drafted_boss(at(240.0, 30.0, 0.0))?);
    add("cup", cup(at(320.0, 30.0, 0.0))?);
    add("tube", tube(at(400.0, 30.0, 0.0))?);

    // Revolving.
    add("pulley", pulley(at(0.0, 120.0, 0.0))?);
    add("vase", vase(at(110.0, 120.0, 0.0))?);
    add("hollow ball", hollow_ball(at(210.0, 120.0, 25.0))?);
    add("elbow", elbow(at(290.0, 100.0, 0.0))?);
    add("knob", knob(at(400.0, 120.0, 0.0))?);

    // Edge treatments, lofts and pipes.
    add(
        "rounded box",
        primitives::rounded_cuboid(&at(0.0, 220.0, 0.0), Vector3::new(40.0, 30.0, 20.0), 6.0)?,
    );
    add(
        "chamfered box",
        primitives::chamfered_cuboid(&at(60.0, 220.0, 0.0), Vector3::new(40.0, 30.0, 20.0), 4.0)?,
    );
    add("rounded plate", rounded_plate(at(120.0, 220.0, 0.0))?);
    add("twisted loft", twisted_loft(at(230.0, 240.0, 0.0))?);
    add("tapered loft", tapered_loft(at(300.0, 240.0, 0.0))?);
    add("bent pipe", bent_pipe(at(350.0, 210.0, 10.0))?);
    add("pipe frame", pipe_frame(at(0.0, 300.0, 5.0))?);

    // Primitives, and a pattern of them.
    add(
        "box",
        primitives::cuboid(&at(0.0, 400.0, 0.0), Vector3::new(30.0, 20.0, 10.0))?,
    );
    add(
        "cylinder",
        primitives::cylinder(&at(60.0, 410.0, 0.0), 10.0, 30.0)?,
    );
    add(
        "cone",
        primitives::cone(&at(100.0, 410.0, 0.0), 12.0, 0.0, 30.0)?,
    );
    add(
        "frustum",
        primitives::cone(&at(140.0, 410.0, 0.0), 14.0, 6.0, 20.0)?,
    );
    add("sphere", primitives::sphere(&at(180.0, 410.0, 15.0), 15.0)?);
    add(
        "torus",
        primitives::torus(&tilted(230.0, 410.0, 15.0), 15.0, 5.0)?,
    );
    let peg = primitives::cylinder(&at(330.0, 410.0, 0.0), 4.0, 15.0)?;
    let ring = circular_pattern(
        &peg,
        Point::new(Vector3::new(310.0, 410.0, 0.0)),
        Vector::new(Vector3::z()),
        6,
        TAU,
    )?;
    for (i, peg) in ring.into_iter().enumerate() {
        add(&format!("peg {}", i + 1), peg);
    }

    println!(
        "{:<14} {:>6} {:>6} {:>9} {:>6} {:>12}",
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
            "{:<14} {:>6} {:>6} {:>9} {:>6} {:>12.1}",
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

/// A plate with filleted corners, two bolt holes and a slot, extruded with its
/// top edges rounded.
///
/// Sketch fillets, circle and slot holes, and a filleted extrusion end (which
/// follows the outline and the holes round).
fn bracket(placement: Placement<World>) -> Result<Body<World>, Box<dyn Error>> {
    let outline =
        Profile::rectangle(Vector2::zeros(), Vector2::new(100.0, 60.0))?.fillet_all(8.0)?;
    let bolt = |x: f64, y: f64| Profile::circle(Vector2::new(x, y), 5.0);
    let slot = Profile::slot(Vector2::new(45.0, 30.0), Vector2::new(75.0, 30.0), 6.0)?;
    let region = Region::new(outline, vec![bolt(15.0, 15.0)?, bolt(15.0, 45.0)?, slot]);
    let ends = ExtrudeEnds {
        bottom: EdgeTreatment::None,
        top: EdgeTreatment::Fillet(2.0),
    };
    Ok(extrude_with(&placement, &region, 8.0, ends)?)
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

/// A rounded hexagon, extruded with its sides drafted in by 8°: planes over
/// its sides and cones round its corners.
fn drafted_boss(placement: Placement<World>) -> Result<Body<World>, Box<dyn Error>> {
    let hexagon = Profile::regular_polygon(Vector2::zeros(), 30.0, 6, 0.0)?.fillet_all(6.0)?;
    Ok(extrude_drafted(
        &placement,
        &hexagon.into(),
        25.0,
        8f64.to_radians(),
    )?)
}

/// A cup: a rounded square hollowed out with 3 mm walls and a 4 mm floor.
fn cup(placement: Placement<World>) -> Result<Body<World>, Box<dyn Error>> {
    let square =
        Profile::rounded_rectangle(Vector2::new(-25.0, -25.0), Vector2::new(25.0, 25.0), 8.0)?;
    Ok(extrude_hollow(&placement, &square, 40.0, 3.0, Some(4.0))?)
}

/// A hexagonal tube: open at both ends.
fn tube(placement: Placement<World>) -> Result<Body<World>, Box<dyn Error>> {
    let hexagon = Profile::regular_polygon(Vector2::zeros(), 25.0, 6, PI / 6.0)?;
    Ok(extrude_hollow(&placement, &hexagon, 50.0, 4.0, None)?)
}

/// A pulley: a wheel with a bore through it and a semicircular groove round
/// its rim, revolved about the placement's z axis.
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

/// A ball with a ring-shaped cavity inside it: a half disc with a circular
/// hole, revolved, so the hole becomes a toroidal void (a second shell).
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

/// A pipe elbow: a hollow circle (an annulus) revolved through a quarter turn.
fn elbow(placement: Placement<World>) -> Result<Body<World>, Box<dyn Error>> {
    let annulus = Region::new(
        Profile::circle(Vector2::new(30.0, 0.0), 10.0)?,
        vec![Profile::circle(Vector2::new(30.0, 0.0), 7.0)?],
    );
    Ok(revolve_by(&axial_sketch(&placement), &annulus, PI / 2.0)?)
}

/// A knob: a filleted profile revolved about a sloping axis drawn beside it,
/// through three quarters of a turn.
fn knob(placement: Placement<World>) -> Result<Body<World>, Box<dyn Error>> {
    let profile =
        Profile::rectangle(Vector2::new(5.0, 0.0), Vector2::new(20.0, 15.0))?.fillet_all(3.0)?;
    let axis = (Vector2::new(0.0, 0.0), Vector2::new(0.2, 1.0));
    Ok(revolve_about(
        &placement,
        &profile.into(),
        axis.0,
        axis.1,
        1.5 * PI,
    )?)
}

/// A plate rounded all over its top: sketch fillets round the corners, so the
/// filleted top edge sweeps cylinders along the sides and tori round them.
fn rounded_plate(placement: Placement<World>) -> Result<Body<World>, Box<dyn Error>> {
    let outline =
        Profile::rectangle(Vector2::zeros(), Vector2::new(80.0, 50.0))?.fillet_all(10.0)?;
    let ends = ExtrudeEnds {
        bottom: EdgeTreatment::chamfer(1.0),
        top: EdgeTreatment::Fillet(4.0),
    };
    Ok(extrude_with(&placement, &outline.into(), 12.0, ends)?)
}

/// A square lofted to a smaller square turned 45°: four twisted ruled faces.
fn twisted_loft(placement: Placement<World>) -> Result<Body<World>, Box<dyn Error>> {
    let square = |half: f64, turn: f64| -> Result<Region, Box<dyn Error>> {
        let corners: Vec<Vector2<f64>> = (0..4)
            .map(|i| {
                let angle = turn + PI / 4.0 + PI / 2.0 * i as f64;
                Vector2::new(angle.cos(), angle.sin()) * half
            })
            .collect();
        Ok(Profile::polygon(&corners)?.into())
    };
    let top = Placement::new(
        placement.point(Vector3::new(0.0, 0.0, 40.0)),
        placement.rotation,
    );
    Ok(loft(
        (&placement, &square(25.0, 0.0)?),
        (&top, &square(15.0, PI / 4.0)?),
    )?)
}

/// A circle lofted to a smaller circle off to one side and tipped over: an
/// oblique, skewed cone as exact rational NURBS.
fn tapered_loft(placement: Placement<World>) -> Result<Body<World>, Box<dyn Error>> {
    let top = Placement::from_axes(
        placement.point(Vector3::new(10.0, 0.0, 40.0)),
        placement.vector(Vector3::new(0.3, 0.0, 1.0)),
        placement.vector(Vector3::x()),
    );
    Ok(loft(
        (&placement, &Profile::circle(Vector2::zeros(), 20.0)?.into()),
        (&top, &Profile::circle(Vector2::zeros(), 8.0)?.into()),
    )?)
}

/// A pipe along a path of lines and bends, both ways.
fn bent_pipe(placement: Placement<World>) -> Result<Body<World>, Box<dyn Error>> {
    let path = Profile::builder(Vector2::zeros())
        .line_to(Vector2::new(30.0, 0.0))
        .tangent_arc_to(Vector2::new(50.0, 20.0))
        .line_to(Vector2::new(50.0, 40.0))
        .arc_to(Vector2::new(70.0, 60.0), Vector2::new(70.0, 40.0), false)
        .path()?;
    Ok(pipe(&placement, &path, 6.0)?)
}

/// A closed pipe round a rectangle, mitred at its corners.
fn pipe_frame(placement: Placement<World>) -> Result<Body<World>, Box<dyn Error>> {
    let frame = Path::from(Profile::rectangle(
        Vector2::zeros(),
        Vector2::new(100.0, 60.0),
    )?);
    Ok(pipe(&placement, &frame, 5.0)?)
}
