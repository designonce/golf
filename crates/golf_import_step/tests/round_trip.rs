//! Bodies and assemblies written by golf_export_step read back as they were.

use core::f64::consts::FRAC_PI_2;

use golf_assembly::Assembly;
use golf_brep::Body;
use golf_color::Color;
use golf_export_step::StepFile;
use golf_export_step::StepOptions;
use golf_geom::Placement;
use golf_import_step::read_step;
use golf_manifold::Point;
use golf_manifold::Vector;
use golf_manifold::World;
use golf_mesh::Tolerance;
use golf_mesh::mesh;
use golf_model::EdgeTreatment;
use golf_model::ExtrudeEnds;
use golf_model::Profile;
use golf_model::extrude_with;
use golf_model::primitives;
use nalgebra::Isometry3;
use nalgebra::Point3;
use nalgebra::Translation3;
use nalgebra::UnitQuaternion;
use nalgebra::Vector2;
use nalgebra::Vector3;

fn tilted() -> Placement<World> {
    Placement::from_axes(
        Point::new(Vector3::new(3.0, -2.0, 1.0)),
        Vector::new(Vector3::new(0.2, 0.3, 1.0)),
        Vector::new(Vector3::x()),
    )
}

fn volume(body: &Body<World>) -> f64 {
    let m = mesh(body, &Tolerance::new(1e-3, 0.3)).unwrap();
    assert!(m.is_watertight(), "{} open edges", m.open_edges().len());
    m.signed_volume()
}

fn write_bodies(bodies: &[(&str, &Body<World>)]) -> String {
    let mut file = StepFile::new(StepOptions::default()).unwrap();
    for (name, body) in bodies {
        file.add_body(name, *body).unwrap();
    }
    file.finish()
}

fn bodies() -> Vec<(&'static str, Body<World>)> {
    let tyre = Profile::circle(Vector2::zeros(), 15.0).unwrap();
    vec![
        (
            "box",
            primitives::cuboid(&tilted(), Vector3::new(4.0, 3.0, 2.0)).unwrap(),
        ),
        (
            "cylinder",
            primitives::cylinder(&tilted(), 2.0, 5.0).unwrap(),
        ),
        ("cone", primitives::cone(&tilted(), 2.0, 0.0, 4.0).unwrap()),
        ("sphere", primitives::sphere(&tilted(), 3.0).unwrap()),
        ("torus", primitives::torus(&tilted(), 5.0, 1.5).unwrap()),
        (
            "rounded box",
            primitives::rounded_cuboid(&tilted(), Vector3::new(8.0, 6.0, 4.0), 1.0).unwrap(),
        ),
        (
            "wheel",
            extrude_with(
                &tilted(),
                &tyre.into(),
                8.0,
                ExtrudeEnds::both(EdgeTreatment::Fillet(2.5)),
            )
            .unwrap(),
        ),
    ]
}

#[test]
fn bodies_read_back_with_their_shape() {
    let bodies = bodies();
    let refs: Vec<(&str, &Body<World>)> = bodies.iter().map(|(n, b)| (*n, b)).collect();
    let text = write_bodies(&refs);
    assert!(text.contains("SEAM_CURVE") && text.contains("PCURVE"));
    let import = read_step(text.as_bytes()).unwrap();
    assert!(import.warnings.is_empty(), "{:?}", import.warnings);
    // The file's pcurves are used; only coedges golf left without one (the
    // primitives' caps) are computed.
    assert!(
        import.heal.kept + import.heal.moved > 4 * import.heal.computed,
        "{:?}",
        import.heal
    );
    assert_eq!(import.heal.seams, 0);
    let read: Vec<_> = import.assembly.parts().collect();
    assert_eq!(read.len(), bodies.len());
    for ((name, original), (_, part)) in bodies.iter().zip(read) {
        assert_eq!(&part.name, name);
        assert_eq!(part.body.faces().len(), original.faces().len(), "{name}");
        assert_eq!(part.body.validate(1e-6), Ok(()), "{name}");
        let (a, b) = (volume(original), volume(&part.body));
        assert!((a - b).abs() < 1e-4 * a, "{name}: {a} then {b}");
    }
}

#[test]
fn units_are_converted_to_millimetres() {
    let cube = primitives::cuboid(&tilted(), Vector3::new(1.0, 2.0, 3.0)).unwrap();
    let text = write_bodies(&[("cube", &cube)]);
    assert!(text.contains("SI_UNIT(.MILLI.,.METRE.)"));
    let in_metres = text.replace("SI_UNIT(.MILLI.,.METRE.)", "SI_UNIT($,.METRE.)");
    let import = read_step(in_metres.as_bytes()).unwrap();
    let (_, part) = import.assembly.parts().next().unwrap();
    assert!((volume(&part.body) - 6e9).abs() < 1e-3 * 6e9);
}

#[test]
fn assemblies_read_back_with_their_structure_and_colours() {
    let origin = Placement::<World>::at(Point::new(Vector3::zeros()));
    let mut asm = Assembly::new("cart");
    let chassis = asm.add_part(
        "chassis",
        primitives::cuboid(&origin, Vector3::new(40.0, 20.0, 5.0)).unwrap(),
    );
    let wheel = asm.add_part(
        "wheel",
        primitives::cylinder(&origin, 6.0, 3.0)
            .unwrap()
            .with_color(Color::rgb(20, 20, 20)),
    );
    let root = asm.root();
    asm.add_instance(root, chassis, Isometry3::identity())
        .unwrap();
    let front = asm
        .add_group(root, "axle", Isometry3::translation(5.0, 0.0, 0.0))
        .unwrap();
    let side = |y: f64| {
        Isometry3::from_parts(
            Translation3::new(0.0, y, 0.0),
            UnitQuaternion::from_scaled_axis(Vector3::x() * FRAC_PI_2),
        )
    };
    asm.add_instance(front, wheel, side(-1.0)).unwrap();
    asm.add_instance(front, wheel, side(24.0)).unwrap();
    let back = asm
        .copy_group(front, root, Isometry3::translation(35.0, 0.0, 0.0))
        .unwrap();
    asm.set_color(back, Some(Color::rgb(200, 0, 0))).unwrap();

    let mut file = StepFile::new(StepOptions::default()).unwrap();
    file.add_assembly(&asm).unwrap();
    let import = read_step(file.finish().as_bytes()).unwrap();
    assert!(import.warnings.is_empty(), "{:?}", import.warnings);
    let read = &import.assembly;

    let (before, after) = (asm.occurrences(), read.occurrences());
    assert_eq!(after.len(), before.len());
    for (b, a) in before.iter().zip(&after) {
        assert_eq!(read.path_name(a.frame), asm.path_name(b.frame));
        let moved = (a.placement * Point3::new(1.0, 2.0, 3.0)
            - b.placement * Point3::new(1.0, 2.0, 3.0))
        .norm();
        assert!(moved < 1e-9, "{}", asm.path_name(b.frame));
        assert_eq!(a.color, b.color, "{}", asm.path_name(b.frame));
    }
    // The plain wheel is one part placed twice; the red one its own.
    let wheels: Vec<usize> = read
        .parts()
        .filter(|(_, p)| p.name == "wheel")
        .map(|(id, _)| read.instance_count(id))
        .collect();
    assert_eq!(wheels, vec![2, 2]);
}
