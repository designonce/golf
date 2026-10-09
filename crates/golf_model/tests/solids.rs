//! Bodies built by each operation: valid, of the expected genus, and meshing to
//! a watertight solid of the right volume.

use core::f64::consts::PI;

use golf_brep::Body;
use golf_frame::Frame;
use golf_frame::FrameTree;
use golf_geom::Placement;
use golf_manifold::Point;
use golf_manifold::Space;
use golf_manifold::Vector;
use golf_manifold::World;
use golf_mesh::Tolerance;
use golf_mesh::mesh;
use golf_model::ModelError;
use golf_model::Profile;
use golf_model::Region;
use golf_model::extrude;
use golf_model::primitives;
use golf_model::revolve;
use nalgebra::Isometry3;
use nalgebra::Vector2;
use nalgebra::Vector3;

fn tilted() -> Placement<World> {
    Placement::from_axes(
        Point::new(Vector3::new(1.0, -2.0, 0.5)),
        Vector::new(Vector3::new(0.3, -0.2, 1.0)),
        Vector::new(Vector3::new(1.0, 1.0, 0.0)),
    )
}

/// Valid, of `genus`, and meshing watertight with `volume` to within 1%.
fn assert_solid<S: Space<3>>(body: &Body<S>, volume: f64, genus: usize) {
    assert_eq!(body.validate(1e-9), Ok(()));
    assert_eq!(body.genus(), Some(genus));
    let mesh = mesh(body, &Tolerance::new(2e-3, 0.5)).unwrap();
    let open = mesh.open_edges();
    assert!(open.is_empty(), "{} open edges", open.len());
    let v = mesh.signed_volume();
    assert!(
        ((v - volume) / volume).abs() < 1e-2,
        "volume {v}, expected {volume}"
    );
}

#[test]
fn cuboid() {
    let body = primitives::cuboid(&tilted(), Vector3::new(1.0, 2.0, 3.0)).unwrap();
    assert_eq!(body.faces().len(), 6);
    assert_solid(&body, 6.0, 0);
}

#[test]
fn cylinder() {
    let body = primitives::cylinder(&tilted(), 1.5, 4.0).unwrap();
    assert_eq!(
        (
            body.vertices().len(),
            body.edges().len(),
            body.faces().len()
        ),
        (2, 3, 3)
    );
    assert_solid(&body, PI * 2.25 * 4.0, 0);
}

#[test]
fn plate_with_a_hole() {
    let outer = Profile::rectangle(Vector2::zeros(), Vector2::new(4.0, 3.0)).unwrap();
    let hole = Profile::circle(Vector2::new(2.0, 1.5), 0.5).unwrap();
    let body = extrude(&tilted(), &Region::new(outer, vec![hole]), 1.0).unwrap();
    assert_solid(&body, 12.0 - PI * 0.25, 1);
}

#[test]
fn convex_and_concave_arcs() {
    // A stadium with a semicircular bite out of one long side.
    let profile = Profile::builder(Vector2::new(0.0, -1.0))
        .line_to(Vector2::new(1.5, -1.0))
        .arc_to(Vector2::new(2.5, -1.0), Vector2::new(2.0, -1.0), false)
        .line_to(Vector2::new(4.0, -1.0))
        .arc_to(Vector2::new(4.0, 1.0), Vector2::new(4.0, 0.0), true)
        .line_to(Vector2::new(0.0, 1.0))
        .arc_to(Vector2::new(0.0, -1.0), Vector2::zeros(), true)
        .close()
        .unwrap();
    let area = 8.0 + PI - PI * 0.25 / 2.0;
    assert!((profile.signed_area() - area).abs() < 1e-12);
    let body = extrude(&tilted(), &profile.into(), 2.0).unwrap();
    assert_solid(&body, area * 2.0, 0);
}

#[test]
fn negative_extrusion_goes_backwards() {
    let square = Profile::rectangle(Vector2::zeros(), Vector2::new(1.0, 1.0)).unwrap();
    let placement = Placement::<World>::at(Point::new(Vector3::zeros()));
    let body = extrude(&placement, &square.into(), -2.0).unwrap();
    assert_solid(&body, 2.0, 0);
    assert!(body.vertices().all(|(_, v)| v.point.coords.z <= 0.0));
    let square = Profile::rectangle(Vector2::zeros(), Vector2::new(1.0, 1.0)).unwrap();
    assert_eq!(
        extrude(&placement, &square.into(), 0.0).unwrap_err(),
        ModelError::ZeroDistance
    );
}

#[test]
fn cone_and_frustum() {
    let cone = primitives::cone(&tilted(), 1.5, 0.0, 3.0).unwrap();
    assert_solid(&cone, PI * 2.25 * 3.0 / 3.0, 0);
    let (r1, r2, h) = (2.0, 1.0, 1.5);
    let frustum = primitives::cone(&tilted(), r1, r2, h).unwrap();
    assert_solid(&frustum, PI * h / 3.0 * (r1 * r1 + r1 * r2 + r2 * r2), 0);
}

#[test]
fn sphere() {
    let body = primitives::sphere(&tilted(), 1.5).unwrap();
    assert_eq!(
        (
            body.vertices().len(),
            body.edges().len(),
            body.faces().len()
        ),
        (2, 1, 1)
    );
    assert_solid(&body, 4.0 / 3.0 * PI * 1.5f64.powi(3), 0);
}

#[test]
fn torus() {
    let body = primitives::torus(&tilted(), 3.0, 1.0).unwrap();
    assert_eq!(
        (
            body.vertices().len(),
            body.edges().len(),
            body.faces().len()
        ),
        (1, 2, 1)
    );
    assert_solid(&body, 2.0 * PI * PI * 3.0, 1);
}

#[test]
fn revolved_rectangle_is_a_tube() {
    let rectangle = Profile::rectangle(Vector2::new(1.0, 0.0), Vector2::new(2.0, 3.0)).unwrap();
    let body = revolve(&tilted(), &rectangle.into()).unwrap();
    assert_solid(&body, PI * (4.0 - 1.0) * 3.0, 1);
}

#[test]
fn revolved_profile_with_a_concave_arc() {
    // A cylinder of radius 2 and height 2 with a half-disc groove of radius
    // 0.5 cut round its middle.
    let profile = Profile::builder(Vector2::zeros())
        .line_to(Vector2::new(2.0, 0.0))
        .line_to(Vector2::new(2.0, 0.5))
        .arc_to(Vector2::new(2.0, 1.5), Vector2::new(2.0, 1.0), false)
        .line_to(Vector2::new(2.0, 2.0))
        .line_to(Vector2::new(0.0, 2.0))
        .close()
        .unwrap();
    // Pappus: the groove's half disc has area π/8, centroid 2/(3π) inside x = 2.
    let groove = 2.0 * PI * (PI / 8.0) * (2.0 - 2.0 / (3.0 * PI));
    let body = revolve(&tilted(), &profile.into()).unwrap();
    assert_solid(&body, PI * 4.0 * 2.0 - groove, 0);
}

#[test]
fn revolved_hole_is_a_void() {
    let outer = Profile::rectangle(Vector2::new(1.0, 0.0), Vector2::new(3.0, 4.0)).unwrap();
    let hole = Profile::circle(Vector2::new(2.0, 2.0), 0.5).unwrap();
    let body = revolve(&tilted(), &Region::new(outer, vec![hole])).unwrap();
    assert_eq!(body.shells().len(), 2);
    // A tube (genus 1) with a toroidal void (genus 1).
    assert_solid(&body, PI * 8.0 * 4.0 - 2.0 * PI * PI * 2.0 * 0.25, 2);
}

#[test]
fn revolve_rejects_profiles_across_the_axis() {
    let across = Profile::rectangle(Vector2::new(-1.0, 0.0), Vector2::new(1.0, 1.0)).unwrap();
    assert!(matches!(
        revolve(&tilted(), &across.into()),
        Err(ModelError::CrossesAxis { .. })
    ));
    let spindle = Profile::circle(Vector2::new(0.5, 0.0), 0.5).unwrap();
    assert!(matches!(
        revolve(&tilted(), &spindle.into()),
        Err(ModelError::SpindleTorus { .. })
    ));
}

#[test]
fn bodies_in_a_frame() {
    let mut tree = FrameTree::new();
    let frame = tree.add(tree.root(), Isometry3::translation(5.0, 0.0, 0.0));
    let placement: Placement<Frame> = Placement::at(frame.point(Vector3::zeros()));
    let body = primitives::torus(&placement, 2.0, 0.5).unwrap();
    assert_eq!(body.tag(), frame);
    assert_eq!(body.validate(1e-9), Ok(()));
    let cuboid = primitives::cuboid(&placement, Vector3::new(1.0, 1.0, 1.0)).unwrap();
    assert!(
        mesh(&cuboid, &Tolerance::new(1e-3, 0.5))
            .unwrap()
            .is_watertight()
    );
}

#[test]
fn solids_export_to_step() {
    use golf_export_step::StepFile;
    use golf_export_step::StepOptions;
    let outer = Profile::rectangle(Vector2::new(1.0, 0.0), Vector2::new(3.0, 4.0)).unwrap();
    let hole = Profile::circle(Vector2::new(2.0, 2.0), 0.5).unwrap();
    let mut file = StepFile::new(StepOptions::default()).unwrap();
    file.add_body("sphere", &primitives::sphere(&tilted(), 1.0).unwrap())
        .unwrap();
    file.add_body("torus", &primitives::torus(&tilted(), 3.0, 1.0).unwrap())
        .unwrap();
    file.add_body("cone", &primitives::cone(&tilted(), 1.0, 0.0, 2.0).unwrap())
        .unwrap();
    file.add_body(
        "void",
        &revolve(&tilted(), &Region::new(outer, vec![hole])).unwrap(),
    )
    .unwrap();
    let text = file.finish();
    let count = |entity: &str| text.matches(&format!("={entity}(")).count();
    assert_eq!(count("MANIFOLD_SOLID_BREP"), 3);
    assert_eq!(count("BREP_WITH_VOIDS"), 1);
    assert_eq!(count("ORIENTED_CLOSED_SHELL"), 1);
    assert_eq!(count("SPHERICAL_SURFACE"), 1);
    assert_eq!(count("TOROIDAL_SURFACE"), 2);
    assert_eq!(count("CONICAL_SURFACE"), 1);
}

#[test]
fn partial_revolution_of_an_off_axis_rectangle() {
    // Pappus: volume = angle · area · centroid distance from the axis.
    let rectangle = Profile::rectangle(Vector2::new(1.0, 0.0), Vector2::new(2.0, 1.0)).unwrap();
    let body = golf_model::revolve_by(&tilted(), &rectangle.into(), PI / 2.0).unwrap();
    assert_eq!(body.faces().len(), 6);
    assert_solid(&body, PI / 2.0 * 1.0 * 1.5, 0);
}

#[test]
fn partial_revolutions_touching_the_axis() {
    // A wedge of sphere, a three-quarter cone, and half a torus.
    let r = 1.5;
    let half_disc = Profile::new(vec![
        golf_model::Segment::Arc { center: Vector2::zeros(), radius: r, start_angle: -PI / 2.0, sweep: PI },
        golf_model::Segment::Line { start: Vector2::new(0.0, r), end: Vector2::new(0.0, -r) },
    ])
    .unwrap();
    let wedge = golf_model::revolve_by(&tilted(), &half_disc.into(), 2.0).unwrap();
    assert_solid(&wedge, 2.0 / (2.0 * PI) * 4.0 / 3.0 * PI * r.powi(3), 0);

    let triangle = Profile::polygon(&[Vector2::zeros(), Vector2::new(1.0, 0.0), Vector2::new(0.0, 2.0)]).unwrap();
    let cone = golf_model::revolve_by(&tilted(), &triangle.into(), 1.5 * PI).unwrap();
    assert_solid(&cone, 0.75 * PI * 1.0 * 2.0 / 3.0, 0);

    let tube = Profile::circle(Vector2::new(3.0, 0.0), 1.0).unwrap();
    let half_torus = golf_model::revolve_by(&tilted(), &tube.into(), PI).unwrap();
    assert_solid(&half_torus, PI * PI * 3.0, 0);
}

#[test]
fn partial_revolution_turns_a_hole_into_a_tunnel() {
    let outer = Profile::rectangle(Vector2::new(1.0, 0.0), Vector2::new(3.0, 4.0)).unwrap();
    let hole = Profile::circle(Vector2::new(2.0, 2.0), 0.5).unwrap();
    let body = golf_model::revolve_by(&tilted(), &Region::new(outer, vec![hole]), PI).unwrap();
    assert_eq!(body.shells().len(), 1);
    assert_solid(&body, PI * (8.0 * 2.0 - PI * 0.25 * 2.0), 1);
}

#[test]
fn revolution_angles_are_checked() {
    let square = Profile::rectangle(Vector2::new(1.0, 0.0), Vector2::new(2.0, 1.0)).unwrap();
    for angle in [0.0, -1.0, 7.0, f64::NAN] {
        let result = golf_model::revolve_by(&tilted(), &square.clone().into(), angle);
        assert!(matches!(result, Err(ModelError::BadAngle(_))), "{angle}");
    }
}
