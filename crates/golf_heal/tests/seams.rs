//! Faces whose loops wrap round their surface, as some files have them, are
//! cut open with seams.

use core::f64::consts::FRAC_PI_2;
use core::f64::consts::PI;
use core::f64::consts::TAU;

use golf_brep::Body;
use golf_brep::Coedge;
use golf_brep::Loop;
use golf_geom::Circle;
use golf_geom::Cylinder;
use golf_geom::Placement;
use golf_geom::Plane;
use golf_geom::Sphere;
use golf_heal::HealOptions;
use golf_heal::heal;
use golf_manifold::Point;
use golf_manifold::Vector;
use golf_manifold::World;
use golf_mesh::Tolerance;
use golf_mesh::mesh;
use nalgebra::Vector3;

fn frame(z: f64) -> Placement<World> {
    Placement::from_axes(
        Point::new(Vector3::new(0.0, 0.0, z)),
        Vector::new(Vector3::z()),
        Vector::new(Vector3::x()),
    )
}

/// A cylinder of radius `r` and height `h` with no seam: its side bounded by
/// two whole circles, with their vertices at angles `a` and `b`.
fn seamless_cylinder(r: f64, h: f64, a: f64, b: f64) -> Body<World> {
    let mut body = Body::new();
    let at = |z: f64, angle: f64| Point::new(Vector3::new(r * angle.cos(), r * angle.sin(), z));
    let (va, vb) = (
        body.add_vertex(at(0.0, a)).unwrap(),
        body.add_vertex(at(h, b)).unwrap(),
    );
    let bottom = body
        .add_edge(Circle::new(frame(0.0), r), (a, a + TAU), va, va)
        .unwrap();
    let top = body
        .add_edge(Circle::new(frame(h), r), (b, b + TAU), vb, vb)
        .unwrap();
    let side = body
        .add_face(
            Cylinder::new(frame(0.0), r),
            true,
            vec![
                Loop::new(vec![Coedge::new(bottom, false)]),
                Loop::new(vec![Coedge::new(top, true)]),
            ],
        )
        .unwrap();
    let floor = body
        .add_face(
            Plane::new(frame(0.0)),
            false,
            vec![Loop::new(vec![Coedge::new(bottom, true)])],
        )
        .unwrap();
    let lid = body
        .add_face(
            Plane::new(frame(h)),
            true,
            vec![Loop::new(vec![Coedge::new(top, false)])],
        )
        .unwrap();
    body.add_shell(vec![side, floor, lid]).unwrap();
    body
}

/// A hemisphere of radius `r`: its curved face bounded by the equator alone.
fn seamless_hemisphere(r: f64) -> Body<World> {
    let mut body = Body::new();
    let v = body
        .add_vertex(Point::new(Vector3::new(r, 0.0, 0.0)))
        .unwrap();
    let equator = body
        .add_edge(Circle::new(frame(0.0), r), (0.0, TAU), v, v)
        .unwrap();
    let dome = body
        .add_face(
            Sphere::new(Point::new(Vector3::zeros()), r),
            true,
            vec![Loop::new(vec![Coedge::new(equator, false)])],
        )
        .unwrap();
    let floor = body
        .add_face(
            Plane::new(frame(0.0)),
            false,
            vec![Loop::new(vec![Coedge::new(equator, true)])],
        )
        .unwrap();
    body.add_shell(vec![dome, floor]).unwrap();
    body
}

fn healed_volume(mut body: Body<World>, seams: usize, vertices_added: usize) -> f64 {
    let vertices = body.vertices().len();
    let report = heal(&mut body, &HealOptions::default());
    assert!(report.is_clean(), "{:?}", report.issues);
    assert_eq!(report.seams, seams);
    assert_eq!(body.vertices().len(), vertices + vertices_added);
    // Computed pcurves are as accurate as healing asks.
    assert_eq!(body.validate(HealOptions::default().tolerance), Ok(()));
    let m = mesh(&body, &Tolerance::new(1e-3, 0.3)).unwrap();
    assert!(m.is_watertight(), "{} open edges", m.open_edges().len());
    m.signed_volume()
}

#[test]
fn a_band_between_aligned_circles_gets_a_seam() {
    let volume = healed_volume(seamless_cylinder(2.0, 5.0, 0.3, 0.3), 1, 0);
    assert!((volume - PI * 4.0 * 5.0).abs() < 1e-3 * volume);
}

#[test]
fn a_band_between_unaligned_circles_splits_one() {
    // The seam leaves the bottom circle's vertex and splits the top circle.
    let volume = healed_volume(seamless_cylinder(2.0, 5.0, 0.3, 2.0), 1, 1);
    assert!((volume - PI * 4.0 * 5.0).abs() < 1e-3 * volume);
}

#[test]
fn a_cap_gets_a_seam_to_its_pole() {
    let volume = healed_volume(seamless_hemisphere(3.0), 1, 1);
    let expected = 2.0 / 3.0 * PI * 27.0;
    assert!(
        (volume - expected).abs() < 1e-3 * expected,
        "{volume} vs {expected}"
    );
    let _ = FRAC_PI_2;
}
