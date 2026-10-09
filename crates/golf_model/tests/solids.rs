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
        golf_model::Segment::Arc {
            center: Vector2::zeros(),
            radius: r,
            start_angle: -PI / 2.0,
            sweep: PI,
        },
        golf_model::Segment::Line {
            start: Vector2::new(0.0, r),
            end: Vector2::new(0.0, -r),
        },
    ])
    .unwrap();
    let wedge = golf_model::revolve_by(&tilted(), &half_disc.into(), 2.0).unwrap();
    assert_solid(&wedge, 2.0 / (2.0 * PI) * 4.0 / 3.0 * PI * r.powi(3), 0);

    let triangle = Profile::polygon(&[
        Vector2::zeros(),
        Vector2::new(1.0, 0.0),
        Vector2::new(0.0, 2.0),
    ])
    .unwrap();
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

#[test]
fn linear_and_circular_patterns() {
    let placement = Placement::<World>::at(Point::new(Vector3::new(5.0, 0.0, 0.0)));
    let peg = primitives::cylinder(&placement, 0.5, 2.0).unwrap();
    let row =
        golf_model::linear_pattern(&peg, Vector::new(Vector3::new(0.0, 2.0, 0.0)), 4).unwrap();
    assert_eq!(row.len(), 4);
    let centre = |b: &Body<World>| {
        let (sum, n) = b
            .vertices()
            .fold((Vector3::zeros(), 0.0), |(s, n), (_, v)| {
                (s + v.point.coords, n + 1.0)
            });
        sum / n
    };
    for (i, b) in row.iter().enumerate() {
        assert_eq!(b.validate(1e-9), Ok(()));
        assert!((centre(b) - centre(&peg) - Vector3::new(0.0, 2.0 * i as f64, 0.0)).norm() < 1e-12);
    }
    let ring = golf_model::circular_pattern(
        &peg,
        Point::new(Vector3::zeros()),
        Vector::new(Vector3::z()),
        6,
        2.0 * PI,
    )
    .unwrap();
    assert_eq!(ring.len(), 6);
    for (i, b) in ring.iter().enumerate() {
        assert_solid(b, PI * 0.25 * 2.0, 0);
        let angle = PI / 3.0 * i as f64;
        let c = centre(&peg);
        let expected = Vector3::new(
            c.x * angle.cos() - c.y * angle.sin(),
            c.x * angle.sin() + c.y * angle.cos(),
            c.z,
        );
        assert!((centre(b) - expected).norm() < 1e-9, "{i}");
    }
}

#[test]
fn loft_between_squares_is_a_frustum() {
    let (a, b, h) = (2.0, 1.0, 3.0);
    let bottom = Profile::rectangle(
        Vector2::new(-a / 2.0, -a / 2.0),
        Vector2::new(a / 2.0, a / 2.0),
    )
    .unwrap();
    let top = Profile::rectangle(
        Vector2::new(-b / 2.0, -b / 2.0),
        Vector2::new(b / 2.0, b / 2.0),
    )
    .unwrap();
    let base = tilted();
    let raised = Placement::new(base.point(Vector3::new(0.0, 0.0, h)), base.rotation);
    let body = golf_model::loft((&base, &bottom.into()), (&raised, &top.into())).unwrap();
    assert_solid(&body, h / 3.0 * (a * a + a * b + b * b), 0);
}

#[test]
fn loft_between_circles_is_exact() {
    // Coaxial: a frustum of a cone. Offset sideways: the same volume
    // (Cavalieri), through a skewed ruled surface.
    let (r1, r2, h) = (2.0, 1.0, 3.0);
    let base = tilted();
    for shift in [0.0, 1.5] {
        let raised = Placement::new(base.point(Vector3::new(shift, 0.0, h)), base.rotation);
        let body = golf_model::loft(
            (
                &base,
                &Profile::circle(Vector2::zeros(), r1).unwrap().into(),
            ),
            (
                &raised,
                &Profile::circle(Vector2::zeros(), r2).unwrap().into(),
            ),
        )
        .unwrap();
        assert_solid(&body, PI * h / 3.0 * (r1 * r1 + r1 * r2 + r2 * r2), 0);
    }
}

#[test]
fn twisted_loft_matches_the_prismatoid_formula() {
    // A square to the same square turned 45°: ruled, so the prismatoid formula
    // V = h/6 (A₀ + 4 A½ + A₁) is exact, the middle section being the polygon
    // through the midpoints of matching corners.
    let h = 2.0;
    let square = |turn: f64| -> Vec<Vector2<f64>> {
        (0..4)
            .map(|i| {
                let angle = turn + PI / 2.0 * i as f64 + PI / 4.0;
                Vector2::new(angle.cos(), angle.sin())
            })
            .collect()
    };
    let (lower, upper) = (square(0.0), square(PI / 4.0));
    let middle: Vec<Vector2<f64>> = lower
        .iter()
        .zip(&upper)
        .map(|(a, b)| (a + b) / 2.0)
        .collect();
    let area = |p: &[Vector2<f64>]| {
        (0..p.len())
            .map(|i| p[i].perp(&p[(i + 1) % p.len()]))
            .sum::<f64>()
            / 2.0
    };
    let volume = h / 6.0 * (area(&lower) + 4.0 * area(&middle) + area(&upper));
    let base = tilted();
    let raised = Placement::new(base.point(Vector3::new(0.0, 0.0, h)), base.rotation);
    let body = golf_model::loft(
        (&base, &Profile::polygon(&lower).unwrap().into()),
        (&raised, &Profile::polygon(&upper).unwrap().into()),
    )
    .unwrap();
    assert_solid(&body, volume, 0);
}

#[test]
fn loft_with_a_hole_and_arcs() {
    // A slot with a round hole, lofted to a smaller copy of itself.
    let shape = |scale: f64| {
        let slot = Profile::builder(Vector2::new(-1.0, -1.0) * scale)
            .line_to(Vector2::new(1.0, -1.0) * scale)
            .arc_to(
                Vector2::new(1.0, 1.0) * scale,
                Vector2::new(1.0, 0.0) * scale,
                true,
            )
            .line_to(Vector2::new(-1.0, 1.0) * scale)
            .arc_to(
                Vector2::new(-1.0, -1.0) * scale,
                Vector2::new(-1.0, 0.0) * scale,
                true,
            )
            .close()
            .unwrap();
        Region::new(
            slot,
            vec![Profile::circle(Vector2::zeros(), 0.5 * scale).unwrap()],
        )
    };
    let (h, s) = (2.0, 0.5);
    let base = tilted();
    let raised = Placement::new(base.point(Vector3::new(0.0, 0.0, h)), base.rotation);
    let body = golf_model::loft((&base, &shape(1.0)), (&raised, &shape(s))).unwrap();
    // Similar sections scaling linearly: V = h · A₀ · (1 + s + s²)/3.
    let area = 4.0 + PI - PI * 0.25;
    assert_solid(&body, h * area * (1.0 + s + s * s) / 3.0, 1);
    let mismatched = golf_model::loft(
        (
            &base,
            &Profile::circle(Vector2::zeros(), 1.0).unwrap().into(),
        ),
        (&raised, &shape(1.0)),
    );
    assert!(matches!(mismatched, Err(ModelError::LoftMismatch)));
}

mod treated {
    use golf_model::EdgeTreatment;
    use golf_model::ExtrudeEnds;
    use golf_model::extrude_drafted;
    use golf_model::extrude_with;

    use super::*;

    /// The area a quarter-round fillet of radius r removes from a right-angled
    /// corner, and its first moment about the outer side (area · centroid
    /// distance in from it).
    fn fillet_section(r: f64) -> (f64, f64) {
        (r * r * (1.0 - PI / 4.0), r.powi(3) * (5.0 / 6.0 - PI / 4.0))
    }

    #[test]
    fn box_with_a_filleted_top_meets_in_mitres() {
        // The rounding swept round the rectangle with mitred corners:
        // removed volume = section area · length of the path its centroid takes.
        let (a, b, h, r) = (4.0, 3.0, 2.0, 0.5);
        let rectangle = Profile::rectangle(Vector2::zeros(), Vector2::new(a, b)).unwrap();
        let ends = ExtrudeEnds {
            bottom: EdgeTreatment::None,
            top: EdgeTreatment::Fillet(r),
        };
        let body = extrude_with(&tilted(), &rectangle.into(), h, ends).unwrap();
        let (area, moment) = fillet_section(r);
        let removed = area * 2.0 * (a + b) - 8.0 * moment;
        assert_solid(&body, a * b * h - removed, 0);
    }

    #[test]
    fn box_with_a_chamfered_top() {
        let (a, b, h, d) = (4.0, 3.0, 2.0, 0.5);
        let rectangle = Profile::rectangle(Vector2::zeros(), Vector2::new(a, b)).unwrap();
        let ends = ExtrudeEnds {
            bottom: EdgeTreatment::None,
            top: EdgeTreatment::chamfer(d),
        };
        let body = extrude_with(&tilted(), &rectangle.into(), h, ends).unwrap();
        let slab = d / 6.0 * (a * b + 4.0 * (a - d) * (b - d) + (a - 2.0 * d) * (b - 2.0 * d));
        assert_solid(&body, a * b * (h - d) + slab, 0);
    }

    #[test]
    fn rounded_cuboid_has_spherical_corners() {
        // Steiner: the inner box grown by r.
        let (size, r) = (Vector3::new(4.0, 3.0, 2.0), 0.5);
        let body = primitives::rounded_cuboid(&tilted(), size, r).unwrap();
        let (p, q, s) = (size.x - 2.0 * r, size.y - 2.0 * r, size.z - 2.0 * r);
        let volume = p * q * s
            + 2.0 * (p * q + q * s + s * p) * r
            + PI * r * r * (p + q + s)
            + 4.0 / 3.0 * PI * r.powi(3);
        assert_eq!(body.faces().len(), 6 + 12 + 8);
        assert_solid(&body, volume, 0);
    }

    #[test]
    fn chamfered_cuboid() {
        // An octagon (the rectangle with chamfered corners) with chamfered ends.
        // Mitred insets of a polygon have area A − Pδ + cδ², c = Σ tan(turn/2).
        let (a, b, h, d) = (4.0, 3.0, 2.0, 0.4);
        let body = primitives::chamfered_cuboid(&tilted(), Vector3::new(a, b, h), d).unwrap();
        let (area, perimeter) = (
            a * b - 2.0 * d * d,
            2.0 * (a + b) - 8.0 * d + 4.0 * 2f64.sqrt() * d,
        );
        let c = 8.0 * (PI / 8.0).tan();
        let inset = |delta: f64| area - perimeter * delta + c * delta * delta;
        let slab = d / 6.0 * (inset(0.0) + 4.0 * inset(d / 2.0) + inset(d));
        assert_solid(&body, area * (h - 2.0 * d) + 2.0 * slab, 0);
    }

    #[test]
    fn puck_with_rounded_edges() {
        // Pappus: each fillet removes its section swept round the circle its
        // centroid is on.
        let (radius, h, r) = (3.0, 2.0, 0.5);
        let disc = Profile::circle(Vector2::zeros(), radius).unwrap();
        let body = extrude_with(
            &tilted(),
            &disc.into(),
            h,
            ExtrudeEnds::both(EdgeTreatment::Fillet(r)),
        )
        .unwrap();
        let (area, moment) = fillet_section(r);
        let removed = 2.0 * 2.0 * PI * (radius * area - moment);
        assert_solid(&body, PI * radius * radius * h - removed, 0);
    }

    #[test]
    fn ring_with_rounded_edges_inside_and_out() {
        // The hole's edge is concave: its fillets sweep tori outside the tube.
        let (outer, inner, h, r) = (3.0, 1.5, 2.0, 0.4);
        let region = Region::new(
            Profile::circle(Vector2::zeros(), outer).unwrap(),
            vec![Profile::circle(Vector2::zeros(), inner).unwrap()],
        );
        let body = extrude_with(
            &tilted(),
            &region,
            h,
            ExtrudeEnds::both(EdgeTreatment::Fillet(r)),
        )
        .unwrap();
        let (area, moment) = fillet_section(r);
        let removed = 2.0 * 2.0 * PI * ((outer * area - moment) + (inner * area + moment));
        assert_solid(&body, PI * (outer * outer - inner * inner) * h - removed, 1);
    }

    #[test]
    fn drafted_extrusions_are_frusta() {
        let (a, h, draft) = (4.0, 2.0, 0.2);
        let square = Profile::rectangle(Vector2::zeros(), Vector2::new(a, a)).unwrap();
        let body = extrude_drafted(&tilted(), &square.into(), h, draft).unwrap();
        let b = a - 2.0 * h * draft.tan();
        assert_solid(&body, h / 3.0 * (a * a + a * b + b * b), 0);

        let (r1, r2) = (2.0, 2.0 + h * 0.3f64.tan());
        let disc = Profile::circle(Vector2::zeros(), r1).unwrap();
        let body = extrude_drafted(&tilted(), &disc.into(), h, -0.3).unwrap();
        assert_solid(&body, PI * h / 3.0 * (r1 * r1 + r1 * r2 + r2 * r2), 0);
    }

    #[test]
    fn treatments_check_their_inputs() {
        let square =
            || Region::from(Profile::rectangle(Vector2::zeros(), Vector2::new(2.0, 2.0)).unwrap());
        let tall = ExtrudeEnds::both(EdgeTreatment::Fillet(0.6));
        assert_eq!(
            extrude_with(&tilted(), &square(), 1.0, tall).unwrap_err(),
            ModelError::TreatmentsTooTall
        );
        let wide = ExtrudeEnds {
            bottom: EdgeTreatment::None,
            top: EdgeTreatment::chamfer(1.5),
        };
        assert!(matches!(
            extrude_with(&tilted(), &square(), 3.0, wide),
            Err(ModelError::InsetCollapses { .. })
        ));
        // A sharp corner between a line and an arc.
        let d = Profile::builder(Vector2::zeros())
            .line_to(Vector2::new(2.0, 0.0))
            .arc_to(Vector2::zeros(), Vector2::new(1.0, 0.0), true)
            .close()
            .unwrap();
        let top = ExtrudeEnds {
            bottom: EdgeTreatment::None,
            top: EdgeTreatment::Fillet(0.1),
        };
        assert!(matches!(
            extrude_with(&tilted(), &d.into(), 1.0, top),
            Err(ModelError::UnsupportedCorner { .. })
        ));
        // A fillet between half the arc's radius and all of it.
        let disc = Region::from(Profile::circle(Vector2::zeros(), 1.0).unwrap());
        let top = ExtrudeEnds {
            bottom: EdgeTreatment::None,
            top: EdgeTreatment::Fillet(0.7),
        };
        assert!(matches!(
            extrude_with(&tilted(), &disc, 2.0, top),
            Err(ModelError::FilletTooLarge { .. })
        ));
    }
}

#[test]
fn revolving_about_any_axis_in_the_sketch() {
    // Pappus again: volume = angle · area · distance from the axis to the
    // centroid.
    let rectangle = || {
        Region::from(Profile::rectangle(Vector2::new(3.0, 0.0), Vector2::new(4.0, 1.0)).unwrap())
    };
    // About the upright line x = 1: the centroid is 2.5 away.
    let body = golf_model::revolve_about(
        &tilted(),
        &rectangle(),
        Vector2::new(1.0, 0.0),
        Vector2::y(),
        2.0 * PI,
    )
    .unwrap();
    assert_solid(&body, 2.0 * PI * 1.0 * 2.5, 1);
    // About a sloping line with the region on its left, through half a turn.
    let (point, direction): (Vector2<f64>, Vector2<f64>) =
        (Vector2::new(5.0, -2.0), Vector2::new(1.0, 2.0));
    let centroid = Vector2::new(3.5, 0.5);
    let distance = (centroid - point).perp(&direction.normalize()).abs();
    let body = golf_model::revolve_about(&tilted(), &rectangle(), point, direction, PI).unwrap();
    assert_solid(&body, PI * 1.0 * distance, 0);
}

mod pipes {
    use golf_model::Path;
    use golf_model::pipe;

    use super::*;

    #[test]
    fn straight_and_bent_pipes() {
        // A tube's volume is its section's area times its centre line's length
        // (Pappus), round bends and mitres alike.
        let r = 0.5;
        let straight = Profile::builder(Vector2::zeros())
            .line_to(Vector2::new(4.0, 0.0))
            .path()
            .unwrap();
        assert_solid(&pipe(&tilted(), &straight, r).unwrap(), PI * r * r * 4.0, 0);

        let bent = Profile::builder(Vector2::zeros())
            .line_to(Vector2::new(3.0, 0.0))
            .arc_to(Vector2::new(5.0, 2.0), Vector2::new(3.0, 2.0), true)
            .line_to(Vector2::new(5.0, 4.0))
            .arc_to(Vector2::new(7.0, 6.0), Vector2::new(7.0, 4.0), false)
            .path()
            .unwrap();
        assert_solid(
            &pipe(&tilted(), &bent, r).unwrap(),
            PI * r * r * bent.length(),
            0,
        );

        let mitred = Profile::builder(Vector2::zeros())
            .line_to(Vector2::new(3.0, 0.0))
            .line_to(Vector2::new(3.0, 3.0))
            .line_to(Vector2::new(6.0, 5.0))
            .path()
            .unwrap();
        assert_solid(
            &pipe(&tilted(), &mitred, r).unwrap(),
            PI * r * r * mitred.length(),
            0,
        );
    }

    #[test]
    fn closed_paths_make_rings() {
        let r = 0.4;
        let stadium =
            Path::from(Profile::slot(Vector2::zeros(), Vector2::new(4.0, 0.0), 2.0).unwrap());
        assert_solid(
            &pipe(&tilted(), &stadium, r).unwrap(),
            PI * r * r * stadium.length(),
            1,
        );
        let square =
            Path::from(Profile::rectangle(Vector2::zeros(), Vector2::new(3.0, 3.0)).unwrap());
        assert_solid(
            &pipe(&tilted(), &square, r).unwrap(),
            PI * r * r * square.length(),
            1,
        );
    }

    #[test]
    fn pipes_too_thick_for_their_path() {
        let tight = Profile::builder(Vector2::zeros())
            .line_to(Vector2::new(1.0, 0.0))
            .arc_to(Vector2::new(1.6, 0.6), Vector2::new(1.0, 0.6), true)
            .path()
            .unwrap();
        assert!(matches!(
            pipe(&tilted(), &tight, 0.7),
            Err(ModelError::PipeTooThick { index: 1 })
        ));
        let short = Profile::builder(Vector2::zeros())
            .line_to(Vector2::new(1.0, 0.0))
            .line_to(Vector2::new(1.0, 0.2))
            .line_to(Vector2::new(3.0, 0.2))
            .path()
            .unwrap();
        assert!(matches!(
            pipe(&tilted(), &short, 0.5),
            Err(ModelError::PipeTooThick { .. })
        ));
    }
}
