//! Curves and surfaces as STEP geometry entities.

use golf_geom::AnyCurve;
use golf_geom::AnySurface;
use golf_geom::Circle;
use golf_geom::Cone;
use golf_geom::Cylinder;
use golf_geom::Ellipse;
use golf_geom::Line;
use golf_geom::NurbsCurve;
use golf_geom::NurbsSurface;
use golf_geom::Placement;
use golf_geom::Plane;
use golf_geom::Sphere;
use golf_geom::Torus;
use golf_manifold::Point;
use golf_manifold::Space;
use golf_nurbs::KnotVector;
use nalgebra::UnitQuaternion;
use nalgebra::Vector3;

use crate::error::StepError;
use crate::writer::Ref;
use crate::writer::Writer;
use crate::writer::logical;
use crate::writer::real;
use crate::writer::reals;
use crate::writer::refs;

/// A curve that can be written as a STEP curve entity.
pub trait StepCurve {
    fn write(&self, w: &mut Writer) -> Result<Ref, StepError>;
}

/// A surface that can be written as a STEP surface entity whose normal is this
/// surface's.
pub trait StepSurface {
    fn write(&self, w: &mut Writer) -> Result<Ref, StepError>;
}

impl<S: Space<3>> StepCurve for Line<S, 3> {
    fn write(&self, w: &mut Writer) -> Result<Ref, StepError> {
        let point = w.point(self.origin.coords)?;
        let direction = w.direction(self.direction.coords)?;
        let magnitude = real(self.direction.coords.norm())?;
        let vector = w.add(format!("VECTOR('',{direction},{magnitude})"));
        Ok(w.add(format!("LINE('',{point},{vector})")))
    }
}

impl<S: Space<3>> StepCurve for Circle<S> {
    fn write(&self, w: &mut Writer) -> Result<Ref, StepError> {
        let position = w.placement(&self.placement)?;
        Ok(w.add(format!("CIRCLE('',{position},{})", real(self.radius)?)))
    }
}

impl<S: Space<3>> StepCurve for Ellipse<S> {
    fn write(&self, w: &mut Writer) -> Result<Ref, StepError> {
        let position = w.placement(&self.placement)?;
        Ok(w.add(format!(
            "ELLIPSE('',{position},{},{})",
            real(self.major_radius)?,
            real(self.minor_radius)?
        )))
    }
}

impl<S: Space<3>> StepCurve for NurbsCurve<S, 3> {
    fn write(&self, w: &mut Writer) -> Result<Ref, StepError> {
        let curve = self.geometry();
        let points = curve
            .control_points()
            .map(|p| w.point(p))
            .collect::<Result<Vec<_>, _>>()?;
        let degree = curve.degree();
        let (multiplicities, knots) = knot_lists(curve.knots())?;
        let points = refs(points);
        let count = curve.control_point_count();
        let closed = logical(same(curve.control_point(0), curve.control_point(count - 1)));
        if !curve.is_rational() {
            return Ok(w.add(format!(
                "B_SPLINE_CURVE_WITH_KNOTS('',{degree},{points},.UNSPECIFIED.,{closed},.F.,{multiplicities},{knots},.UNSPECIFIED.)"
            )));
        }
        let weights = reals((0..count).map(|i| curve.weight(i)))?;
        Ok(w.add(format!(
            "(BOUNDED_CURVE()B_SPLINE_CURVE({degree},{points},.UNSPECIFIED.,{closed},.F.)\
             B_SPLINE_CURVE_WITH_KNOTS({multiplicities},{knots},.UNSPECIFIED.)CURVE()\
             GEOMETRIC_REPRESENTATION_ITEM()RATIONAL_B_SPLINE_CURVE({weights})REPRESENTATION_ITEM(''))"
        )))
    }
}

impl<S: Space<3>> StepCurve for AnyCurve<S> {
    fn write(&self, w: &mut Writer) -> Result<Ref, StepError> {
        match self {
            Self::Line(c) => c.write(w),
            Self::Circle(c) => c.write(w),
            Self::Ellipse(c) => c.write(w),
            Self::Nurbs(c) => c.write(w),
            other => Err(StepError::Unsupported(format!("curve {other:?}"))),
        }
    }
}

impl<S: Space<3>> StepSurface for Plane<S> {
    fn write(&self, w: &mut Writer) -> Result<Ref, StepError> {
        let position = w.placement(&self.placement)?;
        Ok(w.add(format!("PLANE('',{position})")))
    }
}

impl<S: Space<3>> StepSurface for Cylinder<S> {
    fn write(&self, w: &mut Writer) -> Result<Ref, StepError> {
        let position = w.placement(&self.placement)?;
        Ok(w.add(format!(
            "CYLINDRICAL_SURFACE('',{position},{})",
            real(self.radius)?
        )))
    }
}

impl<S: Space<3>> StepSurface for Cone<S> {
    /// STEP cones have a non-negative radius at their placement and a positive
    /// half angle, so this one is placed at its apex (radius 0), its axis
    /// flipped if its half angle is negative. Both keep the surface, its
    /// parametrisation's handedness and its outward normal.
    fn write(&self, w: &mut Writer) -> Result<Ref, StepError> {
        let Some(apex) = self.apex_height() else {
            return Cylinder::new(self.placement, self.radius).write(w);
        };
        let origin: Point<S, 3> = self.placement.point(Vector3::new(0.0, 0.0, apex));
        let mut rotation = self.placement.rotation;
        if self.half_angle < 0.0 {
            // A half turn about local x: z and y reverse, x stays.
            rotation *= UnitQuaternion::from_axis_angle(&Vector3::x_axis(), core::f64::consts::PI);
        }
        let position = w.placement(&Placement::new(origin, rotation))?;
        Ok(w.add(format!(
            "CONICAL_SURFACE('',{position},0.,{})",
            real(self.half_angle.abs())?
        )))
    }
}

impl<S: Space<3>> StepSurface for Sphere<S> {
    fn write(&self, w: &mut Writer) -> Result<Ref, StepError> {
        let position = w.placement(&self.placement)?;
        Ok(w.add(format!(
            "SPHERICAL_SURFACE('',{position},{})",
            real(self.radius)?
        )))
    }
}

impl<S: Space<3>> StepSurface for Torus<S> {
    fn write(&self, w: &mut Writer) -> Result<Ref, StepError> {
        let position = w.placement(&self.placement)?;
        Ok(w.add(format!(
            "TOROIDAL_SURFACE('',{position},{},{})",
            real(self.major_radius)?,
            real(self.minor_radius)?
        )))
    }
}

impl<S: Space<3>> StepSurface for NurbsSurface<S> {
    fn write(&self, w: &mut Writer) -> Result<Ref, StepError> {
        let surface = self.geometry();
        let (count_u, count_v) = surface.control_grid_size();
        let (degree_u, degree_v) = surface.degrees();
        let mut rows = Vec::with_capacity(count_u);
        for i in 0..count_u {
            let row = (0..count_v)
                .map(|j| w.point(surface.control_point(i, j)))
                .collect::<Result<Vec<_>, _>>()?;
            rows.push(refs(row));
        }
        let points = format!("({})", rows.join(","));
        // Closed in u if the first and last rows of control points coincide;
        // readers need it to give a seam edge its two pcurves.
        let u_closed = logical((0..count_v).all(|j| {
            same(
                surface.control_point(0, j),
                surface.control_point(count_u - 1, j),
            )
        }));
        let v_closed = logical((0..count_u).all(|i| {
            same(
                surface.control_point(i, 0),
                surface.control_point(i, count_v - 1),
            )
        }));
        let (multiplicities_u, knots_u) = knot_lists(surface.knots_u())?;
        let (multiplicities_v, knots_v) = knot_lists(surface.knots_v())?;
        if !surface.is_rational() {
            return Ok(w.add(format!(
                "B_SPLINE_SURFACE_WITH_KNOTS('',{degree_u},{degree_v},{points},.UNSPECIFIED.,{u_closed},{v_closed},.F.,\
                 {multiplicities_u},{multiplicities_v},{knots_u},{knots_v},.UNSPECIFIED.)"
            )));
        }
        let weights = (0..count_u)
            .map(|i| reals((0..count_v).map(|j| surface.weight(i, j))))
            .collect::<Result<Vec<_>, _>>()?;
        let weights = format!("({})", weights.join(","));
        Ok(w.add(format!(
            "(BOUNDED_SURFACE()B_SPLINE_SURFACE({degree_u},{degree_v},{points},.UNSPECIFIED.,{u_closed},{v_closed},.F.)\
             B_SPLINE_SURFACE_WITH_KNOTS({multiplicities_u},{multiplicities_v},{knots_u},{knots_v},.UNSPECIFIED.)\
             GEOMETRIC_REPRESENTATION_ITEM()RATIONAL_B_SPLINE_SURFACE({weights})REPRESENTATION_ITEM('')SURFACE())"
        )))
    }
}

impl<S: Space<3>> StepSurface for AnySurface<S> {
    fn write(&self, w: &mut Writer) -> Result<Ref, StepError> {
        match self {
            Self::Plane(s) => s.write(w),
            Self::Sphere(s) => s.write(w),
            Self::Cylinder(s) => s.write(w),
            Self::Cone(s) => s.write(w),
            Self::Torus(s) => s.write(w),
            Self::Nurbs(s) => s.write(w),
            other => Err(StepError::Unsupported(format!("surface {other:?}"))),
        }
    }
}

/// A knot vector as STEP's lists of multiplicities and distinct knots.
fn knot_lists(knots: &KnotVector) -> Result<(String, String), StepError> {
    let (values, multiplicities): (Vec<f64>, Vec<usize>) = knots.distinct().unzip();
    let multiplicities: Vec<String> = multiplicities.iter().map(usize::to_string).collect();
    Ok((format!("({})", multiplicities.join(",")), reals(values)?))
}

/// Whether two control points are the same point, to rounding.
fn same(a: Vector3<f64>, b: Vector3<f64>) -> bool {
    (a - b).norm() <= 1e-12 * (1.0 + a.norm().max(b.norm()))
}

#[cfg(test)]
mod tests {
    use golf_manifold::Vector;
    use golf_manifold::World;

    use super::*;

    fn text(write: impl FnOnce(&mut Writer) -> Result<Ref, StepError>) -> String {
        let mut w = Writer::default();
        write(&mut w).unwrap();
        w.finish("", "")
    }

    #[test]
    fn cone_with_negative_half_angle_is_written_from_its_apex_axis_flipped() {
        let placement = Placement::<World>::at(Point::new(Vector3::zeros()));
        let cone = Cone::new(placement, 2.0, -0.3);
        let out = text(|w| cone.write(w));
        // Apex at height radius / tan(0.3) above the origin, axis pointing down.
        let apex = 2.0 / 0.3f64.tan();
        assert!(
            out.contains(&format!(
                "CARTESIAN_POINT('',(0.0,0.0,{}))",
                real(apex).unwrap()
            )),
            "{out}"
        );
        let axis: Vec<f64> = out
            .lines()
            .find_map(|l| l.strip_prefix("#2=DIRECTION('',("))
            .and_then(|l| l.strip_suffix("));"))
            .expect("the axis is entity #2")
            .split(',')
            .map(|x| x.parse().unwrap())
            .collect();
        assert!(
            (Vector3::from_vec(axis) - Vector3::new(0.0, 0.0, -1.0)).norm() < 1e-12,
            "{out}"
        );
        assert!(
            out.contains(&format!("CONICAL_SURFACE('',#4,0.,{})", real(0.3).unwrap())),
            "{out}"
        );
    }

    #[test]
    fn closed_b_splines_say_so() {
        let full = golf_nurbs::NurbsCurve::circular_arc(
            Vector3::zeros(),
            Vector3::x(),
            Vector3::y(),
            1.0,
            0.0,
            0.0,
        );
        let half = golf_nurbs::NurbsCurve::circular_arc(
            Vector3::zeros(),
            Vector3::x(),
            Vector3::y(),
            1.0,
            0.0,
            3.0,
        );
        let closed = text(|w| NurbsCurve::<World, 3>::new(full).write(w));
        let open = text(|w| NurbsCurve::<World, 3>::new(half).write(w));
        assert!(closed.contains(".UNSPECIFIED.,.T.,.F.)"), "{closed}");
        assert!(open.contains(".UNSPECIFIED.,.F.,.F.)"), "{open}");
    }

    #[test]
    fn lines_keep_their_parameter_scale() {
        let line = Line::<World, 3>::new(
            Point::new(Vector3::new(1.0, 2.0, 3.0)),
            Vector::new(Vector3::new(0.0, 3.0, 4.0)),
        );
        let out = text(|w| line.write(w));
        assert!(out.contains("DIRECTION('',(0.0,0.6,0.8))"), "{out}");
        assert!(out.contains("VECTOR('',#2,5.0)"), "{out}");
    }
}
