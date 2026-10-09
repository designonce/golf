//! Points, placements, curves and surfaces, in millimetres and radians.

use core::f64::consts::TAU;

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
use golf_manifold::Vector;
use golf_manifold::World;
use golf_step::Entity;
use golf_step::Id;
use golf_step::Record;
use golf_step::StepData;
use golf_step::Value;
use nalgebra::SVector;
use nalgebra::Vector3;

use crate::error::At;
use crate::error::Read;
use crate::error::invalid;
use crate::error::unsupported;
use crate::units::Units;

/// Reads geometry from `data` in one representation's units.
#[derive(Clone, Copy)]
pub(crate) struct Geometry<'a> {
    pub data: &'a StepData,
    pub units: Units,
}

/// An edge's curve, and the pcurves a file gives with it, each as the surface
/// it's on and its 2D curve.
pub(crate) struct EdgeCurve {
    pub curve: AnyCurve<World>,
    pub pcurves: Vec<(Id, Id)>,
}

/// A face's surface, and how to take the file's pcurves on it into golf's
/// parameters: each axis scaled, or `None` if they don't carry over.
pub(crate) struct FaceSurface {
    pub surface: AnySurface<World>,
    pub uv_scale: Option<[f64; 2]>,
}

impl Geometry<'_> {
    fn entity(&self, id: Id) -> Read<&Entity> {
        self.data.entity(id).at(id)
    }

    fn record(&self, id: Id, name: &str) -> Read<&Record> {
        self.data.record(id, name).at(id)
    }

    /// A `CARTESIAN_POINT` with `N` coordinates, scaled to millimetres.
    pub fn point<const N: usize>(&self, id: Id) -> Read<SVector<f64, N>> {
        let coords = self.record(id, "CARTESIAN_POINT")?.reals(1).at(id)?;
        if coords.len() < N {
            return Err(invalid(
                id,
                format!("{} coordinates, not {N}", coords.len()),
            ));
        }
        Ok(SVector::from_fn(|i, _| coords[i] * self.units.length))
    }

    /// A `DIRECTION`, normalised.
    pub fn direction<const N: usize>(&self, id: Id) -> Read<SVector<f64, N>> {
        let ratios = self.record(id, "DIRECTION")?.reals(1).at(id)?;
        if ratios.len() < N {
            return Err(invalid(id, format!("{} ratios, not {N}", ratios.len())));
        }
        SVector::<f64, N>::from_fn(|i, _| ratios[i])
            .try_normalize(1e-300)
            .ok_or_else(|| invalid(id, "a zero direction"))
    }

    /// A `VECTOR`: its direction times its magnitude, in millimetres.
    pub fn vector<const N: usize>(&self, id: Id) -> Read<SVector<f64, N>> {
        let record = self.record(id, "VECTOR")?;
        let direction = self.direction::<N>(record.reference(1).at(id)?)?;
        Ok(direction * record.real(2).at(id)? * self.units.length)
    }

    /// An `AXIS2_PLACEMENT_3D`, missing axes taken as Part 42 does.
    pub fn placement(&self, id: Id) -> Read<Placement<World>> {
        let record = self.record(id, "AXIS2_PLACEMENT_3D")?;
        let origin = self.point::<3>(record.reference(1).at(id)?)?;
        let z = match record.optional_reference(2).at(id)? {
            Some(axis) => self.direction::<3>(axis)?,
            None => Vector3::z(),
        };
        let x = match record.optional_reference(3).at(id)? {
            Some(reference) => self.direction::<3>(reference)?,
            None if (z - Vector3::x()).norm() > 1e-12 => Vector3::x(),
            None => Vector3::z(),
        };
        let x = (x - z * x.dot(&z))
            .try_normalize(1e-12)
            .ok_or_else(|| invalid(id, "its reference direction is along its axis"))?;
        Ok(Placement::from_axes(
            Point::new(origin),
            Vector::new(z),
            Vector::new(x),
        ))
    }

    /// An edge's curve, looking through `SURFACE_CURVE`s and the like to their
    /// 3D curve and pcurves.
    pub fn edge_curve(&self, id: Id) -> Read<EdgeCurve> {
        let entity = self.entity(id)?;
        for name in [
            "SURFACE_CURVE",
            "SEAM_CURVE",
            "INTERSECTION_CURVE",
            "BOUNDED_SURFACE_CURVE",
        ] {
            if let Some(record) = entity.record(name) {
                let curve = self.curve(record.reference(1).at(id)?)?;
                let mut pcurves = Vec::new();
                for associated in record.references(2).at(id)? {
                    let Ok(pcurve) = self.record(associated, "PCURVE") else {
                        continue;
                    };
                    let surface = pcurve.reference(1).at(associated)?;
                    // DEFINITIONAL_REPRESENTATION('', (curve), context).
                    let definition = pcurve.reference(2).at(associated)?;
                    let items = self
                        .entity(definition)?
                        .records
                        .iter()
                        .find_map(|r| r.references(1).ok())
                        .unwrap_or_default();
                    if let Some(&curve_2d) = items.first() {
                        pcurves.push((surface, curve_2d));
                    }
                }
                return Ok(EdgeCurve { curve, pcurves });
            }
        }
        Ok(EdgeCurve {
            curve: self.curve(id)?,
            pcurves: Vec::new(),
        })
    }

    /// A 3D curve.
    pub fn curve(&self, id: Id) -> Read<AnyCurve<World>> {
        let entity = self.entity(id)?;
        if let Some(record) = entity.simple() {
            match record.name.as_str() {
                "LINE" => {
                    let origin = self.point::<3>(record.reference(1).at(id)?)?;
                    let direction = self.vector::<3>(record.reference(2).at(id)?)?;
                    return Ok(Line::new(Point::new(origin), Vector::new(direction)).into());
                }
                "CIRCLE" => {
                    let placement = self.placement(record.reference(1).at(id)?)?;
                    let radius = record.real(2).at(id)? * self.units.length;
                    return Ok(Circle::new(placement, radius).into());
                }
                "ELLIPSE" => {
                    let placement = self.placement(record.reference(1).at(id)?)?;
                    let (a, b) = (record.real(2).at(id)?, record.real(3).at(id)?);
                    return Ok(Ellipse::new(
                        placement,
                        a * self.units.length,
                        b * self.units.length,
                    )
                    .into());
                }
                "TRIMMED_CURVE" => return self.curve(record.reference(1).at(id)?),
                "POLYLINE" => {
                    let points = record
                        .references(1)
                        .at(id)?
                        .into_iter()
                        .map(|p| self.point::<3>(p))
                        .collect::<Read<Vec<_>>>()?;
                    return Ok(NurbsCurve::new(polyline(id, points)?).into());
                }
                "SURFACE_CURVE" | "SEAM_CURVE" | "INTERSECTION_CURVE" | "BOUNDED_SURFACE_CURVE" => {
                    return self.curve(record.reference(1).at(id)?);
                }
                _ => {}
            }
        }
        if entity.is("B_SPLINE_CURVE")
            || entity.is("B_SPLINE_CURVE_WITH_KNOTS")
            || entity.is("BEZIER_CURVE")
            || entity.is("QUASI_UNIFORM_CURVE")
            || entity.is("UNIFORM_CURVE")
        {
            return Ok(NurbsCurve::new(self.b_spline_curve::<3>(entity)?).into());
        }
        Err(unsupported(id, format!("the curve {}", entity.name())))
    }

    /// A B-spline curve of any kind, simple or complex (rational).
    pub fn b_spline_curve<const N: usize>(
        &self,
        entity: &Entity,
    ) -> Read<golf_nurbs::NurbsCurve<N>> {
        let id = entity.id;
        // The inherited parameters: in one record for a simple instance, or
        // split across its parts for a complex one.
        let (base, base_at) = match entity.simple() {
            Some(record) => (record, 1),
            None => (entity.expect("B_SPLINE_CURVE").at(id)?, 0),
        };
        let degree = base.integer(base_at).at(id)? as usize;
        let points = base
            .references(base_at + 1)
            .at(id)?
            .into_iter()
            .map(|p| self.point::<N>(p))
            .collect::<Read<Vec<_>>>()?;
        let knots = match (entity.simple(), entity.record("B_SPLINE_CURVE_WITH_KNOTS")) {
            (Some(record), Some(_)) => expand(record.integers(6).at(id)?, record.reals(7).at(id)?),
            (None, Some(record)) => expand(record.integers(0).at(id)?, record.reals(1).at(id)?),
            _ => implicit_knots(entity, degree, points.len()),
        };
        let weights = match entity.record("RATIONAL_B_SPLINE_CURVE") {
            Some(record) => Some(record.reals(0).at(id)?),
            None => None,
        };
        golf_nurbs::NurbsCurve::new(degree, knots, points, weights)
            .map_err(|e| invalid(id, e.to_string()))
    }

    /// A surface, given points on the face it bounds (to size the part of an
    /// unbounded swept surface to make).
    pub fn surface(&self, id: Id, face_points: &[Vector3<f64>]) -> Read<FaceSurface> {
        let entity = self.entity(id)?;
        let length = self.units.length;
        let exact = |surface: AnySurface<World>, scale: [f64; 2]| FaceSurface {
            surface,
            uv_scale: (self.units.angle == 1.0).then_some(scale),
        };
        if let Some(record) = entity.simple() {
            let placement = || self.placement(record.reference(1).at(id)?);
            match record.name.as_str() {
                "PLANE" => return Ok(exact(Plane::new(placement()?).into(), [length, length])),
                "CYLINDRICAL_SURFACE" => {
                    let radius = record.real(2).at(id)? * length;
                    return Ok(exact(
                        Cylinder::new(placement()?, radius).into(),
                        [1.0, length],
                    ));
                }
                "CONICAL_SURFACE" => {
                    let radius = record.real(2).at(id)? * length;
                    let half_angle = record.real(3).at(id)? * self.units.angle;
                    return Ok(exact(
                        Cone::new(placement()?, radius, half_angle).into(),
                        [1.0, length],
                    ));
                }
                "SPHERICAL_SURFACE" => {
                    let radius = record.real(2).at(id)? * length;
                    return Ok(exact(
                        Sphere::from_placement(placement()?, radius).into(),
                        [1.0, 1.0],
                    ));
                }
                "TOROIDAL_SURFACE" | "DEGENERATE_TOROIDAL_SURFACE" => {
                    let (major, minor) = (
                        record.real(2).at(id)? * length,
                        record.real(3).at(id)? * length,
                    );
                    if major > minor {
                        return Ok(exact(
                            Torus::new(placement()?, major, minor).into(),
                            [1.0, 1.0],
                        ));
                    }
                    // A spindle or horn torus: the part of its tube outside
                    // its axis (or, if the file selects it, inside), revolved.
                    let outer = match record.name.as_str() {
                        "DEGENERATE_TOROIDAL_SURFACE" => record.boolean(4).at(id)?,
                        _ => true,
                    };
                    let p = placement()?;
                    let (origin, x, z) = (
                        p.origin.coords,
                        p.vector(Vector3::x()).coords,
                        p.vector(Vector3::z()).coords,
                    );
                    // Where the tube crosses the axis, as an angle round it.
                    let crossing = (-major / minor).clamp(-1.0, 1.0).acos();
                    let (from, to) = match outer {
                        true => (-crossing, crossing),
                        false => (crossing, TAU - crossing),
                    };
                    let tube = golf_nurbs::NurbsCurve::circular_arc(
                        origin + x * major,
                        x,
                        z,
                        minor,
                        from,
                        to,
                    );
                    let surface = revolve(id, &tube, origin, z)?;
                    return Ok(FaceSurface {
                        surface: NurbsSurface::new(surface).into(),
                        uv_scale: None,
                    });
                }
                "RECTANGULAR_TRIMMED_SURFACE" => {
                    return self.surface(record.reference(1).at(id)?, face_points);
                }
                "SURFACE_OF_LINEAR_EXTRUSION" => return self.extrusion(record, id, face_points),
                "SURFACE_OF_REVOLUTION" => return self.revolution(record, id, face_points),
                _ => {}
            }
        }
        if entity.is("B_SPLINE_SURFACE")
            || entity.is("B_SPLINE_SURFACE_WITH_KNOTS")
            || entity.is("BEZIER_SURFACE")
            || entity.is("QUASI_UNIFORM_SURFACE")
            || entity.is("UNIFORM_SURFACE")
        {
            let surface = self.b_spline_surface(entity)?;
            return Ok(exact(NurbsSurface::new(surface).into(), [1.0, 1.0]));
        }
        Err(unsupported(id, format!("the surface {}", entity.name())))
    }

    fn b_spline_surface(&self, entity: &Entity) -> Read<golf_nurbs::NurbsSurface<3>> {
        let id = entity.id;
        let (base, at) = match entity.simple() {
            Some(record) => (record, 1),
            None => (entity.expect("B_SPLINE_SURFACE").at(id)?, 0),
        };
        let (degree_u, degree_v) = (
            base.integer(at).at(id)? as usize,
            base.integer(at + 1).at(id)? as usize,
        );
        let rows = base.list(at + 2).at(id)?;
        let points = rows
            .iter()
            .map(|row| {
                row.as_list()
                    .ok_or_else(|| invalid(id, "control points aren't a list of lists"))?
                    .iter()
                    .map(|p| match p {
                        Value::Ref(p) => self.point::<3>(*p),
                        _ => Err(invalid(id, "a control point isn't a reference")),
                    })
                    .collect::<Read<Vec<_>>>()
            })
            .collect::<Read<Vec<_>>>()?;
        let (count_u, count_v) = (points.len(), points.first().map_or(0, Vec::len));
        let knots = match (
            entity.simple(),
            entity.record("B_SPLINE_SURFACE_WITH_KNOTS"),
        ) {
            (Some(record), Some(_)) => (
                expand(record.integers(8).at(id)?, record.reals(10).at(id)?),
                expand(record.integers(9).at(id)?, record.reals(11).at(id)?),
            ),
            (None, Some(record)) => (
                expand(record.integers(0).at(id)?, record.reals(2).at(id)?),
                expand(record.integers(1).at(id)?, record.reals(3).at(id)?),
            ),
            _ => (
                implicit_knots(entity, degree_u, count_u),
                implicit_knots(entity, degree_v, count_v),
            ),
        };
        let weights = match entity.record("RATIONAL_B_SPLINE_SURFACE") {
            Some(record) => Some(
                record
                    .list(0)
                    .at(id)?
                    .iter()
                    .map(|row| {
                        row.as_list()
                            .ok_or_else(|| invalid(id, "weights aren't a list of lists"))?
                            .iter()
                            .map(|w| {
                                w.as_real()
                                    .ok_or_else(|| invalid(id, "a weight isn't a number"))
                            })
                            .collect::<Read<Vec<_>>>()
                    })
                    .collect::<Read<Vec<_>>>()?,
            ),
            None => None,
        };
        golf_nurbs::NurbsSurface::new((degree_u, knots.0), (degree_v, knots.1), points, weights)
            .map_err(|e| invalid(id, e.to_string()))
    }

    /// The profile of a swept surface as a NURBS curve, and whether its
    /// parameter is the file's. A line is cut to cover `face_points`.
    fn profile(
        &self,
        id: Id,
        face_points: &[Vector3<f64>],
    ) -> Read<(golf_nurbs::NurbsCurve<3>, bool)> {
        match self.curve(id)? {
            AnyCurve::Nurbs(nurbs) => Ok((nurbs.geometry().clone(), true)),
            AnyCurve::Line(line) => {
                let (o, d) = (line.origin.coords, line.direction.coords);
                let (lo, hi) = extent(face_points, |p| (p - o).dot(&d) / d.norm_squared());
                let line = golf_nurbs::NurbsCurve::new(
                    1,
                    vec![lo, lo, hi, hi],
                    vec![o + d * lo, o + d * hi],
                    None,
                )
                .map_err(|e| invalid(id, e.to_string()))?;
                Ok((line, true))
            }
            AnyCurve::Circle(circle) => {
                let p = &circle.placement;
                let (x, y) = (p.vector(Vector3::x()).coords, p.vector(Vector3::y()).coords);
                let arc = golf_nurbs::NurbsCurve::circular_arc(
                    p.origin.coords,
                    x,
                    y,
                    circle.radius,
                    0.0,
                    TAU,
                );
                Ok((arc, false))
            }
            _ => Err(unsupported(id, "this profile for a swept surface")),
        }
    }

    /// `SURFACE_OF_LINEAR_EXTRUSION`: the profile swept along a vector, as a
    /// NURBS surface ruled in v.
    fn extrusion(
        &self,
        record: &Record,
        id: Id,
        face_points: &[Vector3<f64>],
    ) -> Read<FaceSurface> {
        let (profile, kept) = self.profile(record.reference(1).at(id)?, face_points)?;
        let along = self.vector::<3>(record.reference(2).at(id)?)?;
        let start = profile.point(profile.domain().0);
        let (lo, hi) = extent(face_points, |p| {
            (p - start).dot(&along) / along.norm_squared()
        });
        let count = profile.control_point_count();
        let points = (0..count)
            .map(|i| {
                vec![
                    profile.control_point(i) + along * lo,
                    profile.control_point(i) + along * hi,
                ]
            })
            .collect();
        let weights = (0..count).map(|i| vec![profile.weight(i); 2]).collect();
        let surface = golf_nurbs::NurbsSurface::new(
            (profile.degree(), profile.knots().knots().to_vec()),
            (1, vec![lo, lo, hi, hi]),
            points,
            Some(weights),
        )
        .map_err(|e| invalid(id, e.to_string()))?;
        Ok(FaceSurface {
            surface: NurbsSurface::new(surface).into(),
            uv_scale: (kept && self.units.angle == 1.0).then_some([1.0, 1.0]),
        })
    }

    /// `SURFACE_OF_REVOLUTION`: the profile turned a full turn about an axis,
    /// as a rational NURBS surface (u along the profile, v round the axis,
    /// though not at the file's rate, so its pcurves don't carry over).
    fn revolution(
        &self,
        record: &Record,
        id: Id,
        face_points: &[Vector3<f64>],
    ) -> Read<FaceSurface> {
        let (profile, _) = self.profile(record.reference(1).at(id)?, face_points)?;
        let axis_id = record.reference(2).at(id)?;
        let axis = self.record(axis_id, "AXIS1_PLACEMENT")?;
        let origin = self.point::<3>(axis.reference(1).at(axis_id)?)?;
        let z = match axis.optional_reference(2).at(axis_id)? {
            Some(d) => self.direction::<3>(d)?,
            None => Vector3::z(),
        };
        let surface = revolve(id, &profile, origin, z)?;
        Ok(FaceSurface {
            surface: NurbsSurface::new(surface).into(),
            uv_scale: None,
        })
    }

    /// A pcurve's 2D curve, in the file's surface parameters (scaled after),
    /// if it's a kind that keeps its edge's parameter.
    pub fn curve_2d(&self, id: Id) -> Read<Curve2> {
        let entity = self.entity(id)?;
        if let Some(record) = entity.simple() {
            if record.name == "LINE" {
                let origin = self.raw_point::<2>(record.reference(1).at(id)?)?;
                let vector = self.record(record.reference(2).at(id)?, "VECTOR")?;
                let vector_id = record.reference(2).at(id)?;
                let direction = self.direction::<2>(vector.reference(1).at(vector_id)?)?;
                let magnitude = vector.real(2).at(vector_id)?;
                return Ok(Curve2::Line(origin, direction * magnitude));
            }
        }
        if entity.is("B_SPLINE_CURVE")
            || entity.is("B_SPLINE_CURVE_WITH_KNOTS")
            || entity.is("BEZIER_CURVE")
            || entity.is("QUASI_UNIFORM_CURVE")
            || entity.is("UNIFORM_CURVE")
        {
            // Parameters, not lengths: no unit scaling here.
            let raw = Geometry {
                data: self.data,
                units: Units::default(),
            };
            return Ok(Curve2::Nurbs(raw.b_spline_curve::<2>(entity)?));
        }
        Err(unsupported(id, format!("the pcurve {}", entity.name())))
    }

    fn raw_point<const N: usize>(&self, id: Id) -> Read<SVector<f64, N>> {
        Geometry {
            data: self.data,
            units: Units::default(),
        }
        .point::<N>(id)
    }
}

/// A file's 2D curve, before scaling into golf's surface parameters.
pub(crate) enum Curve2 {
    /// `origin + t·direction`.
    Line(nalgebra::Vector2<f64>, nalgebra::Vector2<f64>),
    Nurbs(golf_nurbs::NurbsCurve<2>),
}

/// `profile` turned a full turn about the axis through `origin` along unit
/// `z`, as a rational NURBS surface: u along the profile, v round the axis.
fn revolve(
    id: Id,
    profile: &golf_nurbs::NurbsCurve<3>,
    origin: Vector3<f64>,
    z: Vector3<f64>,
) -> Read<golf_nurbs::NurbsSurface<3>> {
    // A fixed frame about the axis, x towards the profile.
    let count = profile.control_point_count();
    let radial = |p: Vector3<f64>| {
        let d = p - origin;
        d - z * d.dot(&z)
    };
    let x = (0..count)
        .map(|i| radial(profile.control_point(i)))
        .find(|r| r.norm() > 1e-9)
        .and_then(|r| r.try_normalize(1e-300))
        .ok_or_else(|| invalid(id, "its profile lies on its axis"))?;
    let y = z.cross(&x);
    let unit = golf_nurbs::NurbsCurve::<2>::circular_arc(
        nalgebra::Vector2::zeros(),
        nalgebra::Vector2::x(),
        nalgebra::Vector2::y(),
        1.0,
        0.0,
        TAU,
    );
    let arc_count = unit.control_point_count();
    let mut points = Vec::with_capacity(count);
    let mut weights = Vec::with_capacity(count);
    for i in 0..count {
        let d = profile.control_point(i) - origin;
        let height = d.dot(&z);
        let (a, b) = (d.dot(&x), d.dot(&y));
        let centre = origin + z * height;
        // Rotating (a, b) is linear in (cos, sin), so the arc's control
        // points rotate it exactly.
        points.push(
            (0..arc_count)
                .map(|j| {
                    let c = unit.control_point(j);
                    centre + x * (a * c.x - b * c.y) + y * (a * c.y + b * c.x)
                })
                .collect(),
        );
        weights.push(
            (0..arc_count)
                .map(|j| profile.weight(i) * unit.weight(j))
                .collect(),
        );
    }
    golf_nurbs::NurbsSurface::new(
        (profile.degree(), profile.knots().knots().to_vec()),
        (2, unit.knots().knots().to_vec()),
        points,
        Some(weights),
    )
    .map_err(|e| invalid(id, e.to_string()))
}

/// Knots from distinct values and their multiplicities.
fn expand(multiplicities: Vec<i64>, values: Vec<f64>) -> Vec<f64> {
    multiplicities
        .into_iter()
        .zip(values)
        .flat_map(|(m, k)| core::iter::repeat_n(k, m.max(0) as usize))
        .collect()
}

/// The knots of a B-spline whose kind implies them: a Bézier (each piece of
/// `degree` sharing ends), quasi-uniform (clamped, unit spacing) or uniform.
fn implicit_knots(entity: &Entity, degree: usize, count: usize) -> Vec<f64> {
    let total = count + degree + 1;
    if entity.is("BEZIER_CURVE") || entity.is("BEZIER_SURFACE") {
        let pieces = (count - 1) / degree.max(1);
        let mut knots = vec![0.0; degree + 1];
        for i in 1..pieces {
            knots.extend(core::iter::repeat_n(i as f64, degree));
        }
        knots.extend(core::iter::repeat_n(pieces as f64, degree + 1));
        knots
    } else if entity.is("QUASI_UNIFORM_CURVE") || entity.is("QUASI_UNIFORM_SURFACE") {
        let interior = total - 2 * (degree + 1);
        let mut knots = vec![0.0; degree + 1];
        knots.extend((1..=interior).map(|i| i as f64));
        knots.extend(core::iter::repeat_n((interior + 1) as f64, degree + 1));
        knots
    } else {
        (0..total).map(|i| i as f64 - degree as f64).collect()
    }
}

/// Points joined by straight segments, segment `i` over `[i, i + 1]`.
fn polyline(id: Id, points: Vec<Vector3<f64>>) -> Read<golf_nurbs::NurbsCurve<3>> {
    let n = points.len();
    let mut knots = vec![0.0];
    knots.extend((0..n).map(|i| i as f64));
    knots.push((n - 1) as f64);
    golf_nurbs::NurbsCurve::new(1, knots, points, None).map_err(|e| invalid(id, e.to_string()))
}

/// The range of `measure` over `points`, widened by half its span (and at
/// least a little), to cover a face with room to spare.
fn extent(points: &[Vector3<f64>], measure: impl Fn(Vector3<f64>) -> f64) -> (f64, f64) {
    let values: Vec<f64> = points.iter().map(|&p| measure(p)).collect();
    let lo = values.iter().copied().fold(f64::INFINITY, f64::min);
    let hi = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !(lo.is_finite() && hi.is_finite()) {
        return (-1.0, 1.0);
    }
    let margin = ((hi - lo) / 2.0).max(1e-3 * (1.0 + lo.abs().max(hi.abs())));
    (lo - margin, hi + margin)
}
