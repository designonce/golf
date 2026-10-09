//! NURBS curves and surfaces from `golf_nurbs`, placed in a space.

use golf_manifold::Axis;
use golf_manifold::Domain;
use golf_manifold::Embedding;
use golf_manifold::HasT;
use golf_manifold::HasUv;
use golf_manifold::Mapping;
use golf_manifold::Point;
use golf_manifold::ProjectError;
use golf_manifold::Space;
use golf_manifold::T;
use golf_manifold::Uv;
use golf_manifold::newton_project;
use nalgebra::Matrix3x2;
use nalgebra::SMatrix;
use nalgebra::SVector;

/// A NURBS curve in an `N`-dimensional space `S`.
///
/// Its domain is the knot vector's, periodic if the curve is closed.
#[derive(Clone, Debug)]
pub struct NurbsCurve<S: Space<N>, const N: usize> {
    geometry: golf_nurbs::NurbsCurve<N>,
    tag: S::Tag,
    domain: Domain<1>,
}

impl<S: Space<N, Tag = ()>, const N: usize> NurbsCurve<S, N> {
    pub fn new(geometry: golf_nurbs::NurbsCurve<N>) -> Self {
        Self::with_tag(geometry, ())
    }
}

impl<S: Space<N>, const N: usize> NurbsCurve<S, N> {
    /// `geometry` in the instance of `S` that `tag` names.
    pub fn with_tag(geometry: golf_nurbs::NurbsCurve<N>, tag: S::Tag) -> Self {
        let (start, end) = geometry.domain();
        let tolerance = tolerance(geometry.control_points());
        let closed = (geometry.point(start) - geometry.point(end)).norm() <= tolerance;
        let axis = match closed {
            true => Axis::periodic(start, end),
            false => Axis::bounded(start, end),
        };
        Self {
            geometry,
            tag,
            domain: Domain::new([axis]),
        }
    }

    pub fn geometry(&self) -> &golf_nurbs::NurbsCurve<N> {
        &self.geometry
    }

    pub fn tag(&self) -> S::Tag {
        self.tag
    }

    /// Parameters spread through every knot span, for seeding projection.
    fn samples(&self) -> impl Iterator<Item = f64> + '_ {
        let per_span = 2 * self.geometry.degree() + 2;
        spans(self.geometry.knots()).flat_map(move |(a, b)| {
            (0..=per_span).map(move |i| a + (b - a) * i as f64 / per_span as f64)
        })
    }
}

impl<S: Space<N>, const N: usize> HasT for NurbsCurve<S, N> {}

impl<S: Space<N>, const N: usize> Mapping<1, N> for NurbsCurve<S, N> {
    type From = T<NurbsCurve<S, N>>;
    type To = S;

    /// A closed curve's parameter is periodic: outside its knots it wraps
    /// round, as an analytic closed curve's does.
    fn apply(&self, t: Point<Self::From, 1>) -> Point<S, N> {
        let t = self.domain.wrap(t);
        Point::with_tag(self.geometry.point(t.coords.x), self.tag)
    }

    fn jacobian(&self, t: Point<Self::From, 1>) -> SMatrix<f64, N, 1> {
        let t = self.domain.wrap(t);
        self.geometry.point_and_tangent(t.coords.x).1
    }
}

impl<S: Space<N>, const N: usize> Embedding<1, N> for NurbsCurve<S, N> {
    fn domain(&self) -> Domain<1> {
        self.domain
    }

    /// Newton iteration from the hint, or else from the nearest of a few samples
    /// per knot span. Finds the nearest point near that start, which is the
    /// nearest overall unless the curve doubles back closer than the sampling.
    fn project(
        &self,
        p: Point<S, N>,
        hint: Option<Point<Self::From, 1>>,
    ) -> Result<Point<Self::From, 1>, ProjectError> {
        let start = match hint {
            Some(hint) => hint,
            None => nearest(
                self.samples().map(|t| Point::new(SVector::from([t]))),
                |t| (self.apply(t) - p).coords.norm_squared(),
            ),
        };
        newton_project(self, p, start)
    }
}

/// A NURBS surface in a 3D space `S`.
///
/// Its domain is the knot vectors': an axis is periodic if the surface closes up
/// across it, and a bound is singular where its whole edge collapses to a point.
#[derive(Clone, Debug)]
pub struct NurbsSurface<S: Space<3>> {
    geometry: golf_nurbs::NurbsSurface<3>,
    tag: S::Tag,
    domain: Domain<2>,
}

impl<S: Space<3, Tag = ()>> NurbsSurface<S> {
    pub fn new(geometry: golf_nurbs::NurbsSurface<3>) -> Self {
        Self::with_tag(geometry, ())
    }
}

impl<S: Space<3>> NurbsSurface<S> {
    /// `geometry` in the instance of `S` that `tag` names.
    pub fn with_tag(geometry: golf_nurbs::NurbsSurface<3>, tag: S::Tag) -> Self {
        let (count_u, count_v) = geometry.control_grid_size();
        let tolerance = tolerance(
            (0..count_u)
                .flat_map(|i| (0..count_v).map(move |j| (i, j)))
                .map(|(i, j)| geometry.control_point(i, j)),
        );
        let domain = Domain::new([
            surface_axis(&geometry, tolerance),
            surface_axis(&geometry.transposed(), tolerance),
        ]);
        Self {
            geometry,
            tag,
            domain,
        }
    }

    pub fn geometry(&self) -> &golf_nurbs::NurbsSurface<3> {
        &self.geometry
    }

    pub fn tag(&self) -> S::Tag {
        self.tag
    }

    /// A grid of parameters spread through every pair of knot spans, for
    /// seeding projection.
    fn samples(&self) -> impl Iterator<Item = (f64, f64)> + '_ {
        let (p, q) = self.geometry.degrees();
        let along = |knots, per_span: usize| -> Vec<f64> {
            spans(knots)
                .flat_map(|(a, b)| {
                    (0..=per_span).map(move |i| a + (b - a) * i as f64 / per_span as f64)
                })
                .collect()
        };
        let us = along(self.geometry.knots_u(), p + 2);
        let vs = along(self.geometry.knots_v(), q + 2);
        us.into_iter()
            .flat_map(move |u| vs.clone().into_iter().map(move |v| (u, v)))
    }
}

impl<S: Space<3>> HasUv for NurbsSurface<S> {}

impl<S: Space<3>> Mapping<2, 3> for NurbsSurface<S> {
    type From = Uv<NurbsSurface<S>>;
    type To = S;

    /// A parameter round an axis the surface closes across is periodic:
    /// outside its knots it wraps round, as an analytic surface's does.
    fn apply(&self, uv: Point<Self::From, 2>) -> Point<S, 3> {
        let uv = self.domain.wrap(uv);
        Point::with_tag(self.geometry.point(uv.coords.x, uv.coords.y), self.tag)
    }

    fn jacobian(&self, uv: Point<Self::From, 2>) -> Matrix3x2<f64> {
        let uv = self.domain.wrap(uv);
        let (_, du, dv) = self.geometry.point_and_partials(uv.coords.x, uv.coords.y);
        Matrix3x2::from_columns(&[du, dv])
    }
}

impl<S: Space<3>> Embedding<2, 3> for NurbsSurface<S> {
    fn domain(&self) -> Domain<2> {
        self.domain
    }

    /// Newton iteration from the hint, or else from each of the few nearest of
    /// a grid of samples per pair of knot spans, keeping the closest result: a
    /// surface that folds back can put the nearest sample on the wrong fold.
    fn project(
        &self,
        p: Point<S, 3>,
        hint: Option<Point<Self::From, 2>>,
    ) -> Result<Point<Self::From, 2>, ProjectError> {
        const SEEDS: usize = 4;
        if let Some(hint) = hint {
            return newton_project(self, p, hint);
        }
        let distance = |uv: Point<Self::From, 2>| (self.apply(uv) - p).coords.norm_squared();
        let mut seeds: Vec<(f64, Point<Self::From, 2>)> = self
            .samples()
            .map(|(u, v)| Point::new(SVector::from([u, v])))
            .map(|uv| (distance(uv), uv))
            .collect();
        seeds.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut best: Option<(f64, Point<Self::From, 2>)> = None;
        let mut error = ProjectError::Ambiguous;
        for &(_, seed) in seeds.iter().take(SEEDS) {
            match newton_project(self, p, seed) {
                Ok(uv) => {
                    let d = distance(uv);
                    if best.is_none_or(|(b, _)| d < b) {
                        best = Some((d, uv));
                    }
                }
                Err(e) => error = e,
            }
        }
        best.map(|(_, uv)| uv).ok_or(error)
    }
}

/// The `u` axis of `surface`: periodic if its `u` ends meet, singular at an end
/// whose edge (the isocurve along `v`) collapses to a point.
fn surface_axis(surface: &golf_nurbs::NurbsSurface<3>, tolerance: f64) -> Axis {
    let ((start, end), (v0, v1)) = surface.domain();
    let closed = (0..=16).all(|i| {
        let v = v0 + (v1 - v0) * i as f64 / 16.0;
        (surface.point(start, v) - surface.point(end, v)).norm() <= tolerance
    });
    let collapses = |u: f64| {
        let edge = surface.isocurve_u(u);
        let first = edge.control_point(0);
        edge.control_points()
            .all(|q| (q - first).norm() <= tolerance)
    };
    let mut axis = match closed {
        true => Axis::periodic(start, end),
        false => Axis::bounded(start, end),
    };
    axis.singular_at_min = collapses(start);
    axis.singular_at_max = collapses(end);
    axis
}

/// The non-empty knot spans within the domain.
fn spans(knots: &golf_nurbs::KnotVector) -> impl Iterator<Item = (f64, f64)> + '_ {
    let (start, end) = knots.domain();
    let inside: Vec<f64> = knots
        .distinct()
        .map(|(k, _)| k)
        .filter(|&k| start <= k && k <= end)
        .collect();
    (0..inside.len() - 1).map(move |i| (inside[i], inside[i + 1]))
}

/// A distance small next to the control points' extent.
fn tolerance<const N: usize>(points: impl Iterator<Item = SVector<f64, N>>) -> f64 {
    let mut lo = SVector::<f64, N>::repeat(f64::INFINITY);
    let mut hi = SVector::<f64, N>::repeat(f64::NEG_INFINITY);
    for p in points {
        lo = lo.inf(&p);
        hi = hi.sup(&p);
    }
    1e-10 * (1.0 + (hi - lo).norm())
}

/// The candidate with the smallest `distance`.
fn nearest<P: Copy>(candidates: impl Iterator<Item = P>, distance: impl Fn(P) -> f64) -> P {
    candidates
        .map(|c| (distance(c), c))
        .min_by(|a, b| a.0.total_cmp(&b.0))
        .expect("a NURBS domain has at least one span")
        .1
}

#[cfg(test)]
mod tests {
    use core::f64::consts::PI;

    use golf_frame::Frame;
    use golf_frame::FrameTree;
    use golf_manifold::Compose;
    use golf_manifold::Surface;
    use golf_manifold::World;
    use nalgebra::Isometry3;
    use nalgebra::Translation3;
    use nalgebra::UnitQuaternion;
    use nalgebra::Vector2;
    use nalgebra::Vector3;

    use super::*;
    use crate::Circle;
    use crate::Cylinder;
    use crate::Placement;
    use crate::Sphere;
    use crate::testing::assert_jacobian;
    use crate::testing::assert_round_trip;

    fn circle_geometry(radius: f64) -> golf_nurbs::NurbsCurve<3> {
        golf_nurbs::NurbsCurve::circular_arc(
            Vector3::zeros(),
            Vector3::x(),
            Vector3::y(),
            radius,
            0.0,
            0.0,
        )
    }

    fn open_cubic() -> NurbsCurve<World, 3> {
        let points = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 2.0, 0.0),
            Vector3::new(2.0, -1.0, 1.0),
            Vector3::new(4.0, 1.0, 0.0),
            Vector3::new(5.0, 0.0, 2.0),
        ];
        let knots = vec![0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0];
        NurbsCurve::new(golf_nurbs::NurbsCurve::new(3, knots, points, None).unwrap())
    }

    /// The `[0, height]` stretch of a cylinder of `radius` about z, as a NURBS.
    fn cylinder_geometry(radius: f64, height: f64) -> golf_nurbs::NurbsSurface<3> {
        let arc = circle_geometry(radius);
        let lift = Vector3::new(0.0, 0.0, height);
        golf_nurbs::NurbsSurface::new(
            (2, arc.knots().knots().to_vec()),
            (1, vec![0.0, 0.0, 1.0, 1.0]),
            arc.control_points().map(|p| vec![p, p + lift]).collect(),
            Some(
                (0..arc.control_point_count())
                    .map(|i| vec![arc.weight(i); 2])
                    .collect(),
            ),
        )
        .unwrap()
    }

    #[test]
    fn closed_curve_is_periodic() {
        let c = NurbsCurve::<World, 3>::new(circle_geometry(2.0));
        assert!(c.domain().axes[0].periodic);
        assert!(!open_cubic().domain().axes[0].periodic);
        for t in [0.0, 0.1, 0.5, 0.9, 0.999] {
            assert_round_trip(&c, [t]);
        }
        // Away from the ends and the double knots, where the derivative jumps.
        for t in [0.1, 0.4, 0.6, 0.9] {
            assert_jacobian(&c, [t]);
        }
    }

    #[test]
    fn curve_projection_matches_analytic_circle() {
        let nurbs = NurbsCurve::<World, 3>::new(circle_geometry(2.0));
        let circle = Circle::new(Placement::at(Point::new(Vector3::zeros())), 2.0);
        for p in [
            [3.0, 1.0, 0.5],
            [-0.5, 0.2, -2.0],
            [0.0, -4.0, 1.0],
            [-3.0, -0.01, 0.0],
        ] {
            let p = Point::new(Vector3::from(p));
            let a = nurbs.apply(nurbs.project(p, None).unwrap());
            let b = circle.apply(circle.project(p, None).unwrap());
            assert!((a - b).coords.norm() < 1e-9, "{p:?}: {a:?} != {b:?}");
        }
    }

    #[test]
    fn open_curve_projection_clamps_to_ends() {
        let c = open_cubic();
        for t in [0.0, 0.2, 0.5, 0.8, 1.0] {
            assert_round_trip(&c, [t]);
        }
        for t in [0.2, 0.8] {
            assert_jacobian(&c, [t]);
        }
        // Behind the start, along the start tangent.
        let behind = Point::new(Vector3::new(-1.0, -2.0, 0.0));
        assert_eq!(c.project(behind, None).unwrap().coords.x, 0.0);
    }

    #[test]
    fn pcurve_composes_with_its_surface() {
        let sphere = Sphere::<World>::new(Point::new(Vector3::zeros()), 3.0);
        let segment = golf_nurbs::NurbsCurve::new(
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![Vector2::new(0.0, 0.0), Vector2::new(PI / 2.0, 0.5)],
            None,
        )
        .unwrap();
        let pcurve = NurbsCurve::<Uv<Sphere<World>>, 2>::new(segment);
        let edge = Compose::<_, _, 2>(pcurve, sphere);
        for t in [0.0, 0.3, 1.0] {
            let p = edge.apply(Point::new(SVector::from([t])));
            assert!((p.coords.norm() - 3.0).abs() < 1e-12);
        }
    }

    #[test]
    fn closed_surface_matches_analytic_cylinder() {
        let nurbs = NurbsSurface::<World>::new(cylinder_geometry(1.5, 4.0));
        let domain = nurbs.domain();
        assert!(domain.axes[0].periodic && !domain.axes[1].periodic);
        assert!(!domain.is_singular(
            Point::<Uv<NurbsSurface<World>>, 2>::new(Vector2::new(0.0, 0.0)),
            1e-9
        ));
        for at in [[0.1, 0.1], [0.5, 0.5], [0.95, 0.9]] {
            assert_round_trip(&nurbs, at);
        }
        for at in [[0.1, 0.1], [0.6, 0.3], [0.9, 0.8]] {
            assert_jacobian(&nurbs, at);
        }
        let cylinder = Cylinder::new(Placement::at(Point::new(Vector3::zeros())), 1.5);
        for p in [
            [3.0, 1.0, 0.5],
            [-0.5, 0.2, 2.0],
            [0.1, -0.4, 3.9],
            [-2.0, -2.0, 1.0],
        ] {
            let p = Point::new(Vector3::from(p));
            let a = nurbs.apply(nurbs.project(p, None).unwrap());
            let b = cylinder.apply(cylinder.project(p, None).unwrap());
            assert!((a - b).coords.norm() < 1e-9, "{p:?}: {a:?} != {b:?}");
            assert!(
                (nurbs.normal(nurbs.project(p, None).unwrap()).coords
                    - Vector3::new(a.coords.x, a.coords.y, 0.0) / 1.5)
                    .norm()
                    < 1e-9
            );
        }
    }

    #[test]
    fn collapsed_edge_is_singular() {
        // A cone: the rim circle at v = 0, the apex at v = 1.
        let rim = circle_geometry(1.0);
        let apex = Vector3::new(0.0, 0.0, 2.0);
        let cone = golf_nurbs::NurbsSurface::new(
            (2, rim.knots().knots().to_vec()),
            (1, vec![0.0, 0.0, 1.0, 1.0]),
            rim.control_points().map(|p| vec![p, apex]).collect(),
            Some(
                (0..rim.control_point_count())
                    .map(|i| vec![rim.weight(i); 2])
                    .collect(),
            ),
        )
        .unwrap();
        let s = NurbsSurface::<World>::new(cone);
        let v = s.domain().axes[1];
        assert!(!v.singular_at_min && v.singular_at_max);
        assert!(!s.domain().axes[0].singular_at_min);
        // Points above the apex project onto it.
        let uv = s
            .project(Point::new(Vector3::new(0.0, 0.0, 5.0)), None)
            .unwrap();
        assert!((s.apply(uv).coords - apex).norm() < 1e-9);
    }

    #[test]
    fn surface_in_a_frame_carries_its_tag() {
        let mut tree = FrameTree::new();
        let frame = tree.add(
            tree.root(),
            Isometry3::from_parts(
                Translation3::new(1.0, 2.0, 3.0),
                UnitQuaternion::from_scaled_axis(Vector3::new(0.3, 0.0, 0.0)),
            ),
        );
        let s = NurbsSurface::<Frame>::with_tag(cylinder_geometry(1.0, 1.0), frame);
        let uv = Point::new(Vector2::new(0.3, 0.5));
        assert_eq!(s.apply(uv).tag(), frame);
        let p = frame.point(Vector3::new(2.0, 0.0, 0.5));
        let q = s.project(p, None).unwrap();
        assert!((s.apply(q).coords - Vector3::new(1.0, 0.0, 0.5)).norm() < 1e-9);
    }
}
