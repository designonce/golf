//! NURBS curves.

use core::f64::consts::FRAC_PI_2;
use core::f64::consts::TAU;

use nalgebra::SVector;

use crate::error::NurbsError;
use crate::homogeneous::Homogeneous;
use crate::homogeneous::binomial;
use crate::knots::KnotVector;
use crate::knots::{self};

/// A NURBS curve in `D` dimensions.
#[derive(Clone, Debug, PartialEq)]
pub struct NurbsCurve<const D: usize> {
    pub(crate) knots: KnotVector,
    pub(crate) points: Vec<Homogeneous<D>>,
}

impl<const D: usize> NurbsCurve<D> {
    /// A curve of `degree` through `points` with `knots`, optionally `weights`
    /// (one per point; all 1 if `None`).
    pub fn new(
        degree: usize,
        knots: Vec<f64>,
        points: Vec<SVector<f64, D>>,
        weights: Option<Vec<f64>>,
    ) -> Result<Self, NurbsError> {
        let knots = KnotVector::new(degree, knots)?;
        let points = homogenise(&points, weights.as_deref())?;
        if points.len() != knots.basis_count() {
            return Err(NurbsError::ControlPointCount {
                expected: knots.basis_count(),
                found: points.len(),
            });
        }
        Ok(Self { knots, points })
    }

    /// The Bézier curve with `points` as control points, on `[0, 1]`.
    pub fn bezier(points: Vec<SVector<f64, D>>) -> Result<Self, NurbsError> {
        let degree = points.len().saturating_sub(1);
        let knots = KnotVector::bezier(degree)?;
        Self::new(degree, knots.knots().to_vec(), points, None)
    }

    /// The circular arc about `center` in the plane of the orthonormal `x_axis`
    /// and `y_axis`, from angle `start` to `end` (radians, anticlockwise from
    /// `x_axis`), as a rational quadratic on `[0, 1]` (A7.1).
    ///
    /// `end` is taken modulo a full turn after `start`, and a full circle when
    /// they are equal.
    pub fn circular_arc(
        center: SVector<f64, D>,
        x_axis: SVector<f64, D>,
        y_axis: SVector<f64, D>,
        radius: f64,
        start: f64,
        end: f64,
    ) -> Self {
        debug_assert!((x_axis.norm() - 1.0).abs() < 1e-9 && (y_axis.norm() - 1.0).abs() < 1e-9);
        debug_assert!(x_axis.dot(&y_axis).abs() < 1e-9);
        let mut sweep = (end - start).rem_euclid(TAU);
        if sweep == 0.0 {
            sweep = TAU;
        }
        let arcs = (sweep / FRAC_PI_2).ceil().clamp(1.0, 4.0) as usize;
        let step = sweep / arcs as f64;
        let w = (step / 2.0).cos();
        let at = |angle: f64, r: f64| center + (x_axis * angle.cos() + y_axis * angle.sin()) * r;

        let mut points = vec![Homogeneous::new(at(start, radius), 1.0)];
        for i in 1..=arcs {
            let angle = start + step * i as f64;
            // The middle control point is where the end tangents meet.
            points.push(Homogeneous::new(at(angle - step / 2.0, radius / w), w));
            points.push(Homogeneous::new(at(angle, radius), 1.0));
        }
        let mut knots = vec![0.0; 3];
        for i in 1..arcs {
            let k = i as f64 / arcs as f64;
            knots.extend([k, k]);
        }
        knots.extend([1.0; 3]);
        Self {
            knots: KnotVector::from_raw(2, knots),
            points,
        }
    }

    pub fn degree(&self) -> usize {
        self.knots.degree()
    }

    pub fn knots(&self) -> &KnotVector {
        &self.knots
    }

    /// The parameter range the curve is defined over.
    pub fn domain(&self) -> (f64, f64) {
        self.knots.domain()
    }

    pub fn control_point_count(&self) -> usize {
        self.points.len()
    }

    pub fn control_point(&self, i: usize) -> SVector<f64, D> {
        self.points[i].point()
    }

    pub fn weight(&self, i: usize) -> f64 {
        self.points[i].w
    }

    pub fn control_points(&self) -> impl ExactSizeIterator<Item = SVector<f64, D>> + '_ {
        self.points.iter().map(Homogeneous::point)
    }

    /// Whether any weight differs from the others, making the curve rational
    /// rather than polynomial.
    pub fn is_rational(&self) -> bool {
        self.points.iter().any(|p| p.w != self.points[0].w)
    }

    /// The point at `u`, which is clamped into the [`domain`](Self::domain).
    pub fn point(&self, u: f64) -> SVector<f64, D> {
        let span = self.knots.span(u);
        let u = self.clamp(u);
        let basis = self.knots.basis(span, u);
        let first = span - self.degree();
        basis
            .iter()
            .zip(&self.points[first..])
            .fold(Homogeneous::zero(), |acc, (&b, &p)| acc + p * b)
            .point()
    }

    /// The point at `u` and its derivatives: `result[k]` is the `k`th
    /// derivative, for `k` up to `order` (A3.2 and A4.2). `u` is clamped into
    /// the domain; at a knot, derivatives are taken from the span after it
    /// (or before it, at the end of the domain).
    pub fn derivatives(&self, u: f64, order: usize) -> Vec<SVector<f64, D>> {
        let span = self.knots.span(u);
        let u = self.clamp(u);
        let ders = self.knots.basis_derivatives(span, u, order);
        let first = span - self.degree();
        let homogeneous: Vec<Homogeneous<D>> = ders
            .iter()
            .map(|row| {
                row.iter()
                    .zip(&self.points[first..])
                    .fold(Homogeneous::zero(), |acc, (&b, &p)| acc + p * b)
            })
            .collect();
        rational_derivatives(&homogeneous)
    }

    /// The same curve with `u` inserted into the knots `times` more times
    /// (capped at the degree). `u` must be strictly inside the domain.
    pub fn insert_knot(&self, u: f64, times: usize) -> Result<Self, NurbsError> {
        let (knots, points) = knots::insert_knot(&self.knots, &self.points, u, times)?;
        Ok(Self { knots, points })
    }

    /// The two pieces either side of `u`, which must be strictly inside the
    /// domain. Each is clamped at the cut.
    pub fn split(&self, u: f64) -> Result<(Self, Self), NurbsError> {
        let (start, end) = self.domain();
        if !(start < u && u < end) {
            return Err(NurbsError::ParameterOutsideDomain { parameter: u });
        }
        Ok((self.before(u), self.after(u)))
    }

    /// The same curve over the same domain, with its knots clamped so it
    /// starts and ends on control points.
    pub fn clamped(&self) -> Self {
        let (start, end) = self.domain();
        self.after(start).before(end)
    }

    /// The Bézier segments making up the curve, one per non-empty knot span,
    /// each keeping its original parameter range.
    pub fn bezier_segments(&self) -> Vec<Self> {
        let (start, end) = self.domain();
        let mut segments = Vec::new();
        let mut rest = self.clamped();
        for (u, _) in self.knots.distinct().filter(|&(k, _)| start < k && k < end) {
            let (left, right) = rest.split(u).expect("cut is inside the remaining domain");
            segments.push(left);
            rest = right;
        }
        segments.push(rest);
        segments
    }

    /// The same curve traversed the other way, over the same domain.
    pub fn reversed(&self) -> Self {
        let knots = self.knots.knots();
        let (first, last) = (knots[0], knots[knots.len() - 1]);
        let knots = knots.iter().rev().map(|k| first + last - k).collect();
        let points = self.points.iter().rev().copied().collect();
        Self {
            knots: KnotVector::from_raw(self.degree(), knots),
            points,
        }
    }

    /// The same curve with its degree raised by `by` (A5.9).
    pub fn elevate_degree(&self, by: usize) -> Result<Self, NurbsError> {
        if by == 0 {
            return Ok(self.clone());
        }
        let clamped = self.clamped();
        let (knots, points) = elevate(&clamped.knots, &clamped.points, by);
        Ok(Self { knots, points })
    }

    /// Applies `f` to every control point. Only meaningful for affine `f`
    /// (rigid motions, scaling, projections), which NURBS are invariant under.
    pub fn map_points<const E: usize>(
        &self,
        f: impl Fn(SVector<f64, D>) -> SVector<f64, E>,
    ) -> NurbsCurve<E> {
        NurbsCurve {
            knots: self.knots.clone(),
            points: self
                .points
                .iter()
                .map(|h| Homogeneous::new(f(h.point()), h.w))
                .collect(),
        }
    }

    fn clamp(&self, u: f64) -> f64 {
        let (start, end) = self.domain();
        u.clamp(start, end)
    }

    /// The part of the curve after `u` (in the domain), clamped at `u`.
    fn after(&self, u: f64) -> Self {
        let p = self.degree();
        if self.knots.multiplicity(u) > p {
            return self.clone();
        }
        // With `u` repeated `p` times the curve passes through a control point
        // there, so the rest is a curve in its own right.
        let (knots, points) = knots::insert_raw(&self.knots, &self.points, u, p);
        let first = knots
            .knots()
            .iter()
            .position(|&k| k == u)
            .expect("u is a knot");
        let mut new_knots = vec![u];
        new_knots.extend_from_slice(&knots.knots()[first..]);
        let count = new_knots.len() - p - 1;
        Self {
            knots: KnotVector::from_raw(p, new_knots),
            points: points[points.len() - count..].to_vec(),
        }
    }

    /// The part of the curve before `u` (in the domain), clamped at `u`.
    fn before(&self, u: f64) -> Self {
        let p = self.degree();
        if self.knots.multiplicity(u) > p {
            return self.clone();
        }
        let (knots, points) = knots::insert_raw(&self.knots, &self.points, u, p);
        let last = knots
            .knots()
            .iter()
            .rposition(|&k| k == u)
            .expect("u is a knot");
        let mut new_knots = knots.knots()[..=last].to_vec();
        new_knots.push(u);
        let count = new_knots.len() - p - 1;
        Self {
            knots: KnotVector::from_raw(p, new_knots),
            points: points[..count].to_vec(),
        }
    }
}

/// Turns points and optional weights into homogeneous control points.
pub(crate) fn homogenise<const D: usize>(
    points: &[SVector<f64, D>],
    weights: Option<&[f64]>,
) -> Result<Vec<Homogeneous<D>>, NurbsError> {
    match weights {
        None => Ok(points.iter().map(|&p| Homogeneous::new(p, 1.0)).collect()),
        Some(weights) => {
            if weights.len() != points.len() {
                return Err(NurbsError::ControlPointCount {
                    expected: points.len(),
                    found: weights.len(),
                });
            }
            if let Some((index, &weight)) = weights
                .iter()
                .enumerate()
                .find(|(_, w)| !(w.is_finite() && **w > 0.0))
            {
                return Err(NurbsError::InvalidWeight { index, weight });
            }
            Ok(points
                .iter()
                .zip(weights)
                .map(|(&p, &w)| Homogeneous::new(p, w))
                .collect())
        }
    }
}

/// Derivatives of a rational curve from those of its homogeneous form (A4.2).
pub(crate) fn rational_derivatives<const D: usize>(h: &[Homogeneous<D>]) -> Vec<SVector<f64, D>> {
    let mut ck: Vec<SVector<f64, D>> = Vec::with_capacity(h.len());
    for k in 0..h.len() {
        let mut v = h[k].wp;
        for i in 1..=k {
            v -= ck[k - i] * (binomial(k, i) * h[i].w);
        }
        ck.push(v / h[0].w);
    }
    ck
}

/// Degree elevation of a clamped B-spline by `t` (A5.9).
///
/// Kept index-for-index with the book's pseudocode so it can be checked against it.
#[allow(clippy::needless_range_loop, clippy::explicit_counter_loop)]
fn elevate<const D: usize>(
    knots: &KnotVector,
    pw: &[Homogeneous<D>],
    t: usize,
) -> (KnotVector, Vec<Homogeneous<D>>) {
    let p = knots.degree();
    let u = knots.knots();
    let m = u.len() - 1;
    let ph = p + t;
    let ph2 = ph / 2;

    // Coefficients for elevating one Bézier segment.
    let mut bezalfs = vec![vec![0.0; p + 1]; ph + 1];
    bezalfs[0][0] = 1.0;
    bezalfs[ph][p] = 1.0;
    for i in 1..=ph2 {
        let inv = 1.0 / binomial(ph, i);
        for j in i.saturating_sub(t)..=p.min(i) {
            bezalfs[i][j] = inv * binomial(p, j) * binomial(t, i - j);
        }
    }
    for i in ph2 + 1..ph {
        for j in i.saturating_sub(t)..=p.min(i) {
            bezalfs[i][j] = bezalfs[ph - i][p - j];
        }
    }

    let segments = knots.distinct().count() - 1;
    let mut uh = vec![0.0; u.len() + segments * t + t];
    let mut qw = vec![Homogeneous::zero(); pw.len() + segments * t + t];
    let mut bpts: Vec<Homogeneous<D>> = pw[..=p].to_vec();
    let mut ebpts = vec![Homogeneous::zero(); ph + 1];
    let mut next_bpts = vec![Homogeneous::zero(); p.saturating_sub(1).max(1)];
    let mut alfs = vec![0.0; p.saturating_sub(1).max(1)];

    let mut mh = ph;
    let mut kind = ph + 1;
    let mut r: isize = -1;
    let mut a = p;
    let mut b = p + 1;
    let mut cind = 1;
    let mut ua = u[0];
    qw[0] = pw[0];
    for k in uh.iter_mut().take(ph + 1) {
        *k = ua;
    }

    while b < m {
        let i = b;
        while b < m && u[b] == u[b + 1] {
            b += 1;
        }
        let mul = b - i + 1;
        mh += mul + t;
        let ub = u[b];
        let oldr = r;
        r = p as isize - mul as isize;
        let lbz = if oldr > 0 {
            ((oldr + 2) / 2) as usize
        } else {
            1
        };
        let rbz = if r > 0 {
            ph - ((r + 1) / 2) as usize
        } else {
            ph
        };

        // Insert `ub` r times to isolate the current Bézier segment.
        if r > 0 {
            let numer = ub - ua;
            for k in (mul + 1..=p).rev() {
                alfs[k - mul - 1] = numer / (u[a + k] - ua);
            }
            for j in 1..=r as usize {
                let save = r as usize - j;
                let s = mul + j;
                for k in (s..=p).rev() {
                    bpts[k] = bpts[k - 1].lerp(bpts[k], alfs[k - s]);
                }
                next_bpts[save] = bpts[p];
            }
        }

        // Elevate it.
        for i in lbz..=ph {
            ebpts[i] = Homogeneous::zero();
            for j in i.saturating_sub(t)..=p.min(i) {
                ebpts[i] += bpts[j] * bezalfs[i][j];
            }
        }

        // Remove the knot `ua` that the previous insertion added.
        if oldr > 1 {
            let oldr = oldr as usize;
            let mut first = kind - 2;
            let mut last = kind;
            let den = ub - ua;
            let bet = (ub - uh[kind - 1]) / den;
            for tr in 1..oldr {
                let mut i = first;
                let mut j = last;
                let mut kj = j as isize - kind as isize + 1;
                while j - i > tr {
                    if i < cind {
                        let alf = (ub - uh[i]) / (ua - uh[i]);
                        qw[i] = qw[i - 1].lerp(qw[i], alf);
                    }
                    if j >= lbz {
                        let kj_ = kj as usize;
                        if j - tr <= kind - ph + oldr {
                            let gam = (ub - uh[j - tr]) / den;
                            ebpts[kj_] = ebpts[kj_ + 1].lerp(ebpts[kj_], gam);
                        } else {
                            ebpts[kj_] = ebpts[kj_ + 1].lerp(ebpts[kj_], bet);
                        }
                    }
                    i += 1;
                    j -= 1;
                    kj -= 1;
                }
                first -= 1;
                last += 1;
            }
        }

        if a != p {
            for _ in 0..(ph as isize - oldr) as usize {
                uh[kind] = ua;
                kind += 1;
            }
        }
        for j in lbz..=rbz {
            qw[cind] = ebpts[j];
            cind += 1;
        }

        if b < m {
            let r = r.max(0) as usize;
            bpts[..r].copy_from_slice(&next_bpts[..r]);
            for j in r..=p {
                bpts[j] = pw[b - p + j];
            }
            a = b;
            b += 1;
            ua = ub;
        } else {
            for i in 0..=ph {
                uh[kind + i] = ub;
            }
        }
    }

    let nh = mh - ph - 1;
    uh.truncate(mh + 1);
    qw.truncate(nh + 1);
    (KnotVector::from_raw(ph, uh), qw)
}

#[cfg(test)]
mod tests {
    use nalgebra::Vector2;
    use nalgebra::Vector3;

    use super::*;

    fn cubic() -> NurbsCurve<3> {
        NurbsCurve::new(
            3,
            vec![0.0, 0.0, 0.0, 0.0, 0.3, 0.5, 0.5, 0.8, 1.0, 1.0, 1.0, 1.0],
            vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(1.0, 2.0, 0.0),
                Vector3::new(2.0, 2.5, 1.0),
                Vector3::new(3.0, 0.0, 2.0),
                Vector3::new(4.0, -1.0, 1.0),
                Vector3::new(5.0, 1.0, 0.0),
                Vector3::new(6.0, 3.0, -1.0),
                Vector3::new(7.0, 0.0, 0.0),
            ],
            Some(vec![1.0, 0.5, 2.0, 1.0, 1.5, 0.7, 1.0, 1.2]),
        )
        .unwrap()
    }

    fn samples(curve: &NurbsCurve<3>) -> impl Iterator<Item = f64> {
        let (a, b) = curve.domain();
        (0..=200).map(move |i| a + (b - a) * i as f64 / 200.0)
    }

    fn assert_same_shape(a: &NurbsCurve<3>, b: &NurbsCurve<3>) {
        assert_eq!(a.domain(), b.domain());
        for u in samples(a) {
            let (pa, pb) = (a.point(u), b.point(u));
            assert!((pa - pb).norm() < 1e-11, "u {u}: {pa:?} != {pb:?}");
        }
    }

    #[test]
    fn rejects_mismatched_inputs() {
        let knots = vec![0.0, 0.0, 1.0, 1.0];
        let points = vec![Vector2::new(0.0, 0.0), Vector2::new(1.0, 0.0)];
        assert!(matches!(
            NurbsCurve::new(1, knots.clone(), points[..1].to_vec(), None),
            Err(NurbsError::ControlPointCount {
                expected: 2,
                found: 1
            })
        ));
        assert!(matches!(
            NurbsCurve::new(1, knots, points, Some(vec![1.0, -1.0])),
            Err(NurbsError::InvalidWeight { index: 1, .. })
        ));
    }

    #[test]
    fn bezier_interpolates_ends_and_matches_de_casteljau() {
        let pts = vec![
            Vector2::new(0.0, 0.0),
            Vector2::new(1.0, 2.0),
            Vector2::new(3.0, 1.0),
        ];
        let c = NurbsCurve::bezier(pts.clone()).unwrap();
        assert_eq!(c.point(0.0), pts[0]);
        assert_eq!(c.point(1.0), pts[2]);
        let t = 0.3;
        let expected =
            pts[0] * (1.0 - t) * (1.0 - t) + pts[1] * 2.0 * t * (1.0 - t) + pts[2] * t * t;
        assert!((c.point(t) - expected).norm() < 1e-14);
    }

    #[test]
    fn circular_arcs_lie_on_the_circle() {
        let center = Vector3::new(1.0, -2.0, 0.5);
        let (x, y) = (Vector3::new(0.0, 0.0, 1.0), Vector3::new(0.0, 1.0, 0.0));
        for (start, end) in [(0.0, 1.0), (0.5, 3.0), (-1.0, 3.5), (2.0, 2.0)] {
            let arc = NurbsCurve::circular_arc(center, x, y, 2.5, start, end);
            for u in samples(&arc) {
                assert!(((arc.point(u) - center).norm() - 2.5).abs() < 1e-12);
            }
            let at = |angle: f64| center + (x * angle.cos() + y * angle.sin()) * 2.5;
            assert!((arc.point(0.0) - at(start)).norm() < 1e-12);
            assert!((arc.point(1.0) - at(end)).norm() < 1e-12);
        }
        let full = NurbsCurve::circular_arc(center, x, y, 2.5, 0.0, 0.0);
        assert_eq!(full.control_point_count(), 9);
        assert!(full.is_rational());
    }

    #[test]
    fn derivatives_match_finite_differences() {
        let c = cubic();
        let h = 1e-5;
        for u in [0.05, 0.2, 0.4, 0.65, 0.9] {
            let d = c.derivatives(u, 3);
            assert_eq!(d[0], c.point(u));
            let lo = c.derivatives(u - h, 2);
            let hi = c.derivatives(u + h, 2);
            for k in 1..=3 {
                let fd = (hi[k - 1] - lo[k - 1]) / (2.0 * h);
                assert!(
                    (fd - d[k]).norm() < 1e-4 * (1.0 + d[k].norm()),
                    "u {u} k {k}"
                );
            }
        }
    }

    #[test]
    fn knot_insertion_preserves_shape() {
        let c = cubic();
        for (u, times) in [(0.1, 1), (0.3, 2), (0.5, 1), (0.5, 5), (0.77, 3)] {
            let refined = c.insert_knot(u, times).unwrap();
            let s = c.knots().multiplicity(u);
            let added = times.min(3 - s);
            assert_eq!(
                refined.control_point_count(),
                c.control_point_count() + added
            );
            assert_same_shape(&c, &refined);
        }
        assert!(c.insert_knot(0.0, 1).is_err());
        assert!(c.insert_knot(1.0, 1).is_err());
    }

    #[test]
    fn split_pieces_reassemble() {
        let c = cubic();
        for cut in [0.2, 0.3, 0.5, 0.9] {
            let (left, right) = c.split(cut).unwrap();
            assert_eq!(left.domain(), (0.0, cut));
            assert_eq!(right.domain(), (cut, 1.0));
            assert!(left.knots().is_clamped() && right.knots().is_clamped());
            for u in samples(&c) {
                let piece = if u <= cut { &left } else { &right };
                assert!(
                    (piece.point(u) - c.point(u)).norm() < 1e-11,
                    "cut {cut} u {u}"
                );
            }
        }
    }

    #[test]
    fn bezier_segments_cover_each_span() {
        let c = cubic();
        let segments = c.bezier_segments();
        let spans: Vec<(f64, f64)> = segments.iter().map(NurbsCurve::domain).collect();
        assert_eq!(spans, vec![(0.0, 0.3), (0.3, 0.5), (0.5, 0.8), (0.8, 1.0)]);
        for segment in &segments {
            assert_eq!(segment.control_point_count(), 4);
            let (a, b) = segment.domain();
            for i in 0..=20 {
                let u = a + (b - a) * i as f64 / 20.0;
                assert!((segment.point(u) - c.point(u)).norm() < 1e-11);
            }
        }
    }

    #[test]
    fn bezier_segments_of_unclamped_curve() {
        let knots: Vec<f64> = (0..9).map(f64::from).collect();
        let points = (0..6)
            .map(|i| Vector3::new(i as f64, (i * i % 5) as f64, 0.0))
            .collect();
        let c = NurbsCurve::new(2, knots, points, None).unwrap();
        let segments = c.bezier_segments();
        assert_eq!(segments.len(), 4);
        for segment in &segments {
            assert!(segment.knots().is_clamped());
            assert_eq!(segment.control_point_count(), 3);
            let (a, b) = segment.domain();
            for i in 0..=10 {
                let u = a + (b - a) * i as f64 / 10.0;
                assert!((segment.point(u) - c.point(u)).norm() < 1e-11);
            }
        }
    }

    #[test]
    fn reversal_runs_backwards() {
        let c = cubic();
        let r = c.reversed();
        assert_eq!(r.domain(), c.domain());
        for u in samples(&c) {
            assert!((r.point(u) - c.point(1.0 - u)).norm() < 1e-11);
        }
        assert_same_shape(&r.reversed(), &c);
    }

    #[test]
    fn degree_elevation_preserves_shape() {
        let c = cubic();
        for by in 1..=3 {
            let e = c.elevate_degree(by).unwrap();
            assert_eq!(e.degree(), 3 + by);
            // Each of the 4 Bézier segments gains `by` control points, less the
            // continuity the original interior knots kept.
            let distinct_interior = 3;
            assert_eq!(
                e.control_point_count(),
                c.control_point_count() + by * (distinct_interior + 1)
            );
            assert_same_shape(&c, &e);
        }
        let circle =
            NurbsCurve::circular_arc(Vector3::zeros(), Vector3::x(), Vector3::y(), 1.0, 0.0, 0.0);
        let e = circle.elevate_degree(1).unwrap();
        assert_same_shape(&circle, &e);
    }

    #[test]
    fn elevating_an_unclamped_curve_clamps_it_first() {
        let knots: Vec<f64> = (0..9).map(f64::from).collect();
        let points = (0..6)
            .map(|i| Vector3::new(i as f64, (i * i % 5) as f64, 1.0))
            .collect();
        let c =
            NurbsCurve::new(2, knots, points, Some(vec![1.0, 2.0, 0.5, 1.0, 3.0, 1.0])).unwrap();
        let e = c.elevate_degree(2).unwrap();
        assert!(e.knots().is_clamped());
        assert_same_shape(&c, &e);
    }

    #[test]
    fn affine_maps_commute_with_evaluation() {
        let c = cubic();
        let f = |p: Vector3<f64>| Vector2::new(2.0 * p.x - p.z + 1.0, p.y * 0.5);
        let mapped = c.map_points(f);
        for u in samples(&c) {
            assert!((mapped.point(u) - f(c.point(u))).norm() < 1e-11);
        }
    }
}
