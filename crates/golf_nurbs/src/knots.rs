//! Knot vectors and the B-spline basis functions they define.

use crate::error::NurbsError;
use crate::homogeneous::Homogeneous;

/// A validated, non-decreasing knot vector for B-splines of a given degree.
#[derive(Clone, Debug, PartialEq)]
pub struct KnotVector {
    degree: usize,
    knots: Vec<f64>,
}

impl KnotVector {
    /// Checks that `knots` define a valid B-spline basis of `degree`:
    /// finite, non-decreasing, at least `2(degree + 1)` of them, a non-empty
    /// domain, and no knot repeated more than `degree` times inside the domain
    /// (`degree + 1` at its ends), so every curve built on it is continuous.
    pub fn new(degree: usize, knots: Vec<f64>) -> Result<Self, NurbsError> {
        if degree == 0 {
            return Err(NurbsError::ZeroDegree);
        }
        if knots.len() < 2 * (degree + 1) {
            return Err(NurbsError::TooFewKnots {
                degree,
                knots: knots.len(),
            });
        }
        if knots.iter().any(|k| !k.is_finite()) || knots.windows(2).any(|w| w[0] > w[1]) {
            return Err(NurbsError::InvalidKnots);
        }
        let vector = Self { degree, knots };
        let (start, end) = vector.domain();
        if start >= end {
            return Err(NurbsError::EmptyDomain);
        }
        for (knot, multiplicity) in vector.distinct() {
            let interior = start < knot && knot < end;
            let limit = if interior { degree } else { degree + 1 };
            if multiplicity > limit {
                return Err(NurbsError::MultiplicityTooHigh { knot, multiplicity });
            }
        }
        Ok(vector)
    }

    /// Knots for `count` basis functions, clamped at 0 and 1 with uniform
    /// interior knots.
    pub fn clamped_uniform(degree: usize, count: usize) -> Result<Self, NurbsError> {
        let interior = count.saturating_sub(degree + 1);
        let knots = core::iter::repeat_n(0.0, degree + 1)
            .chain((1..=interior).map(|i| i as f64 / (interior + 1) as f64))
            .chain(core::iter::repeat_n(1.0, degree + 1))
            .collect();
        Self::new(degree, knots)
    }

    /// The knots of a single Bézier segment on `[0, 1]`.
    pub fn bezier(degree: usize) -> Result<Self, NurbsError> {
        Self::clamped_uniform(degree, degree + 1)
    }

    pub fn degree(&self) -> usize {
        self.degree
    }

    pub fn knots(&self) -> &[f64] {
        &self.knots
    }

    /// How many basis functions (and so control points) the knots define.
    pub fn basis_count(&self) -> usize {
        self.knots.len() - self.degree - 1
    }

    /// The parameter range `[u_p, u_{n+1}]` over which the basis sums to one.
    pub fn domain(&self) -> (f64, f64) {
        (self.knots[self.degree], self.knots[self.basis_count()])
    }

    /// Whether the first and last knots are each repeated `degree + 1` times,
    /// so curves interpolate their end control points.
    pub fn is_clamped(&self) -> bool {
        let p = self.degree;
        let (start, end) = self.domain();
        self.knots[..=p].iter().all(|&k| k == start)
            && self.knots[self.knots.len() - p - 1..]
                .iter()
                .all(|&k| k == end)
    }

    /// Each distinct knot value with its multiplicity, in order.
    pub fn distinct(&self) -> impl Iterator<Item = (f64, usize)> + '_ {
        self.knots
            .chunk_by(|a, b| a == b)
            .map(|run| (run[0], run.len()))
    }

    /// How many times `u` appears in the knot vector.
    pub fn multiplicity(&self, u: f64) -> usize {
        self.knots.iter().filter(|&&k| k == u).count()
    }

    /// The span index `i` with `u_i <= u < u_{i+1}` holding `u` (A2.1).
    ///
    /// `u` is clamped into the domain; the domain's end belongs to the last
    /// non-empty span.
    pub fn span(&self, u: f64) -> usize {
        let p = self.degree;
        let n = self.basis_count() - 1;
        let (start, end) = self.domain();
        if u >= end {
            let mut i = n;
            while self.knots[i] == self.knots[i + 1] {
                i -= 1;
            }
            return i;
        }
        if u <= start {
            let mut i = p;
            while self.knots[i] == self.knots[i + 1] {
                i += 1;
            }
            return i;
        }
        // The last index in p..=n whose knot is <= u.
        p + self.knots[p..=n].partition_point(|&k| k <= u) - 1
    }

    /// The `degree + 1` basis functions non-zero on `span`, at `u` (A2.2).
    pub fn basis(&self, span: usize, u: f64) -> Vec<f64> {
        let p = self.degree;
        let u_knots = &self.knots;
        let mut n = vec![0.0; p + 1];
        let mut left = vec![0.0; p + 1];
        let mut right = vec![0.0; p + 1];
        n[0] = 1.0;
        for j in 1..=p {
            left[j] = u - u_knots[span + 1 - j];
            right[j] = u_knots[span + j] - u;
            let mut saved = 0.0;
            for r in 0..j {
                let temp = n[r] / (right[r + 1] + left[j - r]);
                n[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            n[j] = saved;
        }
        n
    }

    /// The basis functions non-zero on `span` and their derivatives up to
    /// `order`, at `u` (A2.3). `result[k][j]` is the `k`th derivative of the
    /// `j`th function; orders above the degree are zero.
    pub fn basis_derivatives(&self, span: usize, u: f64, order: usize) -> Vec<Vec<f64>> {
        let p = self.degree;
        let u_knots = &self.knots;
        let mut ders = vec![vec![0.0; p + 1]; order + 1];
        let order = order.min(p);

        // ndu[j][r] (j > r) holds knot differences; ndu[r][j] (r <= j) basis values.
        let mut ndu = vec![vec![0.0; p + 1]; p + 1];
        let mut left = vec![0.0; p + 1];
        let mut right = vec![0.0; p + 1];
        ndu[0][0] = 1.0;
        for j in 1..=p {
            left[j] = u - u_knots[span + 1 - j];
            right[j] = u_knots[span + j] - u;
            let mut saved = 0.0;
            for r in 0..j {
                ndu[j][r] = right[r + 1] + left[j - r];
                let temp = ndu[r][j - 1] / ndu[j][r];
                ndu[r][j] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            ndu[j][j] = saved;
        }
        for j in 0..=p {
            ders[0][j] = ndu[j][p];
        }

        let mut a = [vec![0.0; p + 1], vec![0.0; p + 1]];
        for r in 0..=p {
            let (mut s1, mut s2) = (0, 1);
            a[0][0] = 1.0;
            for k in 1..=order {
                let mut d = 0.0;
                let rk = r as isize - k as isize;
                let pk = p - k;
                if rk >= 0 {
                    a[s2][0] = a[s1][0] / ndu[pk + 1][rk as usize];
                    d = a[s2][0] * ndu[rk as usize][pk];
                }
                let j1 = if rk >= -1 { 1 } else { (-rk) as usize };
                let j2 = if r as isize - 1 <= pk as isize {
                    k - 1
                } else {
                    p - r
                };
                for j in j1..=j2 {
                    let idx = (rk + j as isize) as usize;
                    a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][idx];
                    d += a[s2][j] * ndu[idx][pk];
                }
                if r <= pk {
                    a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                    d += a[s2][k] * ndu[r][pk];
                }
                ders[k][r] = d;
                core::mem::swap(&mut s1, &mut s2);
            }
        }

        let mut factor = p as f64;
        for (k, row) in ders.iter_mut().enumerate().take(order + 1).skip(1) {
            for value in row.iter_mut() {
                *value *= factor;
            }
            factor *= (p - k) as f64;
        }
        ders
    }

    pub(crate) fn from_raw(degree: usize, knots: Vec<f64>) -> Self {
        debug_assert!(
            Self::new(degree, knots.clone()).is_ok(),
            "invalid knots {knots:?}"
        );
        Self { degree, knots }
    }
}

/// Inserts `u` into a B-spline `times` times (A5.1), returning the new knots and
/// control points. `u` must be strictly inside the domain, and `times` is
/// capped so the knot's multiplicity doesn't exceed the degree.
pub(crate) fn insert_knot<const D: usize>(
    knots: &KnotVector,
    points: &[Homogeneous<D>],
    u: f64,
    times: usize,
) -> Result<(KnotVector, Vec<Homogeneous<D>>), NurbsError> {
    let (start, end) = knots.domain();
    if !(start < u && u < end) {
        return Err(NurbsError::ParameterOutsideDomain { parameter: u });
    }
    let (knots, points) = insert_raw(knots, points, u, times);
    Ok((KnotVector::from_raw(knots.degree, knots.knots), points))
}

/// A5.1 with no domain check, so it also serves to clamp at the domain's ends.
/// `u` must lie in the domain; `times` is capped at `degree - multiplicity`.
pub(crate) fn insert_raw<const D: usize>(
    knots: &KnotVector,
    points: &[Homogeneous<D>],
    u: f64,
    times: usize,
) -> (KnotVector, Vec<Homogeneous<D>>) {
    let p = knots.degree;
    let s = knots.multiplicity(u);
    let r = times.min(p.saturating_sub(s));
    if r == 0 {
        return (knots.clone(), points.to_vec());
    }
    let up = &knots.knots;
    // The last knot at or before `u`; unlike `span`, not clamped to the domain,
    // so an unclamped end is inserted next to its own copies.
    let k = up.partition_point(|&x| x <= u) - 1;
    let np = points.len() - 1;

    let mut uq = Vec::with_capacity(up.len() + r);
    uq.extend_from_slice(&up[..=k]);
    uq.extend(core::iter::repeat_n(u, r));
    uq.extend_from_slice(&up[k + 1..]);

    let mut qw = vec![Homogeneous::zero(); np + 1 + r];
    qw[..=k - p].copy_from_slice(&points[..=k - p]);
    qw[k - s + r..].copy_from_slice(&points[k - s..]);
    let mut rw: Vec<_> = points[k - p..=k - s].to_vec();
    let mut l = 0;
    for j in 1..=r {
        l = k - p + j;
        for i in 0..=p - j - s {
            let alpha = (u - up[l + i]) / (up[i + k + 1] - up[l + i]);
            rw[i] = rw[i].lerp(rw[i + 1], alpha);
        }
        qw[l] = rw[0];
        qw[k + r - j - s] = rw[p - j - s];
    }
    if l + 1 < k - s {
        qw[l + 1..k - s].copy_from_slice(&rw[1..k - s - l]);
    }
    // Unvalidated: clamping inserts at the ends, where the multiplicity
    // briefly sits between the interior and end limits.
    (
        KnotVector {
            degree: p,
            knots: uq,
        },
        qw,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cubic() -> KnotVector {
        KnotVector::new(
            3,
            vec![0.0, 0.0, 0.0, 0.0, 0.3, 0.5, 0.5, 0.8, 1.0, 1.0, 1.0, 1.0],
        )
        .unwrap()
    }

    #[test]
    fn rejects_invalid_knots() {
        assert_eq!(
            KnotVector::new(0, vec![0.0, 1.0]),
            Err(NurbsError::ZeroDegree)
        );
        assert!(matches!(
            KnotVector::new(2, vec![0.0; 5]),
            Err(NurbsError::TooFewKnots { .. })
        ));
        assert_eq!(
            KnotVector::new(1, vec![0.0, 0.0, 1.0, 0.5]),
            Err(NurbsError::InvalidKnots)
        );
        assert_eq!(
            KnotVector::new(1, vec![0.0, 0.0, 0.0, 0.0]),
            Err(NurbsError::EmptyDomain)
        );
        assert!(matches!(
            KnotVector::new(2, vec![0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0]),
            Err(NurbsError::MultiplicityTooHigh { .. })
        ));
    }

    #[test]
    fn span_finds_containing_interval() {
        let k = cubic();
        assert_eq!(k.span(0.0), 3);
        assert_eq!(k.span(0.2), 3);
        assert_eq!(k.span(0.3), 4);
        assert_eq!(k.span(0.5), 6);
        assert_eq!(k.span(0.9), 7);
        assert_eq!(k.span(1.0), 7);
        assert_eq!(k.span(7.0), 7);
        assert_eq!(k.span(-1.0), 3);
    }

    #[test]
    fn basis_is_a_partition_of_unity() {
        let k = cubic();
        for i in 0..=50 {
            let u = i as f64 / 50.0;
            let basis = k.basis(k.span(u), u);
            assert!(basis.iter().all(|&b| b >= -1e-15));
            assert!((basis.iter().sum::<f64>() - 1.0).abs() < 1e-14);
        }
    }

    #[test]
    fn basis_derivatives_match_finite_differences() {
        let k = cubic();
        for u in [0.1, 0.35, 0.6, 0.95] {
            let span = k.span(u);
            let ders = k.basis_derivatives(span, u, 4);
            assert_eq!(ders[0], k.basis(span, u));
            let h = 1e-6;
            let lo = k.basis_derivatives(span, u - h, 3);
            let hi = k.basis_derivatives(span, u + h, 3);
            for order in 1..=3 {
                for j in 0..=3 {
                    let fd = (hi[order - 1][j] - lo[order - 1][j]) / (2.0 * h);
                    let scale = 1.0 + ders[order][j].abs();
                    assert!(
                        (fd - ders[order][j]).abs() < 1e-4 * scale,
                        "u {u} order {order} j {j}"
                    );
                }
            }
            assert!(ders[4].iter().all(|&d| d == 0.0));
            // Derivatives of a partition of unity sum to zero.
            for row in &ders[1..] {
                assert!(row.iter().sum::<f64>().abs() < 1e-9);
            }
        }
    }

    #[test]
    fn clamped_uniform_is_clamped() {
        let k = KnotVector::clamped_uniform(2, 5).unwrap();
        assert_eq!(
            k.knots(),
            &[0.0, 0.0, 0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0, 1.0, 1.0]
        );
        assert!(k.is_clamped());
        assert_eq!(k.basis_count(), 5);
        let unclamped = KnotVector::new(2, (0..8).map(f64::from).collect()).unwrap();
        assert!(!unclamped.is_clamped());
        assert_eq!(unclamped.domain(), (2.0, 5.0));
    }
}
