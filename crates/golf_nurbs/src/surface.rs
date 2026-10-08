//! NURBS surfaces.

use nalgebra::SVector;

use crate::curve::NurbsCurve;
use crate::curve::homogenise;
use crate::error::NurbsError;
use crate::homogeneous::Homogeneous;
use crate::homogeneous::binomial;
use crate::knots::KnotVector;
use crate::knots::{self};

/// A tensor-product NURBS surface in `D` dimensions.
///
/// Control points form a grid indexed `(i, j)`, `i` along `u` and `j` along `v`.
#[derive(Clone, Debug, PartialEq)]
pub struct NurbsSurface<const D: usize> {
    knots_u: KnotVector,
    knots_v: KnotVector,
    /// Row-major: `(i, j)` is at `i * count_v + j`.
    points: Vec<Homogeneous<D>>,
}

impl<const D: usize> NurbsSurface<D> {
    /// A surface with `points[i][j]` as control point `(i, j)`, optionally
    /// weighted by `weights[i][j]` (all 1 if `None`).
    pub fn new(
        (degree_u, knots_u): (usize, Vec<f64>),
        (degree_v, knots_v): (usize, Vec<f64>),
        points: Vec<Vec<SVector<f64, D>>>,
        weights: Option<Vec<Vec<f64>>>,
    ) -> Result<Self, NurbsError> {
        let knots_u = KnotVector::new(degree_u, knots_u)?;
        let knots_v = KnotVector::new(degree_v, knots_v)?;
        let (count_u, count_v) = (knots_u.basis_count(), knots_v.basis_count());
        if points.len() != count_u {
            return Err(NurbsError::ControlPointCount {
                expected: count_u,
                found: points.len(),
            });
        }
        if let Some(row) = points.iter().find(|row| row.len() != count_v) {
            return Err(NurbsError::ControlPointCount {
                expected: count_v,
                found: row.len(),
            });
        }
        let flat: Vec<_> = points.into_iter().flatten().collect();
        let weights = match weights {
            Some(weights) => {
                if weights.len() != count_u || weights.iter().any(|row| row.len() != count_v) {
                    return Err(NurbsError::ControlPointCount {
                        expected: count_u * count_v,
                        found: weights.iter().map(Vec::len).sum(),
                    });
                }
                Some(weights.into_iter().flatten().collect::<Vec<_>>())
            }
            None => None,
        };
        let points = homogenise(&flat, weights.as_deref())?;
        Ok(Self {
            knots_u,
            knots_v,
            points,
        })
    }

    pub fn degrees(&self) -> (usize, usize) {
        (self.knots_u.degree(), self.knots_v.degree())
    }

    pub fn knots_u(&self) -> &KnotVector {
        &self.knots_u
    }

    pub fn knots_v(&self) -> &KnotVector {
        &self.knots_v
    }

    /// The parameter ranges in `u` and `v`.
    pub fn domain(&self) -> ((f64, f64), (f64, f64)) {
        (self.knots_u.domain(), self.knots_v.domain())
    }

    /// The control grid's size, `(count_u, count_v)`.
    pub fn control_grid_size(&self) -> (usize, usize) {
        (self.knots_u.basis_count(), self.knots_v.basis_count())
    }

    pub fn control_point(&self, i: usize, j: usize) -> SVector<f64, D> {
        self.at(i, j).point()
    }

    pub fn weight(&self, i: usize, j: usize) -> f64 {
        self.at(i, j).w
    }

    /// Whether any weight differs from the others.
    pub fn is_rational(&self) -> bool {
        self.points.iter().any(|p| p.w != self.points[0].w)
    }

    /// The point at `(u, v)`, each clamped into its domain.
    pub fn point(&self, u: f64, v: f64) -> SVector<f64, D> {
        self.derivatives(u, v, 0)[0][0]
    }

    /// The point at `(u, v)` and its partial derivatives up to total `order`:
    /// `result[k][l]` is `∂^(k+l) S / ∂u^k ∂v^l`, for `k + l <= order`
    /// (A3.6 and A4.4).
    pub fn derivatives(&self, u: f64, v: f64, order: usize) -> Vec<Vec<SVector<f64, D>>> {
        let (p, q) = self.degrees();
        let ((u0, u1), (v0, v1)) = self.domain();
        let (span_u, span_v) = (self.knots_u.span(u), self.knots_v.span(v));
        let (u, v) = (u.clamp(u0, u1), v.clamp(v0, v1));
        let nu = self.knots_u.basis_derivatives(span_u, u, order);
        let nv = self.knots_v.basis_derivatives(span_v, v, order);

        // Homogeneous derivatives, A3.6.
        let mut a = vec![vec![Homogeneous::zero(); order + 1]; order + 1];
        for k in 0..=order {
            // Combine along u first: one homogeneous point per v-index.
            let column: Vec<Homogeneous<D>> = (0..=q)
                .map(|s| {
                    (0..=p).fold(Homogeneous::zero(), |acc, r| {
                        acc + self.at(span_u - p + r, span_v - q + s) * nu[k][r]
                    })
                })
                .collect();
            for l in 0..=order - k {
                a[k][l] = (0..=q).fold(Homogeneous::zero(), |acc, s| acc + column[s] * nv[l][s]);
            }
        }

        // Rational derivatives, A4.4.
        let w = |k: usize, l: usize| a[k][l].w;
        let mut skl = vec![Vec::new(); order + 1];
        for k in 0..=order {
            for l in 0..=order - k {
                let mut value = a[k][l].wp;
                for j in 1..=l {
                    value -= skl[k][l - j] * (binomial(l, j) * w(0, j));
                }
                for i in 1..=k {
                    value -= skl[k - i][l] * (binomial(k, i) * w(i, 0));
                    let mut inner = SVector::<f64, D>::zeros();
                    for j in 1..=l {
                        inner += skl[k - i][l - j] * (binomial(l, j) * w(i, j));
                    }
                    value -= inner * binomial(k, i);
                }
                skl[k].push(value / w(0, 0));
            }
        }
        skl
    }

    /// The curve `v ↦ S(u, v)` at fixed `u`.
    pub fn isocurve_u(&self, u: f64) -> NurbsCurve<D> {
        let p = self.knots_u.degree();
        let (u0, u1) = self.knots_u.domain();
        let span = self.knots_u.span(u);
        let basis = self.knots_u.basis(span, u.clamp(u0, u1));
        let count_v = self.knots_v.basis_count();
        let points = (0..count_v)
            .map(|j| {
                (0..=p).fold(Homogeneous::zero(), |acc, r| {
                    acc + self.at(span - p + r, j) * basis[r]
                })
            })
            .collect();
        NurbsCurve {
            knots: self.knots_v.clone(),
            points,
        }
    }

    /// The curve `u ↦ S(u, v)` at fixed `v`.
    pub fn isocurve_v(&self, v: f64) -> NurbsCurve<D> {
        self.transposed().isocurve_u(v)
    }

    /// The same surface with `u` and `v` swapped.
    pub fn transposed(&self) -> Self {
        let (count_u, count_v) = self.control_grid_size();
        let points = (0..count_v)
            .flat_map(|j| (0..count_u).map(move |i| (i, j)))
            .map(|(i, j)| self.at(i, j))
            .collect();
        Self {
            knots_u: self.knots_v.clone(),
            knots_v: self.knots_u.clone(),
            points,
        }
    }

    /// The same surface with `u` inserted into the `u` knots `times` more times
    /// (capped at the degree). `u` must be strictly inside the `u` domain.
    pub fn insert_knot_u(&self, u: f64, times: usize) -> Result<Self, NurbsError> {
        let (count_u, count_v) = self.control_grid_size();
        let mut knots_u = None;
        let mut columns = Vec::with_capacity(count_v);
        for j in 0..count_v {
            let column: Vec<_> = (0..count_u).map(|i| self.at(i, j)).collect();
            let (new_knots, new_column) = knots::insert_knot(&self.knots_u, &column, u, times)?;
            knots_u = Some(new_knots);
            columns.push(new_column);
        }
        let knots_u = knots_u.expect("a surface has at least one column");
        let new_count_u = knots_u.basis_count();
        let points = (0..new_count_u)
            .flat_map(|i| columns.iter().map(move |column| column[i]))
            .collect();
        Ok(Self {
            knots_u,
            knots_v: self.knots_v.clone(),
            points,
        })
    }

    /// As [`Self::insert_knot_u`], in `v`.
    pub fn insert_knot_v(&self, v: f64, times: usize) -> Result<Self, NurbsError> {
        Ok(self.transposed().insert_knot_u(v, times)?.transposed())
    }

    /// Applies `f` to every control point. Only meaningful for affine `f`.
    pub fn map_points<const E: usize>(
        &self,
        f: impl Fn(SVector<f64, D>) -> SVector<f64, E>,
    ) -> NurbsSurface<E> {
        NurbsSurface {
            knots_u: self.knots_u.clone(),
            knots_v: self.knots_v.clone(),
            points: self
                .points
                .iter()
                .map(|h| Homogeneous::new(f(h.point()), h.w))
                .collect(),
        }
    }

    fn at(&self, i: usize, j: usize) -> Homogeneous<D> {
        self.points[i * self.knots_v.basis_count() + j]
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::Vector3;

    use super::*;

    /// A biquadratic rational surface with interior knots in both directions.
    fn wavy() -> NurbsSurface<3> {
        let points = (0..4)
            .map(|i| {
                (0..5)
                    .map(|j| Vector3::new(i as f64, j as f64, ((i * 3 + j * 7) % 5) as f64 * 0.3))
                    .collect()
            })
            .collect();
        let weights = (0..4)
            .map(|i| {
                (0..5)
                    .map(|j| 1.0 + ((i + 2 * j) % 3) as f64 * 0.4)
                    .collect()
            })
            .collect();
        NurbsSurface::new(
            (2, vec![0.0, 0.0, 0.0, 0.4, 1.0, 1.0, 1.0]),
            (2, vec![0.0, 0.0, 0.0, 0.3, 0.6, 1.0, 1.0, 1.0]),
            points,
            Some(weights),
        )
        .unwrap()
    }

    fn grid() -> impl Iterator<Item = (f64, f64)> {
        (0..=12).flat_map(|i| (0..=12).map(move |j| (i as f64 / 12.0, j as f64 / 12.0)))
    }

    fn assert_same_shape(a: &NurbsSurface<3>, b: &NurbsSurface<3>) {
        for (u, v) in grid() {
            assert!((a.point(u, v) - b.point(u, v)).norm() < 1e-11, "({u}, {v})");
        }
    }

    #[test]
    fn rejects_ragged_grids() {
        let row = vec![Vector3::zeros(); 2];
        let knots = (1, vec![0.0, 0.0, 1.0, 1.0]);
        assert!(
            NurbsSurface::new(
                knots.clone(),
                knots.clone(),
                vec![row.clone(), row.clone()],
                None
            )
            .is_ok()
        );
        assert!(matches!(
            NurbsSurface::new(
                knots.clone(),
                knots.clone(),
                vec![row.clone(), row[..1].to_vec()],
                None
            ),
            Err(NurbsError::ControlPointCount { .. })
        ));
        assert!(matches!(
            NurbsSurface::new(
                knots.clone(),
                knots,
                vec![row.clone(), row],
                Some(vec![vec![1.0; 2]])
            ),
            Err(NurbsError::ControlPointCount { .. })
        ));
    }

    #[test]
    fn bilinear_patch_interpolates() {
        let c = |x: f64, y: f64, z: f64| Vector3::new(x, y, z);
        let s = NurbsSurface::new(
            (1, vec![0.0, 0.0, 1.0, 1.0]),
            (1, vec![0.0, 0.0, 1.0, 1.0]),
            vec![
                vec![c(0.0, 0.0, 0.0), c(0.0, 1.0, 1.0)],
                vec![c(1.0, 0.0, 1.0), c(1.0, 1.0, 0.0)],
            ],
            None,
        )
        .unwrap();
        for (u, v) in grid() {
            let expected = c(u, v, u + v - 2.0 * u * v);
            assert!((s.point(u, v) - expected).norm() < 1e-14);
        }
    }

    #[test]
    fn rational_cylinder_patch_is_exact() {
        let arc =
            NurbsCurve::circular_arc(Vector3::zeros(), Vector3::x(), Vector3::y(), 2.0, 0.0, 2.0);
        let count = arc.control_point_count();
        let points = (0..count)
            .map(|i| {
                vec![
                    arc.control_point(i),
                    arc.control_point(i) + Vector3::new(0.0, 0.0, 3.0),
                ]
            })
            .collect();
        let weights = (0..count).map(|i| vec![arc.weight(i); 2]).collect();
        let s = NurbsSurface::new(
            (2, arc.knots().knots().to_vec()),
            (1, vec![0.0, 0.0, 1.0, 1.0]),
            points,
            Some(weights),
        )
        .unwrap();
        for (u, v) in grid() {
            let p = s.point(u, v);
            assert!((p.xy().norm() - 2.0).abs() < 1e-12);
            assert!((p.z - 3.0 * v).abs() < 1e-12);
        }
    }

    #[test]
    fn derivatives_match_finite_differences() {
        let s = wavy();
        let h = 1e-5;
        for (u, v) in [(0.1, 0.1), (0.35, 0.5), (0.7, 0.2), (0.9, 0.85)] {
            let d = s.derivatives(u, v, 3);
            assert_eq!(d[0][0], s.point(u, v));
            for k in 0..=2 {
                for l in 0..=2 - k {
                    let du = (s.derivatives(u + h, v, 2)[k][l] - s.derivatives(u - h, v, 2)[k][l])
                        / (2.0 * h);
                    let dv = (s.derivatives(u, v + h, 2)[k][l] - s.derivatives(u, v - h, 2)[k][l])
                        / (2.0 * h);
                    let tol = |x: &SVector<f64, 3>| 1e-4 * (1.0 + x.norm());
                    assert!(
                        (du - d[k + 1][l]).norm() < tol(&d[k + 1][l]),
                        "({u},{v}) d{k}{l}/du"
                    );
                    assert!(
                        (dv - d[k][l + 1]).norm() < tol(&d[k][l + 1]),
                        "({u},{v}) d{k}{l}/dv"
                    );
                }
            }
        }
    }

    #[test]
    fn isocurves_lie_on_the_surface() {
        let s = wavy();
        for t in [0.0, 0.25, 0.4, 0.8, 1.0] {
            let (cu, cv) = (s.isocurve_u(t), s.isocurve_v(t));
            for i in 0..=20 {
                let w = i as f64 / 20.0;
                assert!((cu.point(w) - s.point(t, w)).norm() < 1e-12);
                assert!((cv.point(w) - s.point(w, t)).norm() < 1e-12);
            }
        }
    }

    #[test]
    fn knot_insertion_preserves_shape() {
        let s = wavy();
        let refined = s
            .insert_knot_u(0.4, 1)
            .unwrap()
            .insert_knot_u(0.75, 2)
            .unwrap();
        assert_eq!(refined.control_grid_size(), (7, 5));
        assert_same_shape(&s, &refined);
        let refined = s
            .insert_knot_v(0.3, 1)
            .unwrap()
            .insert_knot_v(0.5, 2)
            .unwrap();
        assert_eq!(refined.control_grid_size(), (4, 8));
        assert_same_shape(&s, &refined);
        assert!(s.insert_knot_u(1.0, 1).is_err());
    }

    #[test]
    fn transposing_swaps_parameters() {
        let s = wavy();
        let t = s.transposed();
        for (u, v) in grid() {
            assert!((t.point(v, u) - s.point(u, v)).norm() < 1e-13);
        }
        assert_eq!(t.transposed(), s);
    }
}
