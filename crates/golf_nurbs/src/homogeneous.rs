//! Control points in homogeneous form, `(w·P, w)`.
//!
//! Rational algorithms are the polynomial ones applied to these, which sidesteps
//! needing an `SVector<f64, D + 1>` (not expressible on stable Rust).

use core::ops::Add;
use core::ops::AddAssign;
use core::ops::Mul;
use core::ops::Sub;

use nalgebra::SVector;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Homogeneous<const D: usize> {
    /// The point scaled by its weight.
    pub wp: SVector<f64, D>,
    pub w: f64,
}

impl<const D: usize> Homogeneous<D> {
    pub fn new(point: SVector<f64, D>, weight: f64) -> Self {
        Self {
            wp: point * weight,
            w: weight,
        }
    }

    pub fn zero() -> Self {
        Self {
            wp: SVector::zeros(),
            w: 0.0,
        }
    }

    pub fn point(&self) -> SVector<f64, D> {
        self.wp / self.w
    }

    /// `(1 - t)·self + t·other`.
    pub fn lerp(self, other: Self, t: f64) -> Self {
        self * (1.0 - t) + other * t
    }
}

impl<const D: usize> Add for Homogeneous<D> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self {
            wp: self.wp + rhs.wp,
            w: self.w + rhs.w,
        }
    }
}

impl<const D: usize> AddAssign for Homogeneous<D> {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl<const D: usize> Sub for Homogeneous<D> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self {
            wp: self.wp - rhs.wp,
            w: self.w - rhs.w,
        }
    }
}

impl<const D: usize> Mul<f64> for Homogeneous<D> {
    type Output = Self;
    fn mul(self, k: f64) -> Self {
        Self {
            wp: self.wp * k,
            w: self.w * k,
        }
    }
}

/// `n choose k`, exact for the small values B-spline degrees reach.
pub(crate) fn binomial(n: usize, k: usize) -> f64 {
    (0..k.min(n - k)).fold(1.0, |acc, i| acc * (n - i) as f64 / (i + 1) as f64)
}
