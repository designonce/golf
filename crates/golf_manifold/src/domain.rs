//! The extent of a parameter space: bounds, periodicity and singularities.

use crate::space::Point;
use crate::space::Space;

/// One parameter axis of a [`Domain`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Axis {
    /// Lower bound; may be `-∞`.
    pub min: f64,
    /// Upper bound; may be `+∞`. For a periodic axis, `max - min` is the period.
    pub max: f64,
    /// Whether `min` and `max` are the same place, e.g. longitude on a sphere.
    pub periodic: bool,
    /// Whether the image collapses at `min`: every other parameter there maps to
    /// one point, as longitude does at a sphere's pole.
    pub singular_at_min: bool,
    /// As [`Self::singular_at_min`], at `max`.
    pub singular_at_max: bool,
}

impl Axis {
    /// The whole real line.
    pub const UNBOUNDED: Axis = Axis {
        min: f64::NEG_INFINITY,
        max: f64::INFINITY,
        periodic: false,
        singular_at_min: false,
        singular_at_max: false,
    };

    /// `[min, max]`.
    pub fn bounded(min: f64, max: f64) -> Self {
        Axis {
            min,
            max,
            ..Self::UNBOUNDED
        }
    }

    /// `[min, max)`, with `max` wrapping round to `min`.
    pub fn periodic(min: f64, max: f64) -> Self {
        Axis {
            periodic: true,
            ..Self::bounded(min, max)
        }
    }

    pub fn singular_at_min(self) -> Self {
        Axis {
            singular_at_min: true,
            ..self
        }
    }

    pub fn singular_at_max(self) -> Self {
        Axis {
            singular_at_max: true,
            ..self
        }
    }

    pub fn period(&self) -> Option<f64> {
        self.periodic.then_some(self.max - self.min)
    }

    /// Brings a periodic `x` into `[min, max)`. Other axes are left alone.
    ///
    /// Values already in range come back unchanged, bit for bit.
    pub fn wrap(&self, x: f64) -> f64 {
        match self.period() {
            Some(_) if (self.min..self.max).contains(&x) => x,
            Some(period) => {
                let wrapped = self.min + (x - self.min).rem_euclid(period);
                // Rounding can land a value just below `min` on `max` itself.
                match wrapped < self.max {
                    true => wrapped,
                    false => self.min,
                }
            }
            None => x,
        }
    }

    /// [`Self::wrap`]s a periodic `x`, and clamps any other into `[min, max]`.
    pub fn clamp(&self, x: f64) -> f64 {
        match self.periodic {
            true => self.wrap(x),
            false => x.clamp(self.min, self.max),
        }
    }

    /// The equivalent of `x` nearest `reference`: shifted by whole periods on a
    /// periodic axis, unchanged on any other.
    ///
    /// Keeps a run of parameters continuous across a seam, e.g. along a pcurve.
    pub fn unwrap_near(&self, x: f64, reference: f64) -> f64 {
        match self.period() {
            Some(period) => x + period * ((reference - x) / period).round(),
            None => x,
        }
    }

    pub fn contains(&self, x: f64, tolerance: f64) -> bool {
        self.periodic || (self.min - tolerance <= x && x <= self.max + tolerance)
    }

    /// Whether `x` is within `tolerance` of a singular bound.
    pub fn is_singular(&self, x: f64, tolerance: f64) -> bool {
        (self.singular_at_min && (x - self.min).abs() <= tolerance)
            || (self.singular_at_max && (x - self.max).abs() <= tolerance)
    }
}

/// The extent of an `N`-dimensional parameter space, one [`Axis`] per parameter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Domain<const N: usize> {
    pub axes: [Axis; N],
}

impl<const N: usize> Domain<N> {
    pub fn new(axes: [Axis; N]) -> Self {
        Self { axes }
    }

    /// All of R^N.
    pub fn unbounded() -> Self {
        Self::new([Axis::UNBOUNDED; N])
    }

    /// [`Axis::wrap`] on every axis.
    pub fn wrap<S: Space<N>>(&self, p: Point<S, N>) -> Point<S, N> {
        self.per_axis(p, |axis, x| axis.wrap(x))
    }

    /// [`Axis::clamp`] on every axis.
    pub fn clamp<S: Space<N>>(&self, p: Point<S, N>) -> Point<S, N> {
        self.per_axis(p, |axis, x| axis.clamp(x))
    }

    /// [`Axis::unwrap_near`] on every axis.
    pub fn unwrap_near<S: Space<N>>(&self, p: Point<S, N>, reference: Point<S, N>) -> Point<S, N> {
        let mut coords = p.coords;
        for (i, axis) in self.axes.iter().enumerate() {
            coords[i] = axis.unwrap_near(coords[i], reference.coords[i]);
        }
        Point::with_tag(coords, p.tag())
    }

    pub fn contains<S: Space<N>>(&self, p: Point<S, N>, tolerance: f64) -> bool {
        self.axes
            .iter()
            .zip(p.coords.iter())
            .all(|(axis, &x)| axis.contains(x, tolerance))
    }

    /// Whether `p` is within `tolerance` of a singular bound on any axis.
    pub fn is_singular<S: Space<N>>(&self, p: Point<S, N>, tolerance: f64) -> bool {
        self.axes
            .iter()
            .zip(p.coords.iter())
            .any(|(axis, &x)| axis.is_singular(x, tolerance))
    }

    fn per_axis<S: Space<N>>(&self, p: Point<S, N>, f: impl Fn(&Axis, f64) -> f64) -> Point<S, N> {
        let mut coords = p.coords;
        for (x, axis) in coords.iter_mut().zip(&self.axes) {
            *x = f(axis, *x);
        }
        Point::with_tag(coords, p.tag())
    }
}

#[cfg(test)]
mod tests {
    use core::f64::consts::PI;

    use super::*;

    #[test]
    fn periodic_wrap_is_half_open() {
        let axis = Axis::periodic(-PI, PI);
        assert_eq!(axis.wrap(PI), -PI);
        assert!((axis.wrap(3.0 * PI + 0.5) - (-PI + 0.5)).abs() < 1e-12);
        assert!((axis.wrap(-PI - 0.5) - (PI - 0.5)).abs() < 1e-12);
    }

    #[test]
    fn wrap_leaves_in_range_values_exact_and_stays_half_open() {
        let axis = Axis::periodic(-PI, PI);
        assert_eq!(axis.wrap(0.4), 0.4);
        let tiny_below = -PI - 1e-17;
        assert!(axis.wrap(tiny_below) < PI);
    }

    #[test]
    fn unwrap_near_crosses_the_seam() {
        let axis = Axis::periodic(-PI, PI);
        let x = axis.unwrap_near(-PI + 0.1, PI - 0.1);
        assert!((x - (PI + 0.1)).abs() < 1e-12);
    }

    #[test]
    fn bounded_axis_clamps_and_ignores_periodic_ops() {
        let axis = Axis::bounded(0.0, 1.0);
        assert_eq!(axis.clamp(1.5), 1.0);
        assert_eq!(axis.wrap(1.5), 1.5);
        assert_eq!(axis.unwrap_near(1.5, 10.0), 1.5);
        assert!(axis.contains(1.05, 0.1));
        assert!(!axis.contains(1.2, 0.1));
    }
}
