//! Manifolds: maps from a parameter space that can be (locally) inverted.

use nalgebra::SVector;

use crate::domain::Domain;
use crate::mapping::Mapping;
use crate::space::Point;
use crate::space::Vector;

/// An embedding of an `N`-dimensional parameter space into an `M`-dimensional one.
pub trait Embedding<const N: usize, const M: usize>: Mapping<N, M> {
    /// The parameters `apply` is defined on.
    fn domain(&self) -> Domain<N>;

    /// The parameter of the point on the image nearest `p`, wrapped into the [`domain`].
    ///
    /// `hint` is a nearby parameter. Iterative implementations start from it, and
    /// where several parameters are equally near (a sphere's centre, or longitude
    /// at a pole) the one nearest the hint is chosen.
    ///
    /// [`domain`]: Self::domain
    fn project(
        &self,
        p: Point<Self::To, M>,
        hint: Option<Point<Self::From, N>>,
    ) -> Result<Point<Self::From, N>, ProjectError>;

    /// [`project`](Self::project), failing if `p` is further than `tolerance` from the image.
    fn project_within(
        &self,
        p: Point<Self::To, M>,
        hint: Option<Point<Self::From, N>>,
        tolerance: f64,
    ) -> Result<Point<Self::From, N>, ProjectError> {
        let param = self.project(p, hint)?;
        let distance = (self.apply(param) - p).coords.norm();
        match distance <= tolerance {
            true => Ok(param),
            false => Err(ProjectError::TooFar {
                distance,
                tolerance,
            }),
        }
    }
}

/// Why [`Embedding::project`] found no parameter.
#[derive(Clone, Copy, Debug, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum ProjectError {
    /// Every parameter in some set is equally near, and no hint picked one.
    #[error("projection is ambiguous and no hint was given")]
    Ambiguous,
    /// An iterative solver did not settle.
    #[error("projection did not converge in {iterations} iterations")]
    NoConvergence { iterations: usize },
    /// The nearest point is further away than allowed.
    #[error("point is {distance} from the manifold, more than {tolerance}")]
    TooFar { distance: f64, tolerance: f64 },
}

/// Finds the parameter nearest `p` by Newton iteration from `start`.
///
/// Converges to a local minimum of distance, so `start` must already be close;
/// for an implementation of [`Embedding::project`] that has no closed form.
/// Bounded axes are clamped and periodic ones wrapped after each step.
///
/// Second derivatives are taken by differencing [`Mapping::jacobian`]. Where they
/// don't give a descent direction (near a singularity, or past a focal point) it
/// falls back to a damped Gauss-Newton step.
pub fn newton_project<E, const N: usize, const M: usize>(
    embedding: &E,
    p: Point<E::To, M>,
    start: Point<E::From, N>,
) -> Result<Point<E::From, N>, ProjectError>
where
    E: Embedding<N, M> + ?Sized,
{
    const MAX_ITERATIONS: usize = 64;
    // Relative damping keeps JᵀJ invertible where the image collapses (J loses rank).
    const DAMPING: f64 = 1e-12;

    let domain = embedding.domain();
    let offset = |x: Point<E::From, N>, step: SVector<f64, N>| {
        domain.clamp(Point::with_tag(x.coords + step, x.tag()))
    };
    let objective = |x| (embedding.apply(x) - p).coords.norm_squared() / 2.0;

    let mut x = domain.clamp(start);
    for _ in 0..MAX_ITERATIONS {
        let residual = (embedding.apply(x) - p).coords;
        let j = embedding.jacobian(x);
        let gradient = j.transpose() * residual;
        let jtj = j.transpose() * j;

        let mut hessian = jtj;
        for k in 0..N {
            let h = 1e-6 * (1.0 + x.coords[k].abs());
            let mut e = SVector::<f64, N>::zeros();
            e[k] = h;
            let dj =
                (embedding.jacobian(offset(x, e)) - embedding.jacobian(offset(x, -e))) / (2.0 * h);
            hessian.set_row(
                k,
                &(hessian.row(k) + (dj.transpose() * residual).transpose()),
            );
        }
        let step = match hessian.cholesky() {
            Some(newton) => -newton.solve(&gradient),
            None => {
                let mut damped = jtj;
                let damping = DAMPING * (jtj.trace() + f64::MIN_POSITIVE);
                for i in 0..N {
                    damped[(i, i)] += damping;
                }
                match damped.cholesky() {
                    Some(gauss_newton) => -gauss_newton.solve(&gradient),
                    None => break,
                }
            }
        };

        // Backtrack until the step reduces the distance (Armijo condition).
        let current = objective(x);
        let slope = gradient.dot(&step);
        let mut alpha = 1.0;
        let next = loop {
            let next = offset(x, step * alpha);
            if objective(next) <= current + 1e-4 * alpha * slope {
                break Some(next);
            }
            if alpha < 1e-10 {
                break None;
            }
            alpha /= 2.0;
        };
        // No step reduces the distance: this is the nearest point, to the
        // precision the distance can be computed with.
        let Some(next) = next else {
            return Ok(domain.wrap(x));
        };

        // Converged when the step moves the parameters, or their image, by a
        // rounding error. (Near the solution, the finite-difference Hessian
        // can keep steps a little above the parameters' own rounding.)
        let moved = (next.coords - x.coords).norm();
        let image_moved = (j * (next.coords - x.coords)).norm();
        let scale = 1.0 + p.coords.norm();
        x = next;
        if moved <= 1e-14 * (1.0 + x.coords.norm()) || image_moved <= 1e-13 * scale {
            return Ok(domain.wrap(x));
        }
    }
    Err(ProjectError::NoConvergence {
        iterations: MAX_ITERATIONS,
    })
}

/// A curve in some 3D space.
pub trait Curve: Embedding<1, 3> {
    fn tangent(&self, t: Point<Self::From, 1>) -> Vector<Self::To, 3> {
        Vector::with_tag(self.jacobian(t).column(0).normalize(), self.apply(t).tag())
    }
}
impl<C: Embedding<1, 3>> Curve for C {}

/// A surface in some 3D space.
pub trait Surface: Embedding<2, 3> {
    /// Unit normal, `∂u × ∂v`. Undefined at singular points such as a sphere's poles.
    fn normal(&self, uv: Point<Self::From, 2>) -> Vector<Self::To, 3> {
        let j = self.jacobian(uv);
        Vector::with_tag(
            j.column(0).cross(&j.column(1)).normalize(),
            self.apply(uv).tag(),
        )
    }
}
impl<S: Embedding<2, 3>> Surface for S {}
