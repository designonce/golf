//! Smooth maps between spaces.

use nalgebra::SMatrix;

use crate::space::Point;
use crate::space::Space;
use crate::space::Vector;

/// A differentiable map from an `N`-dimensional space to an `M`-dimensional one.
pub trait Mapping<const N: usize, const M: usize> {
    type From: Space<N>;
    type To: Space<M>;

    fn apply(&self, p: Point<Self::From, N>) -> Point<Self::To, M>;

    /// The Jacobian at `p`; its columns are the partial derivatives of `apply`.
    fn jacobian(&self, p: Point<Self::From, N>) -> SMatrix<f64, M, N>;

    /// Carries a tangent vector at `p` across the map.
    fn push_forward(
        &self,
        p: Point<Self::From, N>,
        v: Vector<Self::From, N>,
    ) -> Vector<Self::To, M> {
        Vector::with_tag(self.jacobian(p) * v.coords, self.apply(p).tag())
    }
}

/// `G ∘ F`, passing through an intermediate space of dimension `B`.
///
/// Only typechecks when `G` starts in the space `F` ends in.
pub struct Compose<F, G, const B: usize>(pub F, pub G);

impl<F, G, const A: usize, const B: usize, const C: usize> Mapping<A, C> for Compose<F, G, B>
where
    F: Mapping<A, B>,
    G: Mapping<B, C, From = F::To>,
{
    type From = F::From;
    type To = G::To;

    fn apply(&self, p: Point<F::From, A>) -> Point<G::To, C> {
        self.1.apply(self.0.apply(p))
    }

    fn jacobian(&self, p: Point<F::From, A>) -> SMatrix<f64, C, A> {
        self.1.jacobian(self.0.apply(p)) * self.0.jacobian(p)
    }
}
