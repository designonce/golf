//! Checks shared by the geometry tests.

use golf_manifold::Embedding;
use golf_manifold::Mapping;
use golf_manifold::Point;
use golf_manifold::Space;
use golf_manifold::newton_project;
use nalgebra::SVector;

fn param<S: Space<N, Tag = ()>, const N: usize>(at: [f64; N]) -> Point<S, N> {
    Point::new(SVector::from(at))
}

/// The jacobian at `at` agrees with central differences of `apply`.
pub(crate) fn assert_jacobian<E, const N: usize, const M: usize>(e: &E, at: [f64; N])
where
    E: Mapping<N, M>,
    E::From: Space<N, Tag = ()>,
{
    let x: Point<E::From, N> = param(at);
    let j = e.jacobian(x);
    let h = 1e-6;
    for k in 0..N {
        let mut d = SVector::<f64, N>::zeros();
        d[k] = h;
        let fd = (e.apply(Point::new(x.coords + d)) - e.apply(Point::new(x.coords - d))).coords
            / (2.0 * h);
        assert!(
            (j.column(k) - fd).norm() < 1e-7 * (1.0 + fd.norm()),
            "column {k} at {at:?}: {} != {fd}",
            j.column(k)
        );
    }
}

/// Projecting `apply(at)` with no hint gives back `at`, up to whole periods.
pub(crate) fn assert_round_trip<E, const N: usize, const M: usize>(e: &E, at: [f64; N])
where
    E: Embedding<N, M>,
    E::From: Space<N, Tag = ()>,
{
    let x: Point<E::From, N> = param(at);
    let q = e.project(e.apply(x), None).unwrap();
    assert!(e.domain().contains(q, 0.0), "{q:?} outside domain");
    let q = e.domain().unwrap_near(q, x);
    assert!(
        (q.coords - x.coords).norm() < 1e-9,
        "{at:?} came back as {q:?}"
    );
}

/// `project(p)` agrees with Newton iteration started at `start`.
pub(crate) fn assert_matches_newton<E, const N: usize, const M: usize>(
    e: &E,
    p: Point<E::To, M>,
    start: [f64; N],
) where
    E: Embedding<N, M>,
    E::From: Space<N, Tag = ()>,
{
    let exact = e.project(p, None).unwrap();
    let newton = newton_project(e, p, param(start)).unwrap();
    let newton = e.domain().unwrap_near(newton, exact);
    assert!(
        (newton.coords - exact.coords).norm() < 1e-8,
        "{p:?}: newton {newton:?} != closed form {exact:?}"
    );
}
