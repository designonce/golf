//! Geometry whose concrete type is only known at runtime.
//!
//! Topology holds curves and surfaces of many kinds side by side, so it stores
//! them as [`AnySurface`], [`AnyCurve`] and [`AnyCurve2`] (a pcurve). Each is an
//! ordinary [`Embedding`] with its own parameter space (`Uv<AnySurface<S>>` and
//! so on): code that needs the concrete shape, such as an analytic intersection,
//! matches on the variant; everything else uses the traits.
//!
//! A typed parameter converts to the erased one by its coordinates, e.g.
//! `Point::new(uv.coords)`, since every variant's parameter space is untagged.
//!
//! Geometry defined outside this crate goes in a `Custom` variant through
//! [`DynEmbedding`], which every suitable [`Embedding`] implements.

use std::sync::Arc;

use nalgebra::SMatrix;
use nalgebra::SVector;

use crate::domain::Domain;
use crate::geom::Circle;
use crate::geom::Cone;
use crate::geom::Cylinder;
use crate::geom::Ellipse;
use crate::geom::Line;
use crate::geom::NurbsCurve;
use crate::geom::NurbsSurface;
use crate::geom::Plane;
use crate::geom::Sphere;
use crate::geom::Torus;
use crate::manifold::Embedding;
use crate::manifold::ProjectError;
use crate::mapping::Mapping;
use crate::space::Point;
use crate::space::Space;
use crate::space::T;
use crate::space::Uv;

/// An object-safe [`Embedding`] of an `N`-dimensional parameter space into `S`,
/// with parameters as bare coordinates.
///
/// Implemented for every `Embedding` whose parameter space is untagged. The
/// methods are prefixed `dyn_` so they don't clash with `Mapping`'s and
/// `Embedding`'s on the same types.
pub trait DynEmbedding<S: Space<M>, const N: usize, const M: usize>: core::fmt::Debug {
    fn dyn_apply(&self, param: SVector<f64, N>) -> Point<S, M>;
    fn dyn_jacobian(&self, param: SVector<f64, N>) -> SMatrix<f64, M, N>;
    fn dyn_domain(&self) -> Domain<N>;
    fn dyn_project(
        &self,
        p: Point<S, M>,
        hint: Option<SVector<f64, N>>,
    ) -> Result<SVector<f64, N>, ProjectError>;
}

impl<E, const N: usize, const M: usize> DynEmbedding<E::To, N, M> for E
where
    E: Embedding<N, M> + core::fmt::Debug,
    E::From: Space<N, Tag = ()>,
{
    fn dyn_apply(&self, param: SVector<f64, N>) -> Point<E::To, M> {
        Mapping::apply(self, Point::new(param))
    }

    fn dyn_jacobian(&self, param: SVector<f64, N>) -> SMatrix<f64, M, N> {
        Mapping::jacobian(self, Point::new(param))
    }

    fn dyn_domain(&self) -> Domain<N> {
        Embedding::domain(self)
    }

    fn dyn_project(
        &self,
        p: Point<E::To, M>,
        hint: Option<SVector<f64, N>>,
    ) -> Result<SVector<f64, N>, ProjectError> {
        Embedding::project(self, p, hint.map(Point::new)).map(|q| q.coords)
    }
}

/// Defines an erased geometry enum: its variants, `From` conversions, parameter
/// space, and `Mapping`/`Embedding` by delegating to the variant.
macro_rules! erased {
    (
        $(#[$meta:meta])*
        $name:ident<S: Space<$m:literal>>, $param:ident<_>, $n:literal,
        { $($(#[$vmeta:meta])* $variant:ident($ty:ty)),* $(,)? }
    ) => {
        $(#[$meta])*
        #[derive(Clone, Debug)]
        #[non_exhaustive]
        pub enum $name<S: Space<$m>> {
            $($(#[$vmeta])* $variant($ty),)*
            /// Any other geometry, e.g. a type defined outside this crate.
            Custom(Arc<dyn DynEmbedding<S, $n, $m> + Send + Sync>),
        }

        $(
            impl<S: Space<$m>> From<$ty> for $name<S> {
                fn from(geometry: $ty) -> Self {
                    Self::$variant(geometry)
                }
            }
        )*

        impl<S: Space<$m>> Space<$n> for $param<$name<S>> {
            type Tag = ();
        }

        impl<S: Space<$m>> Mapping<$n, $m> for $name<S> {
            type From = $param<$name<S>>;
            type To = S;

            fn apply(&self, p: Point<Self::From, $n>) -> Point<S, $m> {
                match self {
                    $(Self::$variant(g) => g.apply(Point::new(p.coords)),)*
                    Self::Custom(g) => g.dyn_apply(p.coords),
                }
            }

            fn jacobian(&self, p: Point<Self::From, $n>) -> SMatrix<f64, $m, $n> {
                match self {
                    $(Self::$variant(g) => g.jacobian(Point::new(p.coords)),)*
                    Self::Custom(g) => g.dyn_jacobian(p.coords),
                }
            }
        }

        impl<S: Space<$m>> Embedding<$n, $m> for $name<S> {
            fn domain(&self) -> Domain<$n> {
                match self {
                    $(Self::$variant(g) => g.domain(),)*
                    Self::Custom(g) => g.dyn_domain(),
                }
            }

            fn project(
                &self,
                p: Point<S, $m>,
                hint: Option<Point<Self::From, $n>>,
            ) -> Result<Point<Self::From, $n>, ProjectError> {
                let q = match self {
                    $(
                        Self::$variant(g) => g
                            .project(p, hint.map(|h| Point::new(h.coords)))
                            .map(|q| q.coords),
                    )*
                    Self::Custom(g) => g.dyn_project(p, hint.map(|h| h.coords)),
                };
                q.map(Point::new)
            }
        }
    };
}

erased! {
    /// Any surface in a 3D space `S`.
    AnySurface<S: Space<3>>, Uv<_>, 2, {
        Plane(Plane<S>),
        Sphere(Sphere<S>),
        Cylinder(Cylinder<S>),
        Cone(Cone<S>),
        Torus(Torus<S>),
        Nurbs(NurbsSurface<S>),
    }
}

erased! {
    /// Any curve in a 3D space `S`.
    AnyCurve<S: Space<3>>, T<_>, 1, {
        Line(Line<S, 3>),
        Circle(Circle<S>),
        Ellipse(Ellipse<S>),
        Nurbs(NurbsCurve<S, 3>),
    }
}

erased! {
    /// Any curve in a 2D space `S`; typically a pcurve, with `S` a surface's
    /// parameter space such as `Uv<AnySurface<World>>`.
    AnyCurve2<S: Space<2>>, T<_>, 1, {
        Line(Line<S, 2>),
        Nurbs(NurbsCurve<S, 2>),
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::Isometry3;
    use nalgebra::Vector2;
    use nalgebra::Vector3;

    use super::*;
    use crate::frame::Frame;
    use crate::frame::FrameTree;
    use crate::geom::Placement;
    use crate::manifold::Curve;
    use crate::manifold::Surface;
    use crate::mapping::Compose;
    use crate::space::Vector;
    use crate::space::World;

    fn placement() -> Placement<World> {
        Placement::from_axes(
            Point::new(Vector3::new(1.0, -1.0, 0.5)),
            Vector::new(Vector3::new(0.2, 0.1, 1.0)),
            Vector::new(Vector3::x()),
        )
    }

    fn surfaces() -> Vec<AnySurface<World>> {
        vec![
            Plane::new(placement()).into(),
            Sphere::new(Point::new(Vector3::new(0.0, 1.0, 0.0)), 2.0).into(),
            Cylinder::new(placement(), 1.5).into(),
            Cone::new(placement(), 1.0, 0.4).into(),
            Torus::new(placement(), 3.0, 1.0).into(),
        ]
    }

    /// Agrees with `typed` at `uv`: same point, jacobian, domain and projection.
    fn assert_agrees<G>(any: &AnySurface<World>, typed: &G, uv: [f64; 2])
    where
        G: Embedding<2, 3, To = World>,
        G::From: Space<2, Tag = ()>,
    {
        let uv = Vector2::from(uv);
        assert_eq!(any.apply(Point::new(uv)), typed.apply(Point::new(uv)));
        assert_eq!(any.jacobian(Point::new(uv)), typed.jacobian(Point::new(uv)));
        assert_eq!(any.domain(), typed.domain());
        let p = Point::new(typed.apply(Point::new(uv)).coords + Vector3::new(0.1, -0.2, 0.3));
        assert_eq!(
            any.project(p, None).map(|q| q.coords),
            typed.project(p, None).map(|q| q.coords)
        );
    }

    #[test]
    fn variants_delegate_to_their_geometry() {
        for (any, uv) in
            surfaces()
                .iter()
                .zip([[0.3, -0.2], [0.4, 0.2], [1.0, 0.5], [-2.0, 0.7], [0.1, 2.0]])
        {
            match any {
                AnySurface::Plane(g) => assert_agrees(any, g, uv),
                AnySurface::Sphere(g) => assert_agrees(any, g, uv),
                AnySurface::Cylinder(g) => assert_agrees(any, g, uv),
                AnySurface::Cone(g) => assert_agrees(any, g, uv),
                AnySurface::Torus(g) => assert_agrees(any, g, uv),
                other => panic!("unexpected {other:?}"),
            }
        }
    }

    #[test]
    fn custom_variant_holds_outside_geometry() {
        // Stands in for a type from another crate: any `Embedding` will do.
        let typed = Cylinder::new(placement(), 0.5);
        let any = AnySurface::Custom(Arc::new(typed.clone()));
        assert_agrees(&any, &typed, [0.7, 1.5]);
        // The blanket `Surface` impl applies to the erased type too.
        let uv = Point::new(Vector2::new(0.7, 1.5));
        assert!(
            (any.normal(uv).coords - typed.normal(Point::new(uv.coords)).coords).norm() < 1e-15
        );
    }

    #[test]
    fn erased_pcurve_composes_with_erased_surface() {
        let surface: AnySurface<World> = Sphere::new(Point::new(Vector3::zeros()), 2.0).into();
        let pcurve: AnyCurve2<Uv<AnySurface<World>>> = Line::new(
            Point::new(Vector2::new(0.0, -0.5)),
            Vector::new(Vector2::new(1.0, 1.0)),
        )
        .into();
        let edge = Compose::<_, _, 2>(pcurve, surface);
        for t in [0.0, 0.25, 1.0] {
            let p = edge.apply(Point::new(SVector::from([t])));
            assert!((p.coords.norm() - 2.0).abs() < 1e-12);
        }
    }

    #[test]
    fn erased_curves_in_a_frame() {
        let mut tree = FrameTree::new();
        let frame = tree.add(tree.root(), Isometry3::identity());
        let circle = Circle::new(Placement::at(frame.point(Vector3::zeros())), 1.0);
        let curves: Vec<AnyCurve<Frame>> = vec![
            circle.into(),
            Line::new(frame.point(Vector3::zeros()), frame.vector(Vector3::x())).into(),
        ];
        for curve in &curves {
            let t = Point::new(SVector::from([0.5]));
            assert_eq!(curve.apply(t).tag(), frame);
            assert_eq!(curve.tangent(t).tag(), frame);
        }
    }
}
