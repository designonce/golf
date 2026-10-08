//! Coordinate spaces and the points/vectors that live in them.

use core::marker::PhantomData;
use core::ops::Add;
use core::ops::Neg;
use core::ops::Sub;

use nalgebra::SVector;

/// A coordinate space of dimension `N`. Implementors are markers.
///
/// The dimension is a trait parameter rather than an associated const because
/// stable Rust can't yet use `S::DIM` as an array length.
///
/// Markers are unit types, so they're required to be `Copy + Debug`; that lets
/// geometry generic over a space derive `Clone` and `Debug`.
pub trait Space<const N: usize>: 'static + Copy + core::fmt::Debug {
    /// Runtime identity of a particular instance of this space, carried by every
    /// point and vector in it. `()` for spaces fully known at compile time.
    type Tag: Copy + PartialEq + core::fmt::Debug;
}

/// Euclidean world space, R3.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct World;
impl Space<3> for World {
    type Tag = ();
}

/// The 2D parameter space of surface type `M`.
pub struct Uv<M>(PhantomData<fn() -> M>);

/// The 1D parameter space of curve type `M`.
pub struct T<M>(PhantomData<fn() -> M>);

// By hand, so they don't require `M: Clone` / `M: Debug`.
macro_rules! impl_marker {
    ($t:ident) => {
        impl<M> Clone for $t<M> {
            fn clone(&self) -> Self {
                *self
            }
        }
        impl<M> Copy for $t<M> {}
        impl<M> core::fmt::Debug for $t<M> {
            fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
                write!(f, "{}<{}>", stringify!($t), core::any::type_name::<M>())
            }
        }
    };
}
impl_marker!(Uv);
impl_marker!(T);

/// Geometry with a 2D parameter space, [`Uv<Self>`](Uv).
///
/// Implementing it makes `Uv<Self>` an untagged [`Space<2>`]; it exists so
/// geometry defined in other crates can have parameter spaces without
/// implementing `Space` for `Uv`, which the orphan rule forbids.
pub trait HasUv: 'static {}

impl<M: HasUv> Space<2> for Uv<M> {
    type Tag = ();
}

/// Geometry with a 1D parameter space, [`T<Self>`](T). See [`HasUv`].
pub trait HasT: 'static {}

impl<M: HasT> Space<1> for T<M> {
    type Tag = ();
}

/// A position in space `S`.
pub struct Point<S: Space<N>, const N: usize> {
    pub coords: SVector<f64, N>,
    tag: S::Tag,
    _space: PhantomData<fn() -> S>,
}

/// A displacement (tangent vector) in space `S`.
pub struct Vector<S: Space<N>, const N: usize> {
    pub coords: SVector<f64, N>,
    tag: S::Tag,
    _space: PhantomData<fn() -> S>,
}

// Implemented by hand so space markers don't need `Clone`/`Copy`/`Debug`.
macro_rules! impl_common {
    ($t:ident) => {
        impl<S: Space<N, Tag = ()>, const N: usize> $t<S, N> {
            pub fn new(coords: SVector<f64, N>) -> Self {
                Self::with_tag(coords, ())
            }
        }
        impl<S: Space<N>, const N: usize> $t<S, N> {
            pub fn with_tag(coords: SVector<f64, N>, tag: S::Tag) -> Self {
                Self {
                    coords,
                    tag,
                    _space: PhantomData,
                }
            }
            pub fn tag(&self) -> S::Tag {
                self.tag
            }
            pub fn scale(self, k: f64) -> Self {
                Self::with_tag(self.coords * k, self.tag)
            }
        }
        impl<S: Space<N>, const N: usize> Clone for $t<S, N> {
            fn clone(&self) -> Self {
                *self
            }
        }
        impl<S: Space<N>, const N: usize> Copy for $t<S, N> {}
        impl<S: Space<N>, const N: usize> PartialEq for $t<S, N> {
            fn eq(&self, other: &Self) -> bool {
                self.tag == other.tag && self.coords == other.coords
            }
        }
        impl<S: Space<N>, const N: usize> core::fmt::Debug for $t<S, N> {
            fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
                f.debug_tuple(stringify!($t))
                    .field(&self.tag)
                    .field(&self.coords.as_slice())
                    .finish()
            }
        }
    };
}
impl_common!(Point);
impl_common!(Vector);

impl<S: Space<N>, const N: usize> Sub for Point<S, N> {
    type Output = Vector<S, N>;
    fn sub(self, rhs: Self) -> Vector<S, N> {
        same_space(self.tag, rhs.tag);
        Vector::with_tag(self.coords - rhs.coords, self.tag)
    }
}

impl<S: Space<N>, const N: usize> Add<Vector<S, N>> for Point<S, N> {
    type Output = Self;
    fn add(self, rhs: Vector<S, N>) -> Self {
        same_space(self.tag, rhs.tag);
        Self::with_tag(self.coords + rhs.coords, self.tag)
    }
}

impl<S: Space<N>, const N: usize> Sub<Vector<S, N>> for Point<S, N> {
    type Output = Self;
    fn sub(self, rhs: Vector<S, N>) -> Self {
        same_space(self.tag, rhs.tag);
        Self::with_tag(self.coords - rhs.coords, self.tag)
    }
}

impl<S: Space<N>, const N: usize> Add for Vector<S, N> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        same_space(self.tag, rhs.tag);
        Self::with_tag(self.coords + rhs.coords, self.tag)
    }
}

impl<S: Space<N>, const N: usize> Sub for Vector<S, N> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        same_space(self.tag, rhs.tag);
        Self::with_tag(self.coords - rhs.coords, self.tag)
    }
}

impl<S: Space<N>, const N: usize> Neg for Vector<S, N> {
    type Output = Self;
    fn neg(self) -> Self {
        Self::with_tag(-self.coords, self.tag)
    }
}

#[track_caller]
fn same_space<T: PartialEq + core::fmt::Debug>(a: T, b: T) {
    assert_eq!(
        a, b,
        "operands are in different instances of the same space"
    );
}
