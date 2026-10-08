//! Indices of entities in a body.

use std::fmt;

macro_rules! id {
    ($name:ident, $prefix:literal, $what:literal) => {
        #[doc = concat!("Index of a ", $what, " in a [`Body`](super::Body).")]
        #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
        pub struct $name(pub(super) u32);

        impl $name {
            pub fn index(self) -> usize {
                self.0 as usize
            }

            pub(super) fn new(index: usize) -> Self {
                Self(u32::try_from(index).expect("too many entities for a u32 index"))
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(f, concat!($prefix, "{}"), self.0)
            }
        }
    };
}

id!(VertexId, "v", "vertex");
id!(EdgeId, "e", "edge");
id!(FaceId, "f", "face");
id!(ShellId, "s", "shell");
