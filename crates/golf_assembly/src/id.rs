use core::fmt;

macro_rules! id {
    ($name:ident, $prefix:literal, $what:literal) => {
        #[doc = concat!("Index of a ", $what, " definition in an [`Assembly`](crate::Assembly).")]
        #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
        pub struct $name(pub(crate) u32);

        impl $name {
            pub fn index(self) -> usize {
                self.0 as usize
            }

            pub(crate) fn new(index: usize) -> Self {
                Self(u32::try_from(index).expect("too many definitions for a u32 index"))
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(f, concat!($prefix, "{}"), self.0)
            }
        }
    };
}

id!(PartId, "p", "part");
id!(GroupId, "g", "group");
