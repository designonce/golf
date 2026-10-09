use core::fmt;

/// Index of a part definition in an [`Assembly`](crate::Assembly).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct PartId(pub(crate) u32);

impl PartId {
    pub fn index(self) -> usize {
        self.0 as usize
    }

    pub(crate) fn new(index: usize) -> Self {
        Self(u32::try_from(index).expect("too many parts for a u32 index"))
    }
}

impl fmt::Display for PartId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "p{}", self.0)
    }
}
