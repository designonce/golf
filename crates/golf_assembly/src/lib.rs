//! Assemblies: parts and sub-assemblies placed relative to one another.
//!
//! An [`Assembly`] holds *definitions*: [parts](Assembly::add_part) (a named
//! body in its own coordinates) and [groups](Assembly::add_group) (named
//! sub-assemblies). A group places parts and other groups with rigid motions,
//! and any definition may be placed many times, in one group or several, so a
//! part used in a hundred places is stored once. The assembly's root group
//! places everything else.
//!
//! Walking the placements from the root gives each [`Occurrence`]: one placed
//! copy of a part, with the path that reached it and its motion into the
//! assembly's coordinates.

mod assembly;
mod error;
#[cfg(feature = "export_step")]
mod export_step;
mod id;

pub use assembly::Assembly;
pub use assembly::Child;
pub use assembly::Group;
pub use assembly::Item;
pub use assembly::Occurrence;
pub use assembly::Part;
pub use error::AssemblyError;
pub use id::GroupId;
pub use id::PartId;
