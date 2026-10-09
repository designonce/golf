//! `use golf::prelude::*;` for the items most modelling code needs.

#[doc(hidden)]
#[cfg(feature = "assembly")]
pub use crate::assembly::prelude::*;
#[doc(hidden)]
pub use crate::brep::prelude::*;
#[doc(hidden)]
pub use crate::color::prelude::*;
#[doc(hidden)]
#[cfg(feature = "export_step")]
pub use crate::export_step::prelude::*;
#[doc(hidden)]
pub use crate::frame::prelude::*;
#[doc(hidden)]
pub use crate::geom::prelude::*;
#[doc(hidden)]
pub use crate::manifold::prelude::*;
#[doc(hidden)]
#[cfg(feature = "mesh")]
pub use crate::mesh::prelude::*;
#[doc(hidden)]
#[cfg(feature = "model")]
pub use crate::model::prelude::*;
#[doc(hidden)]
pub use crate::sketch::prelude::*;
#[doc(hidden)]
#[cfg(feature = "view")]
pub use crate::view::prelude::*;
