//! Reading STEP physical files (ISO 10303-21, "Part 21") into generic
//! entities.
//!
//! A [`StepData`] holds the header and every entity instance by its id, each
//! an [`Entity`]: one [`Record`] (a keyword and its parameters), or several
//! for a complex instance such as a rational B-spline. Nothing about any
//! particular schema is built in, so entities this crate has never heard of
//! read like any other; meaning is given by whoever reads them, through the
//! typed accessors on [`Record`] and [`Value`], whose errors say which
//! instance and parameter were wrong.
//!
//! An instance that can't be parsed is skipped, with a [`ParseIssue`] saying
//! where, rather than failing the file.

mod data;
mod error;
mod parse;
mod string;
mod value;

pub use data::StepData;
pub use error::DecodeError;
pub use error::ParseError;
pub use error::ParseIssue;
pub use value::Entity;
pub use value::Id;
pub use value::Record;
pub use value::Value;
