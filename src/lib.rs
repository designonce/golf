#![cfg_attr(docsrs, feature(doc_cfg))]

//! An open source CAD kernel.
//!
//! Each part of golf is its own crate, re-exported here as a module (`golf_brep`
//! as [`brep`], `golf_model` as [`model`], and so on), with the items most
//! modelling code needs gathered in the [`prelude`]. [`nalgebra`], which golf
//! is written in, is re-exported too.
//!
//!
//! ```
//! use golf::prelude::*;
//! use golf::nalgebra::Vector2;
//! use golf::nalgebra::Vector3;
//!
//! let sketch = Placement::<World>::at(Point::new(Vector3::zeros()));
//! let outline = Profile::rounded_rectangle(Vector2::zeros(), Vector2::new(40.0, 20.0), 4.0)?;
//! let plate: Body<World> = extrude(&sketch, &outline.into(), 5.0)?;
//! let facets = mesh(&plate, &Tolerance::new(0.05, 0.3))?;
//! assert!(facets.is_watertight());
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # Features
//!
//! - `assembly` (default): assemblies of parts placed in frames.
//! - `export_step` (default): writing bodies and assemblies as STEP files.
//! - `import_step` (default): reading STEP files into bodies and assemblies,
//!   healing them (with the Part 21 reader as [`step`] and healing as [`heal`]).
//! - `mesh` (default): triangle meshes of bodies.
//! - `model` (default): sketch-and-extrude style modelling operations.
//! - `view`: a Bevy viewer for meshes. Pulls in Bevy, so off by default.

pub use golf_internal::*;
