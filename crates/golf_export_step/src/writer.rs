//! ISO 10303-21 (Part 21) text: entity instances and value formatting.

use core::fmt;
use core::fmt::Write as _;

use golf_geom::Placement;
use golf_manifold::Space;
use nalgebra::Vector3;

use crate::error::StepError;

/// A reference to a written entity instance, `#n`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Ref(usize);

impl fmt::Display for Ref {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "#{}", self.0)
    }
}

/// Collects the entity instances of a STEP file's data section.
#[derive(Debug, Default)]
pub struct Writer {
    entities: Vec<String>,
}

impl Writer {
    /// Adds an entity instance, e.g. `"CARTESIAN_POINT('',(0.,0.,0.))"`.
    pub fn add(&mut self, instance: impl Into<String>) -> Ref {
        self.entities.push(instance.into());
        Ref(self.entities.len())
    }

    pub fn point(&mut self, p: Vector3<f64>) -> Result<Ref, StepError> {
        Ok(self.add(format!("CARTESIAN_POINT('',{})", reals(p.iter().copied())?)))
    }

    /// A `DIRECTION`, normalised.
    pub fn direction(&mut self, d: Vector3<f64>) -> Result<Ref, StepError> {
        let unit = d
            .try_normalize(f64::MIN_POSITIVE)
            .ok_or(StepError::ZeroDirection)?;
        Ok(self.add(format!("DIRECTION('',{})", reals(unit.iter().copied())?)))
    }

    /// An `AXIS2_PLACEMENT_3D` at the placement's origin, with its local z as
    /// the axis and local x as the reference direction.
    pub fn placement<S: Space<3>>(&mut self, placement: &Placement<S>) -> Result<Ref, StepError> {
        let location = self.point(placement.origin.coords)?;
        let axes = placement.rotation_matrix();
        let axis = self.direction(axes.column(2).into())?;
        let reference = self.direction(axes.column(0).into())?;
        Ok(self.add(format!(
            "AXIS2_PLACEMENT_3D('',{location},{axis},{reference})"
        )))
    }

    /// An `AXIS2_PLACEMENT_3D` for the frame `motion` moves the identity frame
    /// to: at its translation, axis its rotated z, reference its rotated x.
    pub fn motion(&mut self, motion: &nalgebra::Isometry3<f64>) -> Result<Ref, StepError> {
        let location = self.point(motion.translation.vector)?;
        let axis = self.direction(motion.rotation * Vector3::z())?;
        let reference = self.direction(motion.rotation * Vector3::x())?;
        Ok(self.add(format!(
            "AXIS2_PLACEMENT_3D('',{location},{axis},{reference})"
        )))
    }

    /// The whole file.
    pub(crate) fn finish(self, name: &str, time_stamp: &str) -> String {
        let mut out = String::new();
        out.push_str("ISO-10303-21;\nHEADER;\n");
        out.push_str("FILE_DESCRIPTION(('golf'),'2;1');\n");
        let _ = writeln!(
            out,
            "FILE_NAME({},{},(''),(''),'golf','golf','');",
            string(name),
            string(time_stamp)
        );
        out.push_str("FILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\n");
        out.push_str("ENDSEC;\nDATA;\n");
        for (i, instance) in self.entities.iter().enumerate() {
            let _ = writeln!(out, "#{}={instance};", i + 1);
        }
        out.push_str("ENDSEC;\nEND-ISO-10303-21;\n");
        out
    }
}

/// A `REAL`: always with a decimal point, exponent marked `E`.
pub fn real(v: f64) -> Result<String, StepError> {
    if !v.is_finite() {
        return Err(StepError::NonFinite(v));
    }
    let debug = format!("{v:?}");
    let (mantissa, exponent) = match debug.split_once('e') {
        Some((m, e)) => (m, Some(e)),
        None => (debug.as_str(), None),
    };
    let mut out = mantissa.to_string();
    if !out.contains('.') {
        out.push('.');
    }
    if let Some(exponent) = exponent {
        out.push('E');
        out.push_str(exponent);
    }
    Ok(out)
}

/// A list of `REAL`s.
pub fn reals(values: impl IntoIterator<Item = f64>) -> Result<String, StepError> {
    let values: Vec<String> = values.into_iter().map(real).collect::<Result<_, _>>()?;
    Ok(format!("({})", values.join(",")))
}

/// A list of references.
pub fn refs(refs: impl IntoIterator<Item = Ref>) -> String {
    let refs: Vec<String> = refs.into_iter().map(|r| r.to_string()).collect();
    format!("({})", refs.join(","))
}

/// A `LOGICAL` or `BOOLEAN`.
pub fn logical(v: bool) -> &'static str {
    if v { ".T." } else { ".F." }
}

/// A `STRING`, quoted, with `'` and `\` escaped and anything outside printable
/// ASCII encoded as `\X2\…\X0\` (or `\X4\` beyond the basic plane).
pub fn string(s: &str) -> String {
    let mut out = String::from("'");
    let mut run: Option<&str> = None;
    for c in s.chars() {
        let kind = match c {
            ' '..='~' => None,
            '\u{0}'..='\u{FFFF}' => Some("\\X2\\"),
            _ => Some("\\X4\\"),
        };
        if run != kind {
            if run.is_some() {
                out.push_str("\\X0\\");
            }
            if let Some(kind) = kind {
                out.push_str(kind);
            }
            run = kind;
        }
        match (kind, c) {
            (None, '\'') => out.push_str("''"),
            (None, '\\') => out.push_str("\\\\"),
            (None, _) => out.push(c),
            (Some(_), _) if (c as u32) <= 0xFFFF => {
                let _ = write!(out, "{:04X}", c as u32);
            }
            (Some(_), _) => {
                let _ = write!(out, "{:08X}", c as u32);
            }
        }
    }
    if run.is_some() {
        out.push_str("\\X0\\");
    }
    out.push('\'');
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reals_are_part_21() {
        assert_eq!(real(1.0).unwrap(), "1.0");
        assert_eq!(real(-0.25).unwrap(), "-0.25");
        assert_eq!(real(1e-7).unwrap(), "1.E-7");
        assert_eq!(real(2.5e20).unwrap(), "2.5E20");
        assert!(real(f64::NAN).is_err());
    }

    #[test]
    fn strings_are_escaped() {
        assert_eq!(string("it's"), "'it''s'");
        assert_eq!(string("a\\b"), "'a\\\\b'");
        assert_eq!(string("café"), "'caf\\X2\\00E9\\X0\\'");
    }

    #[test]
    fn instances_are_numbered_from_one() {
        let mut w = Writer::default();
        let a = w.point(Vector3::new(1.0, 2.0, 3.0)).unwrap();
        let b = w.direction(Vector3::new(0.0, 0.0, 2.0)).unwrap();
        assert_eq!((a.to_string(), b.to_string()), ("#1".into(), "#2".into()));
        let text = w.finish("part", "");
        assert!(text.contains("#1=CARTESIAN_POINT('',(1.0,2.0,3.0));"));
        assert!(text.contains("#2=DIRECTION('',(0.0,0.0,1.0));"));
        assert!(text.starts_with("ISO-10303-21;") && text.ends_with("END-ISO-10303-21;\n"));
    }
}
