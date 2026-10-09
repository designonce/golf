//! Representation contexts' units, as factors to millimetres and radians.

use golf_step::Entity;
use golf_step::Id;
use golf_step::StepData;

/// Factors taking a representation's lengths to millimetres and its angles
/// to radians.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Units {
    pub length: f64,
    pub angle: f64,
}

impl Default for Units {
    fn default() -> Self {
        Self {
            length: 1.0,
            angle: 1.0,
        }
    }
}

/// The units of representation `representation`'s context, by default where
/// it doesn't say.
pub(crate) fn representation_units(data: &StepData, representation: Id) -> Units {
    let context = data
        .get(representation)
        .and_then(|r| r.records.iter().find_map(|record| record.reference(2).ok()));
    context.map_or_else(Units::default, |c| context_units(data, c))
}

/// The units a `GLOBAL_UNIT_ASSIGNED_CONTEXT` assigns.
pub(crate) fn context_units(data: &StepData, context: Id) -> Units {
    let mut units = Units::default();
    let Some(assigned) = data
        .get(context)
        .and_then(|c| c.record("GLOBAL_UNIT_ASSIGNED_CONTEXT"))
        .and_then(|r| r.references(0).ok())
    else {
        return units;
    };
    for unit in assigned {
        let Some(entity) = data.get(unit) else {
            continue;
        };
        if entity.is("LENGTH_UNIT") {
            if let Some(factor) = factor(data, entity, 1e3, 0) {
                units.length = factor;
            }
        } else if entity.is("PLANE_ANGLE_UNIT") {
            if let Some(factor) = factor(data, entity, 1.0, 0) {
                units.angle = factor;
            }
        }
    }
    units
}

/// A unit's size in the base unit (mm or radians), where `si` is the size of
/// the unprefixed SI unit (a metre is 1000 mm).
fn factor(data: &StepData, unit: &Entity, si: f64, depth: usize) -> Option<f64> {
    if depth > 8 {
        return None;
    }
    if let Some(record) = unit.record("SI_UNIT") {
        let prefix = match record.param(0).ok()?.as_enum() {
            Some(prefix) => si_prefix(prefix)?,
            None => 1.0,
        };
        return Some(prefix * si);
    }
    if let Some(record) = unit.record("CONVERSION_BASED_UNIT") {
        // 'INCH', LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4), #mm).
        let measure = data.get(record.reference(1).ok()?)?;
        let with_unit = measure
            .records
            .iter()
            .find(|r| r.name.ends_with("MEASURE_WITH_UNIT"))?;
        let value = with_unit.real(0).ok()?;
        let base = data.get(with_unit.reference(1).ok()?)?;
        return Some(value * factor(data, base, si, depth + 1)?);
    }
    None
}

fn si_prefix(prefix: &str) -> Option<f64> {
    Some(match prefix {
        "EXA" => 1e18,
        "PETA" => 1e15,
        "TERA" => 1e12,
        "GIGA" => 1e9,
        "MEGA" => 1e6,
        "KILO" => 1e3,
        "HECTO" => 1e2,
        "DECA" => 1e1,
        "DECI" => 1e-1,
        "CENTI" => 1e-2,
        "MILLI" => 1e-3,
        "MICRO" => 1e-6,
        "NANO" => 1e-9,
        "PICO" => 1e-12,
        "FEMTO" => 1e-15,
        "ATTO" => 1e-18,
        _ => return None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn units(lines: &str) -> Units {
        let file = format!(
            "ISO-10303-21;HEADER;ENDSEC;DATA;{lines}#9=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNIT_ASSIGNED_CONTEXT((#1,#2))REPRESENTATION_CONTEXT('',''));ENDSEC;END-ISO-10303-21;"
        );
        context_units(&StepData::parse(file.as_bytes()).unwrap(), Id(9))
    }

    #[test]
    fn si_and_conversion_based_units() {
        let mm = units(
            "#1=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));#2=(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.));",
        );
        assert_eq!(
            mm,
            Units {
                length: 1.0,
                angle: 1.0
            }
        );
        let metres = units(
            "#1=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT($,.METRE.));#2=(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.));",
        );
        assert_eq!(metres.length, 1000.0);
        let inches = units(
            "#3=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));#4=LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4),#3);#1=(CONVERSION_BASED_UNIT('INCH',#4)LENGTH_UNIT()NAMED_UNIT(#5));#5=DIMENSIONAL_EXPONENTS(1.,0.,0.,0.,0.,0.,0.);\
             #6=(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.));#7=PLANE_ANGLE_MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE(0.0174532925),#6);#2=(CONVERSION_BASED_UNIT('DEGREE',#7)NAMED_UNIT(#5)PLANE_ANGLE_UNIT());",
        );
        assert!((inches.length - 25.4).abs() < 1e-12);
        assert!((inches.angle - 0.0174532925).abs() < 1e-12);
    }
}
