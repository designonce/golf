//! Colours from presentation styles.

use std::collections::HashMap;

use golf_color::Color;
use golf_step::Id;
use golf_step::StepData;

/// The colour each styled item (a solid, shell or face) is given.
pub(crate) fn colors(data: &StepData) -> HashMap<Id, Color> {
    let mut colors = HashMap::new();
    // Plain styles first, so overriding ones win.
    for name in ["STYLED_ITEM", "OVER_RIDING_STYLED_ITEM"] {
        for entity in data.all(name) {
            let Some(record) = entity.record(name) else {
                continue;
            };
            let (Ok(styles), Ok(item)) = (record.references(1), record.reference(2)) else {
                continue;
            };
            if let Some(color) = styles.into_iter().find_map(|s| assignment_color(data, s)) {
                colors.insert(item, color);
            }
        }
    }
    colors
}

/// The surface colour a `PRESENTATION_STYLE_ASSIGNMENT` gives, with any
/// transparency.
fn assignment_color(data: &StepData, id: Id) -> Option<Color> {
    let entity = data.get(id)?;
    let styles = entity.records.iter().find_map(|r| r.references(0).ok())?;
    styles.into_iter().find_map(|usage| {
        // SURFACE_STYLE_USAGE(.BOTH., SURFACE_SIDE_STYLE('', (styles...))).
        let usage = data.get(usage)?.record("SURFACE_STYLE_USAGE")?;
        let side = data
            .get(usage.reference(1).ok()?)?
            .record("SURFACE_SIDE_STYLE")?;
        let elements = side.references(1).ok()?;
        let mut color = None;
        let mut transparency = None;
        for element in elements {
            let element = data.get(element)?;
            if let Some(fill) = element.record("SURFACE_STYLE_FILL_AREA") {
                color = color.or_else(|| fill_color(data, fill.reference(0).ok()?));
            } else if let Some(rendering) =
                element.record("SURFACE_STYLE_RENDERING_WITH_PROPERTIES")
            {
                color = color.or_else(|| colour(data, rendering.reference(1).ok()?));
                for property in rendering.references(2).ok()? {
                    if let Some(t) = data
                        .get(property)
                        .and_then(|p| p.record("SURFACE_STYLE_TRANSPARENT"))
                    {
                        transparency = t.real(0).ok();
                    }
                }
            }
        }
        let mut color = color?;
        if let Some(t) = transparency.filter(|&t| t > 0.0) {
            color.a = Some(((1.0 - t).clamp(0.0, 1.0) * 255.0).round() as u8);
        }
        Some(color)
    })
}

/// `FILL_AREA_STYLE('', (FILL_AREA_STYLE_COLOUR('', colour)))`.
fn fill_color(data: &StepData, id: Id) -> Option<Color> {
    let fill = data.get(id)?.record("FILL_AREA_STYLE")?;
    fill.references(1).ok()?.into_iter().find_map(|f| {
        let fill_colour = data.get(f)?.record("FILL_AREA_STYLE_COLOUR")?;
        colour(data, fill_colour.reference(1).ok()?)
    })
}

/// A `COLOUR_RGB`, or a `DRAUGHTING_PRE_DEFINED_COLOUR` by name.
fn colour(data: &StepData, id: Id) -> Option<Color> {
    let entity = data.get(id)?;
    let unit = |x: f64| (x.clamp(0.0, 1.0) * 255.0).round() as u8;
    if let Some(rgb) = entity.record("COLOUR_RGB") {
        let (r, g, b) = (rgb.real(1).ok()?, rgb.real(2).ok()?, rgb.real(3).ok()?);
        return Some(Color::rgb(unit(r), unit(g), unit(b)));
    }
    let named = entity.record("DRAUGHTING_PRE_DEFINED_COLOUR")?;
    let [r, g, b] = match named.string(0).ok()? {
        "red" => [1.0, 0.0, 0.0],
        "green" => [0.0, 1.0, 0.0],
        "blue" => [0.0, 0.0, 1.0],
        "yellow" => [1.0, 1.0, 0.0],
        "magenta" => [1.0, 0.0, 1.0],
        "cyan" => [0.0, 1.0, 1.0],
        "black" => [0.0, 0.0, 0.0],
        "white" => [1.0, 1.0, 1.0],
        _ => return None,
    };
    Some(Color::rgb(unit(r), unit(g), unit(b)))
}
