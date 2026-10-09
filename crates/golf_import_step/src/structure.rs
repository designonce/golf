//! Products and their shapes into an assembly.

use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::collections::HashMap;

use golf_step::Id;
use golf_step::StepData;
use nalgebra::Isometry3;
use nalgebra::Translation3;

use crate::geometry::Geometry;
use crate::units::representation_units;

/// Items that hold a body's shape.
pub(crate) const SOLIDS: [&str; 4] = [
    "MANIFOLD_SOLID_BREP",
    "BREP_WITH_VOIDS",
    "SHELL_BASED_SURFACE_MODEL",
    "FACETED_BREP",
];

/// A product definition: its name, its own shapes, and its components.
#[derive(Debug, Default)]
pub(crate) struct Product {
    pub name: String,
    /// Solid items, each with the representation it's in.
    pub solids: Vec<(Id, Id)>,
    /// Components: their product definition, and its placement in this one.
    pub children: Vec<(Id, Isometry3<f64>)>,
}

/// The product structure: every product definition by id, and the roots
/// (those no other uses).
pub(crate) struct Structure {
    pub products: BTreeMap<Id, Product>,
    pub roots: Vec<Id>,
}

pub(crate) fn structure(data: &StepData) -> Structure {
    let mut products: BTreeMap<Id, Product> = BTreeMap::new();
    for pd in data.all("PRODUCT_DEFINITION") {
        products.insert(
            pd.id,
            Product {
                name: product_name(data, pd.id).unwrap_or_default(),
                ..Product::default()
            },
        );
    }

    // Each product's shape representations, and those joined to them without
    // a transformation (the same coordinates).
    let mut shape_of: HashMap<Id, Vec<Id>> = HashMap::new();
    let mut rep_owner: HashMap<Id, Id> = HashMap::new();
    for sdr in data.all("SHAPE_DEFINITION_REPRESENTATION") {
        let Some(record) = sdr
            .record("SHAPE_DEFINITION_REPRESENTATION")
            .or_else(|| sdr.record("PROPERTY_DEFINITION_REPRESENTATION"))
        else {
            continue;
        };
        let (Ok(pds), Ok(rep)) = (record.reference(0), record.reference(1)) else {
            continue;
        };
        let Some(pd) = data
            .get(pds)
            .and_then(|e| e.records.iter().find_map(|r| r.reference(2).ok()))
        else {
            continue;
        };
        if products.contains_key(&pd) {
            shape_of.entry(pd).or_default().push(rep);
            rep_owner.insert(rep, pd);
        }
    }
    let mut joined: HashMap<Id, Vec<Id>> = HashMap::new();
    let mut relationships = data.all("REPRESENTATION_RELATIONSHIP");
    relationships.extend(data.all("SHAPE_REPRESENTATION_RELATIONSHIP"));
    relationships.sort_by_key(|e| e.id);
    relationships.dedup_by_key(|e| e.id);
    for srr in relationships {
        if srr.is("REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION") {
            continue;
        }
        let Some(record) = srr
            .record("REPRESENTATION_RELATIONSHIP")
            .or_else(|| srr.record("SHAPE_REPRESENTATION_RELATIONSHIP"))
        else {
            continue;
        };
        if let (Ok(a), Ok(b)) = (record.reference(2), record.reference(3)) {
            joined.entry(a).or_default().push(b);
            joined.entry(b).or_default().push(a);
        }
    }
    for (pd, reps) in &shape_of {
        let mut seen: BTreeSet<Id> = BTreeSet::new();
        let mut stack = reps.clone();
        while let Some(rep) = stack.pop() {
            if !seen.insert(rep) {
                continue;
            }
            stack.extend(joined.get(&rep).into_iter().flatten().copied());
            for item in items(data, rep) {
                if data
                    .get(item)
                    .is_some_and(|e| SOLIDS.iter().any(|s| e.is(s)))
                {
                    products
                        .get_mut(pd)
                        .expect("a product")
                        .solids
                        .push((item, rep));
                }
            }
        }
    }

    // Components, placed through their context-dependent shape
    // representations.
    let mut placements: HashMap<Id, Isometry3<f64>> = HashMap::new();
    for cdsr in data.all("CONTEXT_DEPENDENT_SHAPE_REPRESENTATION") {
        let Some(record) = cdsr.record("CONTEXT_DEPENDENT_SHAPE_REPRESENTATION") else {
            continue;
        };
        let (Ok(relation), Ok(pds)) = (record.reference(0), record.reference(1)) else {
            continue;
        };
        let Some(usage) = data
            .get(pds)
            .and_then(|e| e.records.iter().find_map(|r| r.reference(2).ok()))
        else {
            continue;
        };
        if let Some(placement) = relation_placement(data, relation, usage) {
            placements.insert(usage, placement);
        }
    }
    let mut used = BTreeSet::new();
    for name in [
        "NEXT_ASSEMBLY_USAGE_OCCURRENCE",
        "QUANTIFIED_ASSEMBLY_COMPONENT_USAGE",
    ] {
        for nauo in data.all(name) {
            let Some(record) = nauo
                .record(name)
                .or_else(|| nauo.record("PRODUCT_DEFINITION_USAGE"))
            else {
                continue;
            };
            let (Ok(parent), Ok(child)) = (record.reference(3), record.reference(4)) else {
                continue;
            };
            if !(products.contains_key(&parent) && products.contains_key(&child)) {
                continue;
            }
            let placement = placements
                .get(&nauo.id)
                .copied()
                .unwrap_or_else(Isometry3::identity);
            products
                .get_mut(&parent)
                .expect("a product")
                .children
                .push((child, placement));
            used.insert(child);
        }
    }
    let roots = products
        .keys()
        .copied()
        .filter(|id| !used.contains(id))
        .collect();
    Structure { products, roots }
}

/// A representation's items.
pub(crate) fn items(data: &StepData, representation: Id) -> Vec<Id> {
    data.get(representation)
        .and_then(|e| e.records.iter().find_map(|r| r.references(1).ok()))
        .unwrap_or_default()
}

/// A product definition's product's name.
fn product_name(data: &StepData, pd: Id) -> Option<String> {
    let formation = data
        .get(pd)?
        .records
        .iter()
        .find_map(|r| r.reference(2).ok())?;
    let product = data
        .get(formation)?
        .records
        .iter()
        .find_map(|r| r.reference(2).ok())?;
    let record = data.get(product)?.record("PRODUCT")?;
    let name = record
        .string(1)
        .ok()
        .filter(|s| !s.is_empty())
        .or_else(|| record.string(0).ok())?;
    Some(name.to_string())
}

/// The placement of a component in its parent, from the representation
/// relationship placing its shape: its origin item moved onto the parent's.
fn relation_placement(data: &StepData, relation: Id, _usage: Id) -> Option<Isometry3<f64>> {
    let entity = data.get(relation)?;
    let rr = entity.record("REPRESENTATION_RELATIONSHIP")?;
    let (rep_1, rep_2) = (rr.reference(2).ok()?, rr.reference(3).ok()?);
    let transformation = entity
        .record("REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION")?
        .reference(0)
        .ok()?;
    let t = data
        .get(transformation)?
        .record("ITEM_DEFINED_TRANSFORMATION")?;
    let (item_1, item_2) = (t.reference(2).ok()?, t.reference(3).ok()?);
    // Each placement in its own representation's units; rep_1 is the child's.
    let frame = |item: Id, rep: Id| {
        let geometry = Geometry {
            data,
            units: representation_units(data, rep),
        };
        geometry
            .placement(item)
            .ok()
            .map(|p| Isometry3::from_parts(Translation3::from(p.origin.coords), p.rotation))
    };
    let (child, parent) = (frame(item_1, rep_1)?, frame(item_2, rep_2)?);
    Some(parent * child.inverse())
}
