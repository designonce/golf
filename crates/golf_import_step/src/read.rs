use std::collections::HashMap;

use golf_assembly::Assembly;
use golf_assembly::PartId;
use golf_brep::Body;
use golf_frame::FrameId;
use golf_heal::HealOptions;
use golf_heal::HealReport;
use golf_heal::heal;
use golf_manifold::World;
use golf_step::Id;
use golf_step::StepData;
use nalgebra::Isometry3;

use crate::body::BodyBuilder;
use crate::error::ImportError;
use crate::error::ImportWarning;
use crate::error::Problem;
use crate::error::invalid;
use crate::error::problem;
use crate::error::unsupported;
use crate::geometry::Geometry;
use crate::structure::SOLIDS;
use crate::structure::Structure;
use crate::structure::items;
use crate::structure::structure;
use crate::style::colors;
use crate::units::representation_units;

/// How to import.
#[derive(Clone, Debug, PartialEq)]
pub struct ImportOptions {
    /// Healing for each body, or `None` to keep bodies as the file has them
    /// (some pcurves then missing).
    pub heal: Option<HealOptions>,
}

impl Default for ImportOptions {
    fn default() -> Self {
        Self {
            heal: Some(HealOptions::default()),
        }
    }
}

/// A file's contents: its product structure as an assembly of parts, each a
/// body, and everything that was skipped or couldn't be healed.
#[derive(Clone, Debug)]
pub struct StepImport {
    pub assembly: Assembly<World>,
    pub warnings: Vec<ImportWarning>,
    /// What healing did, over every body.
    pub heal: HealReport,
}

/// Reads a STEP file, healing its bodies.
pub fn read_step(bytes: &[u8]) -> Result<StepImport, ImportError> {
    read_step_with(bytes, &ImportOptions::default())
}

/// Reads a STEP file.
pub fn read_step_with(bytes: &[u8], options: &ImportOptions) -> Result<StepImport, ImportError> {
    let data = StepData::parse(bytes)?;
    let mut reader = Reader {
        data: &data,
        options,
        colors: colors(&data),
        parts: HashMap::new(),
        warnings: data
            .issues
            .iter()
            .map(|i| {
                problem(
                    i.id.unwrap_or(Id(0)),
                    Problem::Invalid(format!("line {}: {}", i.line, i.message)),
                )
            })
            .collect(),
        heal: HealReport::default(),
        rep_of: HashMap::new(),
    };
    let structure = structure(&data);
    let assembly = reader.assemble(&structure);
    Ok(StepImport {
        assembly,
        warnings: reader.warnings,
        heal: reader.heal,
    })
}

struct Reader<'a> {
    data: &'a StepData,
    options: &'a ImportOptions,
    colors: HashMap<Id, golf_color::Color>,
    /// Each solid read, as its part (or `None` if it couldn't be).
    parts: HashMap<Id, Option<PartId>>,
    warnings: Vec<ImportWarning>,
    heal: HealReport,
    /// The representation each solid is in, for files without products.
    rep_of: HashMap<Id, Id>,
}

impl Reader<'_> {
    fn assemble(&mut self, structure: &Structure) -> Assembly<World> {
        let named = |id: &Id| structure.products[id].name.clone();
        let mut assembly = match structure.roots.as_slice() {
            [root] => Assembly::new(named(root)),
            _ => Assembly::new("step"),
        };
        let root = assembly.root();
        match structure.roots.as_slice() {
            [single] => self.fill(&mut assembly, structure, *single, root, 0),
            roots => {
                for &r in roots {
                    self.place(&mut assembly, structure, r, root, Isometry3::identity(), 0);
                }
            }
        }
        if structure.products.is_empty() {
            self.orphans(&mut assembly);
        }
        assembly
    }

    /// Places product `pd` in `parent`: a lone solid as an instance of its
    /// part, anything more as a group.
    fn place(
        &mut self,
        assembly: &mut Assembly<World>,
        structure: &Structure,
        pd: Id,
        parent: FrameId,
        placement: Isometry3<f64>,
        depth: usize,
    ) {
        let product = &structure.products[&pd];
        if let ([(solid, rep)], []) = (product.solids.as_slice(), product.children.as_slice()) {
            if let Some(part) = self.part(assembly, *solid, *rep, &product.name) {
                let _ = assembly.add_instance(parent, part, placement);
            }
            return;
        }
        match assembly.add_group(parent, product.name.clone(), placement) {
            Ok(group) => self.fill(assembly, structure, pd, group, depth),
            Err(e) => self.warnings.push(invalid(pd, e.to_string())),
        }
    }

    /// Puts product `pd`'s solids and components into `group`.
    fn fill(
        &mut self,
        assembly: &mut Assembly<World>,
        structure: &Structure,
        pd: Id,
        group: FrameId,
        depth: usize,
    ) {
        if depth > 64 {
            self.warnings
                .push(invalid(pd, "its assembly nests too deep (a cycle?)"));
            return;
        }
        let product = &structure.products[&pd];
        for &(solid, rep) in &product.solids {
            if let Some(part) = self.part(assembly, solid, rep, &product.name) {
                let _ = assembly.add_instance(group, part, Isometry3::identity());
            }
        }
        for &(child, placement) in &product.children {
            self.place(assembly, structure, child, group, placement, depth + 1);
        }
    }

    /// Every solid in a file without product structure, at the root.
    fn orphans(&mut self, assembly: &mut Assembly<World>) {
        let data = self.data;
        for entity in data.entities() {
            for item in items(data, entity.id) {
                self.rep_of.entry(item).or_insert(entity.id);
            }
        }
        let root = assembly.root();
        let mut solids: Vec<Id> = SOLIDS
            .iter()
            .flat_map(|s| data.all(s))
            .map(|e| e.id)
            .collect();
        solids.sort();
        solids.dedup();
        for solid in solids {
            let rep = self.rep_of.get(&solid).copied().unwrap_or(Id(0));
            if let Some(part) = self.part(assembly, solid, rep, "") {
                let _ = assembly.add_instance(root, part, Isometry3::identity());
            }
        }
    }

    /// The part for a solid, read (and healed) on first use.
    fn part(
        &mut self,
        assembly: &mut Assembly<World>,
        solid: Id,
        rep: Id,
        name: &str,
    ) -> Option<PartId> {
        if let Some(&part) = self.parts.get(&solid) {
            return part;
        }
        let part = match self.body(solid, rep) {
            Ok(body) => {
                let name = match name.is_empty() {
                    false => name.to_string(),
                    true => self
                        .data
                        .get(solid)
                        .and_then(|e| e.records[0].string(0).ok().map(str::to_string))
                        .filter(|n| !n.is_empty())
                        .unwrap_or_else(|| "body".to_string()),
                };
                Some(assembly.add_part(name, body))
            }
            Err(warning) => {
                self.warnings.push(warning);
                None
            }
        };
        self.parts.insert(solid, part);
        part
    }

    fn body(&mut self, solid: Id, rep: Id) -> Result<Body<World>, ImportWarning> {
        let data = self.data;
        let geometry = Geometry {
            data,
            units: representation_units(data, rep),
        };
        let entity = data.entity(solid).map_err(|e| problem(solid, e))?;
        let record = &entity.records[0];
        let shells = match record.name.as_str() {
            "MANIFOLD_SOLID_BREP" | "FACETED_BREP" => {
                vec![record.reference(1).map_err(|e| problem(solid, e))?]
            }
            "BREP_WITH_VOIDS" => {
                let mut shells = vec![record.reference(1).map_err(|e| problem(solid, e))?];
                shells.extend(record.references(2).map_err(|e| problem(solid, e))?);
                shells
            }
            "SHELL_BASED_SURFACE_MODEL" => record.references(1).map_err(|e| problem(solid, e))?,
            other => return Err(unsupported(solid, format!("the solid {other}"))),
        };
        let mut builder = BodyBuilder::new(geometry, &self.colors);
        for &shell in &shells {
            if let Err(warning) = builder.add_shell(shell) {
                self.warnings.push(warning);
            }
        }
        let mut color_of = vec![solid];
        color_of.extend(&shells);
        let (mut body, warnings) = builder.finish(&color_of);
        self.warnings.extend(warnings);
        if let Some(options) = &self.options.heal {
            let report = heal(&mut body, options);
            for issue in &report.issues {
                self.warnings
                    .push(problem(solid, Problem::Heal(issue.to_string())));
            }
            self.heal.kept += report.kept;
            self.heal.moved += report.moved;
            self.heal.computed += report.computed;
            self.heal.seams += report.seams;
            self.heal.issues.extend(report.issues);
        }
        Ok(body)
    }
}
