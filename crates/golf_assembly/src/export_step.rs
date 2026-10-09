//! Writing assemblies with `golf_export_step`.

use std::collections::HashMap;

use golf_export_step::Product;
use golf_export_step::StepError;
use golf_export_step::StepFile;
use golf_export_step::StepSource;
use golf_manifold::Space;

use crate::assembly::Assembly;
use crate::assembly::Item;
use crate::id::GroupId;
use crate::id::PartId;

impl<S: Space<3>> Assembly<S>
where
    golf_brep::Body<S>: StepSource,
{
    /// Writes the assembly into `file` as a STEP assembly structure: each part
    /// and group reachable from the root once, as a product, groups placing
    /// their children's products. Returns the root's product.
    pub fn write_step(&self, file: &mut StepFile) -> Result<Product, StepError> {
        let mut written = Written::default();
        self.write_group(file, self.root(), &mut written)
    }

    fn write_group(
        &self,
        file: &mut StepFile,
        group: GroupId,
        written: &mut Written,
    ) -> Result<Product, StepError> {
        if let Some(&product) = written.groups.get(&group) {
            return Ok(product);
        }
        let definition = self.group(group);
        let mut children = Vec::with_capacity(definition.children.len());
        for child in &definition.children {
            let product = match child.item {
                Item::Part(part) => self.write_part(file, part, written)?,
                Item::Group(inner) => self.write_group(file, inner, written)?,
            };
            children.push((product, child.placement));
        }
        let product = file.add_assembly(&definition.name, &children)?;
        written.groups.insert(group, product);
        Ok(product)
    }

    fn write_part(
        &self,
        file: &mut StepFile,
        part: PartId,
        written: &mut Written,
    ) -> Result<Product, StepError> {
        if let Some(&product) = written.parts.get(&part) {
            return Ok(product);
        }
        let definition = self.part(part);
        let product = file.add_body(&definition.name, &definition.body)?;
        written.parts.insert(part, product);
        Ok(product)
    }
}

/// The products already written, so shared definitions are written once.
#[derive(Default)]
struct Written {
    parts: HashMap<PartId, Product>,
    groups: HashMap<GroupId, Product>,
}

#[cfg(test)]
mod tests {
    use golf_export_step::StepOptions;
    use golf_geom::Placement;
    use golf_manifold::Point;
    use golf_manifold::World;
    use golf_model::primitives;
    use nalgebra::Isometry3;
    use nalgebra::Vector3;

    use super::*;

    #[test]
    fn shared_definitions_are_written_once() {
        let origin = Placement::<World>::at(Point::new(Vector3::zeros()));
        let mut asm = Assembly::new("rack");
        let peg = asm.add_part("peg", primitives::cylinder(&origin, 1.0, 5.0).unwrap());
        let row = asm.add_group("row");
        for i in 0..3 {
            asm.place(row, peg, Isometry3::translation(4.0 * i as f64, 0.0, 0.0))
                .unwrap();
        }
        for j in 0..2 {
            asm.place(
                asm.root(),
                row,
                Isometry3::translation(0.0, 4.0 * j as f64, 0.0),
            )
            .unwrap();
        }
        let mut file = StepFile::new(StepOptions::default()).unwrap();
        asm.write_step(&mut file).unwrap();
        let text = file.finish();
        let count = |entity: &str| text.matches(&format!("={entity}(")).count();
        // One solid (the peg), three products (peg, row, rack), five placements.
        assert_eq!(count("MANIFOLD_SOLID_BREP"), 1);
        assert_eq!(count("PRODUCT"), 3);
        assert_eq!(count("NEXT_ASSEMBLY_USAGE_OCCURRENCE"), 3 + 2);
        assert_eq!(count("CONTEXT_DEPENDENT_SHAPE_REPRESENTATION"), 5);
    }
}
