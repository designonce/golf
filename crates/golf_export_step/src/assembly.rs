//! `golf_assembly` assemblies as STEP assembly structures.

use std::collections::HashMap;

use golf_assembly::Assembly;
use golf_assembly::NodeKind;
use golf_assembly::PartId;
use golf_brep::Body;
use golf_color::Color;
use golf_frame::FrameId;
use golf_manifold::Space;

use crate::error::StepError;
use crate::file::Product;
use crate::file::StepFile;
use crate::source::StepSource;

impl StepFile {
    /// Adds `assembly` as a STEP assembly structure: each group as an assembly
    /// product placing its children, and each part as a product written once
    /// per colour it's shown in. Returns the root's product.
    pub fn add_assembly<S: Space<3>>(
        &mut self,
        assembly: &Assembly<S>,
    ) -> Result<Product, StepError>
    where
        Body<S>: StepSource,
    {
        self.add_assembly_group(assembly, assembly.root(), &mut HashMap::new())
    }

    fn add_assembly_group<S: Space<3>>(
        &mut self,
        assembly: &Assembly<S>,
        group: FrameId,
        parts: &mut HashMap<(PartId, Option<Color>), Product>,
    ) -> Result<Product, StepError>
    where
        Body<S>: StepSource,
    {
        let node = assembly.node(group);
        let mut children = Vec::with_capacity(node.children.len());
        for &child in &node.children {
            let product = match assembly.node(child).kind {
                NodeKind::Group => self.add_assembly_group(assembly, child, parts)?,
                NodeKind::Instance(part) => {
                    let color = assembly.color(child);
                    match parts.get(&(part, color)) {
                        Some(&product) => product,
                        None => {
                            let definition = assembly.part(part);
                            // A part shown in another colour than its own is
                            // written again in that colour.
                            let product = match color == definition.body.color() {
                                true => self.add_body(&definition.name, &definition.body)?,
                                false => {
                                    let mut body = definition.body.clone();
                                    body.set_color(color);
                                    self.add_body(&definition.name, &body)?
                                }
                            };
                            parts.insert((part, color), product);
                            product
                        }
                    }
                }
            };
            children.push((product, assembly.frames().placement(child)));
        }
        self.add_group(&node.name, &children)
    }
}

#[cfg(test)]
mod tests {
    use golf_geom::Placement;
    use golf_manifold::Point;
    use golf_manifold::World;
    use golf_model::primitives;
    use nalgebra::Isometry3;
    use nalgebra::Vector3;

    use super::*;
    use crate::StepOptions;

    #[test]
    fn parts_are_written_once_per_colour() {
        let origin = Placement::<World>::at(Point::new(Vector3::zeros()));
        let mut asm = Assembly::new("rack");
        let peg = asm.add_part("peg", primitives::cylinder(&origin, 1.0, 5.0).unwrap());
        let root = asm.root();
        let row = asm.add_group(root, "row", Isometry3::identity()).unwrap();
        for i in 0..3 {
            asm.add_instance(row, peg, Isometry3::translation(4.0 * i as f64, 0.0, 0.0))
                .unwrap();
        }
        let second = asm
            .copy_group(row, root, Isometry3::translation(0.0, 4.0, 0.0))
            .unwrap();
        asm.set_color(second, Some(Color::rgb(0, 0, 255))).unwrap();
        let mut file = StepFile::new(StepOptions::default()).unwrap();
        file.add_assembly(&asm).unwrap();
        let text = file.finish();
        let count = |entity: &str| text.matches(&format!("={entity}(")).count();
        // The peg plain and blue; two rows and the rack.
        assert_eq!(count("MANIFOLD_SOLID_BREP"), 2);
        assert_eq!(count("PRODUCT"), 2 + 3);
        assert_eq!(count("NEXT_ASSEMBLY_USAGE_OCCURRENCE"), 3 + 3 + 2);
        assert_eq!(count("STYLED_ITEM"), 1);
        assert!(text.contains("COLOUR_RGB('',0.0,0.0,1.0)"));
    }
}
