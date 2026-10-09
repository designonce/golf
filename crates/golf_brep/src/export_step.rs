//! Writing bodies with `golf_export_step`.

use golf_export_step::StepEdge;
use golf_export_step::StepFace;
use golf_export_step::StepSource;
use golf_geom::AnyCurve;
use golf_geom::AnySurface;
use golf_manifold::Point;
use golf_manifold::Space;

use crate::Body;
use crate::EdgeId;
use crate::FaceId;
use crate::VertexId;

impl<S: Space<3>> StepSource for Body<S> {
    type Space = S;
    type Curve = AnyCurve<S>;
    type Surface = AnySurface<S>;
    type Vertex = VertexId;
    type Edge = EdgeId;
    type Face = FaceId;

    fn vertices(&self) -> impl Iterator<Item = (VertexId, Point<S, 3>)> + '_ {
        Body::vertices(self).map(|(id, vertex)| (id, vertex.point))
    }

    fn edges(&self) -> impl Iterator<Item = (EdgeId, StepEdge<'_, Self>)> + '_ {
        Body::edges(self).map(|(id, edge)| {
            let data = StepEdge {
                curve: &edge.curve,
                start: edge.start,
                end: edge.end,
            };
            (id, data)
        })
    }

    fn faces(&self) -> impl Iterator<Item = (FaceId, StepFace<'_, Self>)> + '_ {
        Body::faces(self).map(|(id, face)| {
            let loops = face
                .loops
                .iter()
                .map(|l| l.coedges.iter().map(|c| (c.edge, c.reversed)).collect())
                .collect();
            let data = StepFace {
                surface: &face.surface,
                same_sense: face.same_sense,
                loops,
            };
            (id, data)
        })
    }

    fn shells(&self) -> impl Iterator<Item = Vec<FaceId>> + '_ {
        Body::shells(self).map(|(_, shell)| shell.faces.clone())
    }
}

#[cfg(test)]
mod tests {
    use golf_export_step::StepFile;
    use golf_export_step::StepOptions;
    use golf_geom::Placement;
    use golf_manifold::World;
    use nalgebra::Vector3;

    use super::*;
    use crate::cuboid;
    use crate::cylinder;

    fn count(text: &str, entity: &str) -> usize {
        text.matches(&format!("={entity}(")).count()
    }

    #[test]
    fn cuboid_and_cylinder_as_solids() {
        let origin = Placement::<World>::at(Point::new(Vector3::zeros()));
        let mut file = StepFile::new(StepOptions::default()).unwrap();
        file.add_body("box", &cuboid(origin, Vector3::new(1.0, 2.0, 3.0)))
            .unwrap();
        file.add_body("cylinder", &cylinder(origin, 1.0, 2.0))
            .unwrap();
        let text = file.finish();
        assert_eq!(count(&text, "MANIFOLD_SOLID_BREP"), 2);
        assert_eq!(count(&text, "ADVANCED_FACE"), 6 + 3);
        assert_eq!(count(&text, "EDGE_CURVE"), 12 + 3);
        assert_eq!(count(&text, "VERTEX_POINT"), 8 + 2);
        assert_eq!(count(&text, "ORIENTED_EDGE"), 2 * (12 + 3));
        assert_eq!(count(&text, "CYLINDRICAL_SURFACE"), 1);
        assert_eq!(count(&text, "PRODUCT"), 2);
    }
}
