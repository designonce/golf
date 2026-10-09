//! `golf_brep` bodies as STEP sources.

use golf_brep::Body;
use golf_brep::EdgeId;
use golf_brep::FaceId;
use golf_brep::VertexId;
use golf_color::Color;
use golf_geom::AnyCurve;
use golf_geom::AnySurface;
use golf_manifold::Point;
use golf_manifold::Space;

use crate::source::StepEdge;
use crate::source::StepFace;
use crate::source::StepSource;

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

    fn color(&self) -> Option<Color> {
        Body::color(self)
    }

    fn face_color(&self, face: FaceId) -> Option<Color> {
        self.face(face).color
    }
}

#[cfg(test)]
mod tests {
    use golf_geom::Placement;
    use golf_manifold::World;
    use golf_model::primitives::cuboid;
    use golf_model::primitives::cylinder;
    use nalgebra::Vector3;

    use super::*;
    use crate::StepFile;
    use crate::StepOptions;

    fn count(text: &str, entity: &str) -> usize {
        text.matches(&format!("={entity}(")).count()
    }

    #[test]
    fn cuboid_and_cylinder_as_solids() {
        let origin = Placement::<World>::at(Point::new(Vector3::zeros()));
        let mut file = StepFile::new(StepOptions::default()).unwrap();
        file.add_body(
            "box",
            &cuboid(&origin, Vector3::new(1.0, 2.0, 3.0)).unwrap(),
        )
        .unwrap();
        file.add_body("cylinder", &cylinder(&origin, 1.0, 2.0).unwrap())
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
