//! A hierarchy of rigid coordinate frames.

use nalgebra::Isometry3;
use nalgebra::Matrix3;
use nalgebra::Vector3;

use crate::domain::Domain;
use crate::manifold::Embedding;
use crate::manifold::ProjectError;
use crate::mapping::Mapping;
use crate::space::Point;
use crate::space::Space;
use crate::space::Vector;

/// A node in a [`FrameTree`].
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct FrameId(u32);

impl FrameId {
    pub fn point(self, coords: Vector3<f64>) -> Point<Frame, 3> {
        Point::with_tag(coords, self)
    }

    pub fn vector(self, coords: Vector3<f64>) -> Vector<Frame, 3> {
        Vector::with_tag(coords, self)
    }
}

/// Any frame in a [`FrameTree`]; which one is carried at runtime by the
/// [`FrameId`] tag on each point and vector.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Frame;
impl Space<3> for Frame {
    type Tag = FrameId;
}

struct Node {
    parent: Option<FrameId>,
    /// Maps coordinates in this frame to coordinates in `parent`.
    to_parent: Isometry3<f64>,
    depth: u32,
}

/// A tree of frames whose edges are rigid transforms.
pub struct FrameTree {
    nodes: Vec<Node>,
}

impl Default for FrameTree {
    fn default() -> Self {
        Self::new()
    }
}

impl FrameTree {
    /// A tree containing only the root frame.
    pub fn new() -> Self {
        let root = Node {
            parent: None,
            to_parent: Isometry3::identity(),
            depth: 0,
        };
        Self { nodes: vec![root] }
    }

    pub fn root(&self) -> FrameId {
        FrameId(0)
    }

    /// Adds a child of `parent`, placed by `to_parent`.
    pub fn add(&mut self, parent: FrameId, to_parent: Isometry3<f64>) -> FrameId {
        let depth = self.node(parent).depth + 1;
        let id = FrameId(self.nodes.len() as u32);
        self.nodes.push(Node {
            parent: Some(parent),
            to_parent,
            depth,
        });
        id
    }

    /// Moves `frame` (and so everything below it) relative to its parent.
    pub fn set_placement(&mut self, frame: FrameId, to_parent: Isometry3<f64>) {
        self.nodes[frame.0 as usize].to_parent = to_parent;
    }

    /// The rigid transform taking coordinates in `from` to coordinates in `to`.
    ///
    /// The result is a snapshot: it does not follow later `set_placement` calls.
    pub fn transform(&self, from: FrameId, to: FrameId) -> Rigid {
        // Walk both frames up to their lowest common ancestor.
        let (mut a, mut b) = (from, to);
        let (mut up, mut down) = (Isometry3::identity(), Isometry3::identity());
        while a != b {
            if self.node(a).depth >= self.node(b).depth {
                up = self.node(a).to_parent * up;
                a = self.node(a).parent.unwrap();
            } else {
                down = self.node(b).to_parent * down;
                b = self.node(b).parent.unwrap();
            }
        }
        Rigid {
            iso: down.inverse() * up,
            from,
            to,
        }
    }

    fn node(&self, id: FrameId) -> &Node {
        &self.nodes[id.0 as usize]
    }
}

/// A rigid transform between two frames of a [`FrameTree`].
#[derive(Clone, Copy, Debug)]
pub struct Rigid {
    pub iso: Isometry3<f64>,
    pub from: FrameId,
    pub to: FrameId,
}

impl Rigid {
    pub fn inverse(&self) -> Rigid {
        Rigid {
            iso: self.iso.inverse(),
            from: self.to,
            to: self.from,
        }
    }

    /// `other` followed by `self`.
    pub fn after(&self, other: &Rigid) -> Rigid {
        assert_eq!(other.to, self.from, "transforms do not chain");
        Rigid {
            iso: self.iso * other.iso,
            from: other.from,
            to: self.to,
        }
    }
}

impl Mapping<3, 3> for Rigid {
    type From = Frame;
    type To = Frame;

    #[track_caller]
    fn apply(&self, p: Point<Frame, 3>) -> Point<Frame, 3> {
        assert_eq!(
            p.tag(),
            self.from,
            "point is not in this transform's source frame"
        );
        Point::with_tag(self.iso.transform_point(&p.coords.into()).coords, self.to)
    }

    fn jacobian(&self, _p: Point<Frame, 3>) -> Matrix3<f64> {
        self.iso.rotation.to_rotation_matrix().into_inner()
    }
}

impl Embedding<3, 3> for Rigid {
    fn domain(&self) -> Domain<3> {
        Domain::unbounded()
    }

    fn project(
        &self,
        p: Point<Frame, 3>,
        _hint: Option<Point<Frame, 3>>,
    ) -> Result<Point<Frame, 3>, ProjectError> {
        Ok(self.inverse().apply(p))
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::Translation3;
    use nalgebra::UnitQuaternion;

    use super::*;

    fn placement(t: [f64; 3], axis_angle: [f64; 3]) -> Isometry3<f64> {
        Isometry3::from_parts(
            Translation3::new(t[0], t[1], t[2]),
            UnitQuaternion::from_scaled_axis(Vector3::from(axis_angle)),
        )
    }

    #[test]
    fn sibling_transform_goes_through_common_ancestor() {
        let mut tree = FrameTree::new();
        let car = tree.add(tree.root(), placement([10.0, 0.0, 0.0], [0.0, 0.0, 0.5]));
        let wheel = tree.add(car, placement([1.0, 2.0, 0.0], [0.3, 0.0, 0.0]));
        let seat = tree.add(car, placement([0.0, 0.0, 1.0], [0.0, 0.2, 0.0]));

        let p = wheel.point(Vector3::new(0.1, 0.2, 0.3));
        let direct = tree.transform(wheel, seat).apply(p);
        let via_root = tree
            .transform(tree.root(), seat)
            .apply(tree.transform(wheel, tree.root()).apply(p));

        assert_eq!(direct.tag(), seat);
        assert!((direct.coords - via_root.coords).norm() < 1e-12);
    }

    #[test]
    fn project_inverts_apply() {
        let mut tree = FrameTree::new();
        let a = tree.add(tree.root(), placement([1.0, 2.0, 3.0], [0.1, 0.2, 0.3]));
        let t = tree.transform(a, tree.root());
        let p = a.point(Vector3::new(4.0, 5.0, 6.0));
        assert!((t.project(t.apply(p), None).unwrap().coords - p.coords).norm() < 1e-12);
    }

    #[test]
    #[should_panic(expected = "source frame")]
    fn applying_to_wrong_frame_panics() {
        let mut tree = FrameTree::new();
        let a = tree.add(tree.root(), Isometry3::identity());
        let b = tree.add(tree.root(), Isometry3::identity());
        tree.transform(a, b).apply(b.point(Vector3::zeros()));
    }

    #[test]
    #[should_panic(expected = "different instances")]
    fn subtracting_points_in_different_frames_panics() {
        let mut tree = FrameTree::new();
        let a = tree.add(tree.root(), Isometry3::identity());
        let _ = a.point(Vector3::zeros()) - tree.root().point(Vector3::zeros());
    }
}
