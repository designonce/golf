use std::collections::HashMap;

use golf_brep::Body;
use golf_color::Color;
use golf_frame::FrameId;
use golf_frame::FrameTree;
use golf_frame::Rigid;
use golf_manifold::Space;
use nalgebra::Isometry3;

use crate::error::AssemblyError;
use crate::id::PartId;

/// A part definition: a named body in its own coordinates, placed by
/// instances. The body's colour is the part's, unless an instance overrides
/// it.
#[derive(Clone, Debug)]
pub struct Part<S: Space<3>> {
    pub name: String,
    pub body: Body<S>,
}

/// What a frame of the assembly is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NodeKind {
    /// A group (sub-assembly) placing other frames.
    Group,
    /// A placed copy of a part.
    Instance(PartId),
}

/// A frame of the assembly: a group or a part instance.
#[derive(Clone, Debug, PartialEq)]
pub struct Node {
    pub name: String,
    pub kind: NodeKind,
    /// The frames placed in this one, in order. Empty for an instance.
    pub children: Vec<FrameId>,
    /// A colour for this frame and everything below it, overriding parts'
    /// own.
    pub color: Option<Color>,
}

/// One placed part, as found walking the assembly.
#[derive(Clone, Debug, PartialEq)]
pub struct Occurrence {
    pub frame: FrameId,
    pub part: PartId,
    /// The motion from the part's coordinates into the assembly's (the root
    /// frame's).
    pub placement: Isometry3<f64>,
    /// The colour it shows: the nearest override at or above it, else the
    /// part's.
    pub color: Option<Color>,
}

/// Parts placed in a [`FrameTree`] of groups. See the crate docs.
#[derive(Clone, Debug)]
pub struct Assembly<S: Space<3>> {
    parts: Vec<Part<S>>,
    frames: FrameTree,
    nodes: HashMap<FrameId, Node>,
}

impl<S: Space<3>> Assembly<S> {
    /// An empty assembly whose root group is named `name`.
    pub fn new(name: impl Into<String>) -> Self {
        let frames = FrameTree::new();
        let root = Node {
            name: name.into(),
            kind: NodeKind::Group,
            children: Vec::new(),
            color: None,
        };
        let nodes = HashMap::from([(frames.root(), root)]);
        Self {
            parts: Vec::new(),
            frames,
            nodes,
        }
    }

    /// The root group, whose coordinates are the assembly's.
    pub fn root(&self) -> FrameId {
        self.frames.root()
    }

    /// The frames everything is placed in: one per group and per instance.
    pub fn frames(&self) -> &FrameTree {
        &self.frames
    }

    /// The motion taking `from`'s coordinates to `to`'s, for any two frames.
    pub fn transform(&self, from: FrameId, to: FrameId) -> Rigid {
        self.frames.transform(from, to)
    }

    /// Defines a part, to be placed with [`Self::add_instance`].
    pub fn add_part(&mut self, name: impl Into<String>, body: Body<S>) -> PartId {
        self.parts.push(Part {
            name: name.into(),
            body,
        });
        PartId::new(self.parts.len() - 1)
    }

    /// Adds a group in `parent`, placed by `placement` (its coordinates to
    /// `parent`'s).
    pub fn add_group(
        &mut self,
        parent: FrameId,
        name: impl Into<String>,
        placement: Isometry3<f64>,
    ) -> Result<FrameId, AssemblyError> {
        self.add_node(parent, name.into(), NodeKind::Group, placement)
    }

    /// Places `part` in the group `parent` by `placement` (the part's
    /// coordinates to `parent`'s).
    pub fn add_instance(
        &mut self,
        parent: FrameId,
        part: PartId,
        placement: Isometry3<f64>,
    ) -> Result<FrameId, AssemblyError> {
        let name = self
            .parts
            .get(part.index())
            .ok_or(AssemblyError::MissingPart(part))?
            .name
            .clone();
        self.add_node(parent, name, NodeKind::Instance(part), placement)
    }

    /// Copies the group `source`, everything in it included, into `parent` by
    /// `placement`. Parts are shared, not copied.
    pub fn copy_group(
        &mut self,
        source: FrameId,
        parent: FrameId,
        placement: Isometry3<f64>,
    ) -> Result<FrameId, AssemblyError> {
        let node = self.try_node(source)?.clone();
        if node.kind != NodeKind::Group {
            return Err(AssemblyError::NotAGroup(source));
        }
        let copy = self.add_node(parent, node.name, NodeKind::Group, placement)?;
        self.nodes.get_mut(&copy).expect("just added").color = node.color;
        // The children as they were, so copying into the group itself stops.
        for child in node.children {
            let child_node = self.nodes[&child].clone();
            let child_placement = self.frames.placement(child);
            match child_node.kind {
                NodeKind::Group => {
                    self.copy_group(child, copy, child_placement)?;
                }
                NodeKind::Instance(part) => {
                    let instance = self.add_instance(copy, part, child_placement)?;
                    self.nodes.get_mut(&instance).expect("just added").color = child_node.color;
                }
            }
        }
        Ok(copy)
    }

    /// Copies `other` into this assembly, its root becoming a group in
    /// `parent` placed by `placement`.
    pub fn merge(
        &mut self,
        other: &Assembly<S>,
        parent: FrameId,
        placement: Isometry3<f64>,
    ) -> Result<FrameId, AssemblyError> {
        let offset = self.parts.len();
        self.parts.extend(other.parts.iter().cloned());
        let root = other.root();
        let group = self.add_group(parent, other.nodes[&root].name.clone(), placement)?;
        self.nodes.get_mut(&group).expect("just added").color = other.nodes[&root].color;
        self.graft(other, root, group, offset)?;
        Ok(group)
    }

    /// Moves a frame (not the root) relative to its parent.
    pub fn set_placement(
        &mut self,
        frame: FrameId,
        placement: Isometry3<f64>,
    ) -> Result<(), AssemblyError> {
        self.try_node(frame)?;
        if frame == self.root() {
            return Err(AssemblyError::RootPlacement);
        }
        self.frames.set_placement(frame, placement);
        Ok(())
    }

    /// Sets (or clears) the colour of a frame and everything below it.
    pub fn set_color(&mut self, frame: FrameId, color: Option<Color>) -> Result<(), AssemblyError> {
        self.try_node(frame)?;
        self.nodes.get_mut(&frame).expect("checked").color = color;
        Ok(())
    }

    /// The colour a frame shows: its own or the nearest above it, else (for
    /// an instance) its part's.
    pub fn color(&self, frame: FrameId) -> Option<Color> {
        let mut at = Some(frame);
        while let Some(f) = at {
            if let Some(color) = self.nodes.get(&f).and_then(|n| n.color) {
                return Some(color);
            }
            at = self.frames.parent(f);
        }
        match self.nodes.get(&frame)?.kind {
            NodeKind::Instance(part) => self.parts[part.index()].body.color(),
            NodeKind::Group => None,
        }
    }

    pub fn part(&self, id: PartId) -> &Part<S> {
        &self.parts[id.index()]
    }

    pub fn part_mut(&mut self, id: PartId) -> &mut Part<S> {
        &mut self.parts[id.index()]
    }

    pub fn parts(&self) -> impl ExactSizeIterator<Item = (PartId, &Part<S>)> {
        self.parts
            .iter()
            .enumerate()
            .map(|(i, p)| (PartId::new(i), p))
    }

    pub fn node(&self, frame: FrameId) -> &Node {
        &self.nodes[&frame]
    }

    /// Every placed part, depth first from the root in placement order.
    pub fn occurrences(&self) -> Vec<Occurrence> {
        let mut out = Vec::new();
        let mut stack = vec![self.root()];
        while let Some(frame) = stack.pop() {
            let node = &self.nodes[&frame];
            match node.kind {
                NodeKind::Group => stack.extend(node.children.iter().rev()),
                NodeKind::Instance(part) => out.push(Occurrence {
                    frame,
                    part,
                    placement: self.frames.transform(frame, self.root()).iso,
                    color: self.color(frame),
                }),
            }
        }
        out
    }

    /// How many times `part` is placed.
    pub fn instance_count(&self, part: PartId) -> usize {
        self.nodes
            .values()
            .filter(|n| n.kind == NodeKind::Instance(part))
            .count()
    }

    /// A frame's path of names below the root, joined by `/`, with `#n` after
    /// a name shared with siblings.
    pub fn path_name(&self, frame: FrameId) -> String {
        let mut names = Vec::new();
        let mut at = frame;
        while let Some(parent) = self.frames.parent(at) {
            let siblings = &self.nodes[&parent].children;
            let name = &self.nodes[&at].name;
            let same: Vec<FrameId> = siblings
                .iter()
                .copied()
                .filter(|s| &self.nodes[s].name == name)
                .collect();
            names.push(match same.len() {
                1 => name.clone(),
                _ => format!(
                    "{name}#{}",
                    same.iter().position(|&s| s == at).expect("a sibling") + 1
                ),
            });
            at = parent;
        }
        names.reverse();
        names.join("/")
    }

    /// Each placed part's body moved into the assembly's coordinates and given
    /// the colour it shows, with its path name.
    pub fn posed_bodies(&self) -> Result<Vec<(String, Body<S>)>, AssemblyError> {
        self.occurrences()
            .iter()
            .map(|o| {
                let mut body = self.parts[o.part.index()].body.transformed(&o.placement)?;
                body.set_color(o.color);
                Ok((self.path_name(o.frame), body))
            })
            .collect()
    }

    fn add_node(
        &mut self,
        parent: FrameId,
        name: String,
        kind: NodeKind,
        placement: Isometry3<f64>,
    ) -> Result<FrameId, AssemblyError> {
        if self.try_node(parent)?.kind != NodeKind::Group {
            return Err(AssemblyError::NotAGroup(parent));
        }
        let frame = self.frames.add(parent, placement);
        self.nodes.insert(
            frame,
            Node {
                name,
                kind,
                children: Vec::new(),
                color: None,
            },
        );
        self.nodes
            .get_mut(&parent)
            .expect("checked")
            .children
            .push(frame);
        Ok(frame)
    }

    /// Copies `from`'s children in `other` under `to` here, shifting part ids
    /// by `offset`.
    fn graft(
        &mut self,
        other: &Assembly<S>,
        from: FrameId,
        to: FrameId,
        offset: usize,
    ) -> Result<(), AssemblyError> {
        for &child in &other.nodes[&from].children {
            let node = &other.nodes[&child];
            let placement = other.frames.placement(child);
            let kind = match node.kind {
                NodeKind::Group => NodeKind::Group,
                NodeKind::Instance(part) => NodeKind::Instance(PartId::new(part.index() + offset)),
            };
            let copy = self.add_node(to, node.name.clone(), kind, placement)?;
            self.nodes.get_mut(&copy).expect("just added").color = node.color;
            self.graft(other, child, copy, offset)?;
        }
        Ok(())
    }

    fn try_node(&self, frame: FrameId) -> Result<&Node, AssemblyError> {
        self.nodes
            .get(&frame)
            .ok_or(AssemblyError::MissingFrame(frame))
    }
}

#[cfg(test)]
mod tests {
    use core::f64::consts::FRAC_PI_2;

    use golf_geom::Placement;
    use golf_manifold::Mapping;
    use golf_manifold::Point;
    use golf_manifold::World;
    use golf_model::primitives;
    use nalgebra::Point3;
    use nalgebra::Translation3;
    use nalgebra::UnitQuaternion;
    use nalgebra::Vector3;

    use super::*;

    const RED: Color = Color::rgb(200, 30, 30);
    const BLACK: Color = Color::rgb(20, 20, 20);

    /// A cart: a chassis, and an axle group (two wheels) copied to the back.
    struct Cart {
        asm: Assembly<World>,
        chassis: PartId,
        wheel: PartId,
        front: FrameId,
        back: FrameId,
    }

    fn cart() -> Cart {
        let origin = Placement::<World>::at(Point::new(Vector3::zeros()));
        let mut asm = Assembly::new("cart");
        let chassis = asm.add_part(
            "chassis",
            primitives::cuboid(&origin, Vector3::new(40.0, 20.0, 5.0)).unwrap(),
        );
        let wheel = asm.add_part(
            "wheel",
            primitives::cylinder(&origin, 6.0, 3.0)
                .unwrap()
                .with_color(BLACK),
        );
        let root = asm.root();
        asm.add_instance(root, chassis, Isometry3::identity())
            .unwrap();
        let front = asm
            .add_group(root, "axle", Isometry3::translation(5.0, 0.0, 0.0))
            .unwrap();
        let side = |y: f64| {
            Isometry3::from_parts(
                Translation3::new(0.0, y, 0.0),
                UnitQuaternion::from_scaled_axis(Vector3::x() * FRAC_PI_2),
            )
        };
        asm.add_instance(front, wheel, side(-1.0)).unwrap();
        asm.add_instance(front, wheel, side(24.0)).unwrap();
        let back = asm
            .copy_group(front, root, Isometry3::translation(35.0, 0.0, 0.0))
            .unwrap();
        Cart {
            asm,
            chassis,
            wheel,
            front,
            back,
        }
    }

    #[test]
    fn occurrences_compose_frames() {
        let Cart {
            asm,
            chassis,
            wheel,
            back,
            ..
        } = cart();
        let occurrences = asm.occurrences();
        assert_eq!(occurrences.len(), 5);
        assert_eq!(
            (asm.instance_count(wheel), asm.instance_count(chassis)),
            (4, 1)
        );
        let last = &occurrences[4];
        assert_eq!(asm.frames().parent(last.frame), Some(back));
        assert!((last.placement * Point3::origin() - Point3::new(35.0, 24.0, 0.0)).norm() < 1e-12);
        // Wheels turned onto their sides: their axis (local z) along -y.
        assert!(
            (last.placement.rotation * Vector3::z() - Vector3::new(0.0, -1.0, 0.0)).norm() < 1e-12
        );
        assert_eq!(asm.path_name(last.frame), "axle#2/wheel#2");
        assert_eq!(asm.path_name(occurrences[0].frame), "chassis");
    }

    #[test]
    fn frames_relate_any_two_instances() {
        let Cart { mut asm, back, .. } = cart();
        let o = asm.occurrences();
        // A point on one wheel, seen from another: through the frame tree.
        let p = o[1].frame.point(Vector3::new(1.0, 2.0, 3.0));
        let seen = asm.transform(o[1].frame, o[3].frame).apply(p);
        assert_eq!(seen.tag(), o[3].frame);
        let world_a = o[1].placement * Point3::new(1.0, 2.0, 3.0);
        let world_b = o[3].placement * Point3::from(seen.coords);
        assert!((world_a - world_b).norm() < 1e-12);
        // Moving a group moves what's in it.
        asm.set_placement(back, Isometry3::translation(50.0, 0.0, 0.0))
            .unwrap();
        let moved = asm.occurrences();
        assert!(
            (moved[4].placement * Point3::origin() - Point3::new(50.0, 24.0, 0.0)).norm() < 1e-12
        );
    }

    #[test]
    fn colours_come_from_the_nearest_override_or_the_part() {
        let Cart { mut asm, front, .. } = cart();
        let o = asm.occurrences();
        assert_eq!(o[0].color, None);
        assert!(o[1..].iter().all(|w| w.color == Some(BLACK)));
        // Red front axle: both its wheels.
        asm.set_color(front, Some(RED)).unwrap();
        let o = asm.occurrences();
        assert_eq!(
            [o[1].color, o[2].color, o[3].color],
            [Some(RED), Some(RED), Some(BLACK)]
        );
        // An instance's own colour beats its group's.
        asm.set_color(o[2].frame, Some(Color::rgb(0, 200, 0)))
            .unwrap();
        assert_eq!(asm.occurrences()[2].color, Some(Color::rgb(0, 200, 0)));
        let posed = asm.posed_bodies().unwrap();
        assert_eq!(posed[1].1.color(), Some(RED));
        for (name, body) in &posed {
            assert_eq!(body.validate(1e-9), Ok(()), "{name}");
        }
    }

    #[test]
    fn invalid_frames_are_rejected() {
        let Cart { mut asm, wheel, .. } = cart();
        let instance = asm.occurrences()[1].frame;
        let root = asm.root();
        assert_eq!(
            asm.add_instance(instance, wheel, Isometry3::identity()),
            Err(AssemblyError::NotAGroup(instance))
        );
        assert_eq!(
            asm.add_instance(root, PartId(9), Isometry3::identity()),
            Err(AssemblyError::MissingPart(PartId(9)))
        );
        assert_eq!(
            asm.set_placement(root, Isometry3::identity()),
            Err(AssemblyError::RootPlacement)
        );
        assert_eq!(
            asm.copy_group(instance, root, Isometry3::identity()),
            Err(AssemblyError::NotAGroup(instance))
        );
    }

    #[test]
    fn copying_a_group_into_itself_copies_it_once() {
        let Cart { mut asm, front, .. } = cart();
        asm.copy_group(front, front, Isometry3::translation(0.0, 0.0, 10.0))
            .unwrap();
        assert_eq!(asm.node(front).children.len(), 3);
        assert_eq!(asm.occurrences().len(), 7);
    }

    #[test]
    fn merging_assemblies() {
        let Cart { asm: cart, .. } = cart();
        let mut yard = Assembly::<World>::new("yard");
        let root = yard.root();
        yard.merge(&cart, root, Isometry3::identity()).unwrap();
        yard.merge(&cart, root, Isometry3::translation(0.0, 100.0, 0.0))
            .unwrap();
        assert_eq!(yard.occurrences().len(), 10);
        assert_eq!(yard.parts().len(), 4);
        let shifted = &yard.occurrences()[9];
        assert!(
            (shifted.placement * Point3::origin() - Point3::new(35.0, 124.0, 0.0)).norm() < 1e-12
        );
        assert_eq!(yard.path_name(shifted.frame), "cart#2/axle#2/wheel#2");
    }
}
