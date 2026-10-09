use golf_brep::Body;
use golf_manifold::Space;
use nalgebra::Isometry3;

use crate::error::AssemblyError;
use crate::id::GroupId;
use crate::id::PartId;

/// A part definition: a named body in its own coordinates.
#[derive(Clone, Debug)]
pub struct Part<S: Space<3>> {
    pub name: String,
    pub body: Body<S>,
}

/// A group definition: a named sub-assembly placing its children.
#[derive(Clone, Debug)]
pub struct Group {
    pub name: String,
    pub children: Vec<Child>,
}

/// What a group places: a part, or another group.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Item {
    Part(PartId),
    Group(GroupId),
}

impl From<PartId> for Item {
    fn from(id: PartId) -> Self {
        Self::Part(id)
    }
}

impl From<GroupId> for Item {
    fn from(id: GroupId) -> Self {
        Self::Group(id)
    }
}

/// One placement in a group: `item`, moved by `placement` from its own
/// coordinates into the group's.
#[derive(Clone, Debug, PartialEq)]
pub struct Child {
    pub item: Item,
    pub placement: Isometry3<f64>,
}

/// One placed copy of a part, reached from the root.
#[derive(Clone, Debug, PartialEq)]
pub struct Occurrence {
    pub part: PartId,
    /// The groups passed through from the root, each with the index of the
    /// child taken in it; the last child is the part.
    pub path: Vec<(GroupId, usize)>,
    /// The part's motion from its own coordinates into the assembly's.
    pub placement: Isometry3<f64>,
}

/// Parts and groups, placed from a root group. See the crate docs.
#[derive(Clone, Debug)]
pub struct Assembly<S: Space<3>> {
    parts: Vec<Part<S>>,
    groups: Vec<Group>,
}

impl<S: Space<3>> Assembly<S> {
    /// An empty assembly whose root group is named `name`.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            parts: Vec::new(),
            groups: vec![Group {
                name: name.into(),
                children: Vec::new(),
            }],
        }
    }

    /// The group everything is placed from.
    pub fn root(&self) -> GroupId {
        GroupId::new(0)
    }

    /// Defines a part, to be placed with [`Self::place`].
    pub fn add_part(&mut self, name: impl Into<String>, body: Body<S>) -> PartId {
        self.parts.push(Part {
            name: name.into(),
            body,
        });
        PartId::new(self.parts.len() - 1)
    }

    /// Defines an empty group, to be filled and placed with [`Self::place`].
    pub fn add_group(&mut self, name: impl Into<String>) -> GroupId {
        self.groups.push(Group {
            name: name.into(),
            children: Vec::new(),
        });
        GroupId::new(self.groups.len() - 1)
    }

    /// Places `item` in `parent` by `placement`, returning its index among
    /// `parent`'s children. A group may not end up inside itself.
    pub fn place(
        &mut self,
        parent: GroupId,
        item: impl Into<Item>,
        placement: Isometry3<f64>,
    ) -> Result<usize, AssemblyError> {
        let item = item.into();
        self.try_group(parent)?;
        match item {
            Item::Part(part) => {
                self.try_part(part)?;
            }
            Item::Group(child) => {
                self.try_group(child)?;
                if child == parent || self.contains(child, parent) {
                    return Err(AssemblyError::Cycle { parent, child });
                }
            }
        }
        let children = &mut self.groups[parent.index()].children;
        children.push(Child { item, placement });
        Ok(children.len() - 1)
    }

    /// Moves child `index` of `group` to a new placement.
    pub fn set_placement(
        &mut self,
        group: GroupId,
        index: usize,
        placement: Isometry3<f64>,
    ) -> Result<(), AssemblyError> {
        self.try_group(group)?;
        let child = self.groups[group.index()]
            .children
            .get_mut(index)
            .ok_or(AssemblyError::MissingChild { group, index })?;
        child.placement = placement;
        Ok(())
    }

    pub fn part(&self, id: PartId) -> &Part<S> {
        &self.parts[id.index()]
    }

    pub fn part_mut(&mut self, id: PartId) -> &mut Part<S> {
        &mut self.parts[id.index()]
    }

    pub fn group(&self, id: GroupId) -> &Group {
        &self.groups[id.index()]
    }

    pub fn parts(&self) -> impl ExactSizeIterator<Item = (PartId, &Part<S>)> {
        self.parts
            .iter()
            .enumerate()
            .map(|(i, p)| (PartId::new(i), p))
    }

    pub fn groups(&self) -> impl ExactSizeIterator<Item = (GroupId, &Group)> {
        self.groups
            .iter()
            .enumerate()
            .map(|(i, g)| (GroupId::new(i), g))
    }

    /// Every placed copy of a part reachable from the root, depth first in
    /// placement order.
    pub fn occurrences(&self) -> Vec<Occurrence> {
        let mut out = Vec::new();
        self.walk(
            self.root(),
            Isometry3::identity(),
            &mut Vec::new(),
            &mut out,
        );
        out
    }

    /// How many times `part` is placed, counting through every group.
    pub fn instance_count(&self, part: PartId) -> usize {
        self.occurrences().iter().filter(|o| o.part == part).count()
    }

    /// The motion from `from`'s coordinates to `to`'s, for two occurrences of
    /// this assembly.
    pub fn motion_between(&self, from: &Occurrence, to: &Occurrence) -> Isometry3<f64> {
        to.placement.inverse() * from.placement
    }

    /// Each occurrence's body, moved into the assembly's coordinates, with the
    /// occurrence's name (its path of group and part names, joined by `/`).
    pub fn posed_bodies(&self) -> Result<Vec<(String, Body<S>)>, AssemblyError> {
        self.occurrences()
            .iter()
            .map(|o| {
                Ok((
                    self.occurrence_name(o),
                    self.part(o.part).body.transformed(&o.placement)?,
                ))
            })
            .collect()
    }

    /// An occurrence's path of names below the root, ending with its part's,
    /// joined by `/`, with `#n` after any name used more than once in the same
    /// group.
    pub fn occurrence_name(&self, occurrence: &Occurrence) -> String {
        let mut names = Vec::new();
        for &(group, index) in &occurrence.path {
            let children = &self.groups[group.index()].children;
            let item = children[index].item;
            let name = match item {
                Item::Part(p) => &self.parts[p.index()].name,
                Item::Group(g) => &self.groups[g.index()].name,
            };
            let repeats = children.iter().filter(|c| c.item == item).count();
            let ordinal = children[..index].iter().filter(|c| c.item == item).count() + 1;
            names.push(match repeats {
                1 => name.clone(),
                _ => format!("{name}#{ordinal}"),
            });
        }
        names.join("/")
    }

    /// Copies `other`'s parts and groups into this assembly, its root becoming
    /// a group here, which is returned for placing.
    pub fn merge(&mut self, other: Assembly<S>) -> GroupId {
        let (part_offset, group_offset) = (self.parts.len(), self.groups.len());
        self.parts.extend(other.parts);
        self.groups.extend(other.groups.into_iter().map(|mut g| {
            for child in &mut g.children {
                child.item = match child.item {
                    Item::Part(p) => Item::Part(PartId::new(p.index() + part_offset)),
                    Item::Group(g) => Item::Group(GroupId::new(g.index() + group_offset)),
                };
            }
            g
        }));
        GroupId::new(group_offset)
    }

    fn walk(
        &self,
        group: GroupId,
        motion: Isometry3<f64>,
        path: &mut Vec<(GroupId, usize)>,
        out: &mut Vec<Occurrence>,
    ) {
        for (index, child) in self.groups[group.index()].children.iter().enumerate() {
            path.push((group, index));
            let placed = motion * child.placement;
            match child.item {
                Item::Part(part) => out.push(Occurrence {
                    part,
                    path: path.clone(),
                    placement: placed,
                }),
                Item::Group(inner) => self.walk(inner, placed, path, out),
            }
            path.pop();
        }
    }

    /// Whether `group` contains `target` at any depth.
    fn contains(&self, group: GroupId, target: GroupId) -> bool {
        self.groups[group.index()]
            .children
            .iter()
            .any(|c| match c.item {
                Item::Group(g) => g == target || self.contains(g, target),
                Item::Part(_) => false,
            })
    }

    fn try_part(&self, id: PartId) -> Result<&Part<S>, AssemblyError> {
        self.parts
            .get(id.index())
            .ok_or(AssemblyError::MissingPart(id))
    }

    fn try_group(&self, id: GroupId) -> Result<&Group, AssemblyError> {
        self.groups
            .get(id.index())
            .ok_or(AssemblyError::MissingGroup(id))
    }
}

#[cfg(test)]
mod tests {
    use core::f64::consts::FRAC_PI_2;

    use golf_geom::Placement;
    use golf_manifold::Point;
    use golf_manifold::World;
    use golf_model::primitives;
    use nalgebra::Translation3;
    use nalgebra::UnitQuaternion;
    use nalgebra::Vector3;

    use super::*;

    fn at(x: f64, y: f64, z: f64) -> Isometry3<f64> {
        Isometry3::translation(x, y, z)
    }

    /// A cart: a body, and an axle group (axle + two wheels) placed twice.
    fn cart() -> (Assembly<World>, PartId, PartId, GroupId) {
        let origin = Placement::<World>::at(Point::new(Vector3::zeros()));
        let mut cart = Assembly::new("cart");
        let chassis = cart.add_part(
            "chassis",
            primitives::cuboid(&origin, Vector3::new(40.0, 20.0, 5.0)).unwrap(),
        );
        let wheel = cart.add_part("wheel", primitives::cylinder(&origin, 6.0, 3.0).unwrap());
        let axle = cart.add_group("axle");
        let side = |y: f64| {
            Isometry3::from_parts(
                Translation3::new(0.0, y, 0.0),
                UnitQuaternion::from_scaled_axis(Vector3::x() * FRAC_PI_2),
            )
        };
        cart.place(axle, wheel, side(-1.0)).unwrap();
        cart.place(axle, wheel, side(24.0)).unwrap();
        cart.place(cart.root(), chassis, Isometry3::identity())
            .unwrap();
        cart.place(cart.root(), axle, at(5.0, 0.0, 0.0)).unwrap();
        cart.place(cart.root(), axle, at(35.0, 0.0, 0.0)).unwrap();
        (cart, chassis, wheel, axle)
    }

    #[test]
    fn occurrences_compose_placements() {
        let (cart, chassis, wheel, axle) = cart();
        let occurrences = cart.occurrences();
        assert_eq!(occurrences.len(), 5);
        assert_eq!(cart.instance_count(wheel), 4);
        assert_eq!(cart.instance_count(chassis), 1);
        // The second axle's second wheel: axle placement, then wheel placement.
        let last = &occurrences[4];
        assert_eq!(last.path, vec![(cart.root(), 2), (axle, 1)]);
        let origin = last.placement * nalgebra::Point3::origin();
        assert!((origin.coords - Vector3::new(35.0, 24.0, 0.0)).norm() < 1e-12);
        // Wheels turned onto their sides: their axis (local z) is along -y.
        assert!(
            (last.placement.rotation * Vector3::z() - Vector3::new(0.0, -1.0, 0.0)).norm() < 1e-12
        );
        assert_eq!(cart.occurrence_name(last), "axle#2/wheel#2");
        assert_eq!(cart.occurrence_name(&occurrences[0]), "chassis");
    }

    #[test]
    fn motion_between_occurrences() {
        let (cart, ..) = cart();
        let o = cart.occurrences();
        let between = cart.motion_between(&o[1], &o[3]);
        // Wheels on the same side, one axle apart, both turned the same way.
        let world_a = o[1].placement * nalgebra::Point3::new(1.0, 2.0, 3.0);
        let world_b = o[3].placement * (between * nalgebra::Point3::new(1.0, 2.0, 3.0));
        assert!((world_a - world_b).norm() < 1e-12);
    }

    #[test]
    fn posed_bodies_are_moved_copies() {
        let (cart, ..) = cart();
        let posed = cart.posed_bodies().unwrap();
        assert_eq!(posed.len(), 5);
        for (name, body) in &posed {
            assert_eq!(body.validate(1e-9), Ok(()), "{name}");
        }
        let wheel = &posed[4].1;
        let lowest_y = wheel
            .vertices()
            .map(|(_, v)| v.point.coords.y)
            .fold(f64::INFINITY, f64::min);
        assert!(lowest_y > 20.0, "{lowest_y}");
    }

    #[test]
    fn cycles_and_missing_ids_are_rejected() {
        let (mut cart, _, _, axle) = cart();
        let root = cart.root();
        assert_eq!(
            cart.place(axle, root, Isometry3::identity()),
            Err(AssemblyError::Cycle {
                parent: axle,
                child: root
            })
        );
        assert_eq!(
            cart.place(axle, axle, Isometry3::identity()),
            Err(AssemblyError::Cycle {
                parent: axle,
                child: axle
            })
        );
        assert_eq!(
            cart.place(root, PartId(9), Isometry3::identity()),
            Err(AssemblyError::MissingPart(PartId(9)))
        );
        assert!(matches!(
            cart.set_placement(root, 7, Isometry3::identity()),
            Err(AssemblyError::MissingChild { .. })
        ));
    }

    #[test]
    fn merging_assemblies() {
        let (cart, ..) = cart();
        let mut yard = Assembly::<World>::new("yard");
        let first = yard.merge(cart.clone());
        let second = yard.merge(cart);
        yard.place(yard.root(), first, Isometry3::identity())
            .unwrap();
        yard.place(yard.root(), second, at(0.0, 100.0, 0.0))
            .unwrap();
        assert_eq!(yard.occurrences().len(), 10);
        assert_eq!(yard.parts().len(), 4);
        let shifted = &yard.occurrences()[9];
        assert!(
            (shifted.placement * nalgebra::Point3::origin()
                - nalgebra::Point3::new(35.0, 124.0, 0.0))
            .norm()
                < 1e-12
        );
    }
}
