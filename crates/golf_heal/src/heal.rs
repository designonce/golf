use golf_brep::Body;
use golf_brep::EdgeId;
use golf_brep::FaceId;
use golf_brep::FaceUv;
use golf_geom::AnyCurve;
use golf_geom::AnyCurve2;
use golf_geom::AnySurface;
use golf_manifold::Domain;
use golf_manifold::Embedding;
use golf_manifold::Mapping;
use golf_manifold::Point;
use golf_manifold::Space;
use nalgebra::Vector2;
use nalgebra::Vector3;

use crate::error::HealIssue;
use crate::pcurve::fit;
use crate::pcurve::shifted;

/// How closely healed geometry must agree.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HealOptions {
    /// How far a pcurve, mapped through its surface, may stray from its edge.
    pub tolerance: f64,
}

impl Default for HealOptions {
    fn default() -> Self {
        Self { tolerance: 1e-3 }
    }
}

/// What [`heal`] did, and what it couldn't.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HealReport {
    /// Pcurves already present and accurate.
    pub kept: usize,
    /// Accurate pcurves moved by whole periods to join their loop.
    pub moved: usize,
    /// Pcurves computed, where there were none or they were inaccurate.
    pub computed: usize,
    pub issues: Vec<HealIssue>,
}

impl HealReport {
    pub fn is_clean(&self) -> bool {
        self.issues.is_empty()
    }
}

/// Gives every coedge of `body` an accurate pcurve; see the crate docs.
pub fn heal<S: Space<3>>(body: &mut Body<S>, options: &HealOptions) -> HealReport {
    let mut report = HealReport::default();
    let faces: Vec<FaceId> = body.faces().map(|(id, _)| id).collect();
    for face in faces {
        for loop_index in 0..body.face(face).loops.len() {
            heal_loop(body, face, loop_index, options, &mut report);
        }
    }
    report
}

/// A coedge's pcurve, once settled, and its uv at the loop's start and end
/// of it.
struct Settled<S: Space<3>> {
    pcurve: AnyCurve2<FaceUv<S>>,
    start: Vector2<f64>,
    end: Vector2<f64>,
}

fn heal_loop<S: Space<3>>(
    body: &mut Body<S>,
    face: FaceId,
    loop_index: usize,
    options: &HealOptions,
    report: &mut HealReport,
) {
    let tolerance = options.tolerance;
    let surface = body.face(face).surface.clone();
    let same_sense = body.face(face).same_sense;
    let domain = surface.domain();
    let coedges = body.face(face).loops[loop_index].coedges.clone();
    let n = coedges.len();
    // Each coedge's edge, its parameters in loop order, and its existing
    // pcurve if accurate.
    let uses: Vec<Use<S>> = coedges
        .iter()
        .map(|c| {
            let edge = body.edge(c.edge);
            let (a, b) = edge.range;
            let params = if c.reversed { (b, a) } else { (a, b) };
            let accurate = c
                .pcurve
                .clone()
                .filter(|p| max_error(&surface, &edge.curve, p, edge.range) <= tolerance);
            (c.edge, edge.curve.clone(), params, accurate)
        })
        .collect();

    // Start from an accurate pcurve if there is one: it fixes which period
    // the loop is in.
    let first = uses.iter().position(|u| u.3.is_some()).unwrap_or(0);
    let singular = |uv: Vector2<f64>| domain.is_singular(Point::<FaceUv<S>, 2>::new(uv), 1e-9);
    let mut settled: Vec<Option<Settled<S>>> = (0..n).map(|_| None).collect();
    let mut deferred = Vec::new();
    let mut previous_end: Option<Vector2<f64>> = None;
    for k in 0..n {
        let i = (first + k) % n;
        let (edge, curve, (t_start, t_end), accurate) = &uses[i];
        // The second use of a seam: the first's, a period across.
        if accurate.is_none() {
            if let Some(partner) = seam_partner(&uses, &settled, i) {
                if let Some(pcurve) = across_seam(partner, same_sense, &domain) {
                    report.computed += 1;
                    let s = settle(pcurve, (*t_start, *t_end));
                    previous_end = Some(s.end);
                    settled[i] = Some(s);
                    continue;
                }
            }
        }
        // Across a pole, the previous coedge's end says nothing about where
        // this one starts; it's settled from the next one instead.
        let hint = previous_end.filter(|&uv| !singular(uv));
        if previous_end.is_some() && hint.is_none() && accurate.is_none() {
            deferred.push(i);
            previous_end = None;
            continue;
        }
        let result = match accurate {
            Some(pcurve) => Some(align(pcurve, *t_start, hint, &domain, report)),
            None => {
                compute(&surface, curve, (*t_start, *t_end), hint, false, tolerance).inspect(|_| {
                    report.computed += 1;
                })
            }
        };
        match result {
            Some(pcurve) => {
                let s = settle(pcurve, (*t_start, *t_end));
                previous_end = Some(s.end);
                settled[i] = Some(s);
            }
            None => {
                report
                    .issues
                    .push(HealIssue::Projection { face, edge: *edge });
                previous_end = None;
            }
        }
    }
    // Deferred coedges, last first, each seeded from the start of the next.
    for &i in deferred.iter().rev() {
        let (edge, curve, (t_start, t_end), _) = &uses[i];
        if let Some(partner) = seam_partner(&uses, &settled, i) {
            if let Some(pcurve) = across_seam(partner, same_sense, &domain) {
                report.computed += 1;
                settled[i] = Some(settle(pcurve, (*t_start, *t_end)));
                continue;
            }
        }
        let next_start = settled[(i + 1) % n].as_ref().map(|s| s.start);
        match compute(
            &surface,
            curve,
            (*t_start, *t_end),
            next_start,
            true,
            tolerance,
        ) {
            Some(pcurve) => {
                report.computed += 1;
                settled[i] = Some(settle(pcurve, (*t_start, *t_end)));
            }
            None => report
                .issues
                .push(HealIssue::Projection { face, edge: *edge }),
        }
    }

    // A loop must close in parameter space (but for jumps along a pole).
    let ends: Vec<(Vector2<f64>, Vector2<f64>)> =
        settled.iter().flatten().map(|s| (s.start, s.end)).collect();
    if ends.len() == n {
        let opens = (0..n).any(|i| {
            let (end, next) = (ends[i].1, ends[(i + 1) % n].0);
            // A jump of half a period or more: the loop went round the
            // surface instead of enclosing part of it.
            !singular(end)
                && domain.axes.iter().enumerate().any(|(k, axis)| {
                    axis.periodic && (end[k] - next[k]).abs() > (axis.max - axis.min) / 2.0
                })
        });
        if opens {
            report
                .issues
                .push(HealIssue::OpenInParameters { face, loop_index });
        }
    }
    for (i, s) in settled.into_iter().enumerate() {
        if let Some(s) = s {
            body.set_pcurve(face, loop_index, i, Some(s.pcurve));
        }
    }
}

type Use<S> = (
    EdgeId,
    AnyCurve<S>,
    (f64, f64),
    Option<AnyCurve2<FaceUv<S>>>,
);

/// The settled other use, in the same loop, of coedge `i`'s edge: a seam.
fn seam_partner<'a, S: Space<3>>(
    uses: &[Use<S>],
    settled: &'a [Option<Settled<S>>],
    i: usize,
) -> Option<&'a Settled<S>> {
    (0..uses.len())
        .filter(|&j| j != i && uses[j].0 == uses[i].0)
        .find_map(|j| settled[j].as_ref())
}

/// The pcurve of a seam's second use: the first's moved a period across the
/// periodic axis it runs along, to the side the face is on. Walking the first
/// use, the face is on its left in parameter space (on its right if the face
/// is against its surface), and between the two uses.
fn across_seam<S: Space<3>>(
    first: &Settled<S>,
    same_sense: bool,
    domain: &Domain<2>,
) -> Option<AnyCurve2<FaceUv<S>>> {
    let along = first.end - first.start;
    let left = Vector2::new(-along.y, along.x) * if same_sense { 1.0 } else { -1.0 };
    // The periodic axis the seam crosses least along.
    let (k, axis) = domain
        .axes
        .iter()
        .enumerate()
        .filter(|(_, a)| a.periodic)
        .min_by(|(i, _), (j, _)| along[*i].abs().total_cmp(&along[*j].abs()))?;
    let mut shift = Vector2::zeros();
    shift[k] = (axis.max - axis.min) * left[k].signum();
    shifted(&first.pcurve, shift)
}

/// A coedge's pcurve with its uv at the coedge's start and end.
fn settle<S: Space<3>>(pcurve: AnyCurve2<FaceUv<S>>, (t_start, t_end): (f64, f64)) -> Settled<S> {
    let at = |t: f64| pcurve.apply(Point::new([t].into())).coords;
    Settled {
        start: at(t_start),
        end: at(t_end),
        pcurve,
    }
}

/// An accurate pcurve, moved by whole periods to start where the loop is.
fn align<S: Space<3>>(
    pcurve: &AnyCurve2<FaceUv<S>>,
    t_start: f64,
    hint: Option<Vector2<f64>>,
    domain: &Domain<2>,
    report: &mut HealReport,
) -> AnyCurve2<FaceUv<S>> {
    let start = pcurve.apply(Point::new([t_start].into())).coords;
    let shift = hint.map_or(Vector2::zeros(), |h| period_shift(start, h, domain));
    if shift == Vector2::zeros() {
        report.kept += 1;
        return pcurve.clone();
    }
    match shifted(pcurve, shift) {
        Some(moved) => {
            report.moved += 1;
            moved
        }
        None => {
            report.kept += 1;
            pcurve.clone()
        }
    }
}

/// The whole periods taking `uv` nearest `reference`.
fn period_shift(uv: Vector2<f64>, reference: Vector2<f64>, domain: &Domain<2>) -> Vector2<f64> {
    Vector2::from_fn(|k, _| {
        let axis = &domain.axes[k];
        match axis.periodic {
            true => {
                let period = axis.max - axis.min;
                ((reference[k] - uv[k]) / period).round() * period
            }
            false => 0.0,
        }
    })
}

/// A pcurve for the part of `curve` from `t_start` to `t_end` (in loop
/// order), by projecting samples onto `surface`, each from the last; the first
/// from `seed` (the loop's uv there, if known). With `backwards`, `seed` is
/// the uv at `t_end` and sampling runs from there.
fn compute<S: Space<3>>(
    surface: &AnySurface<S>,
    curve: &AnyCurve<S>,
    (t_start, t_end): (f64, f64),
    seed: Option<Vector2<f64>>,
    backwards: bool,
    tolerance: f64,
) -> Option<AnyCurve2<FaceUv<S>>> {
    const INITIAL: usize = 16;
    const MAX_SAMPLES: usize = 20_000;
    let domain = surface.domain();
    let (from, to) = if backwards {
        (t_end, t_start)
    } else {
        (t_start, t_end)
    };
    let point = |t: f64| curve.apply(Point::new([t].into()));
    let project = |t: f64, hint: Option<Vector2<f64>>| -> Option<Vector2<f64>> {
        let hint = hint.map(Point::<FaceUv<S>, 2>::new);
        let uv = surface.project(point(t), hint).ok()?;
        let uv = match hint {
            Some(h) => domain.unwrap_near(uv, h),
            None => uv,
        };
        Some(clamp(&domain, uv.coords)).filter(|uv| uv.iter().all(|x| x.is_finite()))
    };

    // Samples in the order walked, each projected from the one before.
    let mut samples: Vec<(f64, Vector2<f64>)> = Vec::with_capacity(INITIAL + 1);
    let mut hint = seed;
    for i in 0..=INITIAL {
        let t = from + (to - from) * i as f64 / INITIAL as f64;
        let uv = project(t, hint)?;
        samples.push((t, uv));
        hint = Some(uv);
    }
    // Halve any step whose straight uv segment strays from the curve.
    let raise = |uv: Vector2<f64>| surface.apply(Point::new(uv)).coords;
    let mut i = 0;
    while i + 1 < samples.len() && samples.len() < MAX_SAMPLES {
        let ((ta, a), (tb, b)) = (samples[i], samples[i + 1]);
        let tm = (ta + tb) / 2.0;
        let straight = raise(fix_poles(&domain, &[(ta, a), (tb, b)])[0].1.lerp(&b, 0.5));
        let target: Vector3<f64> = point(tm).coords;
        if (straight - target).norm() > tolerance / 2.0
            && (tb - ta).abs() > 1e-12 * (1.0 + ta.abs())
        {
            let uv = project(tm, Some(a))?;
            samples.insert(i + 1, (tm, uv));
        } else {
            i += 1;
        }
    }
    let mut samples = fix_poles(&domain, &samples);
    if samples.first()?.0 > samples.last()?.0 {
        samples.reverse();
    }
    // A line if it follows the edge as well; else the polyline.
    let line = fit::<S>(&[samples[0], *samples.last()?], f64::INFINITY)?;
    let range = (samples[0].0, samples.last()?.0);
    if max_error(surface, curve, &line, range) <= tolerance {
        return Some(line);
    }
    fit(&samples, 0.0)
}

/// `uv` kept within the bounds of the domain's non-periodic axes: just past a
/// pole, say, a surface's normal turns over.
fn clamp(domain: &Domain<2>, uv: Vector2<f64>) -> Vector2<f64> {
    Vector2::from_fn(|k, _| {
        let axis = &domain.axes[k];
        match axis.periodic {
            true => uv[k],
            false => uv[k].clamp(axis.min, axis.max),
        }
    })
}

/// Samples with any at a pole (where one parameter is arbitrary) given the
/// arbitrary parameter of their nearest neighbour off it.
fn fix_poles(domain: &Domain<2>, samples: &[(f64, Vector2<f64>)]) -> Vec<(f64, Vector2<f64>)> {
    let mut fixed = samples.to_vec();
    for (k, axis) in domain.axes.iter().enumerate() {
        let free = 1 - k;
        let at_pole = |uv: Vector2<f64>| axis.is_singular(uv[k], 1e-9 * (1.0 + uv[k].abs()));
        for i in 0..samples.len() {
            if !at_pole(samples[i].1) {
                continue;
            }
            let neighbour = (i + 1..samples.len())
                .chain((0..i).rev())
                .find(|&j| !at_pole(samples[j].1));
            if let Some(j) = neighbour {
                fixed[i].1[free] = samples[j].1[free];
            }
        }
    }
    fixed
}

/// The furthest `pcurve`, through `surface`, strays from `curve` over `range`,
/// at a few points.
fn max_error<S: Space<3>>(
    surface: &AnySurface<S>,
    curve: &AnyCurve<S>,
    pcurve: &AnyCurve2<FaceUv<S>>,
    (a, b): (f64, f64),
) -> f64 {
    (0..=8)
        .map(|i| {
            let t = a + (b - a) * i as f64 / 8.0;
            let on_surface = surface.apply(pcurve.apply(Point::new([t].into()))).coords;
            (on_surface - curve.apply(Point::new([t].into())).coords).norm()
        })
        .fold(
            0.0,
            |m: f64, e| if e.is_nan() { f64::INFINITY } else { m.max(e) },
        )
}

#[cfg(test)]
mod tests {
    use golf_geom::Placement;
    use golf_manifold::World;
    use golf_mesh::Tolerance;
    use golf_mesh::mesh;
    use golf_model::primitives;
    use nalgebra::Vector3;

    use super::*;

    /// Every pcurve stripped, as from a file without them.
    fn stripped(body: &Body<World>) -> Body<World> {
        let mut bare = body.clone();
        let faces: Vec<_> = bare
            .faces()
            .map(|(id, f)| {
                (
                    id,
                    f.loops.iter().map(|l| l.coedges.len()).collect::<Vec<_>>(),
                )
            })
            .collect();
        for (face, loops) in faces {
            for (l, count) in loops.into_iter().enumerate() {
                for c in 0..count {
                    bare.set_pcurve(face, l, c, None);
                }
            }
        }
        bare
    }

    fn origin() -> Placement<World> {
        Placement::from_axes(
            Point::new(Vector3::new(1.0, 2.0, 3.0)),
            golf_manifold::Vector::new(Vector3::new(0.3, -0.2, 1.0)),
            golf_manifold::Vector::new(Vector3::x()),
        )
    }

    #[test]
    fn stripped_pcurves_are_restored() {
        let bodies = [
            primitives::cuboid(&origin(), Vector3::new(4.0, 3.0, 2.0)).unwrap(),
            primitives::cylinder(&origin(), 2.0, 5.0).unwrap(),
            primitives::cone(&origin(), 2.0, 0.0, 4.0).unwrap(),
            primitives::cone(&origin(), 2.0, 1.0, 4.0).unwrap(),
            primitives::sphere(&origin(), 3.0).unwrap(),
            primitives::torus(&origin(), 5.0, 1.5).unwrap(),
            primitives::rounded_cuboid(&origin(), Vector3::new(8.0, 6.0, 4.0), 1.0).unwrap(),
        ];
        let tolerance = Tolerance::new(0.01, 0.3);
        for (index, body) in bodies.into_iter().enumerate() {
            let expected = mesh(&body, &tolerance).unwrap().signed_volume();
            let mut bare = stripped(&body);
            let report = heal(&mut bare, &HealOptions::default());
            assert!(report.is_clean(), "{index}: {:?}", report.issues);
            assert_eq!(report.kept, 0);
            assert_eq!(bare.validate(1e-6), Ok(()));
            let healed = mesh(&bare, &tolerance).unwrap();
            assert!(
                healed.is_watertight(),
                "{index}: {} open edges",
                healed.open_edges().len()
            );
            assert!(
                (healed.signed_volume() - expected).abs() < 1e-3 * expected.abs(),
                "{} vs {expected}",
                healed.signed_volume()
            );
        }
    }

    #[test]
    fn accurate_pcurves_are_kept() {
        // The side's pcurves are kept; only the caps' missing ones computed.
        let mut body = primitives::cylinder(&origin(), 2.0, 5.0).unwrap();
        let missing = body
            .faces()
            .flat_map(|(_, f)| f.loops.iter().flat_map(|l| &l.coedges))
            .filter(|c| c.pcurve.is_none())
            .count();
        let report = heal(&mut body, &HealOptions::default());
        assert!(report.is_clean());
        assert_eq!((report.computed, report.kept), (missing, 4));
    }
}
