//! Reading, healing and meshing a file, away from the window.

use std::path::Path;
use std::path::PathBuf;

use golf_color::Color;
use golf_import_step::read_step;
use golf_manifold::Mapping;
use golf_manifold::Point;
use golf_mesh::Tolerance;
use golf_mesh::mesh;
use golf_view::Scene;
use nalgebra::Vector3;

/// What a load sends back.
pub(crate) enum Loaded {
    /// A file was chosen and is loading.
    Started(PathBuf),
    /// No file was chosen.
    Cancelled,
    Failed(PathBuf, String),
    /// A file's bodies, meshed, and a line saying what came in.
    Scene {
        path: PathBuf,
        scene: Scene,
        summary: String,
    },
}

/// The grey for faces a file doesn't colour.
const PLAIN: Color = Color::rgb(170, 175, 180);

pub(crate) fn load(path: &Path) -> Loaded {
    let failed = |e: String| Loaded::Failed(path.to_path_buf(), e);
    let bytes = match std::fs::read(path) {
        Ok(bytes) => bytes,
        Err(e) => return failed(e.to_string()),
    };
    let import = match read_step(&bytes) {
        Ok(import) => import,
        Err(e) => return failed(e.to_string()),
    };
    let posed = match import.assembly.posed_bodies() {
        Ok(posed) => posed,
        Err(e) => return failed(e.to_string()),
    };
    // One tolerance for the whole file, from its overall size, so small
    // parts beside big ones aren't meshed needlessly finely.
    let size = extent(posed.iter().map(|(_, b)| b));
    let tolerance = Tolerance::new(size * 5e-4, 0.4);
    let mut scene = Scene::new();
    let mut failures = 0;
    for (_, body) in &posed {
        match mesh(body, &tolerance) {
            Ok(m) => {
                scene.add_mesh(&m, |&face| body.face_color(face).unwrap_or(PLAIN));
            }
            Err(_) => failures += 1,
        }
    }
    let mut summary = format!(
        "{} parts, {} placed",
        import.assembly.parts().len(),
        posed.len()
    );
    if failures > 0 {
        summary.push_str(&format!(", {failures} couldn't be meshed"));
    }
    if !import.warnings.is_empty() {
        summary.push_str(&format!(", {} warnings", import.warnings.len()));
    }
    Loaded::Scene {
        path: path.to_path_buf(),
        scene,
        summary,
    }
}

/// The diagonal of the box holding the bodies, from points along their edges
/// (a body's vertices can all coincide, as a torus's can).
fn extent<'a>(bodies: impl Iterator<Item = &'a golf_brep::Body<golf_manifold::World>>) -> f64 {
    let mut lo = Vector3::repeat(f64::INFINITY);
    let mut hi = Vector3::repeat(f64::NEG_INFINITY);
    for body in bodies {
        for (_, edge) in body.edges() {
            for i in 0..=8 {
                let t = edge.range.0 + (edge.range.1 - edge.range.0) * i as f64 / 8.0;
                let p = edge.curve.apply(Point::new([t].into())).coords;
                lo = lo.inf(&p);
                hi = hi.sup(&p);
            }
        }
    }
    let diagonal = (hi - lo).norm();
    if diagonal.is_finite() && diagonal > 0.0 {
        diagonal
    } else {
        1.0
    }
}
