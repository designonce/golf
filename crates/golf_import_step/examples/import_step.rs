//! Imports STEP files and reports what came in: parts, placements, warnings,
//! what healing did, and whether every placed body meshes watertight.
//!
//! ```text
//! cargo run --release -p golf_import_step --example import_step -- part.step [more.step...]
//! cargo run --release -p golf_import_step --example import_step -- --quiet corpus/*.step
//! ```

use std::error::Error;
use std::time::Instant;

use golf_import_step::read_step;
use golf_mesh::Tolerance;
use golf_mesh::mesh;

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let quiet = args.iter().any(|a| a == "--quiet");
    for path in args.iter().filter(|a| !a.starts_with("--")) {
        let bytes = std::fs::read(path)?;
        let start = Instant::now();
        let import = match read_step(&bytes) {
            Ok(import) => import,
            Err(e) => {
                println!("{path}: FAILED {e}");
                continue;
            }
        };
        let read = start.elapsed();
        let asm = &import.assembly;
        let posed = asm.posed_bodies()?;
        let mut faces = 0;
        let mut meshed = 0;
        let mut watertight = 0;
        let mut failures = Vec::new();
        for (name, body) in &posed {
            faces += body.faces().len();
            let size = body
                .vertices()
                .map(|(_, v)| v.point.coords.norm())
                .fold(1.0, f64::max);
            match mesh(body, &Tolerance::new(size * 1e-3, 0.5)) {
                Ok(m) => {
                    meshed += 1;
                    if m.is_watertight() {
                        watertight += 1;
                    }
                }
                Err(e) => failures.push(format!("{name}: {e}")),
            }
        }
        println!(
            "{path}: {} parts, {} placed, {faces} faces in {:.0} ms; {} warnings; pcurves kept {} moved {} computed {}; seams {}; meshed {meshed}/{} ({watertight} watertight)",
            asm.parts().len(),
            posed.len(),
            read.as_secs_f64() * 1e3,
            import.warnings.len(),
            import.heal.kept,
            import.heal.moved,
            import.heal.computed,
            import.heal.seams,
            posed.len(),
        );
        if !quiet {
            for occurrence in asm.occurrences() {
                let at = occurrence.placement * nalgebra::Point3::origin();
                println!(
                    "  {:<40} at ({:8.2}, {:8.2}, {:8.2}) {}",
                    asm.path_name(occurrence.frame),
                    at.x,
                    at.y,
                    at.z,
                    occurrence.color.map_or(String::new(), |c| c.to_string())
                );
            }
        }
        for warning in import.warnings.iter().take(if quiet { 3 } else { 20 }) {
            println!("  warning {warning}");
        }
        for failure in failures.iter().take(if quiet { 3 } else { 20 }) {
            println!("  mesh {failure}");
        }
    }
    Ok(())
}
