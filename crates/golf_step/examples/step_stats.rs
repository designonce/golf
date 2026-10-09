//! Parses STEP files and counts their entities by keyword.
//!
//! ```text
//! cargo run --release -p golf_step --example step_stats -- part.step [more.step...]
//! ```

use std::collections::BTreeMap;
use std::error::Error;
use std::time::Instant;

use golf_step::StepData;

fn main() -> Result<(), Box<dyn Error>> {
    for path in std::env::args().skip(1) {
        let bytes = std::fs::read(&path)?;
        let start = Instant::now();
        let data = StepData::parse(&bytes)?;
        let elapsed = start.elapsed();
        println!(
            "{path}: {} instances in {:.1} ms ({:.0} MB/s), {} skipped",
            data.len(),
            elapsed.as_secs_f64() * 1e3,
            bytes.len() as f64 / 1e6 / elapsed.as_secs_f64(),
            data.issues.len()
        );
        for issue in data.issues.iter().take(5) {
            println!("  line {}: {}", issue.line, issue.message);
        }
        let mut counts: BTreeMap<String, usize> = BTreeMap::new();
        for entity in data.entities() {
            *counts.entry(entity.name()).or_default() += 1;
        }
        let mut counts: Vec<_> = counts.into_iter().collect();
        counts.sort_by_key(|a| std::cmp::Reverse(a.1));
        for (name, count) in counts.iter().take(12) {
            println!("  {count:>8}  {name}");
        }
    }
    Ok(())
}
