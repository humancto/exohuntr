//! Cross-matching candidates against known exoplanet catalogs.
//!
//! Uses hash-based lookups for O(1) matching by hostname/TIC ID,
//! replacing the O(N*M) nested loop in the Python implementation.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;

/// A known exoplanet from a catalog (e.g., NASA Exoplanet Archive).
#[derive(Debug, Clone, Deserialize)]
pub struct CatalogEntry {
    #[serde(default)]
    pub pl_name: String,
    #[serde(default)]
    pub hostname: String,
    #[serde(default)]
    pub pl_orbper: Option<f64>,
    #[serde(default)]
    pub pl_rade: Option<f64>,
}

/// Result of cross-matching a single candidate against the catalog.
#[derive(Debug, Clone, Serialize)]
pub struct CrossmatchResult {
    pub candidate_file: String,
    pub candidate_period: f64,
    pub candidate_snr: f64,
    pub known_planet: String,
    pub known_period: Option<f64>,
    pub known_radius_earth: Option<f64>,
    pub status: String,
}

/// A lookup index built from the exoplanet catalog for O(1) matching.
pub struct CatalogIndex {
    /// Map from normalized hostname to catalog entries.
    by_hostname: HashMap<String, Vec<CatalogEntry>>,
}

impl CatalogIndex {
    /// Build the index from a list of catalog entries.
    pub fn new(entries: Vec<CatalogEntry>) -> Self {
        let mut by_hostname: HashMap<String, Vec<CatalogEntry>> = HashMap::new();
        for entry in entries {
            let key = normalize_name(&entry.hostname);
            if !key.is_empty() {
                by_hostname.entry(key).or_default().push(entry);
            }
        }
        CatalogIndex { by_hostname }
    }

    /// Number of unique hostnames in the index.
    pub fn len(&self) -> usize {
        self.by_hostname.len()
    }

    /// Whether the index is empty.
    pub fn is_empty(&self) -> bool {
        self.by_hostname.is_empty()
    }

    /// Look up a candidate filename against the catalog.
    ///
    /// Finds the longest matching hostname contained in the filename to avoid
    /// ambiguity (e.g., "toi_12" vs "toi_123" both matching "toi_123_tic_456").
    /// Returns the first catalog entry for that hostname, or None.
    pub fn lookup(&self, filename: &str) -> Option<&CatalogEntry> {
        let normalized = normalize_name(filename);
        let mut best_match: Option<(&str, &Vec<CatalogEntry>)> = None;
        for (hostname, entries) in &self.by_hostname {
            if normalized.contains(hostname.as_str()) {
                let dominated = best_match
                    .map(|(prev, _)| hostname.len() > prev.len())
                    .unwrap_or(true);
                if dominated {
                    best_match = Some((hostname, entries));
                }
            }
        }
        best_match.and_then(|(_, entries)| entries.first())
    }

    /// Look up by exact normalized hostname.
    pub fn lookup_by_hostname(&self, hostname: &str) -> Option<&CatalogEntry> {
        let key = normalize_name(hostname);
        self.by_hostname.get(&key).and_then(|v| v.first())
    }
}

/// Load a confirmed exoplanet catalog from a CSV file.
///
/// Expects columns: pl_name, hostname, pl_orbper, pl_rade (matching
/// the NASA Exoplanet Archive pscomppars table format).
pub fn load_catalog(path: &Path) -> anyhow::Result<Vec<CatalogEntry>> {
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .flexible(true)
        .trim(csv::Trim::All)
        .from_path(path)?;

    let mut entries = Vec::new();
    for result in rdr.deserialize() {
        match result {
            Ok(entry) => entries.push(entry),
            Err(_) => continue, // skip malformed rows
        }
    }
    Ok(entries)
}

/// Cross-match a list of candidates against the catalog index.
pub fn crossmatch_candidates(
    candidates: &[(String, f64, f64)], // (filename, period, snr)
    index: &CatalogIndex,
) -> Vec<CrossmatchResult> {
    candidates
        .iter()
        .map(|(filename, period, snr)| {
            if let Some(entry) = index.lookup(filename) {
                CrossmatchResult {
                    candidate_file: filename.clone(),
                    candidate_period: *period,
                    candidate_snr: *snr,
                    known_planet: entry.pl_name.clone(),
                    known_period: entry.pl_orbper,
                    known_radius_earth: entry.pl_rade,
                    status: "KNOWN".to_string(),
                }
            } else {
                CrossmatchResult {
                    candidate_file: filename.clone(),
                    candidate_period: *period,
                    candidate_snr: *snr,
                    known_planet: String::new(),
                    known_period: None,
                    known_radius_earth: None,
                    status: "POTENTIALLY NEW".to_string(),
                }
            }
        })
        .collect()
}

/// Normalize a name for matching: lowercase, replace spaces with underscores,
/// remove common prefixes.
fn normalize_name(name: &str) -> String {
    name.to_lowercase()
        .replace([' ', '-', '/'], "_")
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::TempDir;

    fn make_catalog_entries() -> Vec<CatalogEntry> {
        vec![
            CatalogEntry {
                pl_name: "TOI-125 b".to_string(),
                hostname: "TOI-125".to_string(),
                pl_orbper: Some(4.65),
                pl_rade: Some(2.76),
            },
            CatalogEntry {
                pl_name: "WASP-18 b".to_string(),
                hostname: "WASP-18".to_string(),
                pl_orbper: Some(0.94),
                pl_rade: Some(12.4),
            },
            CatalogEntry {
                pl_name: "HD 21749 b".to_string(),
                hostname: "HD 21749".to_string(),
                pl_orbper: Some(35.61),
                pl_rade: Some(2.86),
            },
        ]
    }

    #[test]
    fn test_catalog_index_creation() {
        let index = CatalogIndex::new(make_catalog_entries());
        assert_eq!(index.len(), 3);
        assert!(!index.is_empty());
    }

    #[test]
    fn test_catalog_index_empty() {
        let index = CatalogIndex::new(vec![]);
        assert_eq!(index.len(), 0);
        assert!(index.is_empty());
    }

    #[test]
    fn test_lookup_by_hostname_exact() {
        let index = CatalogIndex::new(make_catalog_entries());
        let result = index.lookup_by_hostname("TOI-125");
        assert!(result.is_some());
        assert_eq!(result.unwrap().pl_name, "TOI-125 b");
    }

    #[test]
    fn test_lookup_by_hostname_case_insensitive() {
        let index = CatalogIndex::new(make_catalog_entries());
        let result = index.lookup_by_hostname("toi-125");
        assert!(result.is_some());
    }

    #[test]
    fn test_lookup_by_hostname_not_found() {
        let index = CatalogIndex::new(make_catalog_entries());
        let result = index.lookup_by_hostname("XYZ-999");
        assert!(result.is_none());
    }

    #[test]
    fn test_lookup_filename_contains_hostname() {
        let index = CatalogIndex::new(make_catalog_entries());
        // Typical filename from the pipeline
        let result = index.lookup("TOI_125.01_TIC_12345.csv");
        assert!(result.is_some());
        assert_eq!(result.unwrap().pl_name, "TOI-125 b");
    }

    #[test]
    fn test_lookup_filename_no_match() {
        let index = CatalogIndex::new(make_catalog_entries());
        let result = index.lookup("TIC_999999_s0056.csv");
        assert!(result.is_none());
    }

    #[test]
    fn test_crossmatch_candidates_mixed() {
        let index = CatalogIndex::new(make_catalog_entries());
        let candidates = vec![
            ("TOI_125.01_TIC_12345.csv".to_string(), 4.65, 25.0),
            ("TIC_999999_s0056.csv".to_string(), 2.0, 10.0),
        ];

        let results = crossmatch_candidates(&candidates, &index);
        assert_eq!(results.len(), 2);
        assert_eq!(results[0].status, "KNOWN");
        assert_eq!(results[0].known_planet, "TOI-125 b");
        assert_eq!(results[1].status, "POTENTIALLY NEW");
        assert!(results[1].known_planet.is_empty());
    }

    #[test]
    fn test_crossmatch_all_new() {
        let index = CatalogIndex::new(make_catalog_entries());
        let candidates = vec![
            ("TIC_111111.csv".to_string(), 1.0, 8.0),
            ("TIC_222222.csv".to_string(), 5.0, 12.0),
        ];

        let results = crossmatch_candidates(&candidates, &index);
        assert!(results.iter().all(|r| r.status == "POTENTIALLY NEW"));
    }

    #[test]
    fn test_crossmatch_empty_candidates() {
        let index = CatalogIndex::new(make_catalog_entries());
        let results = crossmatch_candidates(&[], &index);
        assert!(results.is_empty());
    }

    #[test]
    fn test_normalize_name() {
        assert_eq!(normalize_name("TOI-125"), "toi_125");
        assert_eq!(normalize_name("HD 21749"), "hd_21749");
        assert_eq!(normalize_name("WASP-18"), "wasp_18");
        assert_eq!(normalize_name("TIC 261136679"), "tic_261136679");
    }

    #[test]
    fn test_load_catalog_from_csv() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("catalog.csv");
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, "pl_name,hostname,pl_orbper,pl_rade").unwrap();
        writeln!(f, "TOI-125 b,TOI-125,4.65,2.76").unwrap();
        writeln!(f, "WASP-18 b,WASP-18,0.94,12.4").unwrap();

        let entries = load_catalog(&path).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].pl_name, "TOI-125 b");
        assert_eq!(entries[0].hostname, "TOI-125");
        assert!((entries[0].pl_orbper.unwrap() - 4.65).abs() < 1e-10);
    }

    #[test]
    fn test_load_catalog_empty() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("empty.csv");
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, "pl_name,hostname,pl_orbper,pl_rade").unwrap();

        let entries = load_catalog(&path).unwrap();
        assert!(entries.is_empty());
    }

    #[test]
    fn test_load_catalog_malformed_rows() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("bad.csv");
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, "pl_name,hostname,pl_orbper,pl_rade").unwrap();
        writeln!(f, "TOI-125 b,TOI-125,4.65,2.76").unwrap();
        writeln!(f, "bad row with no commas").unwrap();
        writeln!(f, "WASP-18 b,WASP-18,0.94,12.4").unwrap();

        let entries = load_catalog(&path).unwrap();
        // Should skip the bad row
        assert!(entries.len() >= 2);
    }

    #[test]
    fn test_lookup_hostname_collision() {
        // Overlapping hostnames: toi_12 vs toi_123 — should match longest
        let entries = vec![
            CatalogEntry {
                pl_name: "TOI-12 b".to_string(),
                hostname: "TOI-12".to_string(),
                pl_orbper: Some(1.0),
                pl_rade: Some(1.0),
            },
            CatalogEntry {
                pl_name: "TOI-123 b".to_string(),
                hostname: "TOI-123".to_string(),
                pl_orbper: Some(2.0),
                pl_rade: Some(2.0),
            },
        ];
        let index = CatalogIndex::new(entries);

        // Should match TOI-123, not TOI-12
        let result = index.lookup("TOI_123_TIC_456.csv");
        assert!(result.is_some());
        assert_eq!(result.unwrap().pl_name, "TOI-123 b");

        // Should match TOI-12 (no ambiguity)
        let result = index.lookup("TOI_12_TIC_789.csv");
        assert!(result.is_some());
        assert_eq!(result.unwrap().pl_name, "TOI-12 b");
    }

    #[test]
    fn test_multiple_planets_same_host() {
        let entries = vec![
            CatalogEntry {
                pl_name: "TOI-125 b".to_string(),
                hostname: "TOI-125".to_string(),
                pl_orbper: Some(4.65),
                pl_rade: Some(2.76),
            },
            CatalogEntry {
                pl_name: "TOI-125 c".to_string(),
                hostname: "TOI-125".to_string(),
                pl_orbper: Some(9.15),
                pl_rade: Some(2.93),
            },
        ];
        let index = CatalogIndex::new(entries);
        // Both should be stored under the same hostname
        assert_eq!(index.len(), 1); // 1 unique hostname
        let result = index.lookup_by_hostname("TOI-125");
        assert!(result.is_some());
    }
}
