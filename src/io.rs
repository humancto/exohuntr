use anyhow::{Context, Result};
use csv::ReaderBuilder;
use std::fs;
use std::path::{Path, PathBuf};

/// A time-series light curve with flux measurements and errors.
#[derive(Debug, Clone)]
pub struct LightCurve {
    pub filename: String,
    pub time: Vec<f64>,
    pub flux: Vec<f64>,
    pub flux_err: Vec<f64>,
}

/// Load a light curve from a CSV file with columns: time, flux, flux_err.
///
/// Skips rows with NaN/Inf values. The flux_err column is optional and
/// defaults to 0.001 if missing.
pub fn load_lightcurve(path: &Path) -> Result<LightCurve> {
    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .flexible(true)
        .trim(csv::Trim::All)
        .from_path(path)?;

    let mut time = Vec::new();
    let mut flux = Vec::new();
    let mut flux_err = Vec::new();

    for result in rdr.records() {
        let record = result?;
        if record.len() < 2 {
            continue;
        }

        let t: f64 = match record.get(0).unwrap_or("").parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let f: f64 = match record.get(1).unwrap_or("").parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let e: f64 = record.get(2).and_then(|s| s.parse().ok()).unwrap_or(0.001);

        if t.is_finite() && f.is_finite() && e.is_finite() {
            time.push(t);
            flux.push(f);
            flux_err.push(e);
        }
    }

    let filename = path
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .to_string();

    Ok(LightCurve {
        filename,
        time,
        flux,
        flux_err,
    })
}

/// Find all CSV files in a directory (non-recursive).
pub fn find_csv_files(dir: &Path) -> Result<Vec<PathBuf>> {
    let mut files = Vec::new();
    for entry in fs::read_dir(dir).context("Cannot read input directory")? {
        let entry = entry?;
        let path = entry.path();
        if path.extension().map_or(false, |e| e == "csv") {
            files.push(path);
        }
    }
    files.sort();
    Ok(files)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::TempDir;

    fn write_csv(dir: &Path, name: &str, content: &str) -> PathBuf {
        let path = dir.join(name);
        let mut f = fs::File::create(&path).unwrap();
        f.write_all(content.as_bytes()).unwrap();
        path
    }

    #[test]
    fn test_load_valid_csv() {
        let dir = TempDir::new().unwrap();
        let path = write_csv(
            dir.path(),
            "test.csv",
            "time,flux,flux_err\n1.0,0.999,0.001\n2.0,0.998,0.002\n3.0,1.001,0.001\n",
        );
        let lc = load_lightcurve(&path).unwrap();
        assert_eq!(lc.time.len(), 3);
        assert_eq!(lc.flux.len(), 3);
        assert_eq!(lc.flux_err.len(), 3);
        assert_eq!(lc.filename, "test.csv");
        assert!((lc.time[0] - 1.0).abs() < 1e-10);
        assert!((lc.flux[1] - 0.998).abs() < 1e-10);
    }

    #[test]
    fn test_load_csv_skips_nan() {
        let dir = TempDir::new().unwrap();
        let path = write_csv(
            dir.path(),
            "nan.csv",
            "time,flux,flux_err\n1.0,0.999,0.001\nNaN,0.998,0.001\n3.0,NaN,0.001\n4.0,1.0,NaN\n5.0,1.001,0.002\n",
        );
        let lc = load_lightcurve(&path).unwrap();
        // NaN rows are skipped, Inf rows too
        assert_eq!(lc.time.len(), 2);
        assert!((lc.time[0] - 1.0).abs() < 1e-10);
        assert!((lc.time[1] - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_load_csv_missing_flux_err() {
        let dir = TempDir::new().unwrap();
        let path = write_csv(
            dir.path(),
            "two_col.csv",
            "time,flux\n1.0,0.999\n2.0,0.998\n",
        );
        let lc = load_lightcurve(&path).unwrap();
        assert_eq!(lc.time.len(), 2);
        // Default flux_err should be 0.001
        assert!((lc.flux_err[0] - 0.001).abs() < 1e-10);
    }

    #[test]
    fn test_load_csv_empty_file() {
        let dir = TempDir::new().unwrap();
        let path = write_csv(dir.path(), "empty.csv", "time,flux,flux_err\n");
        let lc = load_lightcurve(&path).unwrap();
        assert_eq!(lc.time.len(), 0);
    }

    #[test]
    fn test_load_csv_malformed_rows() {
        let dir = TempDir::new().unwrap();
        let path = write_csv(
            dir.path(),
            "bad.csv",
            "time,flux,flux_err\nabc,0.999,0.001\n2.0,xyz,0.001\n3.0,1.0,0.001\n",
        );
        let lc = load_lightcurve(&path).unwrap();
        // Only the last row is valid
        assert_eq!(lc.time.len(), 1);
        assert!((lc.time[0] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_load_csv_infinity_skipped() {
        let dir = TempDir::new().unwrap();
        let path = write_csv(
            dir.path(),
            "inf.csv",
            "time,flux,flux_err\n1.0,inf,0.001\n2.0,1.0,0.001\n",
        );
        let lc = load_lightcurve(&path).unwrap();
        assert_eq!(lc.time.len(), 1);
    }

    #[test]
    fn test_find_csv_files() {
        let dir = TempDir::new().unwrap();
        write_csv(dir.path(), "a.csv", "time,flux\n1,1\n");
        write_csv(dir.path(), "b.csv", "time,flux\n1,1\n");
        write_csv(dir.path(), "c.txt", "not a csv");

        let files = find_csv_files(dir.path()).unwrap();
        assert_eq!(files.len(), 2);
        // Should be sorted
        let names: Vec<String> = files
            .iter()
            .map(|p| p.file_name().unwrap().to_string_lossy().to_string())
            .collect();
        assert_eq!(names, vec!["a.csv", "b.csv"]);
    }

    #[test]
    fn test_find_csv_files_empty_dir() {
        let dir = TempDir::new().unwrap();
        let files = find_csv_files(dir.path()).unwrap();
        assert_eq!(files.len(), 0);
    }
}
