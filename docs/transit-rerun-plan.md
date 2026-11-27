# Transit pipeline rerun plan

Start with a small pilot before regenerating the full results set. Preserve the existing run as historical output and write pilot results to a separate directory. This document defines proposed work; it does not change the pipeline or report a new scientific result.

## 1. Correct the measurements

- **Flux depth:** In [the BLS search](../src/bls.rs), keep the robustly scaled flux for detection. Refit the selected window on the original relative-flux scale (or convert explicitly using the same baseline and scale) before [exporting depth in ppm and estimating radius ratio](../src/main.rs). Track the detection score separately. Test shallow injections with different noise levels so the physical depth does not change simply because the noise scale changes.
- **Midpoint and coverage:** Convert the winning box start to an initial midpoint, then fit the midpoint and propagate its uncertainty. Record the source time system and export an explicitly labeled time standard (BJD_TDB for fitted midpoints). Count distinct events with usable in-transit cadences rather than elapsed periods. Include a light curve with a gap and a transit window that wraps across the phase boundary.
- **Catalog association:** Keep all catalog entries for each TIC. Match each recovered signal against period, epoch when available, and integer aliases; leave ambiguous matches unresolved. Group repeated labels for the same signal without discarding truly distinct periods.
- **Product selection:** Replace the first-result download in [the sector downloader](../python/download_new_sectors.py) with explicit sector and product selection. Preserve product IDs, cadence, flux column and quality masks in a manifest; normalize sectors separately before combining them. Keep known-TOI recovery and any blind-search cohort separate.

## 2. Run and check a bounded pilot

| Case | Purpose | Pass condition |
| --- | --- | --- |
| Synthetic light curves with two injected periods, several depths and realistic gaps | Test units, timing, coverage and the second-pass search | Recover each injected period once; recover depths and midpoints within tolerances declared before running; do not count unobserved events. |
| TOI-125 b/c | Known multi-signal association control | Recover and match the two known periods (about 4.65 and 9.15 days) separately when coverage supports it. Do not associate the 4.65-day peak with TOI-125.04. See the [NASA system overview](https://exoplanetarchive.ipac.caltech.edu/overview/TOI-125). |
| TOI-210 b | Known-signal timing and depth control | Recover the approximately 9.01-day signal and report fitted values with uncertainty, alongside the [external fit](https://arxiv.org/abs/2609.35747). |
| TOI-155.01 | First targeted remeasurement | Use the exact available TESS products to refit period and midpoint. Compare with the [published transit fit](https://arxiv.org/html/2609.37517v2); report uncertainties, usable events, and whether the longer baseline adds information. |

For the multi-signal test, mask or model the first signal before a second search. Inspect the sampling window and integer aliases. A residual that fails to repeat in separate sectors, where coverage permits, stays an unresolved diagnostic. Use pixel/nearby-source checks before interpreting any repeatable residual.

## 3. Record results before expanding

- Save the commit, settings, catalog snapshot dates, exact [MAST](https://archive.stsci.edu/missions-and-data/tess) products, event counts, uncertainties, phased plots, and decisions. Keep a short before/after table for the same targets.
- Run the Rust and Python test suites plus the synthetic and known-signal controls. If any control fails, stop and retain its diagnostics; do not replace the published tables.
- Once the pilot passes, rerun the original cohort with the same input manifest for comparison. Report unique hosts and signals separately. If trying a blind cohort later, select it independently of TOI status and cross-check any residual against the [TESS threshold-crossing-event tables](https://archive.stsci.edu/tess/bulk_downloads/bulk_downloads_tce.html) and the current catalog.
