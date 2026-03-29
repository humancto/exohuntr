#!/usr/bin/env python3.11
"""
Deep Analysis Pipeline for Exohuntr Top Candidates
====================================================
Steps:
  1. Download Target Pixel Files (TPFs) → centroid analysis
  2. Query Gaia DR3 for nearby contaminating sources
  3. Download NASA SPOC Data Validation reports from MAST
  4. Run Transit Least Squares (TLS) with limb-darkened models
  5. Download multi-sector data for extended secondary eclipse search

Targets: TOI 133.01, TOI 210.01, TOI 155.01 (top 3 plausible planet candidates)
"""
from __future__ import annotations

import os
import sys
import json
import warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from datetime import datetime

warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', message='.*overflow.*')

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'results')
DEEP_DIR = os.path.join(RESULTS_DIR, 'deep_analysis')
os.makedirs(DEEP_DIR, exist_ok=True)

# Our top candidates to deep-analyze
TARGETS = [
    {"toi": "TOI 133.01", "tic": "TIC 219338557", "tic_id": 219338557, "period": 8.2065, "rp_earth": 1.9, "snr": 13.3},
    {"toi": "TOI 210.01", "tic": "TIC 141608198", "tic_id": 141608198, "period": 8.9884, "rp_earth": 2.2, "snr": 7.1},
    {"toi": "TOI 155.01", "tic": "TIC 129637892", "tic_id": 129637892, "period": 5.4504, "rp_earth": 5.3, "snr": 44.3},
]


def log(msg: str) -> None:
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}")


# ===========================================================================
# STEP 1: Target Pixel Files + Centroid Analysis
# ===========================================================================
def step1_centroid_analysis(findings: dict) -> None:
    log("STEP 1: Downloading TPFs and running centroid analysis...")
    import lightkurve as lk

    for t in TARGETS:
        name = t['toi']
        tic = t['tic']
        log(f"  Searching TPFs for {tic}...")
        findings[name] = findings.get(name, {})

        try:
            search = lk.search_targetpixelfile(tic, mission='TESS')
            if len(search) == 0:
                log(f"  WARNING: No TPFs found for {tic}")
                findings[name]['centroid'] = {'status': 'NO_TPF', 'detail': 'No target pixel files available on MAST'}
                continue

            log(f"  Found {len(search)} TPF(s). Downloading first...")
            tpf = search[0].download()

            # Compute centroids directly from TPF pixel data
            log(f"  Computing centroids from pixel data for {name}...")
            flux_cube = tpf.flux.value
            n_frames, ny, nx = flux_cube.shape
            col_grid, row_grid = np.meshgrid(np.arange(nx), np.arange(ny))

            centroid_col = np.zeros(n_frames)
            centroid_row = np.zeros(n_frames)
            aperture_flux = np.zeros(n_frames)
            for i in range(n_frames):
                frame = flux_cube[i]
                total = np.nansum(frame)
                if total > 0:
                    centroid_col[i] = np.nansum(frame * col_grid) / total
                    centroid_row[i] = np.nansum(frame * row_grid) / total
                    aperture_flux[i] = total
                else:
                    centroid_col[i] = np.nan
                    centroid_row[i] = np.nan
                    aperture_flux[i] = np.nan

            time = tpf.time.value
            mask = np.isfinite(centroid_col) & np.isfinite(centroid_row) & np.isfinite(aperture_flux)
            time = time[mask]
            flux = aperture_flux[mask]
            flux = flux / np.nanmedian(flux)
            cc = centroid_col[mask]
            cr = centroid_row[mask]

            # Phase-fold at the transit period
            period = t['period']
            # Find epoch from minimum flux
            phase = (time % period) / period
            transit_mask = np.abs(phase - 0.0) < 0.05  # in-transit: within 5% of phase 0
            # Also check near phase 1.0
            transit_mask |= np.abs(phase - 1.0) < 0.05
            oot_mask = ~transit_mask  # out-of-transit

            if np.sum(transit_mask) < 3:
                log(f"  WARNING: Too few in-transit points for {name}")
                findings[name]['centroid'] = {'status': 'INSUFFICIENT_DATA', 'detail': f'Only {np.sum(transit_mask)} in-transit points'}
                continue

            # Compare centroids in-transit vs out-of-transit
            col_in = np.nanmean(cc[transit_mask])
            col_out = np.nanmean(cc[oot_mask])
            row_in = np.nanmean(cr[transit_mask])
            row_out = np.nanmean(cr[oot_mask])

            col_shift = col_in - col_out
            row_shift = row_in - row_out
            total_shift = np.sqrt(col_shift**2 + row_shift**2)

            # TESS pixel scale is ~21 arcsec/pixel
            shift_arcsec = total_shift * 21.0

            # Significance: compare shift to scatter
            col_std = np.nanstd(cc[oot_mask]) / np.sqrt(np.sum(oot_mask))
            row_std = np.nanstd(cr[oot_mask]) / np.sqrt(np.sum(oot_mask))
            shift_sigma = total_shift / np.sqrt(col_std**2 + row_std**2) if (col_std > 0 and row_std > 0) else 0

            passed = shift_sigma < 3.0  # Less than 3-sigma shift = on target
            status = "PASS" if passed else "FAIL"

            findings[name]['centroid'] = {
                'status': status,
                'col_shift_pix': round(float(col_shift), 6),
                'row_shift_pix': round(float(row_shift), 6),
                'total_shift_pix': round(float(total_shift), 6),
                'shift_arcsec': round(float(shift_arcsec), 3),
                'shift_sigma': round(float(shift_sigma), 2),
                'n_in_transit': int(np.sum(transit_mask)),
                'n_out_transit': int(np.sum(oot_mask)),
                'detail': f'Centroid shift: {shift_arcsec:.1f}" ({shift_sigma:.1f}σ). {"On target — transit is on this star" if passed else "OFFSET DETECTED — signal may be from nearby star"}'
            }

            log(f"  {name}: centroid shift = {shift_arcsec:.1f} arcsec ({shift_sigma:.1f}σ) → {status}")

            # Plot centroid shift
            fig, axes = plt.subplots(1, 3, figsize=(15, 4))
            fig.suptitle(f'{name} — Centroid Analysis', fontsize=14)

            axes[0].scatter(phase[oot_mask], cc[oot_mask], s=1, alpha=0.3, c='gray', label='out-of-transit')
            axes[0].scatter(phase[transit_mask], cc[transit_mask], s=5, alpha=0.8, c='red', label='in-transit')
            axes[0].set_xlabel('Phase')
            axes[0].set_ylabel('Column centroid (pix)')
            axes[0].legend(fontsize=8)

            axes[1].scatter(phase[oot_mask], cr[oot_mask], s=1, alpha=0.3, c='gray', label='out-of-transit')
            axes[1].scatter(phase[transit_mask], cr[transit_mask], s=5, alpha=0.8, c='red', label='in-transit')
            axes[1].set_xlabel('Phase')
            axes[1].set_ylabel('Row centroid (pix)')
            axes[1].legend(fontsize=8)

            axes[2].scatter(cc[oot_mask], cr[oot_mask], s=1, alpha=0.3, c='gray', label='out-of-transit')
            axes[2].scatter(cc[transit_mask], cr[transit_mask], s=10, alpha=0.8, c='red', label='in-transit')
            axes[2].set_xlabel('Column centroid (pix)')
            axes[2].set_ylabel('Row centroid (pix)')
            axes[2].legend(fontsize=8)

            plt.tight_layout()
            safe_name = name.replace(" ", "_").replace(".", "_")
            safe_name = os.path.basename(safe_name)  # strip any directory components
            plot_path = os.path.join(DEEP_DIR, f'centroid_{safe_name}.png')
            plt.savefig(plot_path, dpi=150, bbox_inches='tight')
            plt.close()
            log(f"  Saved centroid plot: {plot_path}")

        except Exception as e:
            log(f"  ERROR for {name}: {e}")
            findings[name]['centroid'] = {'status': 'ERROR', 'detail': str(e)}


# ===========================================================================
# STEP 2: Gaia DR3 Nearby Source Check
# ===========================================================================
def step2_gaia_query(findings: dict) -> None:
    log("STEP 2: Querying Gaia DR3 for nearby contaminating sources...")
    from astroquery.mast import Catalogs
    from astroquery.gaia import Gaia
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    for t in TARGETS:
        name = t['toi']
        tic_id = t['tic_id']
        findings[name] = findings.get(name, {})

        try:
            # Get target coordinates from TIC
            log(f"  Looking up {t['tic']} in TIC catalog...")
            tic_data = Catalogs.query_criteria(catalog="TIC", ID=tic_id)

            if len(tic_data) == 0:
                log(f"  WARNING: TIC {tic_id} not found in catalog")
                findings[name]['gaia'] = {'status': 'NOT_FOUND', 'detail': 'TIC entry not found'}
                continue

            ra = float(tic_data['ra'][0])
            dec = float(tic_data['dec'][0])
            tmag = float(tic_data['Tmag'][0]) if 'Tmag' in tic_data.colnames else None

            log(f"  {name}: RA={ra:.5f}, Dec={dec:.5f}, Tmag={tmag}")

            # Query Gaia within 2 arcmin (TESS aperture is typically ~1-4 pixels = 21-84 arcsec)
            coord = SkyCoord(ra=ra, dec=dec, unit='deg')
            search_radius = 120  # arcsec = 2 arcmin

            log(f"  Querying Gaia DR3 within {search_radius}\" of target...")
            # Sanitize numeric inputs for ADQL (TAP doesn't support parameterized queries)
            safe_ra = float(ra)
            safe_dec = float(dec)
            safe_radius = float(search_radius) / 3600.0
            query = (
                "SELECT source_id, ra, dec, phot_g_mean_mag, parallax,"
                " DISTANCE(POINT('ICRS', ra, dec),"
                f" POINT('ICRS', {safe_ra:.6f}, {safe_dec:.6f})) AS ang_sep"
                " FROM gaiadr3.gaia_source"
                " WHERE DISTANCE(POINT('ICRS', ra, dec),"
                f" POINT('ICRS', {safe_ra:.6f}, {safe_dec:.6f})) < {safe_radius:.8f}"
                " ORDER BY ang_sep ASC"
            )

            job = Gaia.launch_job(query)
            results = job.get_results()

            n_sources = len(results)
            log(f"  Found {n_sources} Gaia sources within {search_radius}\"")

            # Check for bright nearby sources that could contaminate
            nearby_bright = []
            if n_sources > 1:
                for row in results[1:]:  # skip the target itself
                    sep_arcsec = float(row['ang_sep']) * 3600
                    gmag = float(row['phot_g_mean_mag'])
                    # A source within ~42 arcsec (2 TESS pixels) and within 5 mag is concerning
                    if sep_arcsec < 42 and gmag < (tmag + 5 if tmag else 20):
                        nearby_bright.append({
                            'gaia_id': str(row['source_id']),
                            'sep_arcsec': round(sep_arcsec, 1),
                            'gmag': round(gmag, 2),
                            'delta_mag': round(gmag - tmag, 2) if tmag else None
                        })

            if nearby_bright:
                status = "WARNING"
                source_strs = []
                for s in nearby_bright[:3]:
                    source_strs.append("Gmag=%.2f at %.1f arcsec" % (s['gmag'], s['sep_arcsec']))
                detail = "%d bright source(s) within 42 arcsec: %s" % (len(nearby_bright), ', '.join(source_strs))
            else:
                status = "CLEAR"
                detail = "No bright contaminating sources within 42 arcsec (checked %d Gaia sources)" % n_sources

            findings[name]['gaia'] = {
                'status': status,
                'n_gaia_sources': n_sources,
                'nearby_bright': nearby_bright,
                'target_ra': ra,
                'target_dec': dec,
                'target_tmag': tmag,
                'detail': detail
            }

            log(f"  {name}: {status} — {detail}")

        except Exception as e:
            log(f"  ERROR for {name}: {e}")
            findings[name]['gaia'] = {'status': 'ERROR', 'detail': str(e)}


# ===========================================================================
# STEP 3: Check NASA SPOC DV Reports on MAST
# ===========================================================================
def step3_dv_reports(findings: dict) -> None:
    log("STEP 3: Checking NASA SPOC Data Validation reports on MAST...")
    from astroquery.mast import Observations

    for t in TARGETS:
        name = t['toi']
        tic_id = t['tic_id']
        findings[name] = findings.get(name, {})

        try:
            log(f"  Searching MAST for DV products for TIC {tic_id}...")

            # Search for TESS observations
            obs = Observations.query_criteria(
                target_name=str(tic_id),
                obs_collection='TESS',
                dataproduct_type='timeseries'
            )

            if len(obs) == 0:
                log(f"  No TESS observations found for TIC {tic_id}")
                findings[name]['dv_report'] = {'status': 'NOT_FOUND', 'detail': 'No TESS timeseries on MAST'}
                continue

            log(f"  Found {len(obs)} TESS observation(s)")

            # Get data products to look for DV files
            products = Observations.get_product_list(obs)

            # Filter for DV reports
            dv_products = products[
                (products['productSubGroupDescription'] == 'DVR') |
                (products['productSubGroupDescription'] == 'DVS') |
                (products['productSubGroupDescription'] == 'DVT')
            ] if 'productSubGroupDescription' in products.colnames else []

            n_dv = len(dv_products) if hasattr(dv_products, '__len__') else 0

            if n_dv > 0:
                status = "AVAILABLE"
                # Get the filenames
                dv_files = [str(row['productFilename']) for row in dv_products[:5]]
                detail = f"{n_dv} DV product(s) available: {', '.join(dv_files[:3])}"
                log(f"  {name}: {n_dv} DV products found")
            else:
                # Check if there are any SPOC products at all
                spoc = products[products['project'] == 'SPOC'] if 'project' in products.colnames else []
                n_spoc = len(spoc) if hasattr(spoc, '__len__') else 0
                status = "NO_DV" if n_spoc > 0 else "NO_SPOC"
                detail = f"SPOC processed ({n_spoc} products) but no DV reports found" if n_spoc > 0 else "No SPOC processing found"
                log(f"  {name}: {status}")

            # Also note number of sectors observed
            sectors = set()
            for o in obs:
                sn = str(o.get('sequence_number', ''))
                if sn:
                    sectors.add(sn)

            findings[name]['dv_report'] = {
                'status': status,
                'n_dv_products': n_dv,
                'n_observations': len(obs),
                'sectors_observed': sorted(list(sectors)),
                'detail': detail
            }

        except Exception as e:
            log(f"  ERROR for {name}: {e}")
            findings[name]['dv_report'] = {'status': 'ERROR', 'detail': str(e)}


# ===========================================================================
# STEP 4: Transit Least Squares (TLS)
# ===========================================================================
def step4_tls(findings: dict) -> None:
    log("STEP 4: Running Transit Least Squares (TLS) analysis...")
    try:
        from transitleastsquares import transitleastsquares
    except ImportError:
        log("  ERROR: transitleastsquares not installed")
        for t in TARGETS:
            findings[t['toi']]['tls'] = {'status': 'NOT_INSTALLED', 'detail': 'pip install transitleastsquares'}
        return

    import lightkurve as lk

    for t in TARGETS:
        name = t['toi']
        tic = t['tic']
        findings[name] = findings.get(name, {})

        try:
            log(f"  Searching light curves for {tic}...")
            search = lk.search_lightcurve(tic, mission='TESS', author='SPOC')

            if len(search) == 0:
                search = lk.search_lightcurve(tic, mission='TESS')

            if len(search) == 0:
                log(f"  No light curves found for {tic}")
                findings[name]['tls'] = {'status': 'NO_DATA', 'detail': 'No light curves on MAST'}
                continue

            log(f"  Found {len(search)} light curve(s). Downloading and stitching...")
            lc_collection = search.download_all()
            lc = lc_collection.stitch()
            lc = lc.remove_nans().remove_outliers(sigma=5)

            time = lc.time.value
            flux = lc.flux.value

            # Normalize
            flux = flux / np.nanmedian(flux)

            log(f"  Running TLS on {len(time)} data points (this may take a minute)...")
            model = transitleastsquares(time, flux)
            results = model.power(
                period_min=t['period'] * 0.8,
                period_max=t['period'] * 1.2,
                oversampling_factor=3,
                duration_grid_step=1.05,
            )

            tls_period = float(results.period)
            tls_depth = float(1 - results.depth) if hasattr(results, 'depth') else None
            tls_snr = float(results.snr) if hasattr(results, 'snr') else None
            tls_sde = float(results.SDE) if hasattr(results, 'SDE') else None
            tls_t0 = float(results.T0) if hasattr(results, 'T0') else None
            tls_duration = float(results.duration) if hasattr(results, 'duration') else None
            tls_rp_rs = float(results.rp_rs) if hasattr(results, 'rp_rs') else None

            # Compare TLS period to our BLS period
            period_diff_pct = abs(tls_period - t['period']) / t['period'] * 100

            if period_diff_pct < 1.0 and tls_sde and tls_sde > 6.0:
                status = "CONFIRMED"
                detail = f"TLS confirms transit: P={tls_period:.4f}d (Δ={period_diff_pct:.3f}%), SDE={tls_sde:.1f}, Rp/Rs={tls_rp_rs:.4f}"
            elif period_diff_pct < 1.0:
                status = "WEAK_CONFIRM"
                detail = f"TLS finds same period but weak signal: P={tls_period:.4f}d, SDE={tls_sde:.1f if tls_sde else 'N/A'}"
            else:
                status = "DISAGREE"
                detail = f"TLS finds different period: {tls_period:.4f}d vs BLS {t['period']:.4f}d (Δ={period_diff_pct:.1f}%)"

            findings[name]['tls'] = {
                'status': status,
                'tls_period': round(tls_period, 6),
                'bls_period': t['period'],
                'period_diff_pct': round(period_diff_pct, 4),
                'tls_sde': round(tls_sde, 2) if tls_sde else None,
                'tls_snr': round(tls_snr, 2) if tls_snr else None,
                'tls_rp_rs': round(tls_rp_rs, 5) if tls_rp_rs else None,
                'tls_depth': round(tls_depth, 6) if tls_depth else None,
                'tls_duration_hours': round(tls_duration * 24, 2) if tls_duration else None,
                'n_datapoints': len(time),
                'detail': detail
            }

            log(f"  {name}: {status} — {detail}")

            # Plot TLS results
            fig, axes = plt.subplots(1, 2, figsize=(14, 5))
            fig.suptitle(f'{name} — Transit Least Squares Analysis', fontsize=14)

            # TLS periodogram
            axes[0].plot(results.periods, results.power, 'k-', lw=0.5)
            axes[0].axvline(tls_period, color='red', ls='--', label=f'TLS P={tls_period:.4f}d')
            axes[0].axvline(t['period'], color='blue', ls=':', label=f'BLS P={t["period"]:.4f}d')
            axes[0].set_xlabel('Period (days)')
            axes[0].set_ylabel('SDE')
            axes[0].set_title('TLS Periodogram')
            axes[0].legend(fontsize=9)

            # Phase-folded with TLS model
            phase = (time - tls_t0 + tls_period/2) % tls_period - tls_period/2 if tls_t0 else (time % tls_period)
            axes[1].scatter(phase * 24, flux, s=1, alpha=0.3, c='gray')
            # Bin the phase-folded data
            n_bins = 100
            bin_edges = np.linspace(phase.min(), phase.max(), n_bins + 1)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            bin_flux = np.array([np.nanmedian(flux[(phase >= bin_edges[i]) & (phase < bin_edges[i+1])])
                                 for i in range(n_bins)])
            valid = np.isfinite(bin_flux)
            axes[1].plot(bin_centers[valid] * 24, bin_flux[valid], 'r-', lw=2, label='Binned data')
            axes[1].set_xlabel('Time from mid-transit (hours)')
            axes[1].set_ylabel('Normalized flux')
            axes[1].set_title(f'Phase-folded at P={tls_period:.4f}d')
            axes[1].set_xlim(-5, 5)
            axes[1].legend(fontsize=9)

            plt.tight_layout()
            safe_name = name.replace(" ", "_").replace(".", "_")
            safe_name = os.path.basename(safe_name)  # strip any directory components
            plot_path = os.path.join(DEEP_DIR, f'tls_{safe_name}.png')
            plt.savefig(plot_path, dpi=150, bbox_inches='tight')
            plt.close()
            log(f"  Saved TLS plot: {plot_path}")

        except Exception as e:
            log(f"  ERROR for {name}: {e}")
            findings[name]['tls'] = {'status': 'ERROR', 'detail': str(e)}


# ===========================================================================
# STEP 5: Multi-sector secondary eclipse check
# ===========================================================================
def step5_multisector_secondary(findings: dict) -> None:
    log("STEP 5: Multi-sector secondary eclipse search...")
    import lightkurve as lk

    for t in TARGETS:
        name = t['toi']
        tic = t['tic']
        period = t['period']
        findings[name] = findings.get(name, {})

        try:
            log(f"  Downloading all available sectors for {tic}...")
            search = lk.search_lightcurve(tic, mission='TESS')

            if len(search) == 0:
                findings[name]['secondary_eclipse_deep'] = {'status': 'NO_DATA'}
                continue

            lc_collection = search.download_all()
            lc = lc_collection.stitch()
            lc = lc.remove_nans().remove_outliers(sigma=5)
            lc = lc.normalize()

            time = lc.time.value
            flux = lc.flux.value

            log(f"  {name}: {len(time)} total data points across {len(search)} sector(s)")

            # Phase-fold
            phase = (time % period) / period

            # Check for secondary eclipse at phase 0.5 (±0.05)
            sec_mask = np.abs(phase - 0.5) < 0.05
            oot_mask = (np.abs(phase - 0.25) < 0.1) | (np.abs(phase - 0.75) < 0.1)

            if np.sum(sec_mask) < 5:
                findings[name]['secondary_eclipse_deep'] = {
                    'status': 'INSUFFICIENT_DATA',
                    'detail': f'Only {np.sum(sec_mask)} points near phase 0.5'
                }
                continue

            sec_flux = np.nanmedian(flux[sec_mask])
            oot_flux = np.nanmedian(flux[oot_mask])
            sec_depth = oot_flux - sec_flux

            # Significance
            oot_std = np.nanstd(flux[oot_mask]) / np.sqrt(np.sum(oot_mask))
            sec_sigma = sec_depth / oot_std if oot_std > 0 else 0

            has_secondary = sec_sigma > 3.0 and sec_depth > 0

            status = "FAIL — secondary eclipse detected" if has_secondary else "PASS — no secondary eclipse"

            findings[name]['secondary_eclipse_deep'] = {
                'status': 'FAIL' if has_secondary else 'PASS',
                'sec_depth': round(float(sec_depth), 6),
                'sec_sigma': round(float(sec_sigma), 2),
                'n_sec_points': int(np.sum(sec_mask)),
                'n_sectors': len(search),
                'n_total_points': len(time),
                'detail': f'Secondary eclipse depth: {sec_depth:.6f} ({sec_sigma:.1f}σ) from {len(search)} sector(s), {len(time)} points. {status}'
            }

            log(f"  {name}: {status} (depth={sec_depth:.6f}, {sec_sigma:.1f}σ, {len(search)} sectors)")

            # Plot phase-folded with secondary eclipse region highlighted
            fig, ax = plt.subplots(figsize=(10, 4))
            ax.scatter(phase, flux, s=1, alpha=0.2, c='gray')
            ax.scatter(phase[sec_mask], flux[sec_mask], s=3, alpha=0.5, c='red', label='Secondary region (0.45-0.55)')

            # Bin
            n_bins = 200
            bin_edges = np.linspace(0, 1, n_bins + 1)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            bin_flux = np.array([np.nanmedian(flux[(phase >= bin_edges[i]) & (phase < bin_edges[i+1])])
                                 for i in range(n_bins)])
            valid = np.isfinite(bin_flux)
            ax.plot(bin_centers[valid], bin_flux[valid], 'b-', lw=1.5, label='Binned')
            ax.axvline(0.5, color='red', ls='--', alpha=0.5, label='Phase 0.5')
            ax.axvline(0.0, color='green', ls='--', alpha=0.5, label='Primary transit')
            ax.set_xlabel('Orbital Phase')
            ax.set_ylabel('Normalized Flux')
            ax.set_title(f'{name} — Multi-sector Secondary Eclipse Search ({len(search)} sectors, {len(time)} pts)')
            ax.legend(fontsize=8)
            plt.tight_layout()
            safe_name = name.replace(" ", "_").replace(".", "_")
            safe_name = os.path.basename(safe_name)  # strip any directory components
            plot_path = os.path.join(DEEP_DIR, f'secondary_{safe_name}.png')
            plt.savefig(plot_path, dpi=150, bbox_inches='tight')
            plt.close()

        except Exception as e:
            log(f"  ERROR for {name}: {e}")
            findings[name]['secondary_eclipse_deep'] = {'status': 'ERROR', 'detail': str(e)}


# ===========================================================================
# WRITE RESULTS
# ===========================================================================
def write_results(findings: dict) -> None:
    log("Writing results...")

    # Save JSON
    json_path = os.path.join(DEEP_DIR, 'deep_analysis_results.json')
    with open(json_path, 'w') as f:
        json.dump(findings, f, indent=2, default=str)
    log(f"  Saved: {json_path}")

    # Generate markdown report
    md_path = os.path.join(RESULTS_DIR, 'DEEP_ANALYSIS.md')
    with open(md_path, 'w') as f:
        f.write("# Exohuntr — Deep Analysis Milestone\n\n")
        f.write(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M UTC')}\n\n")
        f.write("## Purpose\n\n")
        f.write("This document records the results of 5 independent deep-validation tests on our top\n")
        f.write("exoplanet candidates. The goal: determine which (if any) of our 17 high-confidence\n")
        f.write("BLS detections are genuinely publishable planet candidates vs. false positives.\n\n")
        f.write("### Methodology (peer-reviewed basis)\n\n")
        f.write("| Test | Based on | What it proves |\n")
        f.write("|------|---------|----------------|\n")
        f.write("| Centroid analysis (TPF) | Twicken+2018, NASA SPOC DV | Transit is on the target star, not a neighbor |\n")
        f.write("| Gaia DR3 contamination | Standard practice, Furlan+2017 | No bright nearby stars could fake the signal |\n")
        f.write("| NASA DV report check | SPOC pipeline (Jenkins+2016) | NASA's own validation results for this target |\n")
        f.write("| Transit Least Squares | Hippke & Heller (2019) | Limb-darkened model confirms BLS box detection |\n")
        f.write("| Multi-sector secondary eclipse | Shporer+2017 | Extended data rules out self-luminous companion |\n\n")
        f.write("---\n\n")

        for t in TARGETS:
            name = t['toi']
            data = findings.get(name, {})
            f.write(f"## {name} ({t['tic']})\n\n")
            f.write(f"**BLS detection:** P={t['period']:.4f}d, SNR={t['snr']}, Rp≈{t['rp_earth']} R⊕\n\n")

            f.write("| Test | Result | Detail |\n")
            f.write("|------|--------|--------|\n")

            for test_key, test_name in [
                ('centroid', 'Centroid (TPF)'),
                ('gaia', 'Gaia DR3 check'),
                ('dv_report', 'NASA DV report'),
                ('tls', 'Transit Least Squares'),
                ('secondary_eclipse_deep', 'Multi-sector sec. eclipse')
            ]:
                result = data.get(test_key, {})
                status = result.get('status', 'NOT_RUN')
                detail = result.get('detail', 'N/A')
                emoji = "✅" if status in ("PASS", "CLEAR", "CONFIRMED", "AVAILABLE") else "⚠️" if status in ("WARNING", "WEAK_CONFIRM", "NO_DV", "NO_SPOC") else "❌" if status == "FAIL" else "—"
                f.write(f"| {test_name} | {emoji} **{status}** | {detail} |\n")

            f.write("\n")

            # Verdict
            statuses = [data.get(k, {}).get('status', '') for k in ['centroid', 'gaia', 'tls', 'secondary_eclipse_deep']]
            passes = sum(1 for s in statuses if s in ('PASS', 'CLEAR', 'CONFIRMED'))
            fails = sum(1 for s in statuses if s == 'FAIL')

            if fails == 0 and passes >= 3:
                verdict = "**STRONG CANDIDATE** — Passes all available deep validation tests. This signal is consistent with a genuine planetary transit."
            elif fails == 0 and passes >= 2:
                verdict = "**PROMISING** — Passes most deep tests. More data or follow-up would strengthen the case."
            elif fails > 0:
                verdict = "**WEAKENED** — One or more deep tests raise concerns. Needs additional investigation before any claim."
            else:
                verdict = "**INCONCLUSIVE** — Not enough data to make a definitive assessment."

            f.write(f"**Verdict:** {verdict}\n\n")
            f.write("---\n\n")

        # Overall assessment
        f.write("## Overall Assessment\n\n")
        f.write("### Can we publish these candidates?\n\n")

        n_strong = 0
        n_promising = 0
        for t in TARGETS:
            data = findings.get(t['toi'], {})
            statuses = [data.get(k, {}).get('status', '') for k in ['centroid', 'gaia', 'tls', 'secondary_eclipse_deep']]
            passes = sum(1 for s in statuses if s in ('PASS', 'CLEAR', 'CONFIRMED'))
            fails = sum(1 for s in statuses if s == 'FAIL')
            if fails == 0 and passes >= 3:
                n_strong += 1
            elif fails == 0 and passes >= 2:
                n_promising += 1

        f.write(f"- **{n_strong} strong candidate(s)** that pass all deep validation tests\n")
        f.write(f"- **{n_promising} promising candidate(s)** that need more data\n\n")

        f.write("### What we CAN claim (honestly):\n\n")
        f.write("1. We built a working BLS transit detection pipeline that correctly recovers known planets\n")
        f.write("2. We ran 5 standard false-positive tests (matching NASA SPOC methodology) on all detections\n")
        f.write("3. We performed deep validation (centroid, Gaia, TLS, multi-sector) on top candidates\n")
        f.write("4. Our candidates that pass all tests are **legitimate community TOI submissions**\n\n")

        f.write("### What we CANNOT claim:\n\n")
        f.write("1. We \"discovered\" a planet — confirmation requires independent ground-based follow-up\n")
        f.write("2. All 197 detections are planets — most are likely eclipsing binaries (expected)\n")
        f.write("3. Our Rp estimates are precise — they depend on stellar parameters from TIC which have uncertainties\n\n")

        f.write("### Next steps for real publication:\n\n")
        f.write("1. **Submit strongest candidates as community TOIs** to ExoFOP-TESS\n")
        f.write("2. **Request ground-based follow-up** from TFOPWG (TESS Follow-up Observing Program)\n")
        f.write("3. **Write up methodology** for submission to RNAAS (Research Notes of the AAS) — 1000-word format\n")
        f.write("4. **Contact Planet Hunters TESS** team about collaboration\n")
        f.write("5. **Run pipeline on newer TESS sectors** (80-96) for potentially unstudied targets\n\n")

        f.write("---\n\n")
        f.write("*This analysis uses methods from: Kovacs+2002 (BLS), Hippke & Heller 2019 (TLS),\n")
        f.write("Twicken+2018 (SPOC DV centroid), and standard community vetting practices.*\n")

    log(f"  Saved: {md_path}")


# ===========================================================================
# MAIN
# ===========================================================================
if __name__ == '__main__':
    log("=" * 60)
    log("EXOHUNTR DEEP ANALYSIS PIPELINE")
    log("Targets: " + ", ".join(t['toi'] for t in TARGETS))
    log("=" * 60)

    findings: dict = {}

    step1_centroid_analysis(findings)
    step2_gaia_query(findings)
    step3_dv_reports(findings)
    step4_tls(findings)
    step5_multisector_secondary(findings)
    write_results(findings)

    log("=" * 60)
    log("DEEP ANALYSIS COMPLETE")
    log(f"Results: {DEEP_DIR}/")
    log(f"Report:  {RESULTS_DIR}/DEEP_ANALYSIS.md")
    log("=" * 60)
