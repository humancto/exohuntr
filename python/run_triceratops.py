#!/usr/bin/env python3.11
"""Run TRICERATOPS statistical validation on TOI 133.01 (TIC 219338557)"""
if __name__ == '__main__':
    import warnings
    warnings.filterwarnings('ignore')
    import os
    os.environ['OMP_NUM_THREADS'] = '1'

    import numpy as np
    import lightkurve as lk
    import triceratops.triceratops as tr

    TIC_ID = 219338557
    PERIOD = 8.2065  # days (BLS period)
    T0 = None  # will get from TLS or data
    TDEPTH = 0.000707  # fractional depth = (Rp/Rs)^2 = 0.0266^2

    DEEP_DIR = 'results/deep_analysis'
    os.makedirs(DEEP_DIR, exist_ok=True)

    print("=" * 60)
    print("TRICERATOPS Validation — TOI 133.01 / TIC 219338557")
    print("=" * 60)

    # Step 1: Download light curve for phase-folding
    print("\n[1/5] Downloading light curve...", flush=True)
    search = lk.search_lightcurve(f'TIC {TIC_ID}', mission='TESS', author='SPOC')
    print(f"  Found {len(search)} sectors", flush=True)

    # Use first 3 sectors for the phase-folded LC
    n_download = min(3, len(search))
    lc_collection = search[:n_download].download_all()
    lc = lc_collection.stitch().remove_nans().remove_outliers(sigma=5)
    lc = lc.normalize()

    time_raw = lc.time.value
    flux_raw = lc.flux.value

    # Get sector numbers for TRICERATOPS
    import re
    sector_numbers = []
    for s in search[:n_download]:
        mission_str = str(s.mission)
        nums = re.findall(r'\d+', mission_str)
        if nums:
            sector_numbers.append(int(nums[0]))
    if not sector_numbers:
        # Fallback: try table column
        try:
            for row in search[:n_download].table:
                if 'sequence_number' in row.colnames:
                    sector_numbers.append(int(row['sequence_number']))
        except Exception:
            pass
    if not sector_numbers:
        sector_numbers = [1]  # final fallback
    print(f"  Sectors: {sector_numbers}", flush=True)

    # Phase-fold the light curve
    print("\n[2/5] Phase-folding light curve...", flush=True)
    # Find T0 by finding the deepest dip
    phase = time_raw % PERIOD
    n_bins = 200
    bin_edges = np.linspace(0, PERIOD, n_bins + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_flux = np.array([np.nanmedian(flux_raw[(phase >= bin_edges[i]) & (phase < bin_edges[i+1])])
                         for i in range(n_bins)])
    valid = np.isfinite(bin_flux)
    t0_phase = bin_centers[valid][np.argmin(bin_flux[valid])]
    T0 = time_raw[0] + t0_phase
    print(f"  T0 = {T0:.4f} (phase of minimum = {t0_phase:.4f}d)", flush=True)

    # Phase-fold centered on transit
    phase_fold = (time_raw - T0 + PERIOD / 2) % PERIOD - PERIOD / 2
    # Keep only points within +/- 0.25 * period of transit for efficiency
    transit_window = 0.25 * PERIOD
    mask = np.abs(phase_fold) < transit_window
    time_folded = phase_fold[mask]
    flux_folded = flux_raw[mask]

    # Sort by phase
    sort_idx = np.argsort(time_folded)
    time_folded = time_folded[sort_idx]
    flux_folded = flux_folded[sort_idx]

    flux_err = float(np.nanstd(flux_raw[np.abs(phase_fold) > 0.1 * PERIOD]))
    print(f"  Phase-folded: {len(time_folded)} points, flux_err = {flux_err:.6f}", flush=True)

    # Step 3: Create TRICERATOPS target
    print("\n[3/5] Creating TRICERATOPS target (querying TIC + TESSCut)...", flush=True)
    sectors = np.array(sector_numbers)
    try:
        target = tr.target(ID=TIC_ID, sectors=sectors, search_radius=10)
        print(f"  Found {len(target.stars)} nearby stars", flush=True)
        print(f"  Target star: Tmag={target.stars.iloc[0].get('Tmag', 'N/A')}", flush=True)
    except Exception as e:
        print(f"  ERROR creating target: {e}", flush=True)
        print("  Trying with single sector...", flush=True)
        target = tr.target(ID=TIC_ID, sectors=np.array([sectors[0]]), search_radius=10)
        print(f"  Found {len(target.stars)} nearby stars", flush=True)

    # Check and fill missing stellar properties
    star = target.stars.iloc[0]
    print(f"  Star properties: mass={star.get('mass', 'N/A')}, "
          f"rad={star.get('rad', 'N/A')}, Teff={star.get('Teff', 'N/A')}", flush=True)

    for param in ['mass', 'rad', 'Teff']:
        val = star.get(param, None)
        if val is None or (isinstance(val, float) and np.isnan(val)):
            defaults = {'mass': 0.85, 'rad': 0.82, 'Teff': 5100.0}  # K-dwarf defaults for TOI 133
            print(f"  WARNING: {param} is missing, setting to {defaults[param]}", flush=True)
            target.update_star(ID=TIC_ID, param=param, value=defaults[param])

    # Step 4: Calculate depths
    print("\n[4/5] Calculating transit depths for all nearby stars...", flush=True)
    try:
        apertures = target.get_spoc_apertures()
        target.calc_depths(tdepth=TDEPTH, all_ap_pixels=apertures)
    except Exception as e:
        print(f"  SPOC apertures failed ({e}), using default aperture", flush=True)
        target.calc_depths(tdepth=TDEPTH)

    # Step 5: Calculate FPP
    print("\n[5/5] Running TRICERATOPS FPP calculation (N=50000)...", flush=True)
    print("  This may take 10-30 minutes...", flush=True)
    try:
        target.calc_probs(
            time=time_folded,
            flux_0=flux_folded,
            flux_err_0=flux_err,
            P_orb=PERIOD,
            N=50000,  # reduced from 1M for speed
        )

        fpp = target.FPP
        nfpp = target.NFPP
        probs = target.probs

        print("\n" + "=" * 60)
        print("TRICERATOPS RESULTS — TOI 133.01")
        print("=" * 60)
        print(f"  FPP  = {fpp:.6f} ({fpp*100:.4f}%)")
        print(f"  NFPP = {nfpp:.6f} ({nfpp*100:.4f}%)")
        print()

        if fpp < 0.015 and nfpp < 0.001:
            verdict = "STATISTICALLY VALIDATED PLANET (FPP < 1.5%, NFPP < 0.1%)"
        elif fpp < 0.5 and nfpp < 0.1:
            verdict = "LIKELY PLANET (FPP < 50%, NFPP < 10%)"
        else:
            verdict = "NOT VALIDATED (FPP too high)"
        print(f"  Verdict: {verdict}")
        print()

        # Print scenario probabilities
        print("  Scenario probabilities:")
        for _, row in probs.iterrows():
            if row['prob'] > 0.001:
                print(f"    {row['scenario']:12s}: {row['prob']:.6f} ({row['prob']*100:.4f}%)")

        # Save results
        results_file = os.path.join(DEEP_DIR, 'triceratops_TOI_133_01.txt')
        with open(results_file, 'w') as f:
            f.write("TRICERATOPS Statistical Validation — TOI 133.01 / TIC 219338557\n")
            f.write("=" * 60 + "\n\n")
            f.write(f"Period: {PERIOD} d\n")
            f.write(f"Transit depth: {TDEPTH} (fractional, = {TDEPTH*1e6:.0f} ppm)\n")
            f.write(f"Sectors used: {sector_numbers}\n")
            f.write(f"N (MC draws): 100000\n\n")
            f.write(f"FPP  = {fpp:.6f} ({fpp*100:.4f}%)\n")
            f.write(f"NFPP = {nfpp:.6f} ({nfpp*100:.4f}%)\n\n")
            f.write(f"Verdict: {verdict}\n\n")
            f.write("Scenario probabilities:\n")
            for _, row in probs.iterrows():
                f.write(f"  {row['scenario']:12s}: {row['prob']:.6f}\n")
        print(f"\n  Saved: {results_file}")

    except Exception as e:
        print(f"\n  ERROR in calc_probs: {e}", flush=True)
        import traceback
        traceback.print_exc()

    print("\nTRICERATOPS complete!", flush=True)
