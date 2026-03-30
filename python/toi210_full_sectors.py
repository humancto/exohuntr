#!/usr/bin/env python3.11
"""Download all available sectors for TOI 210.01 and run definitive secondary eclipse search"""
if __name__ == '__main__':
    import warnings
    warnings.filterwarnings('ignore')
    import os
    os.environ['OMP_NUM_THREADS'] = '1'

    import numpy as np
    import lightkurve as lk
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy import stats

    TIC_ID = 'TIC 141608198'
    PERIOD = 8.9884  # days
    TOI_NAME = 'TOI 210.01'

    DEEP_DIR = 'results/deep_analysis'
    os.makedirs(DEEP_DIR, exist_ok=True)

    print("=" * 60)
    print(f"Full-Sector Secondary Eclipse — {TOI_NAME} / {TIC_ID}")
    print("=" * 60)

    # Step 1: Find all available sectors
    print("\n[1/4] Searching for all available sectors...", flush=True)
    search = lk.search_lightcurve(TIC_ID, mission='TESS', author='SPOC')
    print(f"  Found {len(search)} sectors", flush=True)

    # Download in batches to avoid memory issues
    batch_size = 10
    all_time = []
    all_flux = []
    n_sectors_downloaded = 0

    print(f"\n[2/4] Downloading in batches of {batch_size}...", flush=True)
    for batch_start in range(0, len(search), batch_size):
        batch_end = min(batch_start + batch_size, len(search))
        print(f"  Downloading sectors {batch_start+1}-{batch_end} of {len(search)}...", flush=True)
        try:
            batch = search[batch_start:batch_end].download_all()
            for lc in batch:
                lc = lc.remove_nans().remove_outliers(sigma=5).normalize()
                t = lc.time.value
                f = lc.flux.value
                mask = np.isfinite(t) & np.isfinite(f)
                all_time.append(t[mask])
                all_flux.append(f[mask])
                n_sectors_downloaded += 1
        except Exception as e:
            print(f"  WARNING: Batch failed: {e}", flush=True)
            # Try individual downloads
            for i in range(batch_start, batch_end):
                try:
                    lc = search[i].download()
                    lc = lc.remove_nans().remove_outliers(sigma=5).normalize()
                    t = lc.time.value
                    f = lc.flux.value
                    mask = np.isfinite(t) & np.isfinite(f)
                    all_time.append(t[mask])
                    all_flux.append(f[mask])
                    n_sectors_downloaded += 1
                except Exception as e2:
                    print(f"    Sector {i+1} failed: {e2}", flush=True)

    time = np.concatenate(all_time)
    flux = np.concatenate(all_flux)
    print(f"\n  Downloaded {n_sectors_downloaded} sectors, {len(time)} total data points", flush=True)

    # Step 3: Secondary eclipse search
    print("\n[3/4] Searching for secondary eclipse at phase 0.5...", flush=True)

    # Phase-fold
    # First find T0 from the primary transit
    phase = time % PERIOD
    n_bins = 500
    bin_edges = np.linspace(0, PERIOD, n_bins + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_flux = np.array([np.nanmedian(flux[(phase >= bin_edges[i]) & (phase < bin_edges[i+1])])
                         for i in range(n_bins)])
    valid = np.isfinite(bin_flux)
    t0_phase = bin_centers[valid][np.argmin(bin_flux[valid])]

    # Phase-fold centered on primary transit
    phase_fold = (time - time[0] - t0_phase + PERIOD / 2) % PERIOD - PERIOD / 2

    # Secondary eclipse region: phase 0.4-0.6 of period (centered on 0.5)
    sec_phase = phase_fold + PERIOD / 2  # shift so secondary is at center
    sec_phase = (sec_phase + PERIOD / 2) % PERIOD - PERIOD / 2

    # Compare in-eclipse vs out-of-eclipse
    eclipse_half_width = 0.02 * PERIOD  # ~2% of period
    in_eclipse = np.abs(sec_phase) < eclipse_half_width
    out_eclipse = (np.abs(sec_phase) > 0.1 * PERIOD) & (np.abs(sec_phase) < 0.4 * PERIOD)

    flux_in = flux[in_eclipse]
    flux_out = flux[out_eclipse]

    if len(flux_in) > 10 and len(flux_out) > 10:
        depth = np.nanmedian(flux_out) - np.nanmedian(flux_in)
        # Bootstrap uncertainty
        n_boot = 10000
        depths_boot = np.zeros(n_boot)
        for b in range(n_boot):
            fi = np.random.choice(flux_in, size=len(flux_in), replace=True)
            fo = np.random.choice(flux_out, size=len(flux_out), replace=True)
            depths_boot[b] = np.nanmedian(fo) - np.nanmedian(fi)
        depth_err = np.std(depths_boot)
        depth_sigma = depth / depth_err if depth_err > 0 else 0

        print(f"  In-eclipse points:  {len(flux_in)}")
        print(f"  Out-eclipse points: {len(flux_out)}")
        print(f"  Secondary depth:    {depth:.8f} ({depth*100:.6f}%)")
        print(f"  Depth uncertainty:  {depth_err:.8f}")
        print(f"  Significance:       {depth_sigma:.1f} sigma")

        if depth_sigma < 3:
            verdict = "NO SECONDARY ECLIPSE DETECTED (< 3 sigma)"
        elif depth_sigma < 5:
            verdict = f"MARGINAL DETECTION ({depth_sigma:.1f} sigma)"
        else:
            if depth * 100 < 0.01:
                verdict = f"SHALLOW DETECTION ({depth_sigma:.1f} sigma, {depth*100:.4f}%) — consistent with planetary thermal emission"
            else:
                verdict = f"SIGNIFICANT DETECTION ({depth_sigma:.1f} sigma, {depth*100:.4f}%) — possible binary"
    else:
        depth = 0
        depth_err = 0
        depth_sigma = 0
        verdict = "INSUFFICIENT DATA"

    print(f"\n  Verdict: {verdict}")

    # Step 4: Plot
    print("\n[4/4] Generating plot...", flush=True)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(f'{TOI_NAME} — Full-Sector Secondary Eclipse Search ({n_sectors_downloaded} sectors, {len(time):,} points)', fontsize=12)

    # Left: full phase-folded LC
    phase_hours = phase_fold * 24
    axes[0].scatter(phase_hours, flux, s=0.1, alpha=0.1, c='gray')
    # Binned
    n_plot_bins = 300
    be = np.linspace(phase_fold.min(), phase_fold.max(), n_plot_bins + 1)
    bc = (be[:-1] + be[1:]) / 2
    bf = np.array([np.nanmedian(flux[(phase_fold >= be[i]) & (phase_fold < be[i+1])]) for i in range(n_plot_bins)])
    v = np.isfinite(bf)
    axes[0].plot(bc[v] * 24, bf[v], 'r-', lw=1.5, label='Binned median')
    axes[0].set_xlabel('Hours from mid-transit')
    axes[0].set_ylabel('Normalized flux')
    axes[0].set_title('Primary transit')
    axes[0].set_xlim(-8, 8)
    axes[0].legend(fontsize=9)

    # Right: secondary eclipse region
    sec_hours = sec_phase * 24
    axes[1].scatter(sec_hours, flux, s=0.1, alpha=0.1, c='gray')
    be2 = np.linspace(sec_phase.min(), sec_phase.max(), n_plot_bins + 1)
    bc2 = (be2[:-1] + be2[1:]) / 2
    bf2 = np.array([np.nanmedian(flux[(sec_phase >= be2[i]) & (sec_phase < be2[i+1])]) for i in range(n_plot_bins)])
    v2 = np.isfinite(bf2)
    axes[1].plot(bc2[v2] * 24, bf2[v2], 'b-', lw=1.5, label='Binned median')
    axes[1].axhline(1.0, color='gray', ls='--', alpha=0.5)
    axes[1].axvspan(-eclipse_half_width * 24, eclipse_half_width * 24, alpha=0.2, color='red', label='Eclipse window')
    axes[1].set_xlabel('Hours from expected secondary')
    axes[1].set_ylabel('Normalized flux')
    axes[1].set_title(f'Secondary eclipse ({depth_sigma:.1f}σ, depth={depth*100:.4f}%)')
    axes[1].set_xlim(-8, 8)
    axes[1].legend(fontsize=9)

    plt.tight_layout()
    plot_path = os.path.join(DEEP_DIR, 'secondary_TOI_210_01_full.png')
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {plot_path}", flush=True)

    # Save results
    results_file = os.path.join(DEEP_DIR, 'toi210_full_secondary.txt')
    with open(results_file, 'w') as f:
        f.write(f"Full-Sector Secondary Eclipse Search — {TOI_NAME} / {TIC_ID}\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Period: {PERIOD} d\n")
        f.write(f"Sectors downloaded: {n_sectors_downloaded}\n")
        f.write(f"Total data points: {len(time)}\n\n")
        f.write(f"Secondary eclipse depth: {depth:.8f} ({depth*100:.6f}%)\n")
        f.write(f"Depth uncertainty: {depth_err:.8f}\n")
        f.write(f"Significance: {depth_sigma:.1f} sigma\n\n")
        f.write(f"Verdict: {verdict}\n")
    print(f"  Saved: {results_file}")

    print(f"\nDone! {n_sectors_downloaded} sectors analyzed.", flush=True)
