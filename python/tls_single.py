#!/usr/bin/env python3.11
"""Run TLS on TOI 210.01 and TOI 155.01 — single-threaded to avoid spawn issues"""
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
    from transitleastsquares import transitleastsquares

    DEEP_DIR = 'results/deep_analysis'
    os.makedirs(DEEP_DIR, exist_ok=True)

    targets = [
        {"toi": "TOI 210.01", "tic": "TIC 141608198", "period": 8.9884},
        {"toi": "TOI 155.01", "tic": "TIC 129637892", "period": 5.4504},
    ]

    for t in targets:
        name = t['toi']
        print(f"\n=== {name} ===", flush=True)
        search = lk.search_lightcurve(t['tic'], mission='TESS', author='SPOC')
        if len(search) == 0:
            search = lk.search_lightcurve(t['tic'], mission='TESS')

        print(f"Found {len(search)} light curves", flush=True)

        # Only 2 sectors
        n_download = min(2, len(search))
        print(f"Downloading {n_download} sectors...", flush=True)
        lc_collection = search[:n_download].download_all()
        lc = lc_collection.stitch()
        lc = lc.remove_nans().remove_outliers(sigma=5)

        time = lc.time.value
        flux = lc.flux.value
        flux = flux / np.nanmedian(flux)

        # Aggressive binning to keep fast
        if len(time) > 50000:
            print(f"Binning {len(time)} points...", flush=True)
            lc_norm = lc.normalize()
            lc_binned = lc_norm.bin(time_bin_size=0.02)  # ~30 min bins
            time = lc_binned.time.value
            flux = lc_binned.flux.value
            mask = np.isfinite(time) & np.isfinite(flux)
            time = time[mask]
            flux = flux[mask]

        print(f"Running TLS on {len(time)} data points (single-threaded)...", flush=True)
        model = transitleastsquares(time, flux)
        results = model.power(
            period_min=t['period'] * 0.9,
            period_max=t['period'] * 1.1,
            oversampling_factor=2,
            duration_grid_step=1.15,
            use_threads=1,
        )

        tls_period = float(results.period)
        tls_sde = float(results.SDE)
        tls_rp_rs = float(results.rp_rs) if hasattr(results, 'rp_rs') else None
        tls_t0 = float(results.T0) if hasattr(results, 'T0') else None
        period_diff_pct = abs(tls_period - t['period']) / t['period'] * 100

        if period_diff_pct < 1.0 and tls_sde > 6.0:
            status = "CONFIRMED"
        elif period_diff_pct < 1.0:
            status = "WEAK_CONFIRM"
        else:
            status = "DISAGREE"

        print(f"Result: {status}", flush=True)
        print(f"  TLS period:  {tls_period:.4f}d", flush=True)
        print(f"  BLS period:  {t['period']:.4f}d", flush=True)
        print(f"  Period diff: {period_diff_pct:.3f}%", flush=True)
        print(f"  SDE:         {tls_sde:.1f}", flush=True)
        if tls_rp_rs:
            print(f"  Rp/Rs:       {tls_rp_rs:.5f}", flush=True)

        # Plot
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        fig.suptitle(f'{name} — Transit Least Squares', fontsize=14)
        axes[0].plot(results.periods, results.power, 'k-', lw=0.5)
        axes[0].axvline(tls_period, color='red', ls='--', label=f'TLS P={tls_period:.4f}d')
        axes[0].axvline(t['period'], color='blue', ls=':', label=f'BLS P={t["period"]:.4f}d')
        axes[0].set_xlabel('Period (days)')
        axes[0].set_ylabel('SDE')
        axes[0].legend(fontsize=9)

        if tls_t0:
            phase = (time - tls_t0 + tls_period/2) % tls_period - tls_period/2
        else:
            phase = (time % tls_period)
        axes[1].scatter(phase * 24, flux, s=1, alpha=0.3, c='gray')
        n_bins = 100
        be = np.linspace(phase.min(), phase.max(), n_bins + 1)
        bc = (be[:-1] + be[1:]) / 2
        bf = np.array([np.nanmedian(flux[(phase >= be[i]) & (phase < be[i+1])]) for i in range(n_bins)])
        v = np.isfinite(bf)
        axes[1].plot(bc[v] * 24, bf[v], 'r-', lw=2, label='Binned')
        axes[1].set_xlabel('Hours from mid-transit')
        axes[1].set_ylabel('Normalized flux')
        axes[1].set_xlim(-5, 5)
        axes[1].legend(fontsize=9)
        plt.tight_layout()
        pp = os.path.join(DEEP_DIR, 'tls_%s.png' % name.replace(" ", "_").replace(".", "_"))
        plt.savefig(pp, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {pp}", flush=True)

    print("\nAll TLS complete!", flush=True)
