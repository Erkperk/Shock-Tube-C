#!/usr/bin/env python3
"""
verify_run.py — 3-panel verification plot for shock tube runs.

Usage:
    python3 verify_run.py <timeseries_A.dat> [timeseries_B.dat]

Panel 1: Initial mass flow (first 10s) — should show ~17 t/s spike
Panel 2: Closed-end pressure vs time (full sim) — p_base in bar vs hours
Panel 3: Time shift (only with two files) — t_A(p) - t_B(p) for p_base
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt


def load_timeseries(fname):
    """Load timeseries .dat file, skipping comment lines."""
    data = np.loadtxt(fname, comments='#')
    # Columns: t, massflow[kg/s], u_out, rho_out, p_out, T_out, totmass,
    #          p_base, p_mid, a_out, Mach_out, poutlet_choked, choked_active
    return data


def time_at_pressure(t, p, p_targets):
    """Interpolate t(p) for monotonically decreasing p_base."""
    results = {}
    for pt in p_targets:
        idx = np.where(p < pt)[0]
        if len(idx) > 0:
            i = idx[0]
            if i > 0:
                # Linear interpolation
                frac = (pt - p[i-1]) / (p[i] - p[i-1])
                results[pt] = t[i-1] + frac * (t[i] - t[i-1])
            else:
                results[pt] = t[0]
        else:
            results[pt] = np.nan
    return results


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    fname_a = sys.argv[1]
    fname_b = sys.argv[2] if len(sys.argv) > 2 else None
    two_files = fname_b is not None

    data_a = load_timeseries(fname_a)
    if two_files:
        data_b = load_timeseries(fname_b)

    t_a = data_a[:, 0]
    mf_a = data_a[:, 1]         # mass flow [kg/s]
    pbase_a = data_a[:, 7] / 1e5  # p_base [bar]

    # Derive output prefix from first filename
    base = os.path.basename(fname_a)
    prefix = base.replace('timeseries_', '').replace('.dat', '')
    outname = f"verify_{prefix}.png"

    n_panels = 3 if two_files else 2
    fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels, 4))
    if n_panels == 2:
        axes = list(axes)
        axes.append(None)

    # Panel 1: Initial mass flow (first 10s)
    ax = axes[0]
    mask = t_a < 10.0
    ax.plot(t_a[mask], mf_a[mask] / 1000.0, 'b-', linewidth=0.8, label=os.path.basename(fname_a))
    if two_files:
        t_b = data_b[:, 0]
        mf_b = data_b[:, 1]
        mask_b = t_b < 10.0
        ax.plot(t_b[mask_b], mf_b[mask_b] / 1000.0, 'r--', linewidth=0.8, label=os.path.basename(fname_b))
        ax.legend(fontsize=7)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Mass flow [t/s]')
    ax.set_title('Initial mass flow (first 10s)')
    ax.grid(True, alpha=0.3)

    # Panel 2: Closed-end pressure vs time
    ax = axes[1]
    ax.plot(t_a / 3600.0, pbase_a, 'b-', linewidth=0.8, label=os.path.basename(fname_a))
    if two_files:
        pbase_b = data_b[:, 7] / 1e5
        ax.plot(t_b / 3600.0, pbase_b, 'r--', linewidth=0.8, label=os.path.basename(fname_b))
        ax.legend(fontsize=7)
    ax.set_xlabel('Time [hours]')
    ax.set_ylabel('Closed-end pressure [bar]')
    ax.set_title('Pressure at closed end')
    ax.grid(True, alpha=0.3)

    # Panel 3: Time shift (only with two files)
    if two_files:
        ax = axes[2]
        # Build t(p) for both, compare
        # Use common pressure range
        p_min = max(pbase_a.min(), pbase_b.min())
        p_max = min(pbase_a[0], pbase_b[0])
        p_targets = np.linspace(p_max - 1.0, p_min + 1.0, 200)
        p_targets = p_targets[p_targets > p_min]
        p_targets = p_targets[p_targets < p_max]

        ta_of_p = time_at_pressure(t_a, pbase_a, p_targets)
        tb_of_p = time_at_pressure(t_b, pbase_b, p_targets)

        shifts = []
        ps = []
        for pt in p_targets:
            if pt in ta_of_p and pt in tb_of_p:
                s = ta_of_p[pt] - tb_of_p[pt]
                if not np.isnan(s):
                    shifts.append(s)
                    ps.append(pt)

        ax.plot(ps, shifts, 'k-', linewidth=0.8)

        # Annotate key reference values
        ref_points = {160.0: 30.6, 50.0: 926.0}
        for p_ref, t_ref in ref_points.items():
            if p_min < p_ref < p_max:
                ax.axhline(y=t_ref, color='gray', linestyle=':', alpha=0.5)
                ax.annotate(f'HLLC ref: {t_ref}s @ {p_ref} bar',
                           xy=(p_ref, t_ref), fontsize=7, color='gray')

        ax.set_xlabel('Pressure [bar]')
        ax.set_ylabel('Time shift A - B [s]')
        ax.set_title('Time shift vs pressure')
        ax.grid(True, alpha=0.3)
        ax.invert_xaxis()

    plt.tight_layout()
    plt.savefig(outname, dpi=150)
    print(f"Saved: {outname}")
    plt.close()


if __name__ == '__main__':
    main()
