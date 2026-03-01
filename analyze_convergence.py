#!/usr/bin/env python3
"""
Grid convergence analysis for NS1A vs NS1B German-end time shift.
Compares multiple resolutions, solvers, and merge vs no-merge.
"""
import numpy as np
import matplotlib.pyplot as plt
import os

def load_if_exists(fname):
    path = fname
    if os.path.exists(path):
        try:
            return np.loadtxt(path)
        except ValueError:
            # File may be incomplete (still writing); read valid lines
            lines = []
            with open(path) as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    parts = line.split()
                    if len(parts) == 13:
                        try:
                            lines.append([float(x) for x in parts])
                        except ValueError:
                            continue
            if lines:
                return np.array(lines)
    return None

def compute_time_shift(data_a, data_b, p_levels):
    """Compute time shift: t_A(p) - t_B(p) for closed-end pressure.
    Returns NaN for pressures below the minimum available in either dataset."""
    t_a = data_a[:, 0]
    t_b = data_b[:, 0]
    p_a = data_a[:, 7] / 1e5  # p_base in bar
    p_b = data_b[:, 7] / 1e5

    # Minimum available pressure (end of simulation)
    p_min = max(p_a.min(), p_b.min())

    # Both are monotonically decreasing; flip for interp
    t_at_p_a = np.interp(p_levels, p_a[::-1], t_a[::-1])
    t_at_p_b = np.interp(p_levels, p_b[::-1], t_b[::-1])
    shift = t_at_p_a - t_at_p_b

    # Mask values below minimum available pressure
    shift[p_levels < p_min + 1.0] = np.nan
    return shift

# Pressure levels for analysis
p_levels = np.linspace(165, 10, 5000)

# --- Load all datasets ---
runs = {}

# Grid convergence: LF + merge at different dx
for dx in [50, 20, 10, 5]:
    if dx == 50:
        fa = 'timeseries_ns1a_lf_lf.dat'
        fb = 'timeseries_ns1b_lf_lf.dat'
    else:
        fa = f'timeseries_ns1a_dx{dx}_lf.dat'
        fb = f'timeseries_ns1b_dx{dx}_lf.dat'
    da, db = load_if_exists(fa), load_if_exists(fb)
    if da is not None and db is not None:
        runs[f'LF merge dx={dx}m'] = (da, db)

# HLLC + merge at dx=20
da = load_if_exists('timeseries_ns1a_hllc_dx20.dat')
db = load_if_exists('timeseries_ns1b_hllc_dx20.dat')
if da is not None and db is not None:
    runs['HLLC merge dx=20m'] = (da, db)

# LF no-merge at dx=20
da = load_if_exists('timeseries_ns1a_dx20_nomerge_lf.dat')
db = load_if_exists('timeseries_ns1b_dx20_nomerge_lf.dat')
if da is not None and db is not None:
    runs['LF no-merge dx=20m'] = (da, db)

# HLLC no-merge at dx=50
da = load_if_exists('timeseries_ns1a_hllc_dx50_nomerge.dat')
db = load_if_exists('timeseries_ns1b_hllc_dx50_nomerge.dat')
if da is not None and db is not None:
    runs['HLLC no-merge dx=50m'] = (da, db)
    pa_min = da[:, 7].min() / 1e5
    pb_min = db[:, 7].min() / 1e5
    print(f"  HLLC no-merge dx=50: NS1A pmin={pa_min:.1f} bar, NS1B pmin={pb_min:.1f} bar")

# HLLC no-merge at dx=20
da = load_if_exists('timeseries_ns1a_hllc_dx20_nomerge.dat')
db = load_if_exists('timeseries_ns1b_hllc_dx20_nomerge.dat')
if da is not None and db is not None:
    pa_min = da[:, 7].min() / 1e5
    pb_min = db[:, 7].min() / 1e5
    print(f"  HLLC no-merge dx=20: NS1A pmin={pa_min:.1f} bar, NS1B pmin={pb_min:.1f} bar")
    if pa_min < 55 and pb_min < 55:
        runs['HLLC no-merge dx=20m'] = (da, db)
    else:
        print(f"  -> Skipping HLLC no-merge dx=20 for full analysis (NS1A only down to {pa_min:.0f} bar)")
        # Still use it for high-pressure analysis
        runs['HLLC no-merge dx=20m (partial)'] = (da, db)

print(f"\nLoaded {len(runs)} run pairs: {list(runs.keys())}")
print()

# --- Compute time shifts ---
shifts = {}
for name, (da, db) in runs.items():
    shifts[name] = compute_time_shift(da, db, p_levels)

# --- Print table at key pressures ---
key_pressures = [160, 155, 150, 145, 140, 130, 120, 110, 100, 90, 80, 70, 60, 50, 40, 30, 20, 15]
header = f"{'P [bar]':>8}"
for name in shifts:
    short = name.replace(' merge ', ' ').replace('m', '').replace(' (partial)', '*')
    header += f"  {short:>20}"
print(header)
print("-" * len(header))

for p in key_pressures:
    idx = np.argmin(np.abs(p_levels - p))
    row = f"{p:8.0f}"
    for name in shifts:
        v = shifts[name][idx]
        if np.isnan(v):
            row += f"  {'---':>20}"
        else:
            row += f"  {v:20.1f}"
    print(row)

# --- Find where shift = 15s for each run ---
print("\n--- Time shift = 15 seconds ---")
for name, dt in shifts.items():
    valid = ~np.isnan(dt)
    if np.any(valid):
        dt_valid = dt.copy()
        dt_valid[~valid] = 1e10  # avoid NaN in argmin
        idx = np.argmin(np.abs(dt_valid - 15.0))
        if not np.isnan(dt[idx]):
            print(f"  {name:35s}: P = {p_levels[idx]:.1f} bar, shift = {dt[idx]:.1f} s")
        else:
            print(f"  {name:35s}: data insufficient")

# --- Shift at 50 bar ---
print("\n--- Time shift at 50 bar ---")
idx_50 = np.argmin(np.abs(p_levels - 50.0))
for name, dt in shifts.items():
    v = dt[idx_50]
    if np.isnan(v):
        print(f"  {name:35s}: data insufficient (run still in progress)")
    else:
        print(f"  {name:35s}: {v:.1f} s")

# --- Convergence summary ---
print("\n" + "="*80)
print("CONVERGENCE SUMMARY")
print("="*80)
print(f"\n{'Configuration':>35}  {'@160bar':>8}  {'@140bar':>8}  {'@100bar':>8}  {'@50bar':>8}  {'@30bar':>8}")
print("-"*85)
idx_160 = np.argmin(np.abs(p_levels - 160.0))
idx_140 = np.argmin(np.abs(p_levels - 140.0))
idx_100 = np.argmin(np.abs(p_levels - 100.0))
idx_50 = np.argmin(np.abs(p_levels - 50.0))
idx_30 = np.argmin(np.abs(p_levels - 30.0))
order = ['LF merge dx=50m', 'LF merge dx=20m', 'LF merge dx=10m', 'LF merge dx=5m',
         'HLLC merge dx=20m', 'LF no-merge dx=20m',
         'HLLC no-merge dx=50m', 'HLLC no-merge dx=20m', 'HLLC no-merge dx=20m (partial)']
for name in order:
    if name in shifts:
        dt = shifts[name]
        vals = []
        for idx in [idx_160, idx_140, idx_100, idx_50, idx_30]:
            v = dt[idx]
            vals.append(f'{v:6.1f} s' if not np.isnan(v) else '   --- ')
        print(f"  {name:35s}  {'  '.join(vals)}")

# Analysis
print("\n--- Analysis ---")
# Collect 50 bar results
vals_50 = {}
for name, dt in shifts.items():
    v = dt[idx_50]
    if not np.isnan(v):
        vals_50[name] = v
if vals_50:
    all_v = list(vals_50.values())
    mean_50 = np.mean(all_v)
    std_50 = np.std(all_v)
    min_50 = min(all_v)
    max_50 = max(all_v)
    print(f"  At 50 bar: mean={mean_50:.0f}s, std={std_50:.0f}s, range={min_50:.0f}-{max_50:.0f}s")

    # HLLC no-merge values (highest accuracy)
    hllc_nomerge = {k: v for k, v in vals_50.items() if 'HLLC no-merge' in k}
    if hllc_nomerge:
        hv = list(hllc_nomerge.values())
        print(f"  HLLC no-merge (best): {', '.join(f'{v:.0f}s' for v in hv)}")

# Collect 160 bar results
vals_160 = {}
for name, dt in shifts.items():
    v = dt[idx_160]
    if not np.isnan(v):
        vals_160[name] = v
if vals_160:
    all_v = list(vals_160.values())
    mean_160 = np.mean(all_v)
    std_160 = np.std(all_v)
    print(f"  At 160 bar: mean={mean_160:.1f}s, std={std_160:.1f}s, range={min(all_v):.1f}-{max(all_v):.1f}s")

# Collect 140 bar results
vals_140 = {}
for name, dt in shifts.items():
    v = dt[idx_140]
    if not np.isnan(v):
        vals_140[name] = v
if vals_140:
    all_v = list(vals_140.values())
    print(f"  At 140 bar: mean={np.mean(all_v):.1f}s, std={np.std(all_v):.1f}s, range={min(all_v):.1f}-{max(all_v):.1f}s")

# --- Plot ---
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

colors = {
    'LF merge dx=50m': '#cccccc',
    'LF merge dx=20m': '#1f77b4',
    'LF merge dx=10m': '#2ca02c',
    'LF merge dx=5m':  '#d62728',
    'HLLC merge dx=20m': '#ff7f0e',
    'LF no-merge dx=20m': '#9467bd',
    'HLLC no-merge dx=50m': '#e377c2',
    'HLLC no-merge dx=20m': '#17becf',
    'HLLC no-merge dx=20m (partial)': '#17becf',
}
lstyles = {
    'LF merge dx=50m': '-',
    'LF merge dx=20m': '-',
    'LF merge dx=10m': '-',
    'LF merge dx=5m':  '-',
    'HLLC merge dx=20m': '--',
    'LF no-merge dx=20m': ':',
    'HLLC no-merge dx=50m': '--',
    'HLLC no-merge dx=20m': '--',
    'HLLC no-merge dx=20m (partial)': ':',
}
lwidths = {
    'HLLC no-merge dx=50m': 2.5,
    'HLLC no-merge dx=20m': 2.5,
}

# Left: full time shift vs pressure
ax = axes[0]
for name, dt in shifts.items():
    c = colors.get(name, 'black')
    ls = lstyles.get(name, '-')
    lw = lwidths.get(name, 1.5)
    ax.plot(p_levels, dt, color=c, ls=ls, lw=lw, label=name)
ax.set_xlabel('Closed-end pressure [bar]')
ax.set_ylabel('Time shift NS1A - NS1B [s]')
ax.set_title('Time shift vs pressure level')
ax.legend(fontsize=7, loc='upper right')
ax.grid(True, alpha=0.3)
ax.set_xlim(165, 10)

# Middle: zoom on high pressure (0-200s shift)
ax = axes[1]
for name, dt in shifts.items():
    c = colors.get(name, 'black')
    ls = lstyles.get(name, '-')
    lw = lwidths.get(name, 1.5)
    ax.plot(p_levels, dt, color=c, ls=ls, lw=lw, label=name)
ax.set_xlabel('Closed-end pressure [bar]')
ax.set_ylabel('Time shift NS1A - NS1B [s]')
ax.set_title('Zoom: P > 100 bar')
ax.legend(fontsize=7)
ax.grid(True, which='major', alpha=0.4)
ax.grid(True, which='minor', alpha=0.2)
ax.set_xlim(165, 100)
ax.set_ylim(0, 200)

# Right: convergence at specific pressures (bar chart)
ax = axes[2]
bar_names = ['LF merge\ndx=50', 'LF merge\ndx=20', 'LF merge\ndx=10', 'LF merge\ndx=5',
             'HLLC merge\ndx=20', 'LF nomerge\ndx=20', 'HLLC nomerge\ndx=50', 'HLLC nomerge\ndx=20']
bar_keys = ['LF merge dx=50m', 'LF merge dx=20m', 'LF merge dx=10m', 'LF merge dx=5m',
            'HLLC merge dx=20m', 'LF no-merge dx=20m', 'HLLC no-merge dx=50m', 'HLLC no-merge dx=20m']
x_pos = np.arange(len(bar_names))
width = 0.35

vals_at_160 = []
vals_at_50 = []
for k in bar_keys:
    if k in shifts:
        v160 = shifts[k][idx_160]
        v50 = shifts[k][idx_50]
        vals_at_160.append(v160 if not np.isnan(v160) else 0)
        vals_at_50.append(v50 if not np.isnan(v50) else 0)
    else:
        vals_at_160.append(0)
        vals_at_50.append(0)

bars1 = ax.bar(x_pos - width/2, vals_at_160, width, label='@ 160 bar', color='#1f77b4', alpha=0.7)
bars2 = ax.bar(x_pos + width/2, vals_at_50, width, label='@ 50 bar', color='#d62728', alpha=0.7)
ax.set_xticks(x_pos)
ax.set_xticklabels(bar_names, fontsize=7)
ax.set_ylabel('Time shift [s]')
ax.set_title('Convergence at key pressures')
ax.legend()
ax.grid(True, axis='y', alpha=0.3)

plt.suptitle('Grid convergence: NS1A vs NS1B time shift', fontsize=13, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig('convergence_time_shift.png', dpi=150)
print("\nSaved: convergence_time_shift.png")
