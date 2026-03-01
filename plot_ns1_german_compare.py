#!/usr/bin/env python3
"""
Compare NS1 string A vs string B German end blowdown.
Produces:
  - ns1_german_outlet_30min.png    (2x2: mass flow, Mach, T, P for first 30 min)
  - ns1_german_pressure_shift.png  (closed-end pressure + time-shift difference)

Usage:
  python3 plot_ns1_german_compare.py
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os, sys

def load_timeseries(pfx, allow_lf=False):
    suffixes = ['', '_lf'] if allow_lf else ['']
    for suf in suffixes:
        fname = f'timeseries_{pfx}{suf}.dat'
        if not os.path.exists(fname):
            continue
        lines = open(fname).readlines()
        data = []
        for l in lines:
            if l.startswith('#'):
                continue
            parts = l.split()
            if len(parts) >= 13:
                try:
                    data.append([float(x) for x in parts[:13]])
                except ValueError:
                    continue
        if len(data) > 10:
            print(f'Loaded {fname}: {len(data)} points, '
                  f't = {data[0][0]:.1f} to {data[-1][0]:.1f} s')
            return np.array(data)
    print(f'WARNING: no data found for prefix {pfx}')
    return None

# --- Load data ---
# Use hires (1s resolution, first 2000s) for outlet/pressure plots
# Fall back to dx1 (full blowdown) if hires not available
def load_prefer(pfx1, pfx2, min_time=1500):
    ds = load_timeseries(pfx1)
    if ds is not None and ds[-1, 0] >= min_time:
        return ds
    ds2 = load_timeseries(pfx2)
    if ds2 is not None and ds2[-1, 0] >= min_time:
        return ds2
    return ds if ds is not None else ds2

ds_a = load_prefer('ns1ag_hires', 'ns1ag_dx1')
ds_b = load_prefer('ns1bg_hires', 'ns1bg_dx1')
# Also load full blowdown data if available (for summary)
ds_a_full = load_prefer('ns1ag_full', 'ns1ag_dx1')
ds_b_full = load_prefer('ns1bg_full', 'ns1bg_dx1')
# Load fine-sampled peak data (0.01s dt, captures tube-mode peak at t~0.02s)
ds_a_peak = load_timeseries('ns1ag_peak')
ds_b_peak = load_timeseries('ns1bg_peak')

if ds_a is None or ds_b is None:
    print('ERROR: need both ns1a and ns1b timeseries files')
    sys.exit(1)

# Column indices
COL_T    = 0
COL_MF   = 1   # mass flow [kg/s]
COL_U    = 2   # outlet velocity
COL_RHO  = 3   # outlet density
COL_P    = 4   # outlet pressure
COL_TOUT = 5   # outlet temperature
COL_MTOT = 6   # total mass
COL_PB   = 7   # closed-end (base) pressure
COL_PMID = 8   # mid-pipe pressure
COL_A    = 9   # sound speed at outlet
COL_MACH = 10  # Mach at outlet
COL_PCHK = 11  # choked outlet pressure
COL_CHK  = 12  # choked flag

runs = [
    (ds_a, 'NS1A (224.0 km)', 'C0'),
    (ds_b, 'NS1B (217.7 km)', 'C3'),
]

# === Figure 1: Outlet quantities for first 30 minutes ===
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('NS1 German end: String A (224 km) vs String B (217.7 km)\n'
             'First 30 minutes — HLLC, dx=1m, tube-start 2s', fontsize=13)

t_max_min = 30.0

for ds, lab, col in runs:
    t_min = ds[:, COL_T] / 60.0
    mask = (t_min >= 0.05) & (t_min <= t_max_min)

    axes[0, 0].plot(t_min[mask], ds[mask, COL_MF] / 1000, color=col, linewidth=1.2, label=lab)
    axes[0, 1].plot(t_min[mask], ds[mask, COL_MACH], color=col, linewidth=1.2, label=lab)
    axes[1, 0].plot(t_min[mask], ds[mask, COL_TOUT], color=col, linewidth=1.2, label=lab)
    axes[1, 1].plot(t_min[mask], ds[mask, COL_P] / 1e5, color=col, linewidth=1.2, label=lab)

for ax, title, ylabel in [
    (axes[0, 0], 'Outlet mass flow', 'Mass flow [t/s]'),
    (axes[0, 1], 'Outlet Mach number', 'Mach number'),
    (axes[1, 0], 'Outlet temperature', 'Temperature [K]'),
    (axes[1, 1], 'Outlet pressure', 'Pressure [bar]'),
]:
    ax.set_xlabel('Time [min]')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
axes[0, 1].axhline(y=1.0, color='k', linestyle='--', alpha=0.3, label='M=1')

plt.tight_layout()
fig.savefig('ns1_german_outlet_30min.png', dpi=150)
print(f'\nSaved: ns1_german_outlet_30min.png')
plt.close(fig)

# === Figure 2: Closed-end pressure + time-shift difference ===

# Find wave arrival times (interpolated to sub-sample accuracy)
def find_arrival(ds, threshold=0.1):
    """Find time when pressure drops by threshold [bar] using interpolation."""
    t = ds[:, COL_T]
    pb = ds[:, COL_PB] / 1e5
    drop = pb[0] - pb
    # Find bracket
    for i in range(1, len(drop)):
        if drop[i] >= threshold:
            # Linear interpolation within this interval
            frac = (threshold - drop[i-1]) / (drop[i] - drop[i-1])
            return t[i-1] + frac * (t[i] - t[i-1])
    return np.nan

t_arr_a = find_arrival(ds_a)
t_arr_b = find_arrival(ds_b)
print(f'NS1A wave arrival (0.1 bar drop): {t_arr_a:.1f} s ({t_arr_a/60:.2f} min)')
print(f'NS1B wave arrival (0.1 bar drop): {t_arr_b:.1f} s ({t_arr_b/60:.2f} min)')
print(f'Arrival time difference: {t_arr_a - t_arr_b:.1f} s')
print(f'CoolProp sound speed at 165 bar, 282 K: 480.8 m/s')
print(f'Expected dt from length difference: {6300/480.8:.1f} s')

# Build interpolated pressure functions for the first 30 min window
# Use data from wave arrival onward
def make_pb_interp(ds):
    t = ds[:, COL_T]
    pb = ds[:, COL_PB] / 1e5
    return interp1d(t, pb, kind='linear', bounds_error=False, fill_value=np.nan)

pb_a_fn = make_pb_interp(ds_a)
pb_b_fn = make_pb_interp(ds_b)

# Compute time shift: for each pressure level in B, find when A reaches the same pressure
# Use the shifted time axis (relative to each wave arrival)
def compute_time_shift(ds_a, ds_b, t_arr_a, t_arr_b, t_window=1800):
    """Compute time shift: how many seconds A lags B at each pressure level."""
    t_a = ds_a[:, COL_T]
    pb_a = ds_a[:, COL_PB] / 1e5
    t_b = ds_b[:, COL_T]
    pb_b = ds_b[:, COL_PB] / 1e5

    # Use B as reference: for each point in B after wave arrival,
    # find when A reaches the same pressure
    mask_b = (t_b >= t_arr_b) & (t_b <= t_arr_b + t_window)
    t_b_rel = t_b[mask_b] - t_arr_b  # time since B's wave arrival
    pb_b_vals = pb_b[mask_b]

    # For A: build inverse function pressure -> time (after arrival)
    mask_a = (t_a >= t_arr_a) & (t_a <= t_arr_a + t_window + 300)
    t_a_after = t_a[mask_a] - t_arr_a
    pb_a_after = pb_a[mask_a]

    # Pressure is monotonically decreasing, so invert
    # Remove any non-monotonic points
    mono_mask = np.diff(pb_a_after) < 0
    mono_mask = np.concatenate([[True], mono_mask])
    t_a_mono = t_a_after[mono_mask]
    pb_a_mono = pb_a_after[mono_mask]

    # Interpolate: given pressure, find time in A
    t_of_p_a = interp1d(pb_a_mono[::-1], t_a_mono[::-1],
                        kind='linear', bounds_error=False, fill_value=np.nan)

    # For each B pressure, find the corresponding time in A
    t_a_at_pb = t_of_p_a(pb_b_vals)

    # Time shift: how much later A reaches same pressure (relative to wave arrivals)
    time_shift = t_a_at_pb - t_b_rel  # positive = A lags B

    return t_b_rel, pb_b_vals, time_shift

t_rel, pb_ref, dt_shift = compute_time_shift(ds_a, ds_b, t_arr_a, t_arr_b)

# === Figure 2: Pressure aligned to wave arrival ===
fig, ax = plt.subplots(1, 1, figsize=(10, 7))
fig.suptitle('NS1 German end: Closed-end pressure aligned to wave arrival\n'
             'String A (224 km) vs String B (217.7 km) — HLLC, dx=1m, 165 bar, λ=0.007',
             fontsize=13)
for ds_r, lab, col in runs:
    pb = ds_r[:, COL_PB] / 1e5
    t_arr = find_arrival(ds_r)
    if not np.isnan(t_arr):
        t_shift = (ds_r[:, COL_T] - t_arr) / 60
        mask = (t_shift >= -0.5) & (t_shift <= 30)
        ax.plot(t_shift[mask], pb[mask], color=col, linewidth=1.5,
                label=f'{lab} (arr={t_arr:.1f}s)')
ax.set_xlabel('Time since wave arrival [min]', fontsize=12)
ax.set_ylabel('Pressure [bar]', fontsize=12)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
plt.tight_layout()
fig.savefig('ns1_german_pressure_aligned.png', dpi=150)
print(f'\nSaved: ns1_german_pressure_aligned.png')
plt.close(fig)

# === Figure 3: Time difference (A lags B) ===
fig, ax = plt.subplots(1, 1, figsize=(10, 7))
fig.suptitle('NS1 German end: Pressure drop time difference (NS1A − NS1B)\n'
             'HLLC, dx=1m, 165 bar, λ=0.007', fontsize=13)

valid = ~np.isnan(dt_shift) & (t_rel <= 1800)
ax.plot(t_rel[valid] / 60, dt_shift[valid], color='C2', linewidth=2.0)

# Reference lines
expected_dt = 6300.0 / 480.8
ax.axhline(y=expected_dt, color='k', linestyle='--', alpha=0.5, linewidth=1.0,
           label=f'ΔL/a(165 bar) = {expected_dt:.1f}s  (CoolProp: a=480.8 m/s)')

# Gridlines at every 1 second of time shift
dt_valid = dt_shift[valid]
ymin = max(0, np.floor(np.nanmin(dt_valid)) - 1)
ymax = np.ceil(np.nanmax(dt_valid)) + 1
ax.set_yticks(np.arange(ymin, ymax + 1, 1))
ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
ax.grid(True, alpha=0.6, which='major')
ax.grid(True, alpha=0.2, which='minor')

ax.set_xlabel('Time since wave arrival at NS1B [min]', fontsize=12)
ax.set_ylabel('Time shift [s]  (NS1A lags NS1B)', fontsize=12)
ax.legend(fontsize=10, loc='upper left')
ax.set_xlim(-0.5, 26)

plt.tight_layout()
fig.savefig('ns1_german_pressure_shift.png', dpi=150)
print(f'\nSaved: ns1_german_pressure_shift.png')
plt.close(fig)

# === Figure 4: Full blowdown (if data available) ===
full_runs = [(ds_a_full if ds_a_full is not None else ds_a, 'NS1A (224.0 km)', 'C0'),
             (ds_b_full if ds_b_full is not None else ds_b, 'NS1B (217.7 km)', 'C3')]
if any(d[:, COL_T].max() > 10000 for d, _, _ in full_runs):
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    for ds_f, lab, col in full_runs:
        ax.plot(ds_f[:, COL_T] / 3600, ds_f[:, COL_PB] / 1e5, color=col, linewidth=1.2, label=lab)
    ax.axhline(y=8.0, color='k', linestyle='--', alpha=0.3, label='8 bar')
    ax.set_xlabel('Time [hours]')
    ax.set_ylabel('Pressure [bar]')
    ax.set_title('NS1 German end: Full blowdown — String A vs String B')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig('ns1_german_blowdown_full.png', dpi=150)
    print(f'\nSaved: ns1_german_blowdown_full.png')
    plt.close(fig)

# === Summary table ===
print('\n=== Summary ===')
hdr = f'{"Case":<20} {"Peak MF [t/s]":>14} {"Wave arr [s]":>14} {"To 8 bar [hr]":>14} {"Final P [bar]":>14}'
print(hdr)
summary_runs = [(ds_a, ds_a_full, ds_a_peak, 'NS1A (224.0 km)', 'C0'),
                (ds_b, ds_b_full, ds_b_peak, 'NS1B (217.7 km)', 'C3')]
for ds, ds_f, ds_pk, lab, col in summary_runs:
    # Peak mass flow: use fine-sampled peak data if available (captures tube-mode peak)
    ds_for_peak = ds_pk if ds_pk is not None else ds
    t_pk = ds_for_peak[:, COL_T]
    mf_pk = ds_for_peak[:, COL_MF] / 1000
    settled = np.where(t_pk > 0.01)[0]
    peak_mf = mf_pk[settled].max() if len(settled) else mf_pk.max()

    # Wave arrival (interpolated)
    t_arr = find_arrival(ds)

    # Time to 8 bar — use full blowdown data if available
    ds_8 = ds if ds_f is None else ds_f
    pb_8 = ds_8[:, COL_PB] / 1e5
    t_8_arr = ds_8[:, COL_T]
    idx_8 = np.where(pb_8 <= 8.01)[0]
    t_8 = t_8_arr[idx_8[0]] / 3600 if len(idx_8) else np.nan

    pb = ds_8[:, COL_PB] / 1e5
    p_final = pb[-1]
    t8_str = f'{t_8:>14.1f}' if not np.isnan(t_8) else f'{"N/A":>14}'
    print(f'{lab:<20} {peak_mf:>14.2f} {t_arr:>14.1f} {t8_str} {p_final:>14.1f}')
