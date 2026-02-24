#!/usr/bin/env python3
"""
Generate all verification plots from a simulation run:
  1. Outlet time histories for the first 0.5 s (mass flow, pressure, Mach, etc.)
  2. First 3 bar of pressure drop at the closed end
  3. Full pressure decline at the closed end

Usage:
  python3 plot_all.py <blowdown_prefix> [--hllc <early_prefix>]

  If --hllc is given, figure 1 uses that HLLC timeseries for accurate early data.
  Otherwise, figure 1 uses the blowdown timeseries.
"""
import numpy as np
import matplotlib.pyplot as plt
import sys

# Parse arguments
prefix = 'ns1ag_full'
hllc_prefix = None
args = sys.argv[1:]
i = 0
while i < len(args):
    if args[i] == '--hllc' and i + 1 < len(args):
        hllc_prefix = args[i + 1]
        i += 2
    else:
        prefix = args[i]
        i += 1

def find_file(base, pfx):
    for suf in ['_lf', '']:
        fname = f'{base}_{pfx}{suf}.dat'
        try:
            open(fname)
            return fname
        except FileNotFoundError:
            pass
    raise FileNotFoundError(f'Cannot find {base}_{pfx}[_lf].dat')

def load_timeseries(fname):
    lines = open(fname).readlines()
    data = []
    for l in lines:
        if l.startswith('#'): continue
        parts = l.split()
        if len(parts) >= 13:
            try:
                data.append([float(x) for x in parts[:13]])
            except ValueError:
                continue
    return np.array(data)

# Load blowdown timeseries
ts_file = find_file('timeseries', prefix)
print(f'Blowdown timeseries: {ts_file}')
ts = load_timeseries(ts_file)

# Load early timeseries (HLLC if available, else blowdown)
if hllc_prefix:
    hllc_file = find_file('timeseries', hllc_prefix)
    print(f'HLLC early timeseries: {hllc_file}')
    ts_early = load_timeseries(hllc_file)
else:
    ts_early = ts

t = ts[:, 0]
mf_kg = ts[:, 1]       # mass flow [kg/s]
p_base = ts[:, 7]      # closed-end pressure [Pa]
mf_ts = mf_kg / 1000   # tonnes/s

te = ts_early[:, 0]
mf_e = ts_early[:, 1] / 1000  # tonnes/s
u_e = ts_early[:, 2]
rho_e = ts_early[:, 3]
p_e = ts_early[:, 4]
T_e = ts_early[:, 5]
a_e = ts_early[:, 9]
mach_e = ts_early[:, 10]

print(f'Blowdown: {len(t)} points, t = {t[0]:.4f} to {t[-1]:.1f} s')
print(f'Early: {len(te)} points, t = {te[0]:.4f} to {te[-1]:.1f} s')
settled_idx = np.where(te > 0.005)[0]
if len(settled_idx) > 0:
    pk_i = settled_idx[mf_e[settled_idx].argmax()]
    print(f'Peak mass flow (early): {mf_e[pk_i]:.2f} t/s at t = {te[pk_i]:.4f} s')
else:
    print(f'Peak mass flow (early): {mf_e.max():.2f} t/s at t = {te[mf_e.argmax()]:.4f} s')
print(f'p_base: {p_base[0]/1e5:.1f} to {p_base[-1]/1e5:.1f} bar')

# ================================================================
# FIGURE 1: Outlet time histories for the first 0.5 s
# ================================================================
mask05 = te <= 0.5
t05 = te[mask05]
if len(t05) < 3:
    print('WARNING: fewer than 3 early points in first 0.5 s')
    mask05 = te <= 2.0
    t05 = te[mask05]

solver = 'HLLC' if hllc_prefix else ('LF' if '_lf' in ts_file else 'HLLC')

fig1, axes = plt.subplots(2, 3, figsize=(16, 9))
fig1.suptitle(f'Outlet quantities — first {t05[-1]:.2f} s ({solver}, dx=1m)', fontsize=14)

# Mass flow (hide ghost cell artifact at t < 0.005s)
ax = axes[0, 0]
mf_plot = mf_e[mask05].copy()
mf_plot[t05 < 0.005] = np.nan
ax.plot(t05, mf_plot, 'b-', linewidth=1.5)
# Skip first-timestep ghost cell artifact when finding peak
settled = np.where(te > 0.005)[0]
peak_i = settled[mf_e[settled].argmax()] if len(settled) > 0 else mf_e.argmax()
ax.axhline(y=mf_e[peak_i], color='r', linestyle=':', alpha=0.5)
ax.annotate(f'peak: {mf_e[peak_i]:.1f} t/s @ t={te[peak_i]:.3f}s',
            xy=(te[peak_i], mf_e[peak_i]), fontsize=9, color='r',
            xytext=(0.5, 0.85), textcoords='axes fraction',
            arrowprops=dict(arrowstyle='->', color='r', alpha=0.5))
ax.set_xlabel('Time [s]')
ax.set_ylabel('Mass flow [t/s]')
ax.set_title('Outlet mass flow')
ax.grid(True, alpha=0.3)

# Gas speed
ax = axes[0, 1]
ax.plot(t05, u_e[mask05], 'b-', linewidth=1.5)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Gas velocity [m/s]')
ax.set_title('Outlet gas speed')
ax.grid(True, alpha=0.3)

# Mach number
ax = axes[0, 2]
ax.plot(t05, mach_e[mask05], 'b-', linewidth=1.5)
ax.axhline(y=1.0, color='r', linestyle=':', alpha=0.6, label='M = 1')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Mach number')
ax.set_title('Outlet Mach number')
ax.legend()
ax.grid(True, alpha=0.3)

# Pressure
ax = axes[1, 0]
ax.plot(t05, p_e[mask05] / 1e5, 'b-', linewidth=1.5)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Pressure [bar]')
ax.set_title('Outlet pressure')
ax.grid(True, alpha=0.3)

# Temperature
ax = axes[1, 1]
ax.plot(t05, T_e[mask05], 'b-', linewidth=1.5)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Temperature [K]')
ax.set_title('Outlet temperature')
ax.grid(True, alpha=0.3)

# Density
ax = axes[1, 2]
ax.plot(t05, rho_e[mask05], 'b-', linewidth=1.5)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Density [kg/m³]')
ax.set_title('Outlet density')
ax.grid(True, alpha=0.3)

plt.tight_layout()
fig1.savefig('outlet_early.png', dpi=150)
print('\nSaved: outlet_early.png')
plt.close(fig1)

# ================================================================
# FIGURE 2: First 3 bar of pressure drop at closed end
# ================================================================
fig2, axes2 = plt.subplots(1, 2, figsize=(14, 5))
fig2.suptitle('Pressure drop at closed end (NS1A German, 224 km)', fontsize=13)

pb = p_base / 1e5  # bar
p0 = pb[0]
drop = p0 - pb

# Panel 1: Pressure vs time (full)
ax = axes2[0]
ax.plot(t / 60, pb, 'b-', linewidth=1)
ax.set_xlabel('Time [min]')
ax.set_ylabel('Closed-end pressure [bar]')
ax.set_title('Pressure at closed end')
ax.grid(True, alpha=0.3)

# Mark wave arrival (0.1 bar drop)
idx_arr = np.where(drop > 0.1)[0]
if len(idx_arr) > 0:
    t_arr = t[idx_arr[0]]
    ax.axvline(x=t_arr / 60, color='r', linestyle=':', alpha=0.5,
               label=f'Wave arrives: {t_arr:.0f} s')
    ax.legend(fontsize=9)

# Panel 2: Zoom on first 3 bar of drop
ax = axes2[1]
if len(idx_arr) > 0:
    t_shift = (t - t[idx_arr[0]]) / 60  # minutes since arrival
    mask = (t_shift >= -0.5) & (t_shift <= 10.0)
    ax.plot(t_shift[mask], drop[mask], 'b-', linewidth=1.2)
    ax.axhline(y=3.0, color='k', linestyle='--', alpha=0.3, label='3 bar')
    ax.set_xlabel('Time since wave arrival [min]')
    ax.set_ylabel('Pressure drop [bar]')
    ax.set_title('First 3 bar of pressure drop')
    ax.set_ylim(-0.2, 5.0)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
fig2.savefig('pbase_first3bar.png', dpi=150)
print('Saved: pbase_first3bar.png')
plt.close(fig2)

# ================================================================
# FIGURE 3: Full pressure decline to ~8 bar
# ================================================================
fig3, axes3 = plt.subplots(1, 2, figsize=(14, 5))
fig3.suptitle('Full blowdown: NS1A German (224 km, 165 bar)', fontsize=13)

# Panel 1: Pressure vs time
ax = axes3[0]
ax.plot(t / 3600, pb, 'b-', linewidth=1)
ax.set_xlabel('Time [hours]')
ax.set_ylabel('Closed-end pressure [bar]')
ax.set_title('Closed-end pressure')
ax.axhline(y=8.0, color='r', linestyle='--', alpha=0.4, label='8 bar')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Panel 2: Mass flow at outlet vs time (skip first-timestep artifact)
ax = axes3[1]
mf_clip = mf_ts.copy()
mf_clip[t < 0.005] = np.nan  # hide ghost cell artifact
ax.plot(t / 3600, mf_clip, 'b-', linewidth=1)
ax.set_xlabel('Time [hours]')
ax.set_ylabel('Mass flow at outlet [t/s]')
ax.set_title('Outlet mass flow')
ax.grid(True, alpha=0.3)

plt.tight_layout()
fig3.savefig('full_blowdown.png', dpi=150)
print('Saved: full_blowdown.png')
plt.close(fig3)

print('\nDone — 3 PNG files generated.')
