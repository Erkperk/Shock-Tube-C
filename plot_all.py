#!/usr/bin/env python3
"""
Generate all verification plots from a single simulation run:
  1. Outlet profiles at t≈0.5s (gas speed, Mach, pressure, temperature, mass flow)
  2. First 3 bar of pressure drop at the closed end
  3. Full pressure decline at the closed end
"""
import numpy as np
import matplotlib.pyplot as plt
import sys, glob

prefix = sys.argv[1] if len(sys.argv) > 1 else 'ns1ag_full'

# Try to find the files with _lf suffix first, then without
def find_file(base):
    for suf in ['_lf', '']:
        fname = f'{base}_{prefix}{suf}.dat'
        try:
            open(fname)
            return fname
        except FileNotFoundError:
            pass
    raise FileNotFoundError(f'Cannot find {base}_{prefix}[_lf].dat')

ts_file = find_file('timeseries')
early_file = find_file('early')

print(f'Timeseries: {ts_file}')
print(f'Early profiles: {early_file}')

# Gas constants
R = 8314.46 / 16.71  # J/(kg·K)
GAMMA = 1.27
dia = 1.153
area = 0.25 * np.pi * dia**2

# ======== Load timeseries ========
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

ts = load_timeseries(ts_file)
t = ts[:, 0]
mf_kg = ts[:, 1]       # mass flow [kg/s]
u_out = ts[:, 2]
rho_out = ts[:, 3]
p_out = ts[:, 4]
T_out = ts[:, 5]
totmass = ts[:, 6]
p_base = ts[:, 7]      # closed-end pressure [Pa]
p_mid = ts[:, 8]
a_out = ts[:, 9]
mach_out = ts[:, 10]

print(f'Timeseries: {len(t)} points, t = {t[0]:.1f} to {t[-1]:.1f} s')
print(f'p_base: {p_base[0]/1e5:.1f} to {p_base[-1]/1e5:.1f} bar')

# ======== Load early profiles ========
def load_early_profiles(fname):
    lines = open(fname).readlines()
    snapshots = []
    current = None
    for line in lines:
        if line.startswith('# t ='):
            t_val = float(line.split('=')[1])
            current = {'t': t_val, 'data': []}
            snapshots.append(current)
        elif not line.startswith('#') and line.strip() and current is not None:
            parts = line.split()
            if len(parts) >= 5:
                current['data'].append([float(x) for x in parts[:6]])
    return snapshots

snaps = load_early_profiles(early_file)
print(f'Early profiles: {len(snaps)} snapshots')
for s in snaps:
    print(f'  t = {s["t"]:.3f} s, {len(s["data"])} cells')

# Find snapshot closest to t=0.5s
target_t = 0.5
best = min(snaps, key=lambda s: abs(s['t'] - target_t))
t_snap = best['t']
data = np.array(best['data'])
x = data[:, 0]
rho_s = data[:, 1]
u_s = data[:, 2]
p_s = data[:, 3]
T_s = data[:, 4]

# Compute derived quantities
a_s = np.sqrt(GAMMA * p_s / rho_s)
mach_s = np.where(a_s > 0, u_s / a_s, 0)
mf_s = rho_s * u_s * area / 1000  # tonnes/s

print(f'\nUsing snapshot at t = {t_snap:.3f} s')
i_out = -2
print(f'Outlet: u={u_s[i_out]:.1f} m/s, M={mach_s[i_out]:.3f}, '
      f'p={p_s[i_out]/1e5:.1f} bar, T={T_s[i_out]:.1f} K, '
      f'mf={mf_s[i_out]:.2f} t/s')

# ================================================================
# FIGURE 1: Outlet profiles at t ≈ 0.5 s
# ================================================================
fig1, axes = plt.subplots(2, 3, figsize=(16, 9))
fig1.suptitle(f'Outlet profiles at t = {t_snap:.2f} s (last 250 m)', fontsize=14)

xlim = (-260, 10)

# Gas speed
ax = axes[0, 0]
ax.plot(x, u_s, 'b-', linewidth=1.2)
ax.set_xlabel('x [m]')
ax.set_ylabel('Gas velocity [m/s]')
ax.set_title('Gas speed')
ax.set_xlim(xlim)
ax.grid(True, alpha=0.3)

# Mach number
ax = axes[0, 1]
ax.plot(x, mach_s, 'b-', linewidth=1.2)
ax.axhline(y=1.0, color='r', linestyle=':', alpha=0.6, label='M = 1')
ax.set_xlabel('x [m]')
ax.set_ylabel('Mach number')
ax.set_title('Mach number')
ax.set_xlim(xlim)
ax.legend()
ax.grid(True, alpha=0.3)

# Pressure
ax = axes[0, 2]
ax.plot(x, p_s / 1e5, 'b-', linewidth=1.2)
ax.set_xlabel('x [m]')
ax.set_ylabel('Pressure [bar]')
ax.set_title('Pressure')
ax.set_xlim(xlim)
ax.grid(True, alpha=0.3)

# Temperature
ax = axes[1, 0]
ax.plot(x, T_s, 'b-', linewidth=1.2)
ax.set_xlabel('x [m]')
ax.set_ylabel('Temperature [K]')
ax.set_title('Temperature')
ax.set_xlim(xlim)
ax.grid(True, alpha=0.3)

# Mass flow
ax = axes[1, 1]
ax.plot(x, mf_s, 'b-', linewidth=1.2)
ax.set_xlabel('x [m]')
ax.set_ylabel('Mass flow [t/s]')
ax.set_title('Mass flow')
ax.set_xlim(xlim)
ax.grid(True, alpha=0.3)

# Density
ax = axes[1, 2]
ax.plot(x, rho_s, 'b-', linewidth=1.2)
ax.set_xlabel('x [m]')
ax.set_ylabel('Density [kg/m³]')
ax.set_title('Density')
ax.set_xlim(xlim)
ax.grid(True, alpha=0.3)

plt.tight_layout()
fig1.savefig('outlet_profiles.png', dpi=150)
print('\nSaved: outlet_profiles.png')
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
fig3.suptitle('Full blowdown: NS1A German (224 km, 165 bar → ambient)', fontsize=13)

# Panel 1: Pressure vs time
ax = axes3[0]
ax.plot(t / 3600, pb, 'b-', linewidth=1)
ax.set_xlabel('Time [hours]')
ax.set_ylabel('Closed-end pressure [bar]')
ax.set_title('Closed-end pressure')
ax.axhline(y=8.0, color='r', linestyle='--', alpha=0.4, label='8 bar')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Panel 2: Mass flow at outlet vs time
ax = axes3[1]
ax.plot(t / 3600, mf_kg / 1000, 'b-', linewidth=1)
ax.set_xlabel('Time [hours]')
ax.set_ylabel('Mass flow at outlet [t/s]')
ax.set_title('Outlet mass flow')
ax.grid(True, alpha=0.3)

plt.tight_layout()
fig3.savefig('full_blowdown.png', dpi=150)
print('Saved: full_blowdown.png')
plt.close(fig3)

print('\nDone — 3 PNG files generated.')
