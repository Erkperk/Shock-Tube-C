#!/usr/bin/env python3
"""
High-resolution plot of the first 2.5 seconds: choked outlet from t=0.
Compares three cases:
  1. Ideal gas, no friction (= analytical Riemann solution)
  2. Real gas, no friction
  3. Real gas, with friction (production)

Ghost cell outlet BC applied from t=0 (no vacuum extension).

Usage:
  python3 plot_ns1_german_first2s.py
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, sys

def load_ts(pfx):
    fname = f'timeseries_{pfx}.dat'
    if not os.path.exists(fname):
        return None
    data = []
    for l in open(fname):
        if l.startswith('#'): continue
        p = l.split()
        if len(p) >= 13:
            data.append([float(x) for x in p[:13]])
    if len(data) > 10:
        print(f'Loaded {fname}: {len(data)} points, dt ~ {(data[-1][0]-data[1][0])/(len(data)-2):.4f} s')
        return np.array(data)
    return None

COL_T, COL_MF, COL_U, COL_RHO, COL_P, COL_TOUT = 0, 1, 2, 3, 4, 5
COL_MTOT, COL_PB, COL_PMID, COL_A, COL_MACH, COL_PCHK, COL_CHK = 6, 7, 8, 9, 10, 11, 12

# Load three cases
ds_fric = load_ts('ns1ag_notv')
ds_nosrc = load_ts('ns1ag_nosrc')
ds_ideal = load_ts('ns1ag_ideal')

if ds_fric is None:
    print("Need at least ns1ag_notv (production run)")
    sys.exit(1)

# Analytical ideal gas solution
RGAS = 8314.46 / 16.38
CV = 1647.0
gam = 1.0 + RGAS / CV
a0 = (gam * RGAS * 282.0)**0.5
a_star = 2 * a0 / (gam + 1)
T_star = a_star**2 / (gam * RGAS)
p_star = 165e5 * (T_star / 282.0)**(gam / (gam - 1))
rho_star = p_star / (RGAS * T_star)
A = 3.14159265 / 4 * 1.153**2
mf_analytical = rho_star * a_star * A / 1000  # t/s

# Build runs list
runs = []
if ds_ideal is not None:
    runs.append((ds_ideal, f'Ideal gas, no friction (analytical: {mf_analytical:.2f} t/s)', 'C2', '--'))
if ds_nosrc is not None:
    runs.append((ds_nosrc, 'Real gas, no friction', 'C1', '-'))
runs.append((ds_fric, 'Real gas, with friction ($\\lambda$=0.007)', 'C0', '-'))

t_plot = 2.55

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('NS1 German end: First 2.5 seconds — choked outlet (M=1) from t=0\n'
             'HLLC, dx=1m, ghost cell BC, 165 bar, 1.153m dia', fontsize=13)

for ds, lab, col, ls in runs:
    t = ds[:, COL_T]
    mask = (t >= 0.005) & (t <= t_plot)
    axes[0, 0].plot(t[mask], ds[mask, COL_MF] / 1000, color=col, linewidth=1.5,
                    linestyle=ls, label=lab)
    axes[0, 1].plot(t[mask], ds[mask, COL_MACH], color=col, linewidth=1.5,
                    linestyle=ls, label=lab)
    axes[1, 0].plot(t[mask], ds[mask, COL_P] / 1e5, color=col, linewidth=1.5,
                    linestyle=ls, label=lab)
    axes[1, 1].plot(t[mask], ds[mask, COL_TOUT], color=col, linewidth=1.5,
                    linestyle=ls, label=lab)

# Analytical reference lines
axes[0, 0].axhline(y=mf_analytical, color='C2', linestyle=':', alpha=0.5, linewidth=0.8)
axes[1, 0].axhline(y=p_star/1e5, color='C2', linestyle=':', alpha=0.5, linewidth=0.8)
axes[1, 1].axhline(y=T_star, color='C2', linestyle=':', alpha=0.5, linewidth=0.8)

# Mass flow
ax = axes[0, 0]
ax.set_ylabel('Mass flow [t/s]')
ax.set_title('Outlet mass flow')
ax.legend(fontsize=9, loc='upper right')
ax.grid(True, alpha=0.3)
# Annotate analytical value
ax.annotate(f'Analytical: {mf_analytical:.2f} t/s',
            xy=(2.2, mf_analytical), fontsize=9, color='C2',
            textcoords='offset points', xytext=(-80, 8))

# Mach
ax = axes[0, 1]
ax.set_ylabel('Mach number')
ax.set_title('Outlet Mach number')
ax.axhline(y=1.0, color='k', linestyle=':', alpha=0.3)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_ylim(0.95, 1.05)

# Pressure
ax = axes[1, 0]
ax.set_ylabel('Pressure [bar]')
ax.set_title('Outlet pressure (sonic)')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Temperature
ax = axes[1, 1]
ax.set_ylabel('Temperature [K]')
ax.set_title('Outlet temperature (sonic)')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

for ax in axes.flat:
    ax.set_xlabel('Time [s]')
    ax.set_xlim(0, t_plot)

plt.tight_layout()
fig.savefig('ns1_german_first2s.png', dpi=150)
print(f'\nSaved: ns1_german_first2s.png')
plt.close(fig)

# Print summary
print('\n=== First 2.5 seconds summary ===')
print(f'Analytical (ideal gas, Riemann invariant): {mf_analytical:.3f} t/s (constant)')
print()
for ds, lab, col, ls in runs:
    t = ds[:, COL_T]
    mf = ds[:, COL_MF] / 1000
    # Settled peak (skip startup transient)
    mask = (t > 0.01) & (t < 0.5)
    i_pk = np.argmax(mf[mask])
    mf_peak = mf[mask][i_pk]
    # Values at key times
    vals = []
    for tc in [0.1, 0.5, 1.0, 2.0, 2.5]:
        i = np.argmin(np.abs(t - tc))
        vals.append(f'{mf[i]:.2f}')
    print(f'{lab}:')
    print(f'  Settled peak: {mf_peak:.3f} t/s at t = {t[mask][i_pk]*1000:.0f} ms')
    print(f'  At t=0.1/0.5/1.0/2.0/2.5s: {"/".join(vals)} t/s')
