#!/usr/bin/env python3
"""Compare HLLC vs Lax-Friedrich for NS1B Russian end (1006.3 km) blowdown."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# Load both datasets (same grid: 657 cells, dx_init=150, dx_min=60m)
hllc = np.loadtxt('timeseries_ns1b_russia_hllc.dat', comments='#')
lf   = np.loadtxt('timeseries_ns1b_russia_lf_lf.dat', comments='#')

t_hllc = hllc[:, 0] / 3600.0
p_hllc = hllc[:, 7] / 1e5
m_hllc = hllc[:, 6]

t_lf = lf[:, 0] / 3600.0
p_lf = lf[:, 7] / 1e5
m_lf = lf[:, 6]

print(f"HLLC: final p_base = {p_hllc[-1]:.2f} bar, mass = {m_hllc[-1]:.0f} kg")
print(f"LF:   final p_base = {p_lf[-1]:.2f} bar, mass = {m_lf[-1]:.0f} kg")
print(f"Difference at 125h: {p_hllc[-1] - p_lf[-1]:.2f} bar")

# Key milestones
for label, t, p in [("HLLC", t_hllc, p_hllc), ("LF", t_lf, p_lf)]:
    for target in [100, 50, 20, 10]:
        idx = np.argmax(p < target)
        if p[idx] < target:
            print(f"  {label}: p < {target:3d} bar at t = {t[idx]:6.1f} h")

# --- Plot 1: p_base comparison ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 12))

ax1.plot(t_hllc, p_hllc, '-', color='#d62728', lw=1.5, label='MUSCL-HLLC (SSP-RK2)')
ax1.plot(t_lf, p_lf, '-', color='#1f77b4', lw=1.5, label='Lax-Friedrich')
ax1.axhline(8.1, color='gray', ls='--', lw=0.8, alpha=0.6, label='Ambient (8.1 bar)')

ax1.set_xlabel('Time [hours]', fontsize=13)
ax1.set_ylabel('Pressure at closed (Russian) end [bar]', fontsize=13)
ax1.set_title('NS1B Russian end (1006.3 km): HLLC vs Lax-Friedrich\n'
              '657 cells, dx_min=60 m, dx_max=2000 m, real gas, '
              r'$\lambda$ = 0.00858',
              fontsize=14)
ax1.legend(fontsize=11, loc='upper right')

ax1.xaxis.set_major_locator(MultipleLocator(10))
ax1.xaxis.set_minor_locator(MultipleLocator(5))
ax1.yaxis.set_major_locator(MultipleLocator(20))
ax1.yaxis.set_minor_locator(MultipleLocator(10))
ax1.grid(True, which='major', alpha=0.4, linewidth=0.8)
ax1.grid(True, which='minor', alpha=0.15, linewidth=0.5)
ax1.set_xlim(0, 130)
ax1.set_ylim(0, 175)

# --- Plot 2: difference ---
# Interpolate LF onto HLLC time grid for difference
p_lf_interp = np.interp(t_hllc, t_lf, p_lf)
dp = p_hllc - p_lf_interp

ax2.plot(t_hllc, dp, '-', color='#2ca02c', lw=1.5)
ax2.axhline(0, color='gray', ls='-', lw=0.5)

ax2.set_xlabel('Time [hours]', fontsize=13)
ax2.set_ylabel('p_base(HLLC) - p_base(LF) [bar]', fontsize=13)
ax2.set_title('Pressure difference: HLLC minus Lax-Friedrich', fontsize=14)

ax2.xaxis.set_major_locator(MultipleLocator(10))
ax2.xaxis.set_minor_locator(MultipleLocator(5))
ax2.yaxis.set_major_locator(MultipleLocator(0.5))
ax2.yaxis.set_minor_locator(MultipleLocator(0.25))
ax2.grid(True, which='major', alpha=0.4, linewidth=0.8)
ax2.grid(True, which='minor', alpha=0.15, linewidth=0.5)
ax2.set_xlim(0, 130)

plt.tight_layout()
plt.savefig('russia_hllc_vs_lf.png', dpi=150, bbox_inches='tight')
print("\nSaved: russia_hllc_vs_lf.png")
