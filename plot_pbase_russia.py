#!/usr/bin/env python3
"""Plot full pressure decline at closed (Russian) end of NS1B (1006.3 km)."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

data = np.loadtxt('timeseries_cal_full_lf_lf.dat', comments='#')
t      = data[:, 0]
p_base = data[:, 7]

# Convert to hours and bar
t_h = t / 3600.0
p_bar = p_base / 1e5

print(f"Initial p_base = {p_bar[0]:.2f} bar")
print(f"Final   p_base = {p_bar[-1]:.2f} bar at t = {t_h[-1]:.1f} hours")

# Key milestones
for target in [150, 100, 50, 20, 10]:
    idx = np.argmax(p_bar < target)
    if p_bar[idx] < target:
        print(f"  p < {target:3d} bar at t = {t_h[idx]:6.1f} h")

fig, ax = plt.subplots(figsize=(14, 7))

ax.plot(t_h, p_bar, '-', color='#1f77b4', lw=1.5)

# Reference lines
ax.axhline(8.1, color='red', ls='--', lw=0.8, alpha=0.6, label='Ambient (8.1 bar)')
ax.axhline(165, color='gray', ls=':', lw=0.7, alpha=0.5, label='Initial (165 bar)')

ax.set_xlabel('Time [hours]', fontsize=13)
ax.set_ylabel('Pressure at closed (Russian) end [bar]', fontsize=13)
ax.set_title('NS1B Russian end (1006.3 km): full pressure decline\n'
             r'Real gas, Lax-Friedrich, $\lambda$ = 0.00858 (Colebrook)',
             fontsize=14)
ax.legend(fontsize=11, loc='upper right')

# Gridlines: 10h major / 5h minor for x; 20 bar major / 10 bar minor for y
ax.xaxis.set_major_locator(MultipleLocator(10))
ax.xaxis.set_minor_locator(MultipleLocator(5))
ax.yaxis.set_major_locator(MultipleLocator(20))
ax.yaxis.set_minor_locator(MultipleLocator(10))
ax.grid(True, which='major', alpha=0.4, linewidth=0.8)
ax.grid(True, which='minor', alpha=0.15, linewidth=0.5)

ax.set_xlim(0, 130)
ax.set_ylim(0, 175)

plt.tight_layout()
plt.savefig('pbase_russia_full.png', dpi=150, bbox_inches='tight')
print("\nSaved: pbase_russia_full.png")
