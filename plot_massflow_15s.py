#!/usr/bin/env python3
"""Plot outlet mass flow during first 15 seconds: NS1 (165 bar) vs NS2 (104 bar).
High-resolution comparison with real-gas choked flow references."""
import numpy as np
import matplotlib.pyplot as plt

# Load simulation data
ns1 = np.loadtxt('timeseries__peak165.dat', comments='#')
ns2 = np.loadtxt('timeseries__peak104.dat', comments='#')

t_ns1, mf_ns1 = ns1[:, 0], ns1[:, 1] / 1000  # t/s
t_ns2, mf_ns2 = ns2[:, 0], ns2[:, 1] / 1000

# Skip startup artifact (first ~10ms) — mass flow should start at settled value
t_cut = 0.01
mask1 = t_ns1 >= t_cut
mask2 = t_ns2 >= t_cut

# Settled peak (average over t = 0.02–0.1s, well past transient)
settled1 = mf_ns1[(t_ns1 > 0.02) & (t_ns1 < 0.1)]
settled2 = mf_ns2[(t_ns2 > 0.02) & (t_ns2 < 0.1)]
peak_ns1 = np.mean(settled1)
peak_ns2 = np.mean(settled2)

print(f'NS1 (165 bar): settled peak = {peak_ns1:.2f} t/s')
print(f'NS2 (104 bar): settled peak = {peak_ns2:.2f} t/s')

# --- Plot ---
fig, ax = plt.subplots(figsize=(16, 9))

ax.plot(t_ns1[mask1], mf_ns1[mask1], color='#1f77b4', linewidth=1.8,
        label=f'NS1 — 165 bar, $\\lambda$ = 0.007')
ax.plot(t_ns2[mask2], mf_ns2[mask2], color='#d62728', linewidth=1.8,
        label=f'NS2 — 104 bar, $\\lambda$ = 0.007')

# Real-gas choked peak reference lines
ax.axhline(peak_ns1, color='#1f77b4', ls='--', lw=1.0, alpha=0.6,
           label=f'Peak (real gas): {peak_ns1:.1f} t/s')
ax.axhline(peak_ns2, color='#d62728', ls='--', lw=1.0, alpha=0.6,
           label=f'Peak (real gas): {peak_ns2:.1f} t/s')

ax.set_xlabel('Time [s]', fontsize=14)
ax.set_ylabel('Mass flow at outlet [t/s]', fontsize=14)
ax.set_title('Choked outlet mass flow — first 15 seconds\n'
             'HLLC, dx = 1 m, ghost cell BC, D = 1.153 m, T = 282 K',
             fontsize=14)

# Gridlines every 1 second and 1 t/s
ax.set_xlim(0, 15)
y_max = int(np.ceil(peak_ns1)) + 1
ax.set_ylim(0, y_max)
ax.xaxis.set_major_locator(plt.MultipleLocator(1))
ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
ax.yaxis.set_major_locator(plt.MultipleLocator(1))
ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
ax.grid(True, which='major', alpha=0.5, linewidth=0.8)
ax.grid(True, which='minor', alpha=0.2, linewidth=0.4)

ax.legend(fontsize=12, loc='upper right')
ax.tick_params(labelsize=12)

plt.tight_layout()
fig.savefig('massflow_15s.png', dpi=200)
print(f'\nSaved: massflow_15s.png')
plt.close(fig)
