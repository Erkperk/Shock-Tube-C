#!/usr/bin/env python3
"""Plot first 2 seconds of outlet mass flow: NS1 (165 bar) vs NS2 (104 bar)."""
import numpy as np
import matplotlib.pyplot as plt

# Load timeseries data (col 0 = time, col 1 = mass flow kg/s)
ns1 = np.loadtxt('timeseries_ns1_peak.dat')
ns2 = np.loadtxt('timeseries_ns2_peak.dat')

# Filter startup artifact (t < 10ms)
ns1 = ns1[ns1[:, 0] > 0.01]
ns2 = ns2[ns2[:, 0] > 0.01]

# Settled peak = analytical choked mass flow (max after transient)
peak_ns1 = np.max(ns1[:50, 1]) / 1000  # first ~50 points after filter
peak_ns2 = np.max(ns2[:50, 1]) / 1000

fig, ax = plt.subplots(figsize=(14, 7))

# Plot mass flow
ax.plot(ns1[:, 0] * 1000, ns1[:, 1] / 1000, 'b-', linewidth=1.8,
        label=f'NS1 — 165 bar, $\\lambda$=0.007')
ax.plot(ns2[:, 0] * 1000, ns2[:, 1] / 1000, 'r-', linewidth=1.8,
        label=f'NS2 — 104 bar, $\\lambda$=0.007')

# Analytical choked mass flow reference lines
ax.axhline(peak_ns1, color='blue', linestyle='--', linewidth=1.0, alpha=0.7,
           label=f'Analytical peak: {peak_ns1:.1f} t/s')
ax.axhline(peak_ns2, color='red', linestyle='--', linewidth=1.0, alpha=0.7,
           label=f'Analytical peak: {peak_ns2:.1f} t/s')

ax.set_xlabel('Time [ms]', fontsize=14)
ax.set_ylabel('Mass flow at outlet [t/s]', fontsize=14)
ax.set_title('Outlet Mass Flow — First 2 Seconds\n(Full Wall Model, HLLC, dx = 1 m, D = 1.153 m)',
             fontsize=15)

ax.set_xlim(0, 2000)

# Major ticks every 50ms, labels every 200ms
ax.set_xticks(np.arange(0, 2050, 50))
ax.set_xticklabels([str(int(x)) if x % 200 == 0 else '' for x in np.arange(0, 2050, 50)])

# Gridlines every 50ms
ax.grid(True, which='major', axis='x', linewidth=0.5, alpha=0.5, color='gray')
ax.grid(True, which='major', axis='y', linewidth=0.5, alpha=0.5, color='gray')

ax.legend(fontsize=12, loc='upper right')
ax.tick_params(axis='both', labelsize=12)

plt.tight_layout()
plt.savefig('peak_massflow_compare.png', dpi=150)
print(f'Saved: peak_massflow_compare.png')
print(f'NS1 settled peak: {peak_ns1:.2f} t/s')
print(f'NS2 settled peak: {peak_ns2:.2f} t/s')
plt.close()
