#!/usr/bin/env python3
"""Plot NS1A vs NS1B German-end pressure decline with tight grid."""
import numpy as np
import matplotlib.pyplot as plt

ns1a = np.loadtxt('timeseries_ns1a_lf_lf.dat')
ns1b = np.loadtxt('timeseries_ns1b_lf_lf.dat')

t_1a_hr = ns1a[:, 0] / 3600.0
t_1b_hr = ns1b[:, 0] / 3600.0
p_base_1a = ns1a[:, 7] / 1e5
p_base_1b = ns1b[:, 7] / 1e5
p_out_1a = ns1a[:, 4] / 1e5
p_out_1b = ns1b[:, 4] / 1e5

fig, axes = plt.subplots(2, 1, figsize=(11, 8))

# --- Top: closed-end pressure ---
ax = axes[0]
ax.plot(t_1a_hr, p_base_1a, '#1f77b4', lw=1.5, label='NS1A (224 km, 8.61 bar ambient)')
ax.plot(t_1b_hr, p_base_1b, '#2ca02c', lw=1.5, label='NS1B (218 km, 8.23 bar ambient)')
ax.set_xlabel('Time [hours]')
ax.set_ylabel('Pressure [bar]')
ax.set_title('Closed-end (Lubmin) pressure')
ax.legend(fontsize=10)
ax.grid(True, which='major', alpha=0.5, linewidth=0.8)
ax.grid(True, which='minor', alpha=0.25, linewidth=0.5)
ax.xaxis.set_major_locator(plt.MultipleLocator(1))
ax.xaxis.set_minor_locator(plt.MultipleLocator(0.25))
ax.yaxis.set_major_locator(plt.MultipleLocator(10))
ax.yaxis.set_minor_locator(plt.MultipleLocator(5))
ax.set_xlim(0, 14)
ax.set_ylim(0, 170)

# --- Bottom: outlet pressure at rupture ---
ax = axes[1]
ax.plot(t_1a_hr, p_out_1a, '#1f77b4', lw=1.5, label='NS1A (224 km)')
ax.plot(t_1b_hr, p_out_1b, '#2ca02c', lw=1.5, label='NS1B (218 km)')
ax.set_xlabel('Time [hours]')
ax.set_ylabel('Pressure [bar]')
ax.set_title('Rupture point (outlet) pressure')
ax.legend(fontsize=10)
ax.grid(True, which='major', alpha=0.5, linewidth=0.8)
ax.grid(True, which='minor', alpha=0.25, linewidth=0.5)
ax.xaxis.set_major_locator(plt.MultipleLocator(1))
ax.xaxis.set_minor_locator(plt.MultipleLocator(0.25))
ax.yaxis.set_major_locator(plt.MultipleLocator(5))
ax.yaxis.set_minor_locator(plt.MultipleLocator(1))
ax.set_xlim(0, 14)
ax.set_ylim(0, 40)

plt.suptitle('NS1A vs NS1B — German end (Lubmin)', fontsize=14, fontweight='bold', y=0.98)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('ns1_german_comparison.png', dpi=150)
plt.show()
print("Saved: ns1_german_comparison.png")
