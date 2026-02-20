#!/usr/bin/env python3
"""Plot closed-end pressure drop for German pipeline ends, zoomed to first 3 bar."""
import numpy as np
import matplotlib.pyplot as plt

def safe_load(fname):
    """Load timeseries, skipping incomplete lines."""
    lines = open(fname).readlines()
    data = []
    for l in lines:
        if l.startswith('#'): continue
        parts = l.split()
        if len(parts) == 13:
            try:
                data.append([float(x) for x in parts])
            except ValueError:
                continue
    return np.array(data)

configs = [
    ('timeseries_ns1ag_ghost_lf.dat', 'NS1 String A (224 km, 165 bar)', 'b'),
    ('timeseries_ns1bg_ghost_lf.dat', 'NS1 String B (217.7 km, 165 bar)', 'r'),
    ('timeseries_ns2ag_ghost_lf.dat', 'NS2 String A (153.6 km, 103 bar)', 'g'),
]

# CoolProp reference sound speeds and predicted wave arrival
coolprop_ref = {
    'NS1A': {'L': 224000, 'a': 478.78, 't_pred': 467.9},
    'NS1B': {'L': 217700, 'a': 478.78, 't_pred': 454.7},
    'NS2A': {'L': 153600, 'a': 424.21, 't_pred': 362.1},
}

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Panel 1: Pressure at closed end vs time (zoomed to first 3 bar)
ax = axes[0]
arrival_times = {}
for fname, label, color in configs:
    d = safe_load(fname)
    t = d[:, 0]
    pb = d[:, 7] / 1e5  # p_base in bar
    ax.plot(t / 60.0, pb, color=color, linewidth=1, label=label.split('(')[0].strip())

    # Find wave arrival
    drop = pb[0] - pb
    idx = np.where(drop > 0.1)[0]
    if len(idx) > 0:
        t_arr = t[idx[0]]
        arrival_times[label] = t_arr
        ax.axvline(x=t_arr/60, color=color, linestyle=':', alpha=0.5)

ax.set_xlabel('Time [min]')
ax.set_ylabel('Closed-end pressure [bar]')
ax.set_title('Pressure at closed end')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Panel 2: First 3 bar of drop — zoom on pressure axis
ax = axes[1]
for fname, label, color in configs:
    d = safe_load(fname)
    t = d[:, 0]
    pb = d[:, 7] / 1e5
    p0 = pb[0]
    short_label = label.split('(')[0].strip()

    # Shift time so t=0 is wave arrival
    drop = p0 - pb
    idx = np.where(drop > 0.05)[0]
    if len(idx) > 0:
        t_shift = (t - t[idx[0]]) / 60.0  # minutes since arrival
        mask = (t_shift >= -0.5) & (t_shift <= 5.0)
        ax.plot(t_shift[mask], (p0 - pb[mask]), color=color, linewidth=1.2, label=short_label)

ax.set_xlabel('Time since wave arrival [min]')
ax.set_ylabel('Pressure drop from initial [bar]')
ax.set_title('First 3 bar of pressure drop')
ax.set_ylim(-0.2, 4.0)
ax.axhline(y=3.0, color='k', linestyle='--', alpha=0.3, label='3 bar')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Panel 3: Wave arrival time comparison with CoolProp
ax = axes[2]
labels = ['NS1A\n224 km', 'NS1B\n217.7 km', 'NS2A\n153.6 km']
fnames = [c[0] for c in configs]
colors = [c[2] for c in configs]
refs = ['NS1A', 'NS1B', 'NS2A']
x_pos = np.arange(3)

sim_arrivals = []
pred_arrivals = []
for i, (fname, ref_key) in enumerate(zip(fnames, refs)):
    d = safe_load(fname)
    t = d[:, 0]
    pb = d[:, 7] / 1e5
    drop = pb[0] - pb
    idx = np.where(drop > 0.1)[0]
    t_sim = t[idx[0]] if len(idx) > 0 else np.nan
    sim_arrivals.append(t_sim)
    pred_arrivals.append(coolprop_ref[ref_key]['t_pred'])

bar_w = 0.35
bars1 = ax.bar(x_pos - bar_w/2, sim_arrivals, bar_w, label='Simulation', color=[c for c in colors], alpha=0.7)
bars2 = ax.bar(x_pos + bar_w/2, pred_arrivals, bar_w, label='CoolProp L/a', color=[c for c in colors], alpha=0.3, edgecolor='k', linewidth=1)

# Add time labels on bars
for bar, val in zip(bars1, sim_arrivals):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
            f'{val:.0f}s', ha='center', va='bottom', fontsize=8)
for bar, val in zip(bars2, pred_arrivals):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
            f'{val:.0f}s', ha='center', va='bottom', fontsize=8)

ax.set_xticks(x_pos)
ax.set_xticklabels(labels, fontsize=9)
ax.set_ylabel('Wave arrival time [s]')
ax.set_title('Wave arrival: simulation vs CoolProp')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig('german_pbase_drop.png', dpi=150)
print("Saved: german_pbase_drop.png")

# Print summary
print("\nWave arrival times:")
print(f"  NS1A (224.0 km): sim={sim_arrivals[0]:.1f}s, CoolProp={pred_arrivals[0]:.1f}s, diff={sim_arrivals[0]-pred_arrivals[0]:+.1f}s ({(sim_arrivals[0]-pred_arrivals[0])/pred_arrivals[0]*100:+.1f}%)")
print(f"  NS1B (217.7 km): sim={sim_arrivals[1]:.1f}s, CoolProp={pred_arrivals[1]:.1f}s, diff={sim_arrivals[1]-pred_arrivals[1]:+.1f}s ({(sim_arrivals[1]-pred_arrivals[1])/pred_arrivals[1]*100:+.1f}%)")
print(f"  NS2A (153.6 km): sim={sim_arrivals[2]:.1f}s, CoolProp={pred_arrivals[2]:.1f}s, diff={sim_arrivals[2]-pred_arrivals[2]:+.1f}s ({(sim_arrivals[2]-pred_arrivals[2])/pred_arrivals[2]*100:+.1f}%)")
