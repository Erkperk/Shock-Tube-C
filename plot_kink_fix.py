#!/usr/bin/env python3
"""Compare old BC (discrete switch) vs new BC (smooth isentropic) around the kink."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def load_ts(fname):
    data = np.loadtxt(fname, comments='#')
    return {
        't': data[:, 0], 'massflow': data[:, 1], 'u_out': data[:, 2],
        'p_out': data[:, 4], 'T_out': data[:, 5],
        'p_base': data[:, 7], 'p_mid': data[:, 8],
        'a_out': data[:, 9], 'Mach_out': data[:, 10],
    }

# Old runs (discrete choked BC)
d_old_fric = load_ts('timeseries_ideal_fric.dat')
# New runs (smooth isentropic BC)
d_new_fric = load_ts('timeseries_ideal_fric_v2.dat')
d_new_nosrc = load_ts('timeseries_ideal_choke_v2.dat')

CV    = 1622.0
RGAS  = 8314.46 / 16.71
GAMMA_I = 1.0 + RGAS / CV
L     = 217700.0
a0    = np.sqrt(GAMMA_I * RGAS * 282.0)
t_arrive = L / a0

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('Outlet BC fix: smooth choked→subsonic transition (NS1B, 217.7 km)', fontsize=13)

c_old = '#d62728'
c_new = '#1f77b4'

# ── Panel 1: Mach number zoom around kink ──
ax = axes[0, 0]
mask_old = (d_old_fric['t'] > 780) & (d_old_fric['t'] < 880)
mask_new = (d_new_fric['t'] > 780) & (d_new_fric['t'] < 880)
ax.plot(d_old_fric['t'][mask_old] / 60, d_old_fric['Mach_out'][mask_old],
        '-', color=c_old, lw=1.5, label='Old BC (discrete switch)')
ax.plot(d_new_fric['t'][mask_new] / 60, d_new_fric['Mach_out'][mask_new],
        '-', color=c_new, lw=1.5, label='New BC (smooth isentropic)')
ax.set_xlabel('Time [min]')
ax.set_ylabel('Outlet Mach number')
ax.set_title('Mach number: zoom on kink region')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# ── Panel 2: p_out zoom around kink ──
ax = axes[0, 1]
ax.plot(d_old_fric['t'][mask_old] / 60, d_old_fric['p_out'][mask_old] / 1e5,
        '-', color=c_old, lw=1.5, label='Old BC')
ax.plot(d_new_fric['t'][mask_new] / 60, d_new_fric['p_out'][mask_new] / 1e5,
        '-', color=c_new, lw=1.5, label='New BC')
ax.axhline(8.1, color='gray', ls=':', lw=0.7, label='PRESSURE_R (8.1 bar)')
ax.set_xlabel('Time [min]')
ax.set_ylabel('Outlet pressure [bar]')
ax.set_title('Outlet pressure: zoom on kink region')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# ── Panel 3: Full Mach history (friction) ──
ax = axes[1, 0]
mask_old2 = d_old_fric['t'] > t_arrive
mask_new2 = d_new_fric['t'] > t_arrive
ax.plot(d_old_fric['t'][mask_old2] / 60, d_old_fric['Mach_out'][mask_old2],
        '-', color=c_old, lw=1, alpha=0.7, label='Old BC')
ax.plot(d_new_fric['t'][mask_new2] / 60, d_new_fric['Mach_out'][mask_new2],
        '-', color=c_new, lw=1, alpha=0.7, label='New BC')
ax.set_xlabel('Time [min]')
ax.set_ylabel('Outlet Mach number')
ax.set_title('Full Mach history (with friction)')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# ── Panel 4: Full p_base decline (friction) ──
ax = axes[1, 1]
ax.plot(d_old_fric['t'] / 60, d_old_fric['p_base'] / 1e5,
        '-', color=c_old, lw=1, alpha=0.7, label='Old BC')
ax.plot(d_new_fric['t'] / 60, d_new_fric['p_base'] / 1e5,
        '-', color=c_new, lw=1, alpha=0.7, label='New BC')
ax.axhline(10, color='gray', ls=':', lw=0.7)
ax.set_xlabel('Time [min]')
ax.set_ylabel('Pressure [bar]')
ax.set_title('Pressure at closed end (p_base)')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_ylim(0, 170)

plt.tight_layout()
plt.savefig('kink_fix_comparison.png', dpi=150, bbox_inches='tight')
print("Saved: kink_fix_comparison.png")

# Also check the max Mach for both cases
print(f"\nOld BC: max Mach = {d_old_fric['Mach_out'][mask_old2].max():.6f}")
print(f"New BC: max Mach = {d_new_fric['Mach_out'][mask_new2].max():.6f}")
print(f"New BC no-src: max Mach = {d_new_nosrc['Mach_out'].max():.6f}")

# Check mass flow continuity around the kink
print("\n--- Mass flow around kink (old) ---")
t = d_old_fric['t']; mf = d_old_fric['massflow']
m = (t > 820) & (t < 824)
for i in np.where(m)[0][::5]:
    print(f"  t={t[i]:.1f}  mf={mf[i]:.2f}")

print("\n--- Mass flow around kink (new) ---")
t = d_new_fric['t']; mf = d_new_fric['massflow']
m = (t > 820) & (t < 824)
for i in np.where(m)[0][::5]:
    print(f"  t={t[i]:.1f}  mf={mf[i]:.2f}")
