#!/usr/bin/env python3
"""
Verify HLLC solver against analytical rarefaction fan solution for ideal gas.

For a choked exit (p_L >> p_R), the analytical solution at time t is a
centered rarefaction fan expanding from x=0 (outlet) into still gas at p_L.

In the rarefaction fan region (x_head < x < 0):
  u(x,t) = 2/(gamma+1) * (a_L + x/t)
  a(x,t) = a_L - (gamma-1)/2 * u(x,t)
  p(x,t) = p_L * (a/a_L)^(2*gamma/(gamma-1))
  rho(x,t) = rho_L * (a/a_L)^(2/(gamma-1))

Fan head: x_head = -a_L * t  (moves leftward at sound speed)
Fan tail: x_tail = u_exit * t - a_exit * t  (if subsonic exit)
  or: x_tail = 0 for choked exit (fan reaches outlet)

At choked exit (M=1):
  a_exit = a_L * 2/(gamma+1)  =>  a_exit/a_L = 2/(gamma+1)
  u_exit = a_exit
"""
import numpy as np
import matplotlib.pyplot as plt
import sys

R = 8314.46 / 16.71   # specific gas constant [J/(kg·K)]
CV = 1622.0            # specific heat at constant volume [J/(kg·K)]
gamma = 1.0 + R / CV   # GAMMA_I = 1.3068

# Load the profile data
fname = 'profiles_ideal_verify.dat'
lines = open(fname).readlines()

# Parse snapshots
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
            current['data'].append([float(x) for x in parts[:5]])

# Use the last snapshot (t≈5s)
snap = snapshots[-1]
t_snap = snap['t']
data = np.array(snap['data'])
x = data[:, 0]
rho_num = data[:, 1]
u_num = data[:, 2]
p_num = data[:, 3]
T_num = data[:, 4]

print(f"Using snapshot at t = {t_snap:.4f} s")

# R already defined above as 8314.46/16.71 = 497.57

# Initial conditions
p_L = 165e5  # Pa
T_L = 282.0  # K
rho_L = p_L / (R * T_L)  # ideal gas
a_L = np.sqrt(gamma * R * T_L)

p_R = 7.61e5
rho_R = p_R / (R * T_L)

print(f"rho_L = {rho_L:.2f} kg/m³")
print(f"a_L = {a_L:.1f} m/s")
print(f"rho_R = {rho_R:.2f} kg/m³")

# Check if exit is choked
# At sonic exit: a_star = a_L * 2/(gamma+1)
a_star = a_L * 2.0 / (gamma + 1.0)
u_star = a_star  # M=1
p_star = p_L * (a_star / a_L) ** (2 * gamma / (gamma - 1))
rho_star = rho_L * (a_star / a_L) ** (2 / (gamma - 1))
print(f"\nSonic exit conditions:")
print(f"  a* = {a_star:.1f} m/s")
print(f"  u* = {u_star:.1f} m/s")
print(f"  p* = {p_star/1e5:.2f} bar")
print(f"  rho* = {rho_star:.2f} kg/m³")
print(f"  p* > p_R = {p_star/1e5:.2f} > {p_R/1e5:.2f}: choked = {p_star > p_R}")

# Analytical solution
x_anal = np.linspace(-4400, 0, 10000)
rho_anal = np.full_like(x_anal, rho_L)
u_anal = np.zeros_like(x_anal)
p_anal = np.full_like(x_anal, p_L)
a_anal = np.full_like(x_anal, a_L)

# Fan head position
x_head = -a_L * t_snap
print(f"\nFan head at x = {x_head:.1f} m (t={t_snap:.2f}s)")

for i, xi in enumerate(x_anal):
    if xi < x_head:
        # Undisturbed left state
        pass
    elif xi <= 0:
        # Inside rarefaction fan
        u_val = 2.0 / (gamma + 1.0) * (a_L + xi / t_snap)
        a_val = a_L - (gamma - 1.0) / 2.0 * u_val
        if a_val < 0:
            a_val = 0
        u_anal[i] = u_val
        a_anal[i] = a_val
        rho_anal[i] = rho_L * (a_val / a_L) ** (2.0 / (gamma - 1.0))
        p_anal[i] = p_L * (a_val / a_L) ** (2.0 * gamma / (gamma - 1.0))

# Compute numerical sound speed and Mach
a_num = np.sqrt(gamma * p_num / rho_num)
mach_num = u_num / a_num
mach_anal = u_anal / np.where(a_anal > 0, a_anal, 1)

# Mass flow (assuming area = pi/4 * dia^2, dia = 1.153m)
dia = 1.153
area = 0.25 * np.pi * dia**2
mf_num = rho_num * u_num * area / 1000  # t/s

# Analytical mass flow at exit
mf_anal_exit = rho_star * u_star * area / 1000
print(f"  Analytical choked mass flow: {mf_anal_exit:.2f} t/s")

# Find numerical outlet values
i_out = -2  # last physical cell (before ghost)
print(f"  Numerical outlet (cell N-2, x={x[i_out]:.2f}m):")
print(f"    u = {u_num[i_out]:.1f} m/s, a = {a_num[i_out]:.1f} m/s, M = {mach_num[i_out]:.4f}")
print(f"    p = {p_num[i_out]/1e5:.2f} bar, rho = {rho_num[i_out]:.2f} kg/m³")
print(f"    mass flow = {mf_num[i_out]:.2f} t/s")

# Plot comparison
fig, axes = plt.subplots(2, 3, figsize=(16, 9))
fig.suptitle(f'Ideal gas verification: HLLC vs analytical (t = {t_snap:.2f} s)', fontsize=14)

# Zoom to last 2500m (where the fan is)
xlim = (-2700, 50)

# Pressure
ax = axes[0, 0]
ax.plot(x, p_num/1e5, 'b-', linewidth=1, label='HLLC')
ax.plot(x_anal, p_anal/1e5, 'r--', linewidth=1, label='Analytical')
ax.set_xlabel('x [m]')
ax.set_ylabel('Pressure [bar]')
ax.set_title('Pressure')
ax.set_xlim(xlim)
ax.legend()
ax.grid(True, alpha=0.3)

# Velocity
ax = axes[0, 1]
ax.plot(x, u_num, 'b-', linewidth=1, label='HLLC')
ax.plot(x_anal, u_anal, 'r--', linewidth=1, label='Analytical')
ax.set_xlabel('x [m]')
ax.set_ylabel('Velocity [m/s]')
ax.set_title('Gas velocity')
ax.set_xlim(xlim)
ax.legend()
ax.grid(True, alpha=0.3)

# Density
ax = axes[0, 2]
ax.plot(x, rho_num, 'b-', linewidth=1, label='HLLC')
ax.plot(x_anal, rho_anal, 'r--', linewidth=1, label='Analytical')
ax.set_xlabel('x [m]')
ax.set_ylabel('Density [kg/m³]')
ax.set_title('Density')
ax.set_xlim(xlim)
ax.legend()
ax.grid(True, alpha=0.3)

# Mach number
ax = axes[1, 0]
ax.plot(x, mach_num, 'b-', linewidth=1, label='HLLC')
ax.plot(x_anal, mach_anal, 'r--', linewidth=1, label='Analytical')
ax.axhline(y=1.0, color='k', linestyle=':', alpha=0.5)
ax.set_xlabel('x [m]')
ax.set_ylabel('Mach number')
ax.set_title('Mach number')
ax.set_xlim(xlim)
ax.legend()
ax.grid(True, alpha=0.3)

# Temperature
T_anal = p_anal / (rho_anal * R)
ax = axes[1, 1]
ax.plot(x, T_num, 'b-', linewidth=1, label='HLLC')
ax.plot(x_anal, T_anal, 'r--', linewidth=1, label='Analytical')
ax.set_xlabel('x [m]')
ax.set_ylabel('Temperature [K]')
ax.set_title('Temperature')
ax.set_xlim(xlim)
ax.legend()
ax.grid(True, alpha=0.3)

# Zoom on last 100m — velocity and Mach
ax = axes[1, 2]
mask_n = x > -100
mask_a = x_anal > -100
ax.plot(x[mask_n], mach_num[mask_n], 'b-', linewidth=1, label='HLLC Mach')
ax.plot(x_anal[mask_a], mach_anal[mask_a], 'r--', linewidth=1, label='Analytical Mach')
ax.axhline(y=1.0, color='k', linestyle=':', alpha=0.5, label='M=1')
ax.set_xlabel('x [m]')
ax.set_ylabel('Mach number')
ax.set_title('Mach near outlet (last 100m)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('ideal_gas_verification.png', dpi=150)
print(f"\nSaved: ideal_gas_verification.png")
plt.close()
