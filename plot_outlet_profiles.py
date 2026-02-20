#!/usr/bin/env python3
"""Plot outlet profiles at t~0.5s: gas speed, Mach, pressure, temperature, mass flow."""
import numpy as np
import matplotlib.pyplot as plt

# Load data — extract t=0.5s snapshot
lines = []
in_block = False
with open('early_outlet_05s.dat') as f:
    for line in f:
        if line.startswith('# t = 5.00'):
            in_block = True
            continue
        if in_block and line.startswith('# t ='):
            break
        if in_block and not line.startswith('#') and line.strip():
            lines.append(line)

data = np.array([[float(v) for v in l.split()] for l in lines])
x = data[:, 0]       # position [m]
rho = data[:, 1]     # density [kg/m³]
u = data[:, 2]       # velocity [m/s]
p = data[:, 3]       # pressure [Pa]
T = data[:, 4]       # temperature [K]
e_int = data[:, 5]   # specific internal energy [J/kg]

# Filter to last 250m
mask = x > -250.0
x_p = x[mask]
rho_p = rho[mask]
u_p = u[mask]
p_p = p[mask]
T_p = T[mask]

# Derived quantities
gamma = 1.3
a_p = np.sqrt(gamma * p_p / rho_p)  # sound speed [m/s]
mach_p = u_p / a_p                    # Mach number
dia = 1.153  # NS1 diameter [m]
area = 0.25 * np.pi * dia**2
mflow_p = rho_p * u_p * area / 1000.0  # mass flow [t/s]

fig, axes = plt.subplots(2, 3, figsize=(15, 8))
fig.suptitle('Outlet profiles at t ≈ 0.5 s (last 250 m)', fontsize=14)

# Gas speed
ax = axes[0, 0]
ax.plot(x_p, u_p, 'b-', linewidth=1)
ax.set_xlabel('Position [m]')
ax.set_ylabel('Gas velocity [m/s]')
ax.set_title('Gas speed')
ax.grid(True, alpha=0.3)

# Mach number
ax = axes[0, 1]
ax.plot(x_p, mach_p, 'r-', linewidth=1)
ax.axhline(y=1.0, color='k', linestyle='--', alpha=0.5, label='M=1')
ax.set_xlabel('Position [m]')
ax.set_ylabel('Mach number')
ax.set_title('Mach number')
ax.legend()
ax.grid(True, alpha=0.3)

# Pressure
ax = axes[0, 2]
ax.plot(x_p, p_p / 1e5, 'g-', linewidth=1)
ax.set_xlabel('Position [m]')
ax.set_ylabel('Pressure [bar]')
ax.set_title('Pressure')
ax.grid(True, alpha=0.3)

# Temperature
ax = axes[1, 0]
ax.plot(x_p, T_p, 'm-', linewidth=1)
ax.set_xlabel('Position [m]')
ax.set_ylabel('Temperature [K]')
ax.set_title('Temperature')
ax.grid(True, alpha=0.3)

# Mass flow
ax = axes[1, 1]
ax.plot(x_p, mflow_p, 'c-', linewidth=1)
ax.set_xlabel('Position [m]')
ax.set_ylabel('Mass flow [t/s]')
ax.set_title('Mass flow rate')
ax.grid(True, alpha=0.3)

# Density
ax = axes[1, 2]
ax.plot(x_p, rho_p, 'orange', linewidth=1)
ax.set_xlabel('Position [m]')
ax.set_ylabel('Density [kg/m³]')
ax.set_title('Density')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('outlet_profiles_05s.png', dpi=150)
print("Saved: outlet_profiles_05s.png")
plt.close()

# Print key values at outlet (last cell before ghost)
i_out = -2  # second to last (last physical cell)
print(f"\nOutlet values at t≈0.5s (x={x[i_out]:.1f}m):")
print(f"  Velocity:    {u[i_out]:.1f} m/s")
print(f"  Sound speed: {np.sqrt(gamma * p[i_out] / rho[i_out]):.1f} m/s")
print(f"  Mach:        {u[i_out] / np.sqrt(gamma * p[i_out] / rho[i_out]):.3f}")
print(f"  Pressure:    {p[i_out]/1e5:.1f} bar")
print(f"  Temperature: {T[i_out]:.1f} K")
print(f"  Density:     {rho[i_out]:.1f} kg/m³")
print(f"  Mass flow:   {rho[i_out] * u[i_out] * area / 1000:.2f} t/s")
