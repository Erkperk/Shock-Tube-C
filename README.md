# Shock Tube — Subsea Gas Pipeline Blowdown Simulator

1D compressible Euler solver for simulating natural gas pipeline rupture and blowdown with real-gas thermodynamics. Developed for analysis of the Nord Stream pipeline incidents (September 2022).

Features:
- Real-gas equation of state for Yamal natural gas mixture (98% CH4)
- Two solver modes: MUSCL-HLLC (2nd-order) and Lax-Friedrich (1st-order)
- Darcy-Weisbach wall friction with Colebrook formula
- Transient wall heat transfer with radial conduction through steel and concrete
- Choked and subsonic outlet boundary conditions
- Adaptive grid merging for efficient long-duration blowdown simulations
- Pre-configured Nord Stream 1 and 2 pipeline presets

## Building

Requires a C99 compiler and the math library:

```
make
```

Or directly:

```
cc -O2 -std=c99 -o shock_tube shock_tube.c -lm
```

## Quick Start

```bash
# Short pipe, 20 seconds, default parameters
./shock_tube --prefix test

# Full NS1 string A blowdown from the Russian end (HLLC solver)
./shock_tube --ns1ar

# NS1 German end, Lax-Friedrich solver, stop when pressure reaches 8 bar
./shock_tube --ns1ag --lf --pstop 8e5

# Custom 1 km pipe at 100 bar
./shock_tube --length 1000 --pL 100e5 --prefix custom

# Print all available options
./shock_tube --help
```

## Physical Model

### Governing Equations

The 1D compressible Euler equations in conservative form:

```
d/dt [rho, rho*u, E] + d/dx [rho*u, rho*u^2+p, (E+p)*u] = [0, S_fric, S_heat]
```

where rho is density, u is velocity, E = rho*e + 0.5*rho*u^2 is total energy per unit volume, and p is pressure.

### Real-Gas Equation of State

The equation of state accounts for non-ideal gas behavior:

```
p = rho * Z(T, p) * R_s * T
```

where Z(T, p) is the compression factor and R_s = R/M is the specific gas constant (M = 16.38 g/mol for Yamal natural gas).

The compression factor uses the Dranchuk-Abou-Kassem (DAK) correlation with 19 coefficients, fitted for Yamal pipeline gas (gamma_g = 0.5654). This analytical formula provides Z as a function of reduced temperature T_r = T/T_pc and reduced pressure p_r = p/p_pc, where pseudo-critical properties are computed from the gas specific gravity.

### Thermodynamic Tables

Internal energy U(T, p) and isochoric heat capacity Cv(T, p) are tabulated from CoolProp for the Yamal gas mixture and interpolated using bicubic Catmull-Rom splines (C1-continuous). The tables cover 200-400 K in temperature and 1-200 bar in pressure.

### Sound Speed

The thermodynamic sound speed for a real gas is computed exactly using partial derivatives of the compression factor:

```
a^2 = (p / rho) * (1 + (rho/Z) * dZ/drho|_s) * (c_p / c_v)
```

The partial derivatives dZ/dT and dZ/dp are evaluated from the DAK formula using numerical central differences.

### Friction

Wall friction uses the Darcy-Weisbach formulation:

```
S_fric = -(lambda / D) * 0.5 * rho * u * |u|
```

where lambda is the Darcy friction factor computed from the Colebrook equation for fully-turbulent flow, or overridden via `--lambda`. The pipe diameter D and roughness epsilon are configurable.

### Wall Heat Transfer

The default wall model solves transient 1D radial heat conduction through the pipe wall (5 steel nodes + 3 concrete coating nodes). The inner boundary uses a convective condition with the Dittus-Boelter correlation for the heat transfer coefficient:

```
Nu = 0.023 * Re^0.8 * Pr^0.4
```

The outer boundary is held at the seawater temperature T0. A simplified steady-state model is available via `--simple-wall`.

## Numerical Methods

### MUSCL-HLLC (default)

Second-order spatial reconstruction using MUSCL (Monotone Upstream-centered Schemes for Conservation Laws) with the Van Leer flux limiter, combined with the HLLC (Harten-Lax-van Leer-Contact) approximate Riemann solver. Time integration uses the SSP-RK2 (Strong Stability Preserving Runge-Kutta, 2nd-order) method of Shu & Osher (1988).

Primitive variables are recovered at cell faces using 3 Newton iterations for robust real-gas inversion of conservative variables.

### Lax-Friedrich (--lf)

First-order Lax-Friedrich flux splitting with forward Euler time integration. More diffusive but faster per time step. Useful for quick parameter sweeps and calibration.

### Grid and Merging

The spatial grid supports:
- **Linear transition**: fine cells (dx_init) at the outlet, linearly ramping to coarser cells (dx_max) over a configurable transition length.
- **Adaptive merging**: neighboring cells are merged when their state is sufficiently uniform (controlled by `--merge-tol`), following a geometric schedule. Conservation is maintained to machine precision (~1e-15 relative mass error).
- **Active front tracking**: for large grids (>200 cells), only cells ahead of the expansion wave front are updated, significantly reducing computational cost.

### CFL Condition

The time step is determined by the CFL condition with CFL = 0.6:

```
dt = CFL * min_i(dx_i / (|u_i| + a_i))
```

## Boundary Conditions

- **Outlet (right end)**: choked flow when the isentropic sonic pressure exceeds the back-pressure, using a ghost cell set to M=1 conditions. Transitions smoothly to subsonic outflow using characteristic-based NRBC (non-reflecting boundary condition) with the Riemann invariant J+ = u + 2a/(gamma-1).
- **Closed end (left)**: reflective wall (u = 0).

## Pipeline Presets

| Preset    | Pipeline     | End     | Length    | Pressure | Endtime   |
|-----------|--------------|---------|-----------|----------|-----------|
| `--ns1ag` | NS1 string A | German  | 224 km    | 164 bar  | 50,000 s  |
| `--ns1ar` | NS1 string A | Russian | 999.8 km  | 164 bar  | 450,000 s |
| `--ns1bg` | NS1 string B | German  | 217.7 km  | 164 bar  | 50,000 s  |
| `--ns1br` | NS1 string B | Russian | 1006.2 km | 164 bar  | 450,000 s |
| `--ns2ag` | NS2 string A | German  | 153.6 km  | 103 bar  | 50,000 s  |
| `--ns2ar` | NS2 string A | Russian | 1076.4 km | 103 bar  | 450,000 s |

All presets enable adaptive merging and use a linear transition grid. NS1 presets use a calibrated Darcy friction factor of lambda = 0.007.

## Calibration

The Darcy friction factor for NS1 (lambda = 0.007) was calibrated against published pressure measurements at the Russian closed end (Portovaya compressor station) following the September 26, 2022 rupture of NS1 string A near Bornholm:

| Checkpoint      | Measured | Simulated (HLLC) |
|-----------------|----------|-------------------|
| Midnight Sep 28 | 82 bar   | 82.2 bar          |
| Midnight Sep 29 | 45 bar   | 44.9 bar          |

The initial pressure of 164 bar was also determined from the measurement data. Full blowdown to atmospheric back-pressure (8 bar) takes approximately 123 hours (5.1 days) in the simulation.

## Output Format

### Timeseries file (timeseries_*.dat)

Recorded at 5000 equally-spaced time points, plus high-resolution early recording if enabled.

| Column | Quantity                    | Unit  |
|--------|-----------------------------|-------|
| 1      | Time                        | s     |
| 2      | Outlet mass flow (rho*u*A)  | kg/s  |
| 3      | Outlet pressure             | Pa    |
| 4      | Outlet temperature          | K     |
| 5      | Outlet Mach number          | -     |
| 6      | Outlet sound speed          | m/s   |
| 7      | Closed-end pressure         | Pa    |
| 8      | Number of active cells      | -     |
| 9      | CFL time step               | s     |
| 10     | Choked outlet pressure      | Pa    |
| 11     | Outlet velocity             | m/s   |
| 12     | Closed-end temperature      | K     |
| 13     | Wall surface temperature    | K     |

### Profile files (profiles_*.dat)

Spatial snapshots at regular intervals, separated by blank lines:

| Column | Quantity    | Unit  |
|--------|-------------|-------|
| 1      | Position x  | m     |
| 2      | Density     | kg/m3 |
| 3      | Velocity    | m/s   |
| 4      | Pressure    | Pa    |
| 5      | Temperature | K     |
| 6      | Sound speed | m/s   |
| 7      | Mach number | -     |

## Thermodynamic Table Generation

The thermodynamic tables in `tables.h` are generated by `generate_tables.py` using CoolProp. The gas mixture is 98% methane / 2% ethane (molar), approximating Yamal pipeline gas:

```bash
python3 generate_tables.py
```

This produces tabulated values of U(T, p) and Cv(T, p) on a regular grid, embedded as C arrays in `tables.h`.

## References

- Dranchuk, P. M. & Abou-Kassem, J. H. (1975). Calculation of Z factors for natural gases using equations of state. *Journal of Canadian Petroleum Technology*, 14(3), 34-36.
- Toro, E. F. (2009). *Riemann Solvers and Numerical Methods for Fluid Dynamics* (3rd ed.). Springer.
- Shu, C.-W. & Osher, S. (1988). Efficient implementation of essentially non-oscillatory shock-capturing schemes. *Journal of Computational Physics*, 77(2), 439-471.
- Thompson, K. W. (1987). Time dependent boundary conditions for hyperbolic systems. *Journal of Computational Physics*, 68(1), 1-24.
- Bell, I. H., Wronski, J., Quoilin, S., & Lemort, V. (2014). Pure and pseudo-pure fluid thermophysical property evaluation and the open-source thermophysical property library CoolProp. *Industrial & Engineering Chemistry Research*, 53(6), 2498-2508.
- Colebrook, C. F. (1939). Turbulent flow in pipes, with particular reference to the transition region between the smooth and rough pipe laws. *Journal of the Institution of Civil Engineers*, 11(4), 133-156.
- Dittus, F. W. & Boelter, L. M. K. (1930). Heat transfer in automobile radiators of the tubular type. *University of California Publications in Engineering*, 2, 443-461.

## License

This code is provided for research and educational purposes.
