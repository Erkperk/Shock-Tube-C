#!/usr/bin/env python3
"""
Choking analysis for ideal gas expansion from a closed-end tube.
Compares: analytical MOC solution vs HLLC (choked / unchoked).
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── Gas and pipe parameters (must match shock_tube.c) ──
CV    = 1622.0
CP    = 2125.0
RGAS  = 8314.46 / 16.71
GAMMA = CP / CV
L     = 217700.0      # NS1B pipe length [m]
P0    = 165e5          # initial pressure [Pa]
P_EXT = 8.1e5          # exterior pressure [Pa]
T0    = 282.0          # initial temperature [K]
AREA  = 1.153**2 / 4.0 * np.pi

rho0 = P0 / (RGAS * T0)
a0   = np.sqrt(GAMMA * RGAS * T0)
gm1  = GAMMA - 1
gp1  = GAMMA + 1

print(f"gamma  = {GAMMA:.5f}")
print(f"a0     = {a0:.2f} m/s")
print(f"rho0   = {rho0:.3f} kg/m3")
print(f"R_spec = {RGAS:.2f} J/(kg K)")
t_arrive = L / a0
print(f"t_arrive = L/a0 = {t_arrive:.1f} s = {t_arrive/60:.1f} min")

# ── Analytical MOC solver (vectorized, high resolution) ──
def moc_solve(Nx=10000, dt_cfl=0.45, t_end=8000.0, choked=True):
    """
    Gridded Method of Characteristics for 1D isentropic ideal gas.
    Closed wall at x=0, outlet at x=L.
    If choked: sonic outlet (M=1) while p_sonic > P_EXT.
    Uses vectorized numpy for speed.
    """
    dx = L / (Nx - 1)
    x = np.linspace(0, L, Nx)

    alpha0 = 2.0 * a0 / gm1  # characteristic scale

    # Riemann invariants: J+ = u + 2a/(gamma-1), J- = u - 2a/(gamma-1)
    Jp = np.full(Nx, alpha0)
    Jm = np.full(Nx, -alpha0)

    def get_state(Jp, Jm):
        u = 0.5 * (Jp + Jm)
        a = gm1 / 4.0 * (Jp - Jm)
        a = np.maximum(a, 1.0)
        p = P0 * (a / a0) ** (2 * GAMMA / gm1)
        rho = rho0 * (a / a0) ** (2 / gm1)
        return u, a, p, rho

    # Time series storage
    ts_t = [0.0]
    ts_pb = [P0]
    ts_po = [P0]
    ts_pm = [P0]
    ts_M  = [0.0]
    ts_mf = [0.0]

    t = 0.0
    step = 0
    out_dt = t_end / 5000  # output every ~1.6s
    next_out = out_dt

    while t < t_end:
        u, a, p, rho = get_state(Jp, Jm)

        # CFL time step
        max_s = np.max(np.abs(u) + a)
        dt = dt_cfl * dx / max_s
        if t + dt > t_end:
            dt = t_end - t

        # Characteristic speeds
        cp = u + a  # C+ speed
        cm = u - a  # C- speed

        # Departure points (semi-Lagrangian)
        xdep_p = x - cp * dt  # where C+ came from
        xdep_m = x - cm * dt  # where C- came from

        # Clamp to domain
        xdep_p = np.clip(xdep_p, 0.0, L)
        xdep_m = np.clip(xdep_m, 0.0, L)

        # Interpolation indices
        idx_p = xdep_p / dx
        idx_m = xdep_m / dx

        j_p = np.clip(np.floor(idx_p).astype(int), 0, Nx - 2)
        j_m = np.clip(np.floor(idx_m).astype(int), 0, Nx - 2)

        f_p = idx_p - j_p
        f_m = idx_m - j_m

        # Linear interpolation for interior
        Jp_new = Jp[j_p] * (1 - f_p) + Jp[j_p + 1] * f_p
        Jm_new = Jm[j_m] * (1 - f_m) + Jm[j_m + 1] * f_m

        # Left BC (x=0): wall, u=0 => J+ = -J-
        Jp_new[0] = -Jm_new[0]

        # Right BC (x=L): outlet
        if choked:
            # Check if sonic condition gives p > P_EXT
            # Sonic: J- = -J+ * (3-gamma)/(gamma+1)
            Jp_out = Jp_new[-1]
            Jm_sonic = -Jp_out * (3.0 - GAMMA) / gp1
            u_s, a_s, p_s, _ = get_state(np.array([Jp_out]), np.array([Jm_sonic]))
            if p_s[0] > P_EXT:
                Jm_new[-1] = Jm_sonic
            else:
                # Pressure too low for sonic exit — use pressure BC
                # Set p = P_EXT at outlet
                a_ext = a0 * (P_EXT / P0) ** (gm1 / (2 * GAMMA))
                # J- from incoming + pressure constraint
                # p = P0*(a/a0)^(2g/(g-1)) => a = a0*(P_EXT/P0)^((g-1)/(2g))
                # u = J+ + J-, a = (g-1)/4*(J+ - J-)
                # a = a_ext => J+ - J- = 4*a_ext/(g-1)
                # J- = J+ - 4*a_ext/(g-1)
                Jm_new[-1] = Jp_out - 4.0 * a_ext / gm1
        else:
            # No choke: determine if flow is supersonic or subsonic
            u_out = 0.5 * (Jp_new[-1] + Jm_new[-2])
            a_out = gm1 / 4.0 * (Jp_new[-1] - Jm_new[-2])
            if u_out > 0 and u_out >= a_out:
                # Supersonic: both characteristics from interior (extrapolate)
                pass  # Jm_new[-1] already set from interior interpolation
            else:
                # Subsonic: J+ from interior, impose p = P_EXT
                a_ext = a0 * (P_EXT / P0) ** (gm1 / (2 * GAMMA))
                Jp_out = Jp_new[-1]
                Jm_new[-1] = Jp_out - 4.0 * a_ext / gm1

        Jp = Jp_new
        Jm = Jm_new
        t += dt
        step += 1

        # Output
        if t >= next_out or t >= t_end:
            u_now, a_now, p_now, rho_now = get_state(Jp, Jm)
            M_out = u_now[-1] / a_now[-1] if a_now[-1] > 0 else 0
            ts_t.append(t)
            ts_pb.append(p_now[0])
            ts_po.append(p_now[-1])
            ts_pm.append(p_now[Nx // 2])
            ts_M.append(M_out)
            ts_mf.append(rho_now[-1] * u_now[-1] * AREA)
            next_out += out_dt

    print(f"  {'choked' if choked else 'unchoked'}: {step} steps, "
          f"{len(ts_t)} output points, final p_base = {ts_pb[-1]/1e5:.2f} bar")
    return (np.array(ts_t), np.array(ts_pb), np.array(ts_po),
            np.array(ts_pm), np.array(ts_M), np.array(ts_mf))


print("\nRunning analytical MOC (Nx=10000)...")
t_moc_c, pb_moc_c, po_moc_c, pm_moc_c, M_moc_c, mf_moc_c = \
    moc_solve(Nx=10000, t_end=8000.0, choked=True)
t_moc_u, pb_moc_u, po_moc_u, pm_moc_u, M_moc_u, mf_moc_u = \
    moc_solve(Nx=10000, t_end=8000.0, choked=False)

# ── Load numerical HLLC results ──
def load_ts(fname):
    data = np.loadtxt(fname, comments='#')
    return {
        't': data[:, 0], 'massflow': data[:, 1], 'u_out': data[:, 2],
        'rho_out': data[:, 3], 'p_out': data[:, 4], 'T_out': data[:, 5],
        'totmass': data[:, 6], 'p_base': data[:, 7], 'p_mid': data[:, 8],
        'a_out': data[:, 9], 'Mach_out': data[:, 10],
    }

print("\nLoading numerical HLLC data...")
d_ch = d_nc = d_rg = None
for name, fname, var in [
    ('ideal_choke',   'timeseries_ideal_choke.dat',   'd_ch'),
    ('ideal_nochoke', 'timeseries_ideal_nochoke.dat',  'd_nc'),
    ('realgas_choke', 'timeseries_realgas_choke.dat',  'd_rg'),
]:
    try:
        d = load_ts(fname)
        print(f"  {name}: {len(d['t'])} pts, t=[{d['t'][0]:.0f}, {d['t'][-1]:.0f}] s")
        if var == 'd_ch': d_ch = d
        elif var == 'd_nc': d_nc = d
        else: d_rg = d
    except Exception as e:
        print(f"  {name}: not available ({e})")

# ── Analysis ──
print("\n" + "="*60)
print("CHOKING ANALYSIS")
print("="*60)

# Simple wave analytics
p_out_sw = P0 * (2 / gp1) ** (2 * GAMMA / gm1)
u_out_sw = 2 * a0 / gp1
rho_out_sw = GAMMA * p_out_sw / u_out_sw**2
mf_sw = rho_out_sw * u_out_sw * AREA

print(f"\nWave arrival at closed end: {t_arrive:.1f} s ({t_arrive/60:.1f} min)")
print(f"Round trip time: {2*t_arrive:.0f} s ({2*t_arrive/60:.1f} min)")
print(f"\nSimple wave (exact, t < {t_arrive:.0f} s):")
print(f"  p_out  = {p_out_sw/1e5:.2f} bar ({100*p_out_sw/P0:.1f}% of p0)")
print(f"  u_out  = {u_out_sw:.1f} m/s = a_out (Mach 1)")
print(f"  rho    = {rho_out_sw:.2f} kg/m3")
print(f"  m_dot  = {mf_sw:.1f} kg/s")

# When does choking end?
print(f"\nChoking cessation (Mach < 1 at outlet):")

# MOC choked
idx = np.where((t_moc_c > t_arrive + 100) & (M_moc_c < 0.995))[0]
if len(idx) > 0:
    t_uc = t_moc_c[idx[0]]
    print(f"  MOC analytical:     t = {t_uc:.0f} s ({t_uc/60:.1f} min), "
          f"{(t_uc - t_arrive)/60:.1f} min after arrival")

# HLLC choked
if d_ch is not None:
    idx = np.where((d_ch['t'] > t_arrive + 100) & (d_ch['Mach_out'] < 0.995))[0]
    if len(idx) > 0:
        t_uc = d_ch['t'][idx[0]]
        print(f"  HLLC ideal choked:  t = {t_uc:.0f} s ({t_uc/60:.1f} min), "
              f"{(t_uc - t_arrive)/60:.1f} min after arrival")
    else:
        print(f"  HLLC ideal choked:  flow remains choked through t={d_ch['t'][-1]:.0f} s")

# HLLC real gas
if d_rg is not None:
    idx = np.where((d_rg['t'] > t_arrive + 100) & (d_rg['Mach_out'] < 0.995))[0]
    if len(idx) > 0:
        t_uc = d_rg['t'][idx[0]]
        print(f"  HLLC real gas:      t = {t_uc:.0f} s ({t_uc/60:.1f} min), "
              f"{(t_uc - t_arrive)/60:.1f} min after arrival")
    else:
        print(f"  HLLC real gas:      still choked at t={d_rg['t'][-1]:.0f} s (run incomplete)")

# ── Plotting ──
fig, axes = plt.subplots(3, 2, figsize=(14, 13))
fig.suptitle('Ideal gas tube blowdown: NS1B (217.7 km), no friction/heat transfer',
             fontsize=13, y=0.995)

colors = {'moc_c': 'k', 'moc_u': '0.5', 'hllc_c': '#1f77b4', 'hllc_u': '#d62728', 'real': '#2ca02c'}

def add_arrival_line(ax):
    ax.axvline(t_arrive / 60, color='gray', ls=':', lw=0.7, alpha=0.6)

# ── Panel 1: p_base full view ──
ax = axes[0, 0]
ax.plot(t_moc_c / 60, pb_moc_c / 1e5, '-', color=colors['moc_c'], lw=1.5,
        label='Analytical (choked)')
ax.plot(t_moc_u / 60, pb_moc_u / 1e5, '--', color=colors['moc_u'], lw=1.5,
        label='Analytical (unchoked)')
if d_ch:
    ax.plot(d_ch['t'] / 60, d_ch['p_base'] / 1e5, '-', color=colors['hllc_c'],
            lw=1, alpha=0.8, label='HLLC (choked)')
if d_nc:
    ax.plot(d_nc['t'] / 60, d_nc['p_base'] / 1e5, '-', color=colors['hllc_u'],
            lw=1, alpha=0.8, label='HLLC (no choke)')
if d_rg:
    ax.plot(d_rg['t'] / 60, d_rg['p_base'] / 1e5, '-', color=colors['real'],
            lw=1, alpha=0.8, label='HLLC real gas')
add_arrival_line(ax)
ax.set_ylabel('Pressure [bar]')
ax.set_title('Pressure at closed end')
ax.legend(fontsize=7, loc='upper right')
ax.grid(True, alpha=0.3)

# ── Panel 2: p_out ──
ax = axes[0, 1]
ax.plot(t_moc_c / 60, po_moc_c / 1e5, '-', color=colors['moc_c'], lw=1.5,
        label='Analytical (choked)')
ax.plot(t_moc_u / 60, po_moc_u / 1e5, '--', color=colors['moc_u'], lw=1.5,
        label='Analytical (unchoked)')
if d_ch:
    ax.plot(d_ch['t'] / 60, d_ch['p_out'] / 1e5, '-', color=colors['hllc_c'],
            lw=1, alpha=0.8, label='HLLC (choked)')
if d_nc:
    ax.plot(d_nc['t'] / 60, d_nc['p_out'] / 1e5, '-', color=colors['hllc_u'],
            lw=1, alpha=0.8, label='HLLC (no choke)')
ax.axhline(p_out_sw / 1e5, color='orange', ls=':', lw=0.8,
           label=f'Simple wave: {p_out_sw/1e5:.1f} bar')
ax.axhline(P_EXT / 1e5, color='brown', ls=':', lw=0.8, label=f'p_ext: {P_EXT/1e5:.1f} bar')
add_arrival_line(ax)
ax.set_ylabel('Pressure [bar]')
ax.set_title('Pressure at outlet')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3)

# ── Panel 3: Mach number ──
ax = axes[1, 0]
ax.plot(t_moc_c / 60, M_moc_c, '-', color=colors['moc_c'], lw=1.5,
        label='Analytical (choked)')
ax.plot(t_moc_u / 60, M_moc_u, '--', color=colors['moc_u'], lw=1.5,
        label='Analytical (unchoked)')
if d_ch:
    ax.plot(d_ch['t'] / 60, d_ch['Mach_out'], '-', color=colors['hllc_c'],
            lw=1, alpha=0.8, label='HLLC (choked)')
if d_nc:
    ax.plot(d_nc['t'] / 60, d_nc['Mach_out'], '-', color=colors['hllc_u'],
            lw=1, alpha=0.8, label='HLLC (no choke)')
if d_rg:
    ax.plot(d_rg['t'] / 60, d_rg['Mach_out'], '-', color=colors['real'],
            lw=1, alpha=0.8, label='HLLC real gas')
ax.axhline(1.0, color='orange', ls=':', lw=0.8, label='M = 1')
add_arrival_line(ax)
ax.set_ylabel('Mach number')
ax.set_title('Outlet Mach number')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3)
ax.set_ylim(-0.05, 1.6)

# ── Panel 4: Mass flow ──
ax = axes[1, 1]
ax.plot(t_moc_c / 60, mf_moc_c, '-', color=colors['moc_c'], lw=1.5,
        label='Analytical (choked)')
ax.plot(t_moc_u / 60, mf_moc_u, '--', color=colors['moc_u'], lw=1.5,
        label='Analytical (unchoked)')
if d_ch:
    ax.plot(d_ch['t'] / 60, d_ch['massflow'], '-', color=colors['hllc_c'],
            lw=1, alpha=0.8, label='HLLC (choked)')
if d_nc:
    ax.plot(d_nc['t'] / 60, d_nc['massflow'], '-', color=colors['hllc_u'],
            lw=1, alpha=0.8, label='HLLC (no choke)')
if d_rg:
    ax.plot(d_rg['t'] / 60, d_rg['massflow'], '-', color=colors['real'],
            lw=1, alpha=0.8, label='HLLC real gas')
ax.axhline(mf_sw, color='orange', ls=':', lw=0.8, label=f'Simple wave: {mf_sw:.0f} kg/s')
add_arrival_line(ax)
ax.set_ylabel('Mass flow [kg/s]')
ax.set_title('Mass flow rate at outlet')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3)

# ── Panel 5: Zoom p_base near arrival ──
ax = axes[2, 0]
ax.plot(t_moc_c / 60, pb_moc_c / 1e5, '-', color=colors['moc_c'], lw=1.5,
        label='Analytical (choked)')
ax.plot(t_moc_u / 60, pb_moc_u / 1e5, '--', color=colors['moc_u'], lw=1.5,
        label='Analytical (unchoked)')
if d_ch:
    ax.plot(d_ch['t'] / 60, d_ch['p_base'] / 1e5, '-', color=colors['hllc_c'],
            lw=1, alpha=0.8, label='HLLC (choked)')
if d_nc:
    ax.plot(d_nc['t'] / 60, d_nc['p_base'] / 1e5, '-', color=colors['hllc_u'],
            lw=1, alpha=0.8, label='HLLC (no choke)')
add_arrival_line(ax)
ax.set_xlim((t_arrive - 30) / 60, (t_arrive + 600) / 60)
ax.set_xlabel('Time [min]')
ax.set_ylabel('Pressure [bar]')
ax.set_title('Zoom: p_base near wave arrival')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3)

# ── Panel 6: Choked vs unchoked difference ──
ax = axes[2, 1]
if d_ch and d_nc:
    t_c = d_ch['t']
    pb_nc_i = np.interp(t_c, d_nc['t'], d_nc['p_base'])
    diff = (d_ch['p_base'] - pb_nc_i)
    ax.plot(t_c / 60, diff / 1e5, '-', color=colors['hllc_c'], lw=1,
            label='p_base(choked) - p_base(unchoked)')
    # Also compare mass flow
    ax2 = ax.twinx()
    mf_nc_i = np.interp(t_c, d_nc['t'], d_nc['massflow'])
    ax2.plot(t_c / 60, d_ch['massflow'] - mf_nc_i, '--', color=colors['hllc_u'],
             lw=1, alpha=0.7, label='massflow diff')
    ax2.set_ylabel('Mass flow difference [kg/s]', color=colors['hllc_u'], fontsize=9)
    ax2.tick_params(axis='y', labelcolor=colors['hllc_u'])
    add_arrival_line(ax)
    ax.set_xlabel('Time [min]')
    ax.set_ylabel('Pressure difference [bar]')
    ax.set_title('Effect of choking BC (HLLC ideal gas)')
    ax.legend(fontsize=7, loc='upper left')
    ax2.legend(fontsize=7, loc='upper right')
    ax.grid(True, alpha=0.3)
else:
    ax.text(0.5, 0.5, 'Data not available', transform=ax.transAxes, ha='center')
    ax.set_xlabel('Time [min]')

for ax in axes[0, :]:
    ax.set_xlabel('')
for ax in axes[1, :]:
    ax.set_xlabel('')
axes[2, 0].set_xlabel('Time [min]')
axes[2, 1].set_xlabel('Time [min]')

plt.tight_layout()
plt.savefig('choking_analysis.png', dpi=150, bbox_inches='tight')
print("\nSaved: choking_analysis.png")

# ── Quantitative comparison table ──
print("\n" + "="*70)
print(f"{'Quantity':<35} {'Analytical':>12} {'HLLC choke':>12} {'HLLC nochk':>12}")
print("-"*70)

# p_base at final time
print(f"{'p_base at t=8000 [bar]':<35}", end="")
print(f" {pb_moc_c[-1]/1e5:>12.2f}", end="")
if d_ch: print(f" {d_ch['p_base'][-1]/1e5:>12.2f}", end="")
if d_nc: print(f" {d_nc['p_base'][-1]/1e5:>12.2f}", end="")
print()

# p_out during simple wave
print(f"{'p_out in simple wave [bar]':<35}", end="")
print(f" {p_out_sw/1e5:>12.2f}", end="")
if d_ch:
    mask = (d_ch['t'] > 50) & (d_ch['t'] < t_arrive - 50)
    if mask.any():
        print(f" {np.mean(d_ch['p_out'][mask])/1e5:>12.2f}", end="")
if d_nc:
    mask = (d_nc['t'] > 50) & (d_nc['t'] < t_arrive - 50)
    if mask.any():
        print(f" {np.mean(d_nc['p_out'][mask])/1e5:>12.2f}", end="")
print()

# Mass flow during simple wave
print(f"{'m_dot simple wave [kg/s]':<35}", end="")
print(f" {mf_sw:>12.1f}", end="")
if d_ch:
    mask = (d_ch['t'] > 50) & (d_ch['t'] < t_arrive - 50)
    if mask.any():
        print(f" {np.mean(d_ch['massflow'][mask]):>12.1f}", end="")
if d_nc:
    mask = (d_nc['t'] > 50) & (d_nc['t'] < t_arrive - 50)
    if mask.any():
        print(f" {np.mean(d_nc['massflow'][mask]):>12.1f}", end="")
print()
