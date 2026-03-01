#!/usr/bin/env python3
"""Sensitivity analysis: water temperature and back-pressure effects on NS1A Russian blowdown."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def load_ts(fname):
    data = []
    for l in open(fname):
        if l.startswith('#'): continue
        p = l.split()
        if len(p) >= 13:
            try: data.append([float(x) for x in p[:13]])
            except: pass
    return np.array(data)

runs = [
    ('timeseries_ns1ar_cal.dat',     'Baseline (282 K, 8.61 bar)', 'C0', '-'),
    ('timeseries_ns1ar_T279.dat',    'Cold water (279 K, 8.61 bar)', 'C3', '-'),
    ('timeseries_ns1ar_vacuum.dat',  'Vacuum outlet (282 K, 0 bar)', 'C1', '--'),
    ('timeseries_ns1ar_pR3.dat',     'Low back-pressure (282 K, 3 bar)', 'C2', '--'),
]

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('NS1A Russian end sensitivity: water temperature & back-pressure', fontsize=14)

loaded = []
for fname, lab, col, ls in runs:
    ds = load_ts(fname)
    loaded.append((fname, lab, col, ls, ds))

# Panel 1: Full blowdown pressure
ax = axes[0, 0]
for fname, lab, col, ls, ds in loaded:
    ax.plot(ds[:, 0] / 86400, ds[:, 7] / 1e5, color=col, linestyle=ls, linewidth=1.2, label=lab)
ax.axhline(y=8.0, color='k', linestyle=':', alpha=0.3)
ax.set_xlabel('Time [days]')
ax.set_ylabel('Closed-end pressure [bar]')
ax.set_title('Full blowdown')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Panel 2: Zoom on late blowdown (last 2 days, p < 30 bar)
ax = axes[0, 1]
for fname, lab, col, ls, ds in loaded:
    mask = ds[:, 7] / 1e5 < 30
    if mask.any():
        ax.plot(ds[mask, 0] / 86400, ds[mask, 7] / 1e5, color=col, linestyle=ls, linewidth=1.2, label=lab)
ax.axhline(y=8.0, color='k', linestyle=':', alpha=0.3, label='8 bar (ambient)')
ax.axhline(y=3.0, color='gray', linestyle=':', alpha=0.3, label='3 bar (venturi)')
ax.set_xlabel('Time [days]')
ax.set_ylabel('Pressure [bar]')
ax.set_title('Late blowdown (p < 30 bar)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Panel 3: Mass flow
ax = axes[1, 0]
for fname, lab, col, ls, ds in loaded:
    mf = ds[:, 1] / 1000  # t/s
    mf_plot = mf.copy()
    mf_plot[ds[:, 0] < 5] = np.nan
    ax.plot(ds[:, 0] / 86400, mf_plot, color=col, linestyle=ls, linewidth=1, label=lab)
ax.set_xlabel('Time [days]')
ax.set_ylabel('Mass flow [t/s]')
ax.set_title('Outlet mass flow')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Panel 4: Summary table
ax = axes[1, 1]
ax.axis('off')

# Explosion time offsets
t_sep28 = 111370
t_sep29 = 197770
t_sep30 = 284170
t_oct1  = 370570

rows = []
for fname, lab, col, ls, ds in loaded:
    t = ds[:, 0]
    pb = ds[:, 7] / 1e5

    def p_at(t_target):
        idx = np.argmin(np.abs(t - t_target))
        return f'{pb[idx]:.1f}'

    idx8 = np.where(pb < 8.0)[0]
    t8 = f'{t[idx8[0]]/3600:.0f}' if len(idx8) > 0 else '>125'

    # Cumulative mass released
    dt = np.diff(t)
    mf = ds[:, 1]
    cum = np.sum(0.5 * (mf[:-1] + mf[1:]) * dt) / 1e6  # kt

    short_lab = lab.split('(')[0].strip()
    rows.append([short_lab, p_at(t_sep28), p_at(t_sep29), p_at(t_sep30), t8, f'{cum:.0f}'])

table = ax.table(cellText=rows,
                 colLabels=['Case', 'Sep 28\n[bar]', 'Sep 29\n[bar]', 'Sep 30\n[bar]',
                            'To 8 bar\n[hr]', 'Mass\n[kt]'],
                 loc='center', cellLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1.2, 1.8)

plt.tight_layout()
fig.savefig('sensitivity_russian.png', dpi=150)
print('Saved: sensitivity_russian.png')
