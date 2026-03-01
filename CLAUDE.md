# Shock Tube Project Instructions

## Peak Mass Flow Recording
- `--tube-start` uses the ghost cell outlet BC (choked/subsonic NRBC) from t=0 — NO vacuum extension
- Flow immediately establishes choked conditions (M=1) at the outlet from the first timestep
- Peak mass flow ~20.3 t/s (165 bar, 1.153m dia, λ=0.007) settles by t≈15ms after a brief startup transient
- The peak is independent of pipe length — German and Russian ends give identical peaks at same initial conditions
- When computing peak from timeseries, skip the first ~10ms startup transient (use t > 0.01s)
- Tube-start auto-enables early recording (tube_start+3s at 0.01s dt) to capture the peak
- Timeseries sampling at dt=endtime/5000 (e.g. 10s for 50000s runs) MISSES the peak
- Any reported peak below ~17 t/s from coarse timeseries is a sampling artifact

## Wave Arrival Recording
- Tube-start auto-enables early recording for only tube_start+3s (e.g. 5s at 0.01s dt)
- This leaves a ~450s GAP with no data between t=5s and the first coarse sample
- For German presets, the wave arrives at the closed end at t~450-475s — right in this gap
- To capture wave arrival at high resolution, manually specify `--early-end` beyond the expected arrival time
- Example: `--early-end 1900 --early-dt 1` gives 1s resolution for the first ~30 minutes
- This overrides the auto early recording, so the 0.01s peak resolution is lost — but the peak is still captured at the first timestep (t~0.001s)

## Sound Speed Verification
- ALWAYS verify sound speed with CoolProp when computing expected wave arrival times
- CoolProp sound speed for methane at 165 bar, 282 K: **480.8 m/s** (NOT 346 m/s)
- Sound speed varies with pressure: 481 (165 bar), 437 (120), 421 (80), 424 (40), 434 (10 bar)
- Use `/opt/homebrew/bin/python3` with CoolProp: `CP.PropsSI('A', 'T', T, 'P', P, 'Methane')`
- Expected wave arrival: L / a(initial). NS1A German (224 km): 466s, NS1B German (217.7 km): 454s

## Cell Merging
- Cell merging must not produce visible artifacts in timeseries plots (mass flow, pressure)
- If a sudden jump/drop appears in mass flow coinciding with a merge event, this is a bug that must be investigated and fixed
- Merge conservation should be at machine precision (~1e-15 relative mass error)

## Simulation Runs
- Full blowdown runs: use `--tube-start 2` (enables early recording + delayed coarsening, ghost cell BC from t=0)
- German pipeline presets: `--ns1ag`, `--ns1bg`, `--ns2ag` (50000s endtime)
- Russian pipeline presets: `--ns1ar`, `--ns1br`, `--ns2ar` (450000s endtime)
- Use `--pstop 8e5` to stop when closed-end pressure reaches 8 bar
- Wall model options: default (transient 1D), `--simple-wall` (steady-state), `--noheat` (adiabatic)
