#!/bin/bash
# Regression tests for shock_tube solver (HLLC only).
# Tests key physical quantities against reference values:
#   1. Peak mass flow at choked outlet (HLLC, dx=1m, 2s)
#   2. Outlet Mach number = 1.0 when choked
#   3. Rarefaction wave arrival at closed end (HLLC, dx=1m, 700s)
#   4. Time for initial 1.7 bar drop after arrival
#
# Usage: ./regression_test.sh [--quick]
#   --quick: skip the long wave-arrival test (tests 1-2 only)

set -e
PASS=0
FAIL=0

check() {
    local name="$1" val="$2" ref="$3" tol="$4"
    local diff=$(python3 -c "print(abs($val - $ref))")
    local ok=$(python3 -c "print('PASS' if abs($val - $ref) <= $tol else 'FAIL')")
    if [ "$ok" = "PASS" ]; then
        printf "  %-30s %8.2f  (ref=%6.1f ±%.1f)  PASS\n" "$name" "$val" "$ref" "$tol"
        PASS=$((PASS + 1))
    else
        printf "  %-30s %8.2f  (ref=%6.1f ±%.1f)  FAIL (diff=%.2f)\n" "$name" "$val" "$ref" "$tol" "$diff"
        FAIL=$((FAIL + 1))
    fi
}

echo "=== Shock Tube Regression Tests ==="
echo ""

# --- Test 1-2: Peak mass flow and outlet Mach (HLLC, dx=1m, 2s) ---
echo "Test 1-2: Peak mass flow and outlet Mach (HLLC, dx=1m, 2s)..."
./shock_tube --ns1ag --dx-init 1 --endtime 2 --early-end 1 --early-dt 0.01 --prefix _regtest_mf > /dev/null 2>&1

RESULTS=$(python3 -c "
data = []
for l in open('timeseries__regtest_mf.dat'):
    if l.startswith('#'): continue
    p = l.split()
    if len(p) >= 13:
        try: data.append([float(p[0]), float(p[1])/1000, float(p[10]), int(float(p[12]))])
        except: pass
# Skip first 0.01s (ghost cell artifact before flow establishes)
settled = [(t,mf,M,ch) for t,mf,M,ch in data if t > 0.01]
peak_mf = max(settled, key=lambda x: x[1])
# Average Mach during choked steady state (t=0.05 to 0.5s)
mach_vals = [M for t,mf,M,ch in settled if 0.05 < t < 0.5 and ch == 1]
avg_mach = sum(mach_vals)/len(mach_vals) if mach_vals else 0
print(f'{peak_mf[1]:.3f} {avg_mach:.4f}')
")

MF_PEAK=$(echo $RESULTS | cut -d' ' -f1)
MACH_OUT=$(echo $RESULTS | cut -d' ' -f2)

check "Peak mass flow [t/s]" "$MF_PEAK" 20.4 0.5
check "Outlet Mach (choked)" "$MACH_OUT" 1.0 0.01

# --- Test 5: Outlet stability (no NaN, no oscillation in first 1s) ---
echo ""
echo "Test 5: Outlet stability (HLLC, dx=1m, 1s)..."
./shock_tube --ns1ag --dx-init 1 --endtime 1 --no-merge --early-end 0 --early-dt 0.001 --prefix _regtest_stab > /dev/null 2>&1

STAB_RESULTS=$(python3 -c "
import sys
data = []
for l in open('timeseries__regtest_stab.dat'):
    if l.startswith('#'): continue
    p = l.split()
    if len(p) >= 13:
        try: data.append([float(x) for x in p[:13]])
        except: pass

# Check 1: No NaN in mass flow or total mass
nan_count = sum(1 for r in data if r[1] != r[1] or r[6] != r[6])

# Check 2: No wild oscillations in mass flow after t=0.02s
#   The mass flow should vary by less than 2000 kg/s between consecutive steps
settled = [(r[0], r[1]) for r in data if r[0] > 0.02]
max_jump = 0
for i in range(1, len(settled)):
    jump = abs(settled[i][1] - settled[i-1][1])
    if jump > max_jump:
        max_jump = jump

# Check 3: Sound speed at outlet stays reasonable (< 500 m/s)
max_a = max(r[9] for r in data if r[0] > 0.01)

print(f'{nan_count} {max_jump:.1f} {max_a:.1f}')
")

NAN_COUNT=$(echo $STAB_RESULTS | cut -d' ' -f1)
MAX_JUMP=$(echo $STAB_RESULTS | cut -d' ' -f2)
MAX_A=$(echo $STAB_RESULTS | cut -d' ' -f3)

check "NaN count (outlet)" "$NAN_COUNT" 0 0
check "Max mass flow jump [kg/s]" "$MAX_JUMP" 200.0 1800.0
check "Max outlet sound speed [m/s]" "$MAX_A" 300.0 100.0

if [ "$1" = "--quick" ]; then
    echo ""
    echo "=== Quick mode: skipping wave arrival test ==="
else
    # --- Test 3 & 4: Wave arrival and 1.7 bar drop (HLLC, dx=1m, 700s) ---
    echo ""
    echo "Test 3-4: Wave arrival and pressure drop (HLLC, dx=1m, 700s)..."
    ./shock_tube --ns1ag --dx-init 1 --early-end 1 --early-dt 0.1 --endtime 700 --prefix _regtest_wave > /dev/null 2>&1

    RESULTS=$(python3 -c "
import numpy as np
data = []
for l in open('timeseries__regtest_wave.dat'):
    if l.startswith('#'): continue
    p = l.split()
    if len(p) >= 13:
        try: data.append([float(x) for x in p[:13]])
        except: pass
d = np.array(data)
t, pb = d[:,0], d[:,7]/1e5
drop = pb[0] - pb
idx01 = np.where(drop > 0.1)[0]
t_arr = t[idx01[0]] if len(idx01) > 0 else -1
idx17 = np.where(drop > 1.7)[0]
dt17 = t[idx17[0]] - t_arr if len(idx17) > 0 and t_arr > 0 else -1
print(f'{t_arr:.1f} {dt17:.1f}')
")

    T_ARR=$(echo $RESULTS | cut -d' ' -f1)
    DT17=$(echo $RESULTS | cut -d' ' -f2)

    check "Wave arrival time [s]" "$T_ARR" 466.0 10.0
    check "1.7 bar drop time [s]" "$DT17" 65.0 20.0
fi

# --- Summary ---
echo ""
echo "=== Results: $PASS passed, $FAIL failed ==="
rm -f timeseries__regtest_*.dat profiles__regtest_*.dat early__regtest_*.dat early_timeseries__regtest_*.dat

if [ $FAIL -gt 0 ]; then
    exit 1
fi
