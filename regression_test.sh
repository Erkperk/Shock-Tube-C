#!/bin/bash
# Regression tests for shock_tube solver.
# Tests three key physical quantities against reference values:
#   1. Peak mass flow at choked outlet (HLLC, dx=1m, 2s)
#   2. Rarefaction wave arrival at closed end (LF, dx=1m, 700s)
#   3. Time for initial 1.7 bar drop after arrival
#
# Usage: ./regression_test.sh [--quick]
#   --quick: skip the long wave-arrival test (test 1 only)

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

# --- Test 1: Peak mass flow (HLLC, dx=1m, 2s) ---
echo "Test 1: Peak mass flow (HLLC, dx=1m, 2s)..."
./shock_tube --ns1ag --dx-init 1 --endtime 2 --nprofiles 2 --prefix _regtest_mf > /dev/null 2>&1

MF_PEAK=$(python3 -c "
data = []
for l in open('timeseries__regtest_mf.dat'):
    if l.startswith('#'): continue
    p = l.split()
    if len(p) >= 13:
        try: data.append(float(p[1])/1000)
        except: pass
print(f'{max(data):.3f}')
")
check "Peak mass flow [t/s]" "$MF_PEAK" 21.0 0.5

if [ "$1" = "--quick" ]; then
    echo ""
    echo "=== Quick mode: skipping wave arrival test ==="
else
    # --- Test 2 & 3: Wave arrival and 1.7 bar drop (LF, dx=1m, 700s) ---
    echo ""
    echo "Test 2-3: Wave arrival and pressure drop (LF, dx=1m, 700s)..."
    ./shock_tube --ns1ag --lf --dx-init 1 --early-end 1 --early-dt 0.1 --endtime 700 --prefix _regtest_wave > /dev/null 2>&1

    RESULTS=$(python3 -c "
import numpy as np
data = []
for l in open('timeseries__regtest_wave_lf.dat'):
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

    check "Wave arrival time [s]" "$T_ARR" 461.0 10.0
    check "1.7 bar drop time [s]" "$DT17" 20.0 5.0
fi

# --- Summary ---
echo ""
echo "=== Results: $PASS passed, $FAIL failed ==="
rm -f timeseries__regtest_*.dat profiles__regtest_*.dat early__regtest_*.dat

if [ $FAIL -gt 0 ]; then
    exit 1
fi
