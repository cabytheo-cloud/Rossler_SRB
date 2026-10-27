#!/bin/sh
# Rebuild and rerun the computer-assisted proofs of "Existence of SRB measures for some Rossler systems".
#   CAPD programs and negative tests: always (needs capd-config on PATH).
#   Kv programs:  only if KV points to the directory containing the kv/ header folder.
#   Lean:         only if LEAN_MATHLIB points to a Mathlib checkout at the commit in lean/README.md.
#   FD check:     only if RUN_FD=1 (non-rigorous, a few minutes).
set -e
cd "$(dirname "$0")"
ROOT=$(pwd)
OUT="$ROOT/rerun"
mkdir -p "$OUT"
ok=0; bad=0

run() {   # run <label> <logfile> <cmd...>
    label=$1; log=$2; shift 2
    "$@" > "$log" 2>&1 || true
    if grep -q "ALL CHECKS PASSED: true\|failed: 0" "$log"; then
        echo "  PASS  $label  ($(grep -c '^CHECK' "$log") checks)"; ok=$((ok+1))
    else
        echo "  FAIL  $label  (see $log)"; bad=$((bad+1))
    fi
}

echo "== CAPD (expected: 40 checks each, all true)"
command -v capd-config >/dev/null || { echo "capd-config not found on PATH"; exit 1; }
cd "$ROOT/capd"
for p in tangency_cap tangency_57_n6 tangency_57_n16b tangency_57_n16a; do
    g++ -O2 $p.cpp $(capd-config --cflags --libs) -o $p
done
run "J_0 (Prop. 6.1)" "$OUT/capd_J0.txt" ./tangency_cap
run "J_1 (Prop. 5.1)" "$OUT/capd_J1.txt" ./tangency_57_n6
run "J_2 (Prop. 5.1)" "$OUT/capd_J2.txt" ./tangency_57_n16b
run "J_3 (Prop. 5.1)" "$OUT/capd_J3.txt" ./tangency_57_n16a

echo "== CAPD negative tests (each must fail exactly at (T3) = check (C))"
for spec in "tangency_cap 1e-8 1e-9 3e-8" "tangency_57_n6 1e-8 1e-9 3e-8" "tangency_57_n16b 1e-10 1e-10 3e-10" "tangency_57_n16a 1e-11 1e-11 3e-11"; do
    set -- $spec
    for sg in "" "-"; do
        log="$OUT/neg_$1_${sg}shift.txt"
        ./$1 $2 $3 ${sg}$4 > "$log" 2>&1 || true
        nfail=$(grep '^CHECK' "$log" | grep -vc ': true' || true)
        which=$(grep '^CHECK' "$log" | grep -v ': true' | sed 's/  */ /g' | cut -c1-50)
        echo "  $1 shift ${sg}$4: $nfail failing check(s): $which"
    done
done

g++ -O2 "$ROOT/trap/trap2.cpp" $(capd-config --cflags --libs) -o "$ROOT/trap/trap2"
"$ROOT/trap/trap2" > "$OUT/trap2.txt" 2>&1 && echo "  trapping region: see $OUT/trap2.txt"

if [ "$RUN_FD" = "1" ]; then
    echo "== finite-difference comparison (non-rigorous)"
    g++ -O2 "$ROOT/capd/tests/test_fd.cpp" $(capd-config --cflags --libs) -o "$ROOT/capd/tests/test_fd"
    "$ROOT/capd/tests/test_fd" > "$OUT/fd_J0.txt" 2>&1 && tail -1 "$OUT/fd_J0.txt"
fi

if [ -n "$KV" ]; then
    echo "== Kv (expected: 38 checks each, all true; 7-13 minutes each)"
    cd "$ROOT/kv"
    for p in tangency_kv tangency_kv_n6 tangency_kv_n16a; do
        g++ -std=c++17 -O3 -DNDEBUG -DKV_FASTROUND -frounding-math -ffp-contract=off -I"$KV" $p.cpp -o $p
    done
    run "Kv J_0" "$OUT/kv_J0.txt" ./tangency_kv
    run "Kv J_1" "$OUT/kv_J1.txt" ./tangency_kv_n6
    run "Kv J_3" "$OUT/kv_J3.txt" ./tangency_kv_n16a
else
    echo "== Kv skipped (set KV=/path/to/kv to run it)"
fi

if [ -n "$LEAN_MATHLIB" ]; then
    echo "== Lean"
    mkdir -p "$LEAN_MATHLIB/RosslerTangency"
    cp "$ROOT/lean/Criterion.lean" "$LEAN_MATHLIB/RosslerTangency/"
    (cd "$LEAN_MATHLIB" && lake env lean RosslerTangency/Criterion.lean) > "$OUT/lean.txt" 2>&1 \
        && echo "  Lean file checks; axioms:" && grep -A3 "depends on axioms" "$OUT/lean.txt"
else
    echo "== Lean skipped (set LEAN_MATHLIB=/path/to/mathlib4 to run it)"
fi

echo "== $ok passed, $bad failed"
[ $bad -eq 0 ]
