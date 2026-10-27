#!/bin/sh
# Build the independent Kv-based verifications (J_0, J_1, J_3).
#   Kv 0.4.62 (header-only, http://verifiedby.me/kv/) and Boost uBLAS headers.
#   KV must point to the directory that contains kv/ (e.g. KV=$HOME/kv).
# -DKV_FASTROUND: Kv's SSE2 rounding-mode switch (x86-64). -frounding-math and -ffp-contract=off keep
# the compiler from constant-folding across rounding-mode changes or fusing multiply-adds. No -ffast-math.
set -e
: "${KV:?set KV=/path/to/kv}"
CXX=${CXX:-g++}
cd "$(dirname "$0")"
for p in tangency_kv tangency_kv_n6 tangency_kv_n16a; do
    $CXX -std=c++17 -O3 -DNDEBUG -DKV_FASTROUND -frounding-math -ffp-contract=off -I"$KV" $p.cpp -o $p
    echo "built ./$p"
done
# Run: ./tangency_kv_n16a            (optional args: Delta  taylor_order  c_shift)
