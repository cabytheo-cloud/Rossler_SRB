# Existence of SRB measures for some Rössler systems: code and logs

Code and output archive (release v1.2) for the paper

> T. Caby, *Existence of SRB measures for some Rössler systems* (2026).

The paper (PDF and LaTeX source) is in `paper/`. Repository: https://github.com/cabytheo-cloud/Rossler_SRB

## What is verified

System: x' = −(y+z), y' = x + a y, z' = b + z(x − c) with a = b = 0.2; P_c is the Poincaré map on {x = 0, x' > 0}.

The paper's Main Theorem: the set of c within 6·10⁻⁵ of 5.7 for which the flow has a Hénon-like strange attractor with an SRB measure has positive Lebesgue measure. Its computer-assisted part is Proposition 5.1. The programs verify, by interval arithmetic, the finitely many inequalities of Section 4 of the paper, which imply that P_c has a quadratic homoclinic tangency of its period-one saddle that unfolds generically for some c in each interval below. They also check the numerical **statement** of the propositions as printed in the paper (boxes Q_K and Q, range of λ_u, bound on |λ_u λ_s|, distance of the tangency point; for Kv also the exact endpoints of J and the absence of resonances λ_s^ℓ λ_u^m = 1 of order at most 10³, which is the non-resonance hypothesis of Mora–Viana used in Remark 3.2 of the paper).

| Interval | Paper | CAPD program (checks) | Kv program (checks) |
|---|---|---|---|
| J_1 ≈ 5.69941999 ± 1e-8 (n = 6, N = 15) | Prop. 5.1 | `capd/tangency_57_n6.cpp` (40) | `kv/tangency_kv_n6.cpp` (38) |
| J_2 ≈ 5.6995796079 ± 1e-10 (n = 16, N = 26) | Prop. 5.1 | `capd/tangency_57_n16b.cpp` (40) | — |
| J_3 ≈ 5.70005464311 ± 1e-11 (n = 16, N = 25; closest to 5.7; used for the Main Theorem) | Prop. 5.1 | `capd/tangency_57_n16a.cpp` (40) | `kv/tangency_kv_n16a.cpp` (38) |
| J_0 ≈ 4.50698934631 ± 1e-8 (band merging, N = 15) | Prop. 6.1 | `capd/tangency_cap.cpp` (40) | `kv/tangency_kv.cpp` (38) |

All checks print `true` (logs in `capd/logs/` and `kv/logs/`).

**Negative tests:** each J shifted by ±3·ε_c, so that it no longer contains the tangency (the statement checks are then skipped). In every run exactly one check fails, (T3), at c_L for the positive shift and at c_R for the negative one (`capd/neg/`, `kv/logs/neg_*`).

In the program output, the tangency conditions (T1)–(T4) of the paper are labelled (A), (B), (C), (D): (A) ∂²g/∂s² < 0, (B) sign change of ∂g/∂s, (C) g(·, c_L) < 0 < g(s_m, c_R), (D) ∂g/∂c > 0.

## Layout

```
capd/        CAPD verification programs (main proof)
  tangency_cap.cpp                 template; parameters and statement of Proposition 6.1 (J_0)
  make_variant.py, params_*.json   generate tangency_57_*.cpp (substitutes ONLY numerical parameters and the statement)
  tangency_57_{n6,n16b,n16a}.cpp   programs for J_1, J_2, J_3 (Proposition 5.1)
  logs/, neg/                      full outputs and negative tests (J_0 ... J_3)
  tests/test_fd.cpp, tests/fd_*.txt  non-rigorous finite-difference comparison of the Poincare-map derivatives
  tests/test_hess.cpp, test_sem.cpp  small checks of CAPD conventions
  PROOF_J0.md                      detailed write-up of the J_0 run
kv/          independent re-implementation with the Kv library (fixed-angle Poincare map, no return-time correction)
  tangency_kv.cpp + *.hpp          program for J_0
  make_kv_variant.py               generates tangency_kv_{n6,n16a}.cpp (uses ../capd/params_*.json)
  logs/, sanity/                   outputs, negative tests, finite-difference sanity checks
trap/        trapping region P(W) ⊂ int W at c = 5.7 (Remark 8.1)
lean/        Lean 4 formalisation of the tangency criterion (Lemma 4.5, Remark 4.6)
numerics/    non-rigorous exploration and figure scripts (NOT part of the proof)
paper/       the paper
reproduce.sh rebuild and rerun everything that is installed
```

`python3 capd/make_variant.py tangency_57_n16a capd/params_n16a.json` (and likewise for the other variants and for `kv/make_kv_variant.py`) regenerates the shipped programs byte for byte.

## Requirements

- **CAPD::DynSys** 6.1.0 (github.com/CAPDGroup/CAPD; logs produced with commit `2f06098`), `g++` 13.3; `capd-config` on `PATH`.
- **Kv** 0.4.62 (http://verifiedby.me/kv/, header-only) and Boost uBLAS headers; set `KV` to the directory that contains the `kv/` header folder.
- **Lean** v4.35.0, Mathlib commit `aa35a481e5f0cdb3b04a899ead15d6ca7305a624` (see `lean/README.md`).
- Logs were produced on x86-64 Linux.
- CAPD (GNU GPL v3) and Kv (MIT) are third-party libraries and are not included here.

## Reproduce

```sh
./reproduce.sh     # CAPD programs and negative tests (seconds each); Kv too if KV is set (7-13 min each)
```

or by hand:

```sh
cd capd
g++ -O2 tangency_57_n16a.cpp $(capd-config --cflags --libs) -o tangency_57_n16a
./tangency_57_n16a                       # optional args: eps_c  Delta  c_shift
./tangency_57_n16a 1e-11 1e-11 3e-11     # negative test: only (T3) must fail
```

## Audit notes (Section 7.3 of the paper)

- **Derivatives of the Poincaré map, not of the flow.** CAPD's `computeDP` includes the return-time corrections at first and second order; c is a state variable of the extended field, so the same applies to ∂_c; P^N is computed up to the N-th crossing. `capd/tests/test_fd.cpp` compares these derivatives with central finite differences of P^N (non-rigorous). For J_0 (N = 15), all entries larger than 0.1 agree to within 9·10⁻⁴ relative (most to 10⁻⁵–10⁻⁸). For J_3 (N = 25), at the point used in the proof, the first derivatives along e_u and c and all second derivatives except those in the direction (e_s, e_s) agree to within 6·10⁻³; at a point off the fold the second derivatives vary too fast for double-precision finite differences. The independent Kv implementation, which needs no return-time correction, is the main cross-check.
- **E⁻¹** is enclosed by the interval adjugate formula; the programs abort if the interval determinant contains 0.
- **Decisions.** Every `CHECK` is decided by an interval comparison on endpoints; decimal constants of the statement are converted to intervals that enclose them exactly.
- **Enclosure boxes.** Derivative enclosures are taken on the hull of 𝒦_J + E·B, where 𝒦_J is the Krawczyk enclosure of the fixed point for all c ∈ J — not around an approximate fixed point.


## Provenance

The code was written with AI assistance (Claude Opus 5.5 extra, Anthropic). The Kv implementation is independent of the CAPD implementation in library, method and code.
