# Independent verification with the Kv library

> **Note (release v1.2).** Each program now also checks the non-resonance condition on the multipliers of the
> saddle used for the linearisation hypothesis of Mora--Viana (Remark 3.2 of the paper): it encloses
> rho = log(1/|lambda_s|)/log|lambda_u| for all c in J and certifies that no resonance
> lambda_s = lambda_s^i lambda_u^j or lambda_u = lambda_s^i lambda_u^j of order i+j <= 10^3 occurs (it also prints the
> first order at which a resonance is not excluded: 1200 for J_1, 1229030 for J_3, 45389 for J_0). Each program now performs **38 checks**
> (32 for the hypotheses, 6 for the statement), all true; all other enclosures are unchanged bit for bit
> (logs in `logs/`). The new check is part of the statement block, which is skipped in the negative tests, so the
> negative-test logs of v1.1 remain valid.

> **Note (release v1.1).** In the paper *Existence of SRB measures for some Rössler systems* (T. Caby), the
> band-merging tangency verified by `tangency_kv.cpp` is **Proposition 6.1**; the tangencies near c = 5.7
> (**Proposition 5.1**, on which the Main Theorem rests) are re-verified by the generated variants
> `tangency_kv_n6.cpp` (J_1) and `tangency_kv_n16a.cpp` (J_3). Since v1.1 every program also checks the numerical *statement* of the theorem as
> printed in the paper (exact endpoints of J, boxes Q_K and Q, range of λ_u, bound on |det|): each program now
> performs **37 checks** (32 for the hypotheses, 5 for the statement), all true. The detailed description below
> was written for the previous version (35 checks for J_0, where the band-merging result was called "Theorem A");
> the hypothesis checks and all enclosures are unchanged.


This directory contains a second, independent computer-assisted verification of the hypotheses
of Theorem A of `rossler_tangency.tex` (Section 4, Lemmas 4.1–4.5): a quadratic, generically
unfolding homoclinic tangency of the Poincaré map P_c of the Rössler system

    x' = -(y+z),  y' = x + a y,  z' = b + z (x - c),   a = b = 0.2,

on Σ = {x = 0, x' > 0} (coordinates (y,z)), for c ∈ J = [c_L, c_R].

It uses **Kv 0.4.62** (M. Kashiwagi, header-only, git `c7f8f23`, Boost 1.83 uBLAS) instead of
CAPD. The formula layer is written from the paper. No CAPD code is used, and nothing from the
CAPD program is read at run time.

**Verdict: all 35 checks print `true`** (`proof_output.txt`). Runtime: 416 s wall-clock on one
core of an Intel Xeon at 2.1 GHz, with a niced background build sharing
the 2 cores. A rebuild after comment-only edits reproduced every enclosure bit for bit.

Breakdown of the runtime:
- (a): about 50 s;
- the W^u and W^s boxes: about 20 s each (one revolution, 15-dimensional jet);
- (C): 36 s and 13 s;
- (B): 2 × 33 s;
- (A): 131 s;
- (D): 121 s.

## Files

| File | Content |
|---|---|
| `tangency_kv.cpp` | Main program: all checks of Lemmas 4.1–4.5 (formula layer) |
| `rossler_jet.hpp` | Rössler field in polar coordinates, reparametrised by the polar angle, with forward-mode first/second-order jets (variational equations) |
| `mv_driver.hpp` | Long integration of affine sets: the mean-value algorithm of `kv::ode_maffine`, with a different step-size strategy |
| `poincare_jet.hpp` | N-fold Poincaré map one revolution at a time; exact power-of-two rescaling of the jet components |
| `build.sh` | Build script (flags explained inside) |
| `proof_output.txt` | Full output of the proof run (default parameters) |
| `sanity/` | Non-rigorous finite-difference check (`fd_check.py`), the previous agent's θ' exploration (`explore_theta.py`), and the logs of the negative and robustness runs |
| `dev/` | Exploratory experiments and timing logs. **Not part of the proof**, not maintained (some sources no longer compile against the final headers). `dev/superseded/` holds the previous agent's `rossler_polar.hpp` and `test_timing.cpp`, replaced by `rossler_jet.hpp` |

Build and run:

    KV=/path/to/kv ./build.sh  # KV = directory containing the kv/ header folder
    ./tangency_kv              # optional args: Delta  TaylorOrder  c-shift (the shift is for negative tests only)

Compile flags: `-std=c++17 -O3 -DNDEBUG -DKV_FASTROUND -frounding-math -ffp-contract=off`.
- `-O3 -DNDEBUG` follow the Kv/uBLAS usage.
- `KV_FASTROUND` is Kv's documented SSE2 rounding-mode switch. On x86-64 all double arithmetic is SSE2.
- `-frounding-math` and `-ffp-contract=off` stop constant folding across rounding-mode changes and FMA contraction.
- No `-ffast-math`, no `-march=native`.

## Method

1. **Poincaré map as a fixed-time map.**
   - Write x = r cos θ, y = r sin θ and use t = θ + π/2 as the independent variable:
     - dr/dt = r(a r S² − C z)/D
     - dz/dt = r(b + z(rC − c))/D
     - where C = sin t, S = −cos t and D = r(1 + aCS) + Sz = r θ'.
   - At t ≡ 0 (mod 2π) we have x = 0, y = −r, and x' = rθ'. The other zeros of x (t ≡ π) have x' = −rθ' < 0.
   - So as long as θ' > 0, P_c is exactly the time-2π map in the coordinates (y,z) = (−r, z), and P_c^N is the time-2πN map.
     There is no section-crossing computation, and no return-time correction in the derivatives.
   - θ' > 0 is enforced rigorously:
     - Kv's power-series division evaluates the range of D over every Taylor step, and throws if the range contains 0. A throw makes the step, and then the proof, fail.
     - D = r − z > 0 at t = 0 is checked for every initial set.
     - Hence D > 0 along every enclosed orbit.
   - Numerically θ' ≥ 0.88 on all relevant orbits.

2. **Jets.**
   - The state carries derivatives of (r,z) with respect to up to 3 parameters p_k of the initial condition: (y0,z0) = X + Σ p_k d_k, and c = c0 + p_CSLOT.
   - It optionally carries second derivatives (p0,p0), (p0,p1), (p1,p1).
   - The right-hand side is the field evaluated in forward-mode jet arithmetic (exact chain rule), i.e. exactly the first/second variational equations.
   - Seeds are chosen per run:
     - unit vectors (Krawczyk, dp/dc);
     - e_u, e_s, and p̂' with dc/dp = 1 (local maps);
     - d1 = e_u + w' e_s with **w' an interval** (tangency).
   - So D_N d1 and H_N[d1,d1] are computed directly as ∂_sΓ and the curve part of ∂²_sΓ.
   - The c-slot gives D_N p̂' + ∂_c P^N as one directional derivative, so these two terms are never enclosed separately and cannot cancel.

3. **Sets and long integration.**
   - All sets are Kv affine forms. The noise symbols are s, c, the w-range, the remainder box of p_c, etc.
   - The parametrisation p_c = p̂ + p̂'(c − c_mid) + R keeps the correlation between the base point and the parameter.
   - Each Taylor step uses the mean-value form of `kv::ode_maffine`: x ← φ(mid) + Dφ(hull)·(x − mid) in affine arithmetic.
     - The verified one-step enclosures (`kv::ode` for the centre, and `kv::ode` on `kv::autodif` for the Jacobian over the hull) are Kv's.
     - Correctness of the combination is the mean value theorem, row by row; see the comment in `mv_driver.hpp`.
   - Two engineering changes were needed. Neither affects rigour:
     - **(i) Step size.** Kv's own driver lets the hull integration choose the step, and its heuristic makes tiny steps once the set is stretched. Here the centre point proposes the step, and it is shortened adaptively until the hull integration verifies.
     - **(ii) Rescaling.** Kv's step-size control makes tiny steps when some component is large (step ∝ 1/size, tested). Jet components grow like 1.92^k, and second derivatives reach about 2·10⁷. The stored jet components are therefore rescaled by exact powers of two before each step. The field gets the matching factors: dc/dp → 2^{−E}, and λ = 2^{E_i+E_j−F} in the second-order product rule. The stored jet then solves exactly the variational equations of the rescaled parameters. Components are only scaled down.
   - With these changes, 15 returns of a 9-dimensional jet over S×J take about 2 min.

4. **Formula layer** (written from Section 4 of the paper):
   - Krawczyk (Lemma 4.1);
   - p_c' = −(DP−I)⁻¹∂_cP;
   - χ(±1);
   - the graph-transform hypotheses (U1)–(U6) and (S1)–(S7). The bounds L, K, M, L_s, K_s, M_s are found by iteration with safety factor 1.05 and then verified;
   - Lemma 4.2(c);
   - the chain rule of Lemma 4.4;
   - conditions (A)–(D) of Lemma 4.5, with Remark 4.6 for (C).

## Parameters

Same as the CAPD proof:
- c_mid = 4.50698934631197, ε_c = 10⁻⁸. c_L and c_R are the same doubles as in the paper; this is checked.
- N = 15, s_m = 5.3682155586·10⁻⁴, Δ = 10⁻⁹.
- δ = 6·10⁻⁴, ρ = 10⁻⁵, r_u = r_s = 10⁻⁶.
- E, p̂, p̂' as given.
- Krawczyk radii: (3·10⁻⁸, 10⁻⁹) for J, and 10⁻¹¹ at c_L, c_mid, c_R.
- Taylor order 20.

Δ = 10⁻⁷ was also run as a robustness test; see below.

## What was verified (35 checks, all true)

Every number below is a rigorous enclosure from `proof_output.txt`.

**(a) Fixed point, saddle (Lemma 4.1):**
- Krawczyk K ⊂ int Q holds for all c ∈ J, and for thin c = c_L, c_mid, c_R. K_L, K_mid, K_R ⊂ Q_J, so the four are on the same branch.
- p_c ∈ [−6.9386987435099954, −6.9386987188730753] × [0.036063488559252142, 0.036063488765782005] for all c ∈ J.
- Theorem A(1): the paper's box Q lies inside our Krawczyk box (uniqueness) and contains K_J (existence).
- dp/dc ∈ [−1.2318423479, −1.2318423290] × [−0.00664781948, −0.00664781908].
- Df_c(0): tr ∈ [−1.92147092397, −1.92147089409], det ∈ [3.43502288·10⁻¹¹, 3.43502983·10⁻¹¹].
  - χ(1) ∈ [2.92147, 2.92148] > 0 and χ(−1) ∈ [−0.92148, −0.92147] < 0.
  - |det| < 1, and also |det| < 4.6·10⁻⁸ (Theorem A(2)).
  - The det enclosure is strictly positive, so λ_s = det/λ_u ∈ (−1, 0). This fixes the sign of λ_s; the CAPD enclosure of det contained 0.

**(b) Lemma 4.2 on B_u = [−6·10⁻⁴, 6·10⁻⁴] × [−10⁻⁵, 10⁻⁵], c ∈ J:**
- (U1) mig Φ₁ ≥ 1.92098.
- (U2), both parts: Lδ ≤ ρ, and mag a₂₁ δ + mag a₂₂ ρ ≤ ρ.
- (U3)–(U5) with L = 3.07·10⁻⁷, K = 1.26·10⁻⁴, M = 1.07·10⁻⁶.
- (U6) θ = 1.25·10⁻⁸.
- Part (c): t₀ = 1.24·10⁻¹¹ on J, and 6.67·10⁻¹² at c_L and at c_R. Consistency |t| ≤ L holds.

**(c) Lemma 4.3 on B_s = [−10⁻⁶, 10⁻⁶]², c ∈ J:**
- (S1)–(S7) with L_s = 9.99·10⁻⁷, K_s = 0.969, M_s = 1.43·10⁻⁶.
- mig D ≥ 1.92147.

**(d) Tangency, N = 15, S = s_m ± 10⁻⁹ (exact endpoints):**

| | quantity | enclosure | CAPD enclosure |
|---|---|---|---|
| (C) | g(s_m, c_L) | [−1.72065041·10⁻⁸, −1.72056735·10⁻⁸] | [−1.7531, −1.6882]·10⁻⁸ |
| (C) | ∂_s g(s_m, c_L) | [2.48·10⁻⁵, 3.15·10⁻⁵] | [−5.9·10⁻⁵, 1.2·10⁻⁴] |
| (C) | sup_S g(·,c_L) | ≤ −1.72056735·10⁻⁸ | ≤ −1.688·10⁻⁸ |
| (C) | g(s_m, c_R) | [1.72042500·10⁻⁸, 1.72050059·10⁻⁸] | [1.6881, 1.7529]·10⁻⁸ |
| (B) | ∂_s g(s_m−Δ, J) | [0.021701, 0.022356] | [0.01908, 0.02497] |
| (B) | ∂_s g(s_m+Δ, J) | [−0.022947, −0.022292] | [−0.02557, −0.01967] |
| (A) | ∂²_s g on S×J | [−22323945.2, −22323828.6] | [−2.23388, −2.23090]·10⁷ |
| (D) | ∂_c g on S×J | [1.7205026, 1.7205712] | [1.4164, 2.0247] |
| chart | V on S×J | [6.7·10⁻¹², 2.6·10⁻¹¹] (\|V\| ≤ r_s) | \|V\| ≤ 9.92·10⁻⁹ |

Further details:
- |V| ≤ r_s was also checked at the S-endpoints × J and at (s_m, c_L) and (s_m, c_R).
- S ⊂ (0, δ] was checked.
- The chain rule of Lemma 4.4 is evaluated as follows:
  - ∂_sΓ = D_N d1 and ∂²_sΓ = H_N[d1,d1] + w'' D_N e_s, with d1 = e_u + w' e_s, |w'| ≤ t₀ + K|s| and |w''| ≤ K;
  - (U_c,V_c) = E⁻¹(G₂ − p̂') + (E⁻¹D_N e_u) δp₁ + (E⁻¹D_N e_s)(δp₂ + ∂_c w) − δp, where G₂ = D_N p̂' + ∂_cP^N and δp = E⁻¹(p_c' − p̂'). This is algebraically identical to the formula of Lemma 4.4;
  - g = U − σ(V) is bounded by U + [−L_s|V|, L_s|V|] (σ(0) = 0);
  - ∂_s g = U_s − σ'V_s;
  - ∂²_s g = U_ss − σ''V_s² − σ'V_ss;
  - ∂_c g = U_c − ∂_cσ − σ'V_c.

The enclosures are about 10× (B) to 10⁴× (D) narrower than CAPD's. They agree with the double-precision values in PROOF.md: dG/dc ≈ 1.7205, g'' ≈ −2.2324·10⁷, λ_u ≈ −1.92147, λ_s ≈ −1.79·10⁻¹¹ = det/λ_u.

## Sanity checks (non-rigorous) and negative tests

- `sanity/fd_check.py` computes Df and D²f of the local map at the fixed point by finite differences. It uses SciPy event location in **Cartesian** coordinates, so it is independent of the polar formulation. It gives:
  - a₁₁ = −1.9214709;
  - h¹ = (0.788627, 0.048455, −1.772731);
  - h²₁₁ = −4.397·10⁻⁴.
  All agree with the W^s-box enclosures to finite-difference accuracy, which checks the second-order jet including mixed terms and the y = −r sign conventions.
- Shift J by +3·10⁻⁸ (`sanity/neg_shift_plus3e-8.txt`): only (C) "g(s,c_L) < 0" fails, with g(s_m,c_L) ≈ +3.44·10⁻⁸. This is as predicted by ∂_c g ≈ 1.72.
- Shift J by −3·10⁻⁸ (`sanity/neg_shift_minus3e-8.txt`), Δ = 10⁻¹² (`sanity/neg_delta_1e-12.txt`), and the robustness run with Δ = 10⁻⁷ (`sanity/robust_delta_1e-7.txt`): see the section "Additional runs" at the end of this file.

## Differences from the CAPD proof

**Library and integrator:**
- Kv: interval Taylor method with verified remainder, plus affine arithmetic and the mean-value form.
- CAPD: C²-Lohner Poincaré maps of the 4-D extended field.

**Poincaré map:**
- Here: the polar-angle reparametrisation, so P_c is a fixed-time map with no crossing detection.
- CAPD: section crossing with return-time corrections.

**Derivatives:**
- Here: directional jets seeded along e_u, e_s and p̂'. The tangency run (A) needs only the single directional second derivative H_N[d1,d1], with w' entering as an interval in the seed.
- CAPD: full Hessians in (y,z), contracted with E afterwards.

**c-derivative:**
- Here: one directional derivative along (p̂', 1), which avoids the cancellation between D_N p_c' and ∂_cP^N.
- The paper (remark "robustness and sanity checks", Section 5) reports that with Δ = 10⁻⁷ the CAPD enclosure of ∂_c g becomes too wide to decide (D), because of the loss of correlation between D_N and ∂_cP^N across S. Here Δ = 10⁻⁷ is a robustness test; see below.

**Fixed-point branch check:**
- Here: K_{L,mid,R} ⊂ Q_J, which is what uniqueness in Q_J needs.
- CAPD checked ⊂ K_J. Here K_J is so tight that K_R touches its boundary to within about 10⁻¹⁵, so ⊂ K_J would be a fragile check.

**Bounds obtained:** they are much smaller (table above). The conclusions are the same.

## What is NOT verified / trust base

**Trusted:**
- Kv 0.4.62:
  - interval arithmetic and the rounding control (`KV_FASTROUND`);
  - power-series arithmetic, including its domain mode;
  - `kv::ode`, the Taylor method with Picard verification;
  - the `kv::autodif` variational ODE (default `ODE_AUTODIF_NEW = 1` path);
  - affine arithmetic and `epsilon_reduce2`.
- The compiler (g++ 13.3) and the FPU.
- My own non-Kv code. It is short and its arguments are given in the comments:
  - the mean-value driver;
  - the power-of-two rescaling;
  - the jet arithmetic;
  - the formula layer.
  It has not been independently audited.

**Not re-proved:**
- Lemmas 4.1–4.5 themselves: graph-transform arguments, the Arzelà–Ascoli limit, analyticity, and the criterion. Only their hypotheses are checked.
- The deduction of Theorem B (Mora–Viana, Benedicks–Viana, suspension).
- The "location of the tangency point" remark in the proof of Theorem A.
- Uniqueness of c*, which holds only for the tracked branch, as in the paper.

**Non-rigorous inputs** (they only choose where to look): c_mid, s_m, E, p̂, p̂', and the approximate inverse in Krawczyk.

**Not checked:** anything about c = 5.7.

## Additional runs

All logs are in `sanity/`. Same program (sources identical to the final ones up to comments) and same parameters unless stated.

| run | arguments | result |
|---|---|---|
| J shifted by +3·10⁻⁸ (the approximate fixed point is shifted with it) | `1e-9 20 3e-8` | only (C) `g(s,c_L) < 0` fails: g(s_m,c_L) ∈ [3.44096, 3.44104]·10⁻⁸ > 0 (32 of 33 checks true; the "paper doubles" and "Theorem A(1)" checks are skipped for shifted J) |
| J shifted by −3·10⁻⁸ | `1e-9 20 -3e-8` | only (C) `g(s_m,c_R) > 0` fails (32 of 33 true) |
| Δ = 10⁻¹² | `1e-12` | both (B) checks fail: g_s(s_m−Δ, J) ∈ [−6.0·10⁻⁴, 5.4·10⁻⁵] contains 0, because the true endpoint slope ≈ 2.2·10⁻⁵ is below the variation of g_s over J (33 of 35 true) |
| Δ = 10⁻⁷ (robustness) | `1e-7` | **all 35 checks true** (361 s): g_s(s_m∓Δ, J) ≈ ±2.232, ∂²_s g ∈ [−22329850, −22317924], ∂_c g ∈ [1.71728, 1.72380], \|V\| ≤ 3.8·10⁻⁸ |

The negative tests fail exactly where the paper predicts. The Δ = 10⁻⁷ run gives the same
conclusions as the main run with 100 times more margin in (B). The paper reports that this Δ is
undecidable for (D) with the CAPD set-up.
