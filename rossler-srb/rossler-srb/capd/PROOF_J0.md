# A generic homoclinic tangency in the Rössler family (computer-assisted proof)

*This write-up concerns the band-merging tangency near c ≈ 4.507, which is **Proposition 6.1** of the paper "Existence of SRB measures for some Rössler systems" (T. Caby). The main result of the paper (the Main Theorem, via Proposition 5.1: three tangencies near c = 5.7) is proved with the variants `tangency_57_*.cpp` of the same program; see the paper and `README.md`.*

## Statement

Rössler system: x' = -(y+z), y' = x + a y, z' = b + z(x - c), with **a = b = 0.2**.
P_c is the Poincaré map on Σ = {x = 0, x' > 0}, written in the coordinates (y, z).

**Proposition 6.1 of the paper (computer-assisted).** Let J = [c_L, c_R], where c_L and c_R are the binary doubles
c_L = 4.50698933631197018456759906257502734661102294921875 and
c_R = 4.50698935631197006301817964413203299045562744140625 (that is, 4.50698934631197 ± 10⁻⁸ after rounding).
There is a c* ∈ (c_L, c_R) with the following properties.
1. For every c ∈ J, P_c has a unique fixed point p_c in the box Q_K = [-6.9386987611, -6.9386987012] × [0.0360634877, 0.0360634896], and p_c lies in Q_0 = [-6.9386987439, -6.9386987185] × [0.0360634882, 0.0360634892].
2. p_c is a dissipative flip saddle: λ_u ∈ [-1.92148, -1.92146], -1 < λ_s < 0, |λ_u λ_s| < 5·10⁻⁸.
3. W^u(p_{c*}) and W^s(p_{c*}) have a **quadratic tangency that unfolds generically** in c, in the sense of Palis–Takens / Mora–Viana (normal form with α ≠ 0, v ≠ 0); the tangency point lies within 2·10⁻⁸ of p_{c*}.

The tangency is between the branch Γ(s,c) = P_c^15(W^u_loc) for s ∈ S and W^s_loc(p_c), where S and the charts are as below. c* is unique for this branch. It is not claimed to be the only homoclinic tangency in J: secondary tangencies from longer iterates generically accumulate at c*.

**Proposition 3.5 of the paper (for this tangency).** Combine Proposition 6.1 with
- Mora–Viana (Acta Math. 1993): Hénon-like strange attractors for a positive-measure set of parameters;
- Benedicks–Young (1993) for Hénon maps, and for the Hénon-like families arising from Mora–Viana: Benedicks–Viana (2001), Wang–Young (2001). These attractors carry a unique SRB measure, and almost every point of the topological basin is generic for it (no holes);
- the suspension argument: the Poincaré-map attractor lifts to a flow attractor, and its SRB measure lifts to one for the flow.

Then (relying on Romero 1995 / Gonchenko–Shilnikov–Turaev 2008 for the linearisation hypothesis of Mora–Viana §2) there is a set E ⊂ ℝ of positive Lebesgue measure accumulating at c* ≈ 4.507 such that, for every c ∈ E, the Rössler flow with a = b = 0.2 has a (small, Hénon-like) strange attractor supporting an SRB measure whose basin has positive volume.

*Remark.* The saddle has negative multipliers. Apply Mora–Viana to P_c², which has the same saddle, the same manifolds and the same tangency. If Λ is an attractor of P_c² with SRB measure ν, then Λ ∪ P_c(Λ) is an attractor of P_c with SRB measure (ν + P_c*ν)/2.

## What is NOT proved
- Nothing about c = 5.7. Numerically (not proved), c* ≈ 4.507 is the band-merging (2→1 band) crisis of the attractor.
- The set E is not explicit. Mora–Viana gives positive measure, not a specific parameter.
- The attractors of Proposition 3.5 are the small Hénon-like attractors near the tangency orbit, not necessarily the "big" numerically observed attractor.

## Proof of Proposition 6.1: structure

**Notation.** E = [e_u e_s] is a fixed approximate eigenbasis at c_mid. The local map is
  f_c(u,v) = E⁻¹(P_c(p_c + E(u,v)) - p_c),  so f_c(0,0) = 0.
All derivative enclosures are taken over boxes × J and computed with CAPD's C²-Lohner Poincaré map of the 4D extended field (c' = 0).

**1. Fixed point.** A Krawczyk test on Q × J (and on thin boxes at c_L, c_mid, c_R) gives a unique fixed point p_c. The derivative is enclosed via p_c' = -(DP - I)⁻¹ ∂_c P.

**2. Saddle.** A0 = E⁻¹ DP(Q,J) E. Let χ(λ) = λ² - tr λ + det. Then χ(1) > 0 and χ(-1) < 0, so the multipliers satisfy λ_u < -1 < λ_s < 1. Also |det| < 1.

**3. W^u_loc as a graph** v = w_c(u), |u| ≤ δ_u = 6·10⁻⁴.
- Class 𝒲: w(0) = 0, |w| ≤ ρ, |w'| ≤ L_u, |w''| ≤ K_u, |∂_c w| ≤ M_u.
- Graph transform: ŵ(φ(u)) = ψ(u), where φ = f1(u, w(u)) and ψ = f2(u, w(u)). Its derivatives are
  - ŵ' = ψ'/φ'
  - ŵ'' = (ψ''φ' - ψ'φ'')/φ'³
  - ŵ_c = ψ_c - ŵ'φ_c
- Verified:
  - min|φ'| > 1;
  - Γ(𝒲) ⊂ 𝒲;
  - C⁰-contraction θ = |a21|κ + |a22|(L_u κ + 1) ≈ 4·10⁻⁴, where κ = (|a12|/m)/(1 - |a12| L_u / m) and m = min|a11|.
- Limit argument:
  - Iterate the graph transform from w ≡ 0.
  - 𝒲 is closed in C⁰: Lipschitz bounds on w' and on c ↦ w survive uniform limits (Arzelà–Ascoli).
  - So the unique fixed point w* lies in 𝒲.
  - The true local manifold is analytic in (u,c), so these Lipschitz bounds become pointwise bounds on |w''| and |∂_c w|.
- Backward orbits along the graph contract to p_c, so graph ⊂ W^u(p_c).
- The tangent at p_c is the unstable eigenvector (1, t), with t = a21/(a11 + a12 t - a22). This gives |w'(0)| ≤ 1.1·10⁻⁸ on J and ≤ 6.8·10⁻¹² at c_L, c_R. Hence |w(s)| ≤ |w'(0)||s| + K_u s²/2.
- Bounds obtained: L_u = 2.0·10⁻⁴, K_u = 9.4·10⁻⁴, M_u = 1.1·10⁻³.

**4. W^s_loc as a graph** u = σ_c(v), |v| ≤ r_s = 10⁻⁶ (with |u| ≤ r_u = 10⁻⁶).
- Graph transform: σ̂(v) is the unique u with f1(u,v) = σ(f2(u,v)). Its derivatives are
  - σ̂' = (σ' a22 - a12)/(a11 - σ' a21)
  - σ̂'' = (σ'' F2'² + σ' Q2 - Q1)/(a11 - σ' a21)
  - σ̂_c = (σ' f2c + σ_c - f1c)/(a11 - σ' a21)
- Verified:
  - invariance of the class;
  - contraction factor 1/min|a11 - σ' a21| < 1;
  - v-image stays in [-r_s, r_s];
  - root stays in [-r_u, r_u];
  - the restricted map on the graph contracts, so graph ⊂ W^s(p_c).
- Limit argument: as in §3, iterating from σ ≡ 0.
- Bounds obtained: L_s = 4.3·10⁻⁵, K_s = 0.97, M_s = 2.2·10⁻⁴.

**5. Tangency function.** Let N = 15 and X0(s,c) = p_c + E(s, w_c(s)) ∈ W^u_loc. Define
  - Γ(s,c) = P_c^N(X0(s,c)) ∈ W^u(p_c);
  - (U, V) = E⁻¹(Γ - p_c);
  - g(s,c) = U - σ_c(V).

{g = 0} is exactly "Γ ∈ W^s_loc(p_c)". The derivatives use the C² jet of P^N on an affine set (Lohner doubleton) that contains X0(S×J):
- g_s = U_s - σ'V_s
- g_ss = U_ss - σ''V_s² - σ'V_ss
- g_c = U_c - σ_c - σ'V_c, with ∂_c Γ - p_c' = (DP^N - I)p_c' + DP^N e_s ∂_c w + ∂_c P^N.

Verified on S = s_m + [-Δ, Δ], with s_m = 5.3682155586·10⁻⁴ and Δ = 10⁻⁹ (as doubles). The endpoints of S are evaluated exactly, as s_m plus a thin offset ±Δ:
- **(A)** g_ss ∈ [-2.2339·10⁷, -2.2309·10⁷] < 0 on S × J.
- **(B)** g_s(s_m - Δ, c) ∈ [0.0190, 0.0250] > 0 and g_s(s_m + Δ, c) ∈ [-0.0256, -0.0196] < 0 for all c ∈ J. So there is a unique critical point s*(c) ∈ S.
- **(C)** sup_S g(·, c_L) ≤ g(s_m) + g_s²/(2 min|g_ss|) ≤ -1.688·10⁻⁸ < 0 and g(s_m, c_R) ≥ 1.688·10⁻⁸ > 0.
- **(D)** g_c ∈ [1.416, 2.025] > 0 on S × J.

The invariant manifolds depend analytically on c, so g is C² jointly in (s,c). By the IFT, s*(c) is C¹. Let G(c) = g(s*(c), c). Then G is C¹ with G' = g_c(s*(c), c) > 0, and G(c_L) < 0 < G(c_R). So there is a unique c* with G(c*) = 0. At (s*(c*), c*) we have g = g_s = 0, g_ss ≠ 0 and g_c ≠ 0. This is a quadratic tangency, unfolding generically. ∎

## Consistency with independent double-precision numerics (Python/SciPy)
| Quantity | Python estimate | Enclosure in the proof |
|---|---|---|
| c₀ | 4.50698934631 | contains c* |
| dG/dc | 1.7205 | [1.416, 2.025] |
| g'' | -2.2324·10⁷ | [-2.2339·10⁷, -2.2309·10⁷] |
| λ_u | -1.92147 | [-1.9214717, -1.9214702] |
| λ_s | -1.79·10⁻¹¹ | inside the enclosure |

Negative tests behave as expected:
- shifting J by ±3·10⁻⁸ makes (C) fail, with the sign predicted by dG/dc;
- shrinking S to ±10⁻¹² makes (B) fail.

## Reproduce
- CAPD::DynSys from github.com/CAPDGroup/CAPD (built here from commit 2f06098 with g++ 13.3).
- Build: `g++ -O2 tangency_cap.cpp $(capd-config --cflags --libs) -o tangency_cap`
- Run: `./tangency_cap` (defaults ε_c = 10⁻⁸, Δ = 10⁻⁹; optional arguments ε_c, Δ, and a c-shift for negative tests). It takes about 4 s and prints 40 CHECK lines, all true: 35 for the hypotheses, and 5 that check the numerical statement of Proposition 6.1 as printed in the paper (box Q_K inside the Krawczyk box, K_J inside Q, the range of λ_u, |det| < 5·10⁻⁸, and the distance of the tangency point to p).

## Independent audit
A separate referee pass re-derived every graph-transform and chain-rule formula. It also checked CAPD's Hessian convention (Taylor coefficients) and D(P^n) against finite differences. It found no fatal issue. Its fixes are included:
- exact S endpoints;
- J taken as the exact double interval;
- KL, KM, KR ⊂ KJ (same fixed-point branch);
- positivity of two denominators;
- the uniqueness claim restricted to the tracked branch;
- the limit argument spelled out.

Remark on λ_s: the det enclosure contains 0, so the CAPD computation does not fix the sign of λ_s. Liouville (orientation preservation of P) gives det DP > 0, so λ_s < 0; this is used to apply Mora–Viana to P². The independent Kv enclosure of det is strictly positive.

## Trust base
- Correctness of CAPD's interval arithmetic and C² Poincaré-map enclosures.
- The compiler and FPU rounding mode.
- The cited theorems (Mora–Viana, Benedicks–Young, Benedicks–Viana).

The non-rigorous inputs (c_mid, s_m, the eigenbasis E, approximate fixed points) only choose where to look. They do not enter the logic.
