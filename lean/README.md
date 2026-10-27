# Lean 4 formalisation of the tangency criterion (Lemma 4.5) and Taylor bound (Remark 4.6)

Source: `paper/rossler_tangency.tex`, Section 4.6 ("The tangency criterion"),
`\label{lem:criterion}` and `\label{rem:taylor}`.

File: `mathlib4/RosslerTangency/Criterion.lean` (346 lines, namespace `RosslerTangency`).
It contains no `sorry`, `admit`, `native_decide` or `axiom` declarations, and it compiles with no
errors or warnings.

## Toolchain

* Lean `v4.35.0` (`leanprover/lean4:v4.35.0`, commit `eb6cafb7290e8bd5f2937522d8d804f58f02bafc`), Lake 5.0.0.
* Mathlib `master` at commit `aa35a481e5f0cdb3b04a899ead15d6ca7305a624` (2026-10-15).
* Imports: `Mathlib.Analysis.Calculus.Deriv.MeanValue`, `Mathlib.Topology.Order.Compact`,
  `Mathlib.Topology.Order.IntermediateValue`, `Mathlib.Topology.Order.ProjIcc`.

## How to check

```sh
export PATH=$HOME/.elan/bin:$PATH
cd mathlib4
lake build Mathlib.Analysis.Calculus.Deriv.MeanValue Mathlib.Topology.Order.Compact \
           Mathlib.Topology.Order.IntermediateValue Mathlib.Topology.Order.ProjIcc
lake env lean RosslerTangency/Criterion.lean
grep -nE "sorry|admit|native_decide|axiom " RosslerTangency/Criterion.lean   # prints nothing
```

The last lines of the file print the axioms. Expected output:

```
'RosslerTangency.tangency_criterion' depends on axioms: [propext, Classical.choice, Quot.sound]
'RosslerTangency.tangency_criterion'' depends on axioms: [propext, Classical.choice, Quot.sound]
'RosslerTangency.taylor_bound' depends on axioms: [propext, Classical.choice, Quot.sound]
'RosslerTangency.taylor_quadratic_bound' depends on axioms: [propext, Classical.choice, Quot.sound]
```

Build notes: the Mathlib `.olean` cache (`lake exe cache get`) is blocked by the network policy,
so the required part of Mathlib (about 2100 modules) was compiled from source with `nice`. That
took roughly 2–2.5 hours of wall-clock time on this machine, including one restart. The
build-target list contained `Mathlib.Topology.Algebra.Order.Compact`, which no longer exists at
this commit; that one target failed, and everything else built. The file is not registered as a `lean_lib` in Mathlib's lakefile, so it
is checked with `lake env lean` rather than `lake build`.

## Statements

Notation: `sL, sR` = `s₋, s₊`; `g s c`, `gs s c`, `gss s c`, `gc s c` = `g(s,c)`, `∂ₛg`, `∂ₛ²g`,
`∂_c g`; `maxFun g sL sR c = sSup ((fun s => g s c) '' Icc sL sR)` is the paper's `G(c)`.

### Lemma 4.5 (local version, all hypotheses only on the rectangle)

```lean
theorem tangency_criterion {g gs gss gc : ℝ → ℝ → ℝ} {sL sR cL cR : ℝ}
    (hg : ∀ s ∈ Icc sL sR, ∀ c ∈ Icc cL cR, HasDerivAt (fun s => g s c) (gs s c) s)
    (hgs : ∀ s ∈ Icc sL sR, ∀ c ∈ Icc cL cR, HasDerivAt (fun s => gs s c) (gss s c) s)
    (hgc : ∀ s ∈ Icc sL sR, ∀ c ∈ Icc cL cR, HasDerivAt (fun c => g s c) (gc s c) c)
    (hcont : ContinuousOn (fun p : ℝ × ℝ => g p.1 p.2) (Icc sL sR ×ˢ Icc cL cR))
    (hs : sL < sR) (hc : cL < cR)
    (hA : ∀ s ∈ Icc sL sR, ∀ c ∈ Icc cL cR, gss s c < 0)
    (hB : ∀ c ∈ Icc cL cR, 0 < gs sL c ∧ gs sR c < 0)
    (hC : (∀ s ∈ Icc sL sR, g s cL < 0) ∧ ∃ s₀ ∈ Icc sL sR, 0 < g s₀ cR)
    (hD : ∀ s ∈ Icc sL sR, ∀ c ∈ Icc cL cR, 0 < gc s c) :
    ∃ cstar ∈ Ioo cL cR, ∃ sstar ∈ Ioo sL sR,
      g sstar cstar = 0 ∧ gs sstar cstar = 0 ∧ gss sstar cstar < 0 ∧ 0 < gc sstar cstar ∧
      (∀ c ∈ Icc cL cR, ∃! s, s ∈ Icc sL sR ∧ gs s c = 0) ∧
      maxFun g sL sR cstar = 0 ∧ (∀ c ∈ Icc cL cR, maxFun g sL sR c = 0 → c = cstar) ∧
      (∀ c ∈ Icc cL cR, ∀ s ∈ Icc sL sR, g s c = 0 → gs s c = 0 → c = cstar ∧ s = sstar)
```

`tangency_criterion'` is the same statement with the derivatives assumed at every point
(`∀ s c, HasDerivAt …`) and `Continuous fun p : ℝ × ℝ => g p.1 p.2`. It follows directly from
`tangency_criterion`.

### Remark 4.6

```lean
theorem taylor_bound {h h' h'' : ℝ → ℝ} {sL sR smid γ : ℝ}
    (hh : ∀ s ∈ Icc sL sR, HasDerivAt h (h' s) s)
    (hh' : ∀ s ∈ Icc sL sR, HasDerivAt h' (h'' s) s)
    (hγ : 0 < γ) (hbound : ∀ s ∈ Icc sL sR, h'' s ≤ -γ) (hsm : smid ∈ Icc sL sR) :
    ∀ s ∈ Icc sL sR, h s ≤ h smid + h' smid ^ 2 / (2 * γ)
```

The file also proves the intermediate quadratic bound `taylor_quadratic_bound`:
`h s ≤ h smid + h' smid * (s - smid) - γ / 2 * (s - smid) ^ 2`. It does not need `γ > 0`.

## Proof outline and differences from the paper

* The proof does not use the implicit function theorem. For each `c`, (A) and the mean value
  theorem make `gs (·, c)` strictly decreasing. With (B) and the intermediate value theorem,
  `gs (·, c)` therefore has a unique zero `m(c) ∈ (sL, sR)`, and `m(c)` maximises `g(·, c)`
  (`crit_point`).
* `G = maxFun g sL sR` is continuous on `[cL, cR]` (`continuousOn_maxFun`). This uses
  `IsCompact.continuous_sSup` after extending `g` from the rectangle to `ℝ × ℝ` with `projIcc`.
* `G` is strictly increasing on `[cL, cR]` (`maxFun_strictMonoOn`). For `c₁ < c₂`:
  `G c₁ = g(m(c₁), c₁) < g(m(c₁), c₂) ≤ G c₂`, using (D) and the mean value theorem in `c`.
  This replaces the paper's step "`G` is `C¹` with `G' > 0`".
* (C) gives `G cL < 0 < G cR`. The intermediate value theorem then gives `c⋆ ∈ (cL, cR)` with
  `G c⋆ = 0`, and strict monotonicity makes `c⋆` unique.
* The paper also claims that `s_*` and `G` are `C¹`. That regularity is not formalised, because
  the conclusions above do not need it. Every other conclusion of Lemma 4.5 is proved, plus a
  stronger uniqueness statement: `(s⋆, c⋆)` is the only point of the rectangle where
  `g = ∂ₛg = 0`.
* The formal hypotheses are weaker than the paper's `C²` assumption. The formal theorem needs
  only pointwise derivatives `∂ₛg`, `∂ₛ²g`, `∂_c g` on the rectangle and joint continuity of `g`
  there. It does not need continuity of the derivatives.
