import Mathlib.Analysis.Calculus.Deriv.MeanValue
import Mathlib.Topology.Order.Compact
import Mathlib.Topology.Order.IntermediateValue
import Mathlib.Topology.Order.ProjIcc

/-!
# The tangency criterion

Formalisation of the purely analytic part of Section 4.6 of `rossler_tangency.tex`:

* `RosslerTangency.tangency_criterion` (Lemma 4.5, "tangency criterion"): under the sign
  conditions (A)–(D) on a rectangle `[sL, sR] × [cL, cR]` there is a parameter
  `c⋆ ∈ (cL, cR)` and a point `s⋆ ∈ (sL, sR)` with
  `g = ∂ₛg = 0`, `∂ₛ²g < 0`, `∂_c g > 0` at `(s⋆, c⋆)`; for every `c` the critical point of
  `g(·, c)` in `[sL, sR]` is unique; and `c⋆` is the unique `c ∈ [cL, cR]` with
  `max_{s ∈ [sL,sR]} g(s,c) = 0`.
* `RosslerTangency.taylor_bound` (Remark 4.6): if `h'' ≤ -γ < 0` on `[sL, sR]` and
  `smid ∈ [sL, sR]`, then `h s ≤ h smid + (h' smid)^2 / (2γ)` on `[sL, sR]`.

Notation: `sL, sR` stand for `s₋, s₊` of the paper; `g s c`, `gs s c`, `gss s c`, `gc s c` stand for
`g(s,c)`, `∂ₛg(s,c)`, `∂ₛ²g(s,c)`, `∂_c g(s,c)`.

The proof avoids the implicit function theorem: the maximum function
`G c = sSup (g · c '' [sL, sR])` is continuous (compactness) and strictly increasing
(by (D) and the mean value theorem), so the intermediate value theorem gives `G c⋆ = 0`;
by (A) and (B), `∂ₛg(·, c)` is strictly decreasing and changes sign on `[sL, sR]`, so its unique
zero is interior and is the maximiser of `g(·, c)`.
-/

open Set

namespace RosslerTangency

/-! ### One-variable lemmas, proved from the mean value theorem -/

section OneVariable

variable {f f' : ℝ → ℝ} {a b : ℝ}

/-- Mean value theorem on a subinterval of `[a, b]`, with the intermediate point in `[a, b]`. -/
theorem exists_mvt_Icc (hf : ∀ x ∈ Icc a b, HasDerivAt f (f' x) x) {x y : ℝ}
    (hx : x ∈ Icc a b) (hy : y ∈ Icc a b) (hxy : x < y) :
    ∃ ξ ∈ Icc a b, f y - f x = f' ξ * (y - x) := by
  obtain ⟨ξ, hξ, hξeq⟩ := exists_hasDerivAt_eq_slope f f' hxy
    (fun z hz => (hf z ⟨hx.1.trans hz.1, hz.2.trans hy.2⟩).continuousAt.continuousWithinAt)
    (fun z hz => hf z ⟨hx.1.trans hz.1.le, hz.2.le.trans hy.2⟩)
  refine ⟨ξ, ⟨hx.1.trans hξ.1.le, hξ.2.le.trans hy.2⟩, ?_⟩
  rw [hξeq, div_mul_cancel₀ _ (sub_ne_zero.mpr hxy.ne')]

theorem strictMonoOn_Icc_of_hasDerivAt (hf : ∀ x ∈ Icc a b, HasDerivAt f (f' x) x)
    (hpos : ∀ x ∈ Icc a b, 0 < f' x) : StrictMonoOn f (Icc a b) := by
  intro x hx y hy hxy
  obtain ⟨ξ, hξ, h⟩ := exists_mvt_Icc hf hx hy hxy
  have := mul_pos (hpos ξ hξ) (sub_pos.mpr hxy)
  linarith

theorem strictAntiOn_Icc_of_hasDerivAt (hf : ∀ x ∈ Icc a b, HasDerivAt f (f' x) x)
    (hneg : ∀ x ∈ Icc a b, f' x < 0) : StrictAntiOn f (Icc a b) := by
  intro x hx y hy hxy
  obtain ⟨ξ, hξ, h⟩ := exists_mvt_Icc hf hx hy hxy
  have := mul_pos (neg_pos.mpr (hneg ξ hξ)) (sub_pos.mpr hxy)
  linarith

theorem monotoneOn_Icc_of_hasDerivAt (hf : ∀ x ∈ Icc a b, HasDerivAt f (f' x) x)
    (hnn : ∀ x ∈ Icc a b, 0 ≤ f' x) : MonotoneOn f (Icc a b) := by
  intro x hx y hy hxy
  rcases eq_or_lt_of_le hxy with rfl | hlt
  · exact le_rfl
  obtain ⟨ξ, hξ, h⟩ := exists_mvt_Icc hf hx hy hlt
  have := mul_nonneg (hnn ξ hξ) (sub_pos.mpr hlt).le
  linarith

theorem antitoneOn_Icc_of_hasDerivAt (hf : ∀ x ∈ Icc a b, HasDerivAt f (f' x) x)
    (hnp : ∀ x ∈ Icc a b, f' x ≤ 0) : AntitoneOn f (Icc a b) := by
  intro x hx y hy hxy
  rcases eq_or_lt_of_le hxy with rfl | hlt
  · exact le_rfl
  obtain ⟨ξ, hξ, h⟩ := exists_mvt_Icc hf hx hy hlt
  have := mul_nonneg (neg_nonneg.mpr (hnp ξ hξ)) (sub_pos.mpr hlt).le
  linarith

/-- A function whose derivative is `≥ 0` on `[a, m]` and `≤ 0` on `[m, b]` attains its maximum
over `[a, b]` at `m`. -/
theorem le_of_hasDerivAt_sign {m : ℝ} (hf : ∀ x ∈ Icc a b, HasDerivAt f (f' x) x)
    (hm : m ∈ Icc a b) (hl : ∀ x ∈ Icc a m, 0 ≤ f' x) (hr : ∀ x ∈ Icc m b, f' x ≤ 0) :
    ∀ x ∈ Icc a b, f x ≤ f m := by
  intro x hx
  rcases le_total x m with hxm | hmx
  · have hmono : MonotoneOn f (Icc a m) :=
      monotoneOn_Icc_of_hasDerivAt (fun z hz => hf z ⟨hz.1, hz.2.trans hm.2⟩) hl
    exact hmono ⟨hx.1, hxm⟩ ⟨hm.1, le_rfl⟩ hxm
  · have hanti : AntitoneOn f (Icc m b) :=
      antitoneOn_Icc_of_hasDerivAt (fun z hz => hf z ⟨hm.1.trans hz.1, hz.2⟩) hr
    exact hanti ⟨le_rfl, hm.2⟩ ⟨hmx, hx.2⟩ hmx

/-- One-variable core of the criterion. If `f'' < 0` on `[a, b]` and `f'(a) > 0 > f'(b)`, then
`f'` has a unique zero `m` in `[a, b]`; it lies in `(a, b)` and maximises `f` over `[a, b]`. -/
theorem exists_unique_crit_of_concave {f'' : ℝ → ℝ} (hab : a ≤ b)
    (hf : ∀ x ∈ Icc a b, HasDerivAt f (f' x) x) (hf' : ∀ x ∈ Icc a b, HasDerivAt f' (f'' x) x)
    (hneg : ∀ x ∈ Icc a b, f'' x < 0) (ha : 0 < f' a) (hb : f' b < 0) :
    ∃ m ∈ Ioo a b, f' m = 0 ∧ (∀ x ∈ Icc a b, f' x = 0 → x = m) ∧
      ∀ x ∈ Icc a b, f x ≤ f m := by
  have hanti : StrictAntiOn f' (Icc a b) := strictAntiOn_Icc_of_hasDerivAt hf' hneg
  have hcont : ContinuousOn f' (Icc a b) := fun x hx =>
    (hf' x hx).continuousAt.continuousWithinAt
  obtain ⟨m, hm, hm0⟩ : (0 : ℝ) ∈ f' '' Icc a b :=
    intermediate_value_Icc' hab hcont ⟨hb.le, ha.le⟩
  have hma : a < m := lt_of_le_of_ne hm.1 (by rintro rfl; linarith)
  have hmb : m < b := lt_of_le_of_ne hm.2 (by rintro rfl; linarith)
  refine ⟨m, ⟨hma, hmb⟩, hm0, fun x hx hx0 => hanti.injOn hx hm (hx0.trans hm0.symm), ?_⟩
  refine le_of_hasDerivAt_sign hf hm (fun x hx => ?_) (fun x hx => ?_)
  · have hxI : x ∈ Icc a b := ⟨hx.1, hx.2.trans hm.2⟩
    rw [← hm0]
    exact hanti.antitoneOn hxI hm hx.2
  · have hxI : x ∈ Icc a b := ⟨hm.1.trans hx.1, hx.2⟩
    rw [← hm0]
    exact hanti.antitoneOn hm hxI hx.1

end OneVariable

/-! ### Remark 4.6: the Taylor bound -/

/-- **Remark 4.6.** If `h'' ≤ -γ < 0` on `[sL, sR]` and `smid ∈ [sL, sR]`, then for all
`s ∈ [sL, sR]`, `h s ≤ h smid + (h' smid)^2 / (2γ)`. (We also record the intermediate quadratic
Taylor bound `h s ≤ h smid + h' smid (s - smid) - γ/2 (s - smid)^2`.) -/
theorem taylor_quadratic_bound {h h' h'' : ℝ → ℝ} {sL sR smid γ : ℝ}
    (hh : ∀ s ∈ Icc sL sR, HasDerivAt h (h' s) s)
    (hh' : ∀ s ∈ Icc sL sR, HasDerivAt h' (h'' s) s)
    (hbound : ∀ s ∈ Icc sL sR, h'' s ≤ -γ) (hsm : smid ∈ Icc sL sR) :
    ∀ s ∈ Icc sL sR,
      h s ≤ h smid + h' smid * (s - smid) - γ / 2 * (s - smid) ^ 2 := by
  have hlin : ∀ s : ℝ, HasDerivAt (fun t : ℝ => t - smid) 1 s := fun s =>
    (hasDerivAt_id' (x := s)).sub_const smid
  -- `D t = h t - h smid - h' smid (t - smid) + γ/2 (t - smid)^2` has derivative
  -- `E t = h' t - h' smid + γ (t - smid)`, and `E' = h'' + γ ≤ 0`.
  have hD : ∀ s ∈ Icc sL sR, HasDerivAt
      (fun t => h t - h smid - h' smid * (t - smid) + γ / 2 * ((t - smid) * (t - smid)))
      (h' s - h' smid + γ * (s - smid)) s := by
    intro s hs
    have := (((hh s hs).sub_const (h smid)).sub ((hlin s).const_mul (h' smid))).add
      (((hlin s).mul (hlin s)).const_mul (γ / 2))
    exact this.congr_deriv (by ring)
  have hE : ∀ s ∈ Icc sL sR, HasDerivAt (fun t => h' t - h' smid + γ * (t - smid))
      (h'' s + γ) s := by
    intro s hs
    have := ((hh' s hs).sub_const (h' smid)).add ((hlin s).const_mul γ)
    exact this.congr_deriv (by ring)
  have hEanti : AntitoneOn (fun t => h' t - h' smid + γ * (t - smid)) (Icc sL sR) :=
    antitoneOn_Icc_of_hasDerivAt hE (fun s hs => by linarith [hbound s hs])
  have hDmax := le_of_hasDerivAt_sign hD hsm
    (fun x hx => by
      have := hEanti ⟨hx.1, hx.2.trans hsm.2⟩ hsm hx.2
      simp only [sub_self, mul_zero, add_zero] at this
      exact this)
    (fun x hx => by
      have := hEanti hsm ⟨hsm.1.trans hx.1, hx.2⟩ hx.1
      simp only [sub_self, mul_zero, add_zero] at this
      exact this)
  intro s hs
  have h1 := hDmax s hs
  simp only [sub_self, mul_zero, add_zero] at h1
  nlinarith [h1]

/-- **Remark 4.6 (Taylor bound).** If `h` has derivative `h'`, `h'` has derivative `h''`,
`h'' ≤ -γ` with `γ > 0` on `[sL, sR]`, and `smid ∈ [sL, sR]`, then
`h s ≤ h smid + (h' smid)^2 / (2γ)` for all `s ∈ [sL, sR]`. -/
theorem taylor_bound {h h' h'' : ℝ → ℝ} {sL sR smid γ : ℝ}
    (hh : ∀ s ∈ Icc sL sR, HasDerivAt h (h' s) s)
    (hh' : ∀ s ∈ Icc sL sR, HasDerivAt h' (h'' s) s)
    (hγ : 0 < γ) (hbound : ∀ s ∈ Icc sL sR, h'' s ≤ -γ) (hsm : smid ∈ Icc sL sR) :
    ∀ s ∈ Icc sL sR, h s ≤ h smid + h' smid ^ 2 / (2 * γ) := by
  intro s hs
  have h1 := taylor_quadratic_bound hh hh' hbound hsm s hs
  have h2 : h' smid * (s - smid) - γ / 2 * (s - smid) ^ 2 ≤ h' smid ^ 2 / (2 * γ) := by
    rw [le_div_iff₀ (by linarith)]
    nlinarith [sq_nonneg (h' smid - γ * (s - smid))]
  linarith

/-! ### Lemma 4.5: the tangency criterion -/

/-- The maximum function `G c = max_{s ∈ [sL, sR]} g(s, c)`, written as a supremum. -/
noncomputable def maxFun (g : ℝ → ℝ → ℝ) (sL sR c : ℝ) : ℝ :=
  sSup ((fun s => g s c) '' Icc sL sR)

theorem maxFun_eq_of_le {g : ℝ → ℝ → ℝ} {sL sR c m : ℝ} (hm : m ∈ Icc sL sR)
    (hmax : ∀ x ∈ Icc sL sR, g x c ≤ g m c) : maxFun g sL sR c = g m c :=
  IsGreatest.csSup_eq ⟨⟨m, hm, rfl⟩, by rintro y ⟨x, hx, rfl⟩; exact hmax x hx⟩

/-- `G` is continuous when `g` is jointly continuous. -/
theorem continuous_maxFun {g : ℝ → ℝ → ℝ} (sL sR : ℝ)
    (hcont : Continuous fun p : ℝ × ℝ => g p.1 p.2) : Continuous (maxFun g sL sR) := by
  have : Continuous fun p : ℝ × ℝ => g p.2 p.1 := hcont.comp continuous_swap
  exact isCompact_Icc.continuous_sSup (f := fun c s => g s c) this

/-- `G` is continuous on `[cL, cR]` when `g` is jointly continuous on the rectangle
`[sL, sR] × [cL, cR]` (proved by extending `g` continuously to `ℝ × ℝ` via `projIcc`). -/
theorem continuousOn_maxFun {g : ℝ → ℝ → ℝ} {sL sR cL cR : ℝ} (hs : sL ≤ sR) (hc : cL ≤ cR)
    (hcont : ContinuousOn (fun p : ℝ × ℝ => g p.1 p.2) (Icc sL sR ×ˢ Icc cL cR)) :
    ContinuousOn (maxFun g sL sR) (Icc cL cR) := by
  set g' : ℝ → ℝ → ℝ := fun s c => g (projIcc sL sR hs s) (projIcc cL cR hc c) with hg'def
  have hg' : Continuous fun p : ℝ × ℝ => g' p.1 p.2 := by
    have hpr : Continuous fun p : ℝ × ℝ =>
        ((projIcc sL sR hs p.1 : ℝ), (projIcc cL cR hc p.2 : ℝ)) :=
      ((continuous_subtype_val.comp continuous_projIcc).comp continuous_fst).prodMk
        ((continuous_subtype_val.comp continuous_projIcc).comp continuous_snd)
    exact hcont.comp_continuous hpr (fun p => ⟨(projIcc sL sR hs p.1).2, (projIcc cL cR hc p.2).2⟩)
  have heq : EqOn (maxFun g' sL sR) (maxFun g sL sR) (Icc cL cR) := by
    intro c hc'
    unfold maxFun
    congr 1
    refine image_congr fun s hs' => ?_
    simp only [hg'def, projIcc_of_mem _ hs', projIcc_of_mem _ hc']
  exact (continuous_maxFun sL sR hg').continuousOn.congr heq.symm

section Criterion

variable {g gs gss gc : ℝ → ℝ → ℝ} {sL sR cL cR : ℝ}

/-- For every `c ∈ [cL, cR]`, `∂ₛg(·, c)` has a unique zero `m` in `[sL, sR]`; it lies in
`(sL, sR)` and maximises `g(·, c)` on `[sL, sR]`, so that `G c = g m c`. -/
theorem crit_point
    (hg : ∀ s ∈ Icc sL sR, ∀ c ∈ Icc cL cR, HasDerivAt (fun s => g s c) (gs s c) s)
    (hgs : ∀ s ∈ Icc sL sR, ∀ c ∈ Icc cL cR, HasDerivAt (fun s => gs s c) (gss s c) s)
    (hs : sL < sR)
    (hA : ∀ s ∈ Icc sL sR, ∀ c ∈ Icc cL cR, gss s c < 0)
    (hB : ∀ c ∈ Icc cL cR, 0 < gs sL c ∧ gs sR c < 0) {c : ℝ} (hc : c ∈ Icc cL cR) :
    ∃ m ∈ Ioo sL sR, gs m c = 0 ∧ (∀ x ∈ Icc sL sR, gs x c = 0 → x = m) ∧
      (∀ x ∈ Icc sL sR, g x c ≤ g m c) ∧ maxFun g sL sR c = g m c := by
  obtain ⟨m, hm, hm0, huniq, hmax⟩ :=
    exists_unique_crit_of_concave (f := fun s => g s c) (f' := fun s => gs s c)
      (f'' := fun s => gss s c) hs.le (fun s hs' => hg s hs' c hc) (fun s hs' => hgs s hs' c hc)
      (fun s hs' => hA s hs' c hc) (hB c hc).1 (hB c hc).2
  exact ⟨m, hm, hm0, huniq, hmax, maxFun_eq_of_le (Ioo_subset_Icc_self hm) hmax⟩

/-- `G` is strictly increasing on `[cL, cR]`. -/
theorem maxFun_strictMonoOn
    (hg : ∀ s ∈ Icc sL sR, ∀ c ∈ Icc cL cR, HasDerivAt (fun s => g s c) (gs s c) s)
    (hgs : ∀ s ∈ Icc sL sR, ∀ c ∈ Icc cL cR, HasDerivAt (fun s => gs s c) (gss s c) s)
    (hgc : ∀ s ∈ Icc sL sR, ∀ c ∈ Icc cL cR, HasDerivAt (fun c => g s c) (gc s c) c)
    (hs : sL < sR)
    (hA : ∀ s ∈ Icc sL sR, ∀ c ∈ Icc cL cR, gss s c < 0)
    (hB : ∀ c ∈ Icc cL cR, 0 < gs sL c ∧ gs sR c < 0)
    (hD : ∀ s ∈ Icc sL sR, ∀ c ∈ Icc cL cR, 0 < gc s c) :
    StrictMonoOn (maxFun g sL sR) (Icc cL cR) := by
  intro c₁ hc₁ c₂ hc₂ h12
  obtain ⟨m₁, hm₁, -, -, -, hG₁⟩ := crit_point hg hgs hs hA hB hc₁
  obtain ⟨m₂, -, -, -, hmax₂, hG₂⟩ := crit_point hg hgs hs hA hB hc₂
  have hm₁' : m₁ ∈ Icc sL sR := Ioo_subset_Icc_self hm₁
  -- `c ↦ g m₁ c` is strictly increasing on `[cL, cR]` by (D).
  have hmono : StrictMonoOn (fun c => g m₁ c) (Icc cL cR) :=
    strictMonoOn_Icc_of_hasDerivAt (fun c hc => hgc m₁ hm₁' c hc) (fun c hc => hD m₁ hm₁' c hc)
  rw [hG₁, hG₂]
  exact (hmono hc₁ hc₂ h12).trans_le (hmax₂ m₁ hm₁')

/-- **Lemma 4.5 (tangency criterion).** All hypotheses are local to the closed rectangle
`[sL, sR] × [cL, cR]`: the partial derivatives are only required at its points (as two-sided
derivatives of the functions `g`, `gs` on `ℝ × ℝ`), and `g` is only required to be jointly
continuous on it.

Conclusion: there are `c⋆ ∈ (cL, cR)` and `s⋆ ∈ (sL, sR)` with `g = ∂ₛg = 0`, `∂ₛ²g < 0`,
`∂_c g > 0` at `(s⋆, c⋆)`; for each `c ∈ [cL, cR]` the zero of `∂ₛg(·, c)` in `[sL, sR]` is
unique; `c⋆` is the unique `c ∈ [cL, cR]` with `G c = max_{[sL,sR]} g(·, c) = 0`; and `(s⋆, c⋆)` is
the unique point of the rectangle where `g = ∂ₛg = 0`. -/
theorem tangency_criterion
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
      -- the tangency at `(s⋆, c⋆)`
      g sstar cstar = 0 ∧ gs sstar cstar = 0 ∧ gss sstar cstar < 0 ∧ 0 < gc sstar cstar ∧
      -- uniqueness of the critical point of `g(·, c)` in `[sL, sR]`, for every `c`
      (∀ c ∈ Icc cL cR, ∃! s, s ∈ Icc sL sR ∧ gs s c = 0) ∧
      -- `c⋆` is the unique zero in `[cL, cR]` of `G c = sSup (g · c '' [sL, sR])`
      maxFun g sL sR cstar = 0 ∧ (∀ c ∈ Icc cL cR, maxFun g sL sR c = 0 → c = cstar) ∧
      -- `(s⋆, c⋆)` is the unique point of the rectangle where `g = ∂ₛg = 0`
      (∀ c ∈ Icc cL cR, ∀ s ∈ Icc sL sR, g s c = 0 → gs s c = 0 → c = cstar ∧ s = sstar) := by
  have hcL : cL ∈ Icc cL cR := left_mem_Icc.mpr hc.le
  have hcR : cR ∈ Icc cL cR := right_mem_Icc.mpr hc.le
  have hGmono : StrictMonoOn (maxFun g sL sR) (Icc cL cR) :=
    maxFun_strictMonoOn hg hgs hgc hs hA hB hD
  have hGcont : ContinuousOn (maxFun g sL sR) (Icc cL cR) :=
    continuousOn_maxFun hs.le hc.le hcont
  -- `G cL < 0 < G cR`
  have hGL : maxFun g sL sR cL < 0 := by
    obtain ⟨m, hm, -, -, -, hGm⟩ := crit_point hg hgs hs hA hB hcL
    rw [hGm]
    exact hC.1 m (Ioo_subset_Icc_self hm)
  have hGR : 0 < maxFun g sL sR cR := by
    obtain ⟨m, -, -, -, hmax, hGm⟩ := crit_point hg hgs hs hA hB hcR
    obtain ⟨s₀, hs₀, hpos⟩ := hC.2
    rw [hGm]
    exact hpos.trans_le (hmax s₀ hs₀)
  -- intermediate value theorem
  obtain ⟨cstar, hcstar, hGc⟩ : (0 : ℝ) ∈ maxFun g sL sR '' Ioo cL cR :=
    intermediate_value_Ioo hc.le hGcont ⟨hGL, hGR⟩
  have hcstar' : cstar ∈ Icc cL cR := Ioo_subset_Icc_self hcstar
  obtain ⟨sstar, hsstar, hgs0, huniq, -, hGs⟩ := crit_point hg hgs hs hA hB hcstar'
  have hsstar' : sstar ∈ Icc sL sR := Ioo_subset_Icc_self hsstar
  have hGinj : ∀ c ∈ Icc cL cR, maxFun g sL sR c = 0 → c = cstar := fun c hc0 hGc0 =>
    hGmono.injOn hc0 hcstar' (hGc0.trans hGc.symm)
  refine ⟨cstar, hcstar, sstar, hsstar, hGs.symm.trans hGc, hgs0,
    hA sstar hsstar' cstar hcstar', hD sstar hsstar' cstar hcstar', ?_, hGc, hGinj, ?_⟩
  · intro c hc0
    obtain ⟨m, hm, hm0, hmuniq, -, -⟩ := crit_point hg hgs hs hA hB hc0
    exact ⟨m, ⟨Ioo_subset_Icc_self hm, hm0⟩, fun x hx => hmuniq x hx.1 hx.2⟩
  · intro c hc0 s hs0 hg0 hgs0'
    obtain ⟨m, -, -, hmuniq, -, hGm⟩ := crit_point hg hgs hs hA hB hc0
    have hsm : s = m := hmuniq s hs0 hgs0'
    have hcc : c = cstar := hGinj c hc0 (by rw [hGm, ← hsm]; exact hg0)
    subst hcc
    exact ⟨rfl, huniq s hs0 hgs0'⟩

/-- **Lemma 4.5**, with the hypotheses exactly as in the informal statement: all partial
derivatives exist everywhere. -/
theorem tangency_criterion'
    (hg : ∀ s c, HasDerivAt (fun s => g s c) (gs s c) s)
    (hgs : ∀ s c, HasDerivAt (fun s => gs s c) (gss s c) s)
    (hgc : ∀ s c, HasDerivAt (fun c => g s c) (gc s c) c)
    (hcont : Continuous fun p : ℝ × ℝ => g p.1 p.2)
    (hs : sL < sR) (hc : cL < cR)
    (hA : ∀ s ∈ Icc sL sR, ∀ c ∈ Icc cL cR, gss s c < 0)
    (hB : ∀ c ∈ Icc cL cR, 0 < gs sL c ∧ gs sR c < 0)
    (hC : (∀ s ∈ Icc sL sR, g s cL < 0) ∧ ∃ s₀ ∈ Icc sL sR, 0 < g s₀ cR)
    (hD : ∀ s ∈ Icc sL sR, ∀ c ∈ Icc cL cR, 0 < gc s c) :
    ∃ cstar ∈ Ioo cL cR, ∃ sstar ∈ Ioo sL sR,
      g sstar cstar = 0 ∧ gs sstar cstar = 0 ∧ gss sstar cstar < 0 ∧ 0 < gc sstar cstar ∧
      (∀ c ∈ Icc cL cR, ∃! s, s ∈ Icc sL sR ∧ gs s c = 0) ∧
      maxFun g sL sR cstar = 0 ∧ (∀ c ∈ Icc cL cR, maxFun g sL sR c = 0 → c = cstar) ∧
      (∀ c ∈ Icc cL cR, ∀ s ∈ Icc sL sR, g s c = 0 → gs s c = 0 → c = cstar ∧ s = sstar) :=
  tangency_criterion (fun s _ c _ => hg s c) (fun s _ c _ => hgs s c) (fun s _ c _ => hgc s c)
    hcont.continuousOn hs hc hA hB hC hD

end Criterion

end RosslerTangency

#print axioms RosslerTangency.tangency_criterion
#print axioms RosslerTangency.tangency_criterion'
#print axioms RosslerTangency.taylor_bound
#print axioms RosslerTangency.taylor_quadratic_bound
