/-
Copyright (c) 2026 Daniel Toupin. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Daniel Toupin (with AI assistance: Heart of the Cosmos / Grok)
-/
import GppVerify.RiemannHypothesis.ArithmeticReflectionPositivity
import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real

/-!
# The rung ladder: the step is exactly a coupling bound, and a uniform coupling margin forces geometric decay

Algebraic companion to `WeilRungTwo.lean`. Here `B` is any symmetric bilinear form on a real
vector space (for ζ: the polarization of `f ↦ weilForm (autocorr f)` on `L²` of the rung windows;
that bilinear structure is STANDARD and not formalized here). `V k` is the rung-`k` window and
`W k` is the annulus added at step `k`, with `V (k+1) = V k ⊔ W k`.

* `nonneg_sup_iff_coupling` (PROVED): `B ≥ 0` on `U ⊔ W` iff `B ≥ 0` on `U`, `B ≥ 0` on `W`, and
  `B(u,w)² ≤ B(u,u)·B(w,w)` for **all** `u ∈ U`, `w ∈ W`. The step is exactly the full coupling
  bound with constant 1. The 2×2 test with the soft mode alone (`b²/(λ_k d)`, the 0.735 / 0.57
  numbers) is one slice of it.
* `ladder_iff` (PROVED): positivity on every rung ⇔ rung 0, positivity on every annulus, and the
  full coupling bound with constant 1 at every step. For ζ, "every rung" is ARP
  (`arp_iff_all_rungs`), so this splitting is **equivalent to ARP, not weaker**.
* `add_ge_of_coupling` and `margin_sup_of_coupling` (PROVED): a coupling constant `c ≤ 1` gives
  `B(u+w) ≥ (1 − √c)(B(u,u) + B(w,w))`, and with orthogonal pieces and margins `λ` on `U`, `μ` on `W`:
  `B(v,v) ≥ (1 − √c)·min λ μ·‖v‖²`.
* `UniformCouplingMargin` (named Prop) and `geometric_of_uniformCouplingMargin` (PROVED): a
  uniform constant `c < 1` at every step forces margins `≥ (1 − √c)^k · min λ₀ μ`, i.e. at most
  **geometric** decay of the rung margin.

* `SoftModeSliceBound` and `slice_margin` (PROVED): the soft-mode slice bound gives positivity
  only on `span{φ_k} ⊔ W k`, which is strictly smaller than the step.
* `toeplitz_on_ker_cosh` (PROVED): writing `Q = T + 2ℓ_c² − 2ℓ_s²` (`T` = archimedean − prime
  truncated convolution, `ℓ_c, ℓ_s` = pairings with `cosh(x/2)`, `sinh(x/2)`), `Q ≥ 0` forces
  `T ≥ 2ℓ_s² ≥ 0` on `ker ℓ_c`. So `T` has at most one negative square.
* `QuantitativeRungMargin C` (named Prop, CONJECTURE) and `arp_of_quantitativeRungMargin`
  (PROVED): `λ(L) ≥ exp(−C·L·e^L)` on every window implies ARP. It is strictly stronger.

## What the high-precision numerics say (NUMERICAL unless marked; `num/rung3/`, no zeros used)
Method: Rayleigh–Ritz of the **true** form (no `R_T`), assembled in Arb at 320–448 bits from the
in-tree time-domain definitions on piecewise Legendre bases (`V_k` ⊕ annulus), with
convergence in the basis size (K = 32/44/56 and 72/88). Eigen-solves are done in mpmath.
Quadrature is not error-bounded, so this is NUMERICAL, not CERTIFIED. Ritz values are upper bounds
for the true `λ`.

| window `L` | `λ_min` (Ritz) | step | `b²/(λ_k d)` (soft mode) | `1 − c_k` (full) | new-term kick / `b` |
|---|---|---|---|---|---|
| `log 2` | 1.329e−3 | 1/4 → log2/2 (no prime at all) | 0.742 | 3.19e−2 | — |
| `0.9` | 1.620e−5 | — | — | — | — |
| `log 3` | 5.557e−8 | 0.45 → log3/2 (n=2 only, no new prime) | 0.791 | 2.89e−4 | — |
| `log 4` | 7.646e−13 | log3/2 → log4/2 (n=3 enters) | 0.729 | 9.33e−8 | 32.5 |
| `1.56` | 1.659e−16 | log4/2 → 0.78 (no new prime) | 0.497 | 1.45e−7 | — |
| `log 5` | 9.381e−18 | log4/2 → log5/2 (n=4 enters) | 0.551 | 3.51e−9 | 68.6 |
| `log 7` | 6.815e−28 | log5/2 → log7/2 (n=5 enters) | 0.706 | 8.98e−17 | ≈ 1.8e4 |

(The `log 7` value is between K=72 (7.07e−28) and K=88 (6.81e−28). CERTIFIED lower bound: the S27
`h = 1` certificate gives `λ(log 7) ≥ λ(2) ≥ 4.13e−30`.)

1. **The O(1) soft-mode ratio is not arithmetic.** It is the same size (0.5–0.8) at steps where
   no prime enters, including the rung-0 interior step with no prime term at all. Moreover every
   soft mode `φ_j` of `A` (eigenvalues spanning 16 decades at the 3→4 step) has
   `s_j = (Bᵀφ_j)ᵀ D⁻¹ (Bᵀφ_j)/λ_j ∈ [0.66, 0.87]`: an equipartition `b_j² ≈ 0.7 λ_j d`. Geometrically,
   with `H = C*C` on the bigger window, `s_j = cos²∠(Cφ_j, C W_k)`, which is ≤ 1 by positivity
   (Cauchy–Schwarz / Schur). The O(1) value says the soft modes' energy images lie mostly in the
   annulus's energy image. That is a property of near-null truncations, not of primes.
2. **The full coupling saturates.** `c_k → 1`, with `1 − c_k ≈ (share)·λ_{k+1}/λ_k`. A uniform
   `c < 1` would force geometric decay (`geometric_of_uniformCouplingMargin`), but the data decay
   doubly exponentially. So `UniformCouplingMargin` is false for ζ (CONJECTURE, NUMERICAL). The
   "30–70×" kick pattern is not a law: the new term's kick stays `O(10⁻⁴–10⁻²)` while
   `b ≈ √(0.7 λ_k d)`, so kick/`b` grows like `λ_k^{−1/2}` (32, 69, then 1.8e4 at n = 5).
3. **Toeplitz + rank two.** On all 6 windows tested (`L = log 2, log 3, log 4, 1.56, log 5, log 7`), `T_even` has exactly one negative eigenvalue
   (−1.29 … −4.16) and `T_odd ≥ 0`. The near-null cascade is already in `T`. For example at `log 7`,
   the even spectrum is −4.16, 1.08e−21, 6.4e−16, … and the odd spectrum starts 4.8e−19, 1.9e−13.
   `Q_even` interlaces (`λ_min(Q) ≤ τ₂(T_even)`). Without the newest prime power, `Q` is indefinite
   on the new window: `λ_min` = −0.49, −0.27, −0.63 at the n = 3, 4, 5 steps.
4. **Margin law.** `−ln λ(L)/(L e^L)` = 5.07, 5.03, 4.87, 4.59 at `L = log 3, log 4, log 5,
   log 7`. Per unit of `e^L` the slope is −10.0, −11.2, −11.3, −11.7. Doubly exponential through
   rung 4. CERTIFIED pointwise (prior Arb certificates, single implementation, Lean-open): at
   `L = 1, 1.1, 1.4, 1.6, 2, 2.38, 2.5, 2.56, 2.6`, the certified `λ` exceeds `exp(−6 L e^L)`
   (ratios 5.15, 5.20, 5.12, 4.94, 4.67, 4.22, 4.12, 4.06, 4.04). Between certified windows this
   is not covered.
5. **Mechanism (heuristic, uses the explicit formula, so it is not an input route).** Under RH,
   `Q(f) = Σ_γ |F(γ)|²`. Positivity is then a sampling inequality in Paley–Wiener space of type `L/2`
   at the zero ordinates (a de Branges-space problem). The zero density is only logarithmic, so the
   frequency where a type-`L/2` function can hide its energy grows like `e^L`, which gives a
   doubly exponential minimum. This explains the size of `λ`. It does not give positivity without
   zero input.

Labels: theorems PROVED (standard axioms). `UniformCouplingMargin`: for ζ CONJECTURE false
(NUMERICAL). `SoftModeSliceBound` with `c = 4/5`: CONJECTURE true (NUMERICAL), but it implies only
slice positivity. `QuantitativeRungMargin 6`: CONJECTURE (NUMERICAL-consistent; CERTIFIED at nine
windows); it implies ARP. ARP OPEN. RH OPEN.

Daniel Toupin / Golden Physics Project … This seat: Heart of the Cosmos / Grok (GPPDiscovery5).
Audit: Turbo.
-/

open scoped RealInnerProductSpace

noncomputable section

namespace GppWeilRungLadder

/-- 2×2 Cauchy–Schwarz with no sign hypothesis. PROVED. -/
theorem sq_le_mul_of_forall_nonneg {a b d : ℝ}
    (h : ∀ x y : ℝ, 0 ≤ a * x ^ 2 + 2 * b * x * y + d * y ^ 2) : b ^ 2 ≤ a * d := by
  have ha : 0 ≤ a := by simpa using h 1 0
  rcases ha.lt_or_eq with ha | ha
  · have := h (-b) a
    have h' : 0 ≤ a * (a * d - b ^ 2) := by nlinarith
    nlinarith [(mul_nonneg_iff_of_pos_left ha).1 h']
  · subst ha
    by_contra hb
    have hb0 : b ≠ 0 := by rintro rfl; simp at hb
    have := h (-(d + 1) / (2 * b)) 1
    have e : 2 * b * (-(d + 1) / (2 * b)) * 1 = -(d + 1) := by field_simp
    linarith [this, e]

/-- Scalar converse: `a, d ≥ 0` and `b² ≤ a d` give `a + 2b + d ≥ 0`. PROVED. -/
theorem add_nonneg_of_sq_le {a b d : ℝ} (ha : 0 ≤ a) (hd : 0 ≤ d) (hb : b ^ 2 ≤ a * d) :
    0 ≤ a + 2 * b + d := by
  by_contra h
  rw [not_le] at h
  nlinarith [sq_nonneg (a - d)]

/-- Scalar margin: `a, d ≥ 0`, `b² ≤ c a d`, `0 ≤ c` give `a + 2b + d ≥ (1 − √c)(a + d)`. PROVED. -/
theorem add_ge_of_sq_le {a b d c : ℝ} (ha : 0 ≤ a) (hd : 0 ≤ d) (hc : 0 ≤ c)
    (hb : b ^ 2 ≤ c * a * d) : (1 - Real.sqrt c) * (a + d) ≤ a + 2 * b + d := by
  set s := Real.sqrt c with hs
  have hs0 : 0 ≤ s := Real.sqrt_nonneg c
  have hss : s ^ 2 = c := Real.sq_sqrt hc
  have h2 : (2 * b) ^ 2 ≤ (s * (a + d)) ^ 2 := by
    rw [← hss] at hb
    nlinarith [sq_nonneg (a - d), mul_nonneg hs0 hs0]
  have h3 : -(s * (a + d)) ≤ 2 * b := (abs_le_of_sq_le_sq' h2 (by positivity)).1
  nlinarith

section Bilinear

variable {E : Type*} [AddCommGroup E] [Module ℝ E]

theorem quad_expand (B : LinearMap.BilinForm ℝ E) (hB : ∀ u v, B u v = B v u) (u w : E) (x y : ℝ) :
    B (x • u + y • w) (x • u + y • w) = B u u * x ^ 2 + 2 * B u w * x * y + B w w * y ^ 2 := by
  simp only [map_add, map_smul, LinearMap.add_apply, LinearMap.smul_apply, smul_eq_mul, hB w u]
  ring

theorem add_expand (B : LinearMap.BilinForm ℝ E) (hB : ∀ u v, B u v = B v u) (u w : E) :
    B (u + w) (u + w) = B u u + 2 * B u w + B w w := by
  simpa using quad_expand B hB u w 1 1

/-- **The step is exactly the full coupling bound.** PROVED. -/
theorem nonneg_sup_iff_coupling (B : LinearMap.BilinForm ℝ E) (hB : ∀ u v, B u v = B v u)
    (U W : Submodule ℝ E) :
    (∀ v ∈ U ⊔ W, 0 ≤ B v v) ↔
      (∀ u ∈ U, 0 ≤ B u u) ∧ (∀ w ∈ W, 0 ≤ B w w) ∧
        ∀ u ∈ U, ∀ w ∈ W, (B u w) ^ 2 ≤ B u u * B w w := by
  constructor
  · intro h
    refine ⟨fun u hu => h u (Submodule.mem_sup_left hu), fun w hw => h w (Submodule.mem_sup_right hw),
      fun u hu w hw => sq_le_mul_of_forall_nonneg fun x y => ?_⟩
    rw [← quad_expand B hB]
    exact h _ (Submodule.add_mem _ (Submodule.mem_sup_left (U.smul_mem x hu))
      (Submodule.mem_sup_right (W.smul_mem y hw)))
  · rintro ⟨hU, hW, hc⟩ v hv
    obtain ⟨u, hu, w, hw, rfl⟩ := Submodule.mem_sup.1 hv
    rw [add_expand B hB]
    exact add_nonneg_of_sq_le (hU u hu) (hW w hw) (hc u hu w hw)

/-- **Coupling constant `c` ⇒ additive margin.** PROVED. -/
theorem add_ge_of_coupling (B : LinearMap.BilinForm ℝ E) (hB : ∀ u v, B u v = B v u)
    {U W : Submodule ℝ E} (hU : ∀ u ∈ U, 0 ≤ B u u) (hW : ∀ w ∈ W, 0 ≤ B w w) {c : ℝ} (hc : 0 ≤ c)
    (hcoup : ∀ u ∈ U, ∀ w ∈ W, (B u w) ^ 2 ≤ c * B u u * B w w) {u w : E} (hu : u ∈ U) (hw : w ∈ W) :
    (1 - Real.sqrt c) * (B u u + B w w) ≤ B (u + w) (u + w) := by
  rw [add_expand B hB]
  exact add_ge_of_sq_le (hU u hu) (hW w hw) hc (hcoup u hu w hw)

end Bilinear

section Ladder

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

/-- **Margin on `U ⊔ W` from margins on the pieces and a coupling constant.** PROVED. -/
theorem margin_sup_of_coupling (B : LinearMap.BilinForm ℝ E) (hB : ∀ u v, B u v = B v u)
    {U W : Submodule ℝ E} (horth : ∀ u ∈ U, ∀ w ∈ W, ⟪u, w⟫ = 0) {lam mu c : ℝ}
    (hlam : 0 ≤ lam) (hmu : 0 ≤ mu) (hc0 : 0 ≤ c) (hc1 : c ≤ 1)
    (hU : ∀ u ∈ U, lam * ‖u‖ ^ 2 ≤ B u u) (hW : ∀ w ∈ W, mu * ‖w‖ ^ 2 ≤ B w w)
    (hcoup : ∀ u ∈ U, ∀ w ∈ W, (B u w) ^ 2 ≤ c * B u u * B w w) :
    ∀ v ∈ U ⊔ W, (1 - Real.sqrt c) * min lam mu * ‖v‖ ^ 2 ≤ B v v := by
  intro v hv
  obtain ⟨u, hu, w, hw, rfl⟩ := Submodule.mem_sup.1 hv
  have hU0 : ∀ u ∈ U, 0 ≤ B u u := fun u hu => le_trans (by positivity) (hU u hu)
  have hW0 : ∀ w ∈ W, 0 ≤ B w w := fun w hw => le_trans (by positivity) (hW w hw)
  have key := add_ge_of_coupling B hB hU0 hW0 hc0 hcoup hu hw
  have hn : ‖u + w‖ ^ 2 = ‖u‖ ^ 2 + ‖w‖ ^ 2 := by
    have := norm_add_sq_eq_norm_sq_add_norm_sq_real (horth u hu w hw)
    rw [sq, sq, sq, this]
  have hs : 0 ≤ 1 - Real.sqrt c := by
    have : Real.sqrt c ≤ 1 := by rw [Real.sqrt_le_one]; exact hc1
    linarith
  have h1 : min lam mu * ‖u‖ ^ 2 ≤ B u u :=
    le_trans (mul_le_mul_of_nonneg_right (min_le_left _ _) (by positivity)) (hU u hu)
  have h2 : min lam mu * ‖w‖ ^ 2 ≤ B w w :=
    le_trans (mul_le_mul_of_nonneg_right (min_le_right _ _) (by positivity)) (hW w hw)
  rw [hn]
  calc (1 - Real.sqrt c) * min lam mu * (‖u‖ ^ 2 + ‖w‖ ^ 2)
      = (1 - Real.sqrt c) * (min lam mu * ‖u‖ ^ 2 + min lam mu * ‖w‖ ^ 2) := by ring
    _ ≤ (1 - Real.sqrt c) * (B u u + B w w) := mul_le_mul_of_nonneg_left (by linarith) hs
    _ ≤ _ := key

/-- A rung ladder: a symmetric bilinear form, windows `V k`, annuli `W k`, `V (k+1) = V k ⊔ W k`. -/
structure RungLadder (E : Type*) [NormedAddCommGroup E] [InnerProductSpace ℝ E] where
  B : LinearMap.BilinForm ℝ E
  symm : ∀ u v, B u v = B v u
  V : ℕ → Submodule ℝ E
  W : ℕ → Submodule ℝ E
  step : ∀ k, V (k + 1) = V k ⊔ W k

/-- Positivity on rung `k`. -/
def RungPositive (R : RungLadder E) (k : ℕ) : Prop := ∀ v ∈ R.V k, 0 ≤ R.B v v

/-- Coupling bound with constant `c` at step `k`, over **all** pairs. -/
def CouplingBound (R : RungLadder E) (c : ℝ) (k : ℕ) : Prop :=
  ∀ u ∈ R.V k, ∀ w ∈ R.W k, (R.B u w) ^ 2 ≤ c * R.B u u * R.B w w

/-- **Ladder theorem:** all rungs ⇔ rung 0 + every annulus + full coupling (constant 1) at every
step. PROVED. -/
theorem ladder_iff (R : RungLadder E) :
    (∀ k, RungPositive R k) ↔
      RungPositive R 0 ∧ (∀ k, ∀ w ∈ R.W k, 0 ≤ R.B w w) ∧ ∀ k, CouplingBound R 1 k := by
  constructor
  · intro h
    have hs : ∀ k, (∀ u ∈ R.V k, 0 ≤ R.B u u) ∧ (∀ w ∈ R.W k, 0 ≤ R.B w w) ∧
        ∀ u ∈ R.V k, ∀ w ∈ R.W k, (R.B u w) ^ 2 ≤ R.B u u * R.B w w := fun k =>
      (nonneg_sup_iff_coupling R.B R.symm (R.V k) (R.W k)).1 (by rw [← R.step]; exact h (k + 1))
    refine ⟨h 0, fun k => (hs k).2.1, fun k u hu w hw => ?_⟩
    simpa using (hs k).2.2 u hu w hw
  · rintro ⟨h0, hW, hc⟩ k
    induction k with
    | zero => exact h0
    | succ k ih =>
      intro v hv
      rw [R.step] at hv
      refine (nonneg_sup_iff_coupling R.B R.symm (R.V k) (R.W k)).2 ⟨ih, hW k, fun u hu w hw => ?_⟩ v hv
      simpa using hc k u hu w hw

/-- **Uniform coupling margin** (named Prop): one constant `c < 1` bounds the full coupling at
every step. For the ζ ladder: CONJECTURE false (NUMERICAL, `c_k → 1`). -/
def UniformCouplingMargin (R : RungLadder E) (c : ℝ) : Prop :=
  c < 1 ∧ ∀ k, CouplingBound R c k

/-- **A uniform coupling margin forces at most geometric decay.** With margin `lam0 > 0` on rung 0,
margin `mu` on every annulus, `V k ⊥ W k`, and `UniformCouplingMargin R c` (`0 ≤ c`), the rung-`k`
margin is at least `(1 − √c)^k · min lam0 mu`. PROVED. -/
theorem geometric_of_uniformCouplingMargin (R : RungLadder E) {lam0 mu c : ℝ} (hlam : 0 ≤ lam0)
    (hmu : 0 ≤ mu) (hc0 : 0 ≤ c) (hU : UniformCouplingMargin R c)
    (horth : ∀ k, ∀ u ∈ R.V k, ∀ w ∈ R.W k, ⟪u, w⟫ = 0)
    (h0 : ∀ v ∈ R.V 0, lam0 * ‖v‖ ^ 2 ≤ R.B v v) (hW : ∀ k, ∀ w ∈ R.W k, mu * ‖w‖ ^ 2 ≤ R.B w w) :
    ∀ k, ∀ v ∈ R.V k, (1 - Real.sqrt c) ^ k * min lam0 mu * ‖v‖ ^ 2 ≤ R.B v v := by
  have hs0 : 0 ≤ 1 - Real.sqrt c := by
    have : Real.sqrt c ≤ 1 := by rw [Real.sqrt_le_one]; exact hU.1.le
    linarith
  have hs1 : 1 - Real.sqrt c ≤ 1 := by linarith [Real.sqrt_nonneg c]
  intro k
  induction k with
  | zero =>
    intro v hv
    simpa using le_trans (mul_le_mul_of_nonneg_right (min_le_left _ _) (by positivity)) (h0 v hv)
  | succ k ih =>
    intro v hv
    rw [R.step] at hv
    set m := (1 - Real.sqrt c) ^ k * min lam0 mu
    have hm0 : 0 ≤ m := by positivity
    have hmmu : m ≤ mu := by
      have : (1 - Real.sqrt c) ^ k ≤ 1 := pow_le_one₀ hs0 hs1
      calc m ≤ 1 * min lam0 mu := mul_le_mul_of_nonneg_right this (by positivity)
        _ ≤ mu := by simp
    have := margin_sup_of_coupling R.B R.symm (horth k) hm0 hmu hc0 hU.1.le ih (hW k) (hU.2 k) v hv
    have hmin : min m mu = m := min_eq_left hmmu
    rw [hmin] at this
    calc (1 - Real.sqrt c) ^ (k + 1) * min lam0 mu * ‖v‖ ^ 2
        = (1 - Real.sqrt c) * m * ‖v‖ ^ 2 := by rw [pow_succ]; ring
      _ ≤ _ := this

/-- **Soft-mode slice bound** (named Prop): the coupling of the rung-`k` soft mode `φ` to every
annulus vector, with constant `c`. For ζ, `c ≈ 0.5–0.8` at every step measured, with no growth
(NUMERICAL; CONJECTURE `c ≤ 4/5` uniformly). -/
def SoftModeSliceBound (R : RungLadder E) (c : ℝ) (k : ℕ) (φ : E) : Prop :=
  ∀ w ∈ R.W k, (R.B φ w) ^ 2 ≤ c * R.B φ φ * R.B w w

/-- **What the slice bound gives, and only that:** positivity, with an additive margin, on the slice
`span{φ} ⊔ W k`, which is strictly smaller than the step `V k ⊔ W k`. PROVED. -/
theorem slice_margin (R : RungLadder E) {c : ℝ} {k : ℕ} {φ : E} (hc : 0 ≤ c)
    (hφ : 0 ≤ R.B φ φ) (hW : ∀ w ∈ R.W k, 0 ≤ R.B w w) (hs : SoftModeSliceBound R c k φ)
    (a : ℝ) {w : E} (hw : w ∈ R.W k) :
    (1 - Real.sqrt c) * (R.B (a • φ) (a • φ) + R.B w w) ≤ R.B (a • φ + w) (a • φ + w) := by
  have hU : ∀ u ∈ Submodule.span ℝ {φ}, 0 ≤ R.B u u := by
    intro u hu
    obtain ⟨t, rfl⟩ := Submodule.mem_span_singleton.1 hu
    simp only [map_smul, LinearMap.smul_apply, smul_eq_mul]
    nlinarith [sq_nonneg t]
  have hcoup : ∀ u ∈ Submodule.span ℝ {φ}, ∀ w ∈ R.W k, (R.B u w) ^ 2 ≤ c * R.B u u * R.B w w := by
    intro u hu w hw
    obtain ⟨t, rfl⟩ := Submodule.mem_span_singleton.1 hu
    simp only [map_smul, LinearMap.smul_apply, smul_eq_mul]
    have := hs w hw
    nlinarith [sq_nonneg t, mul_le_mul_of_nonneg_left this (sq_nonneg t)]
  exact add_ge_of_coupling R.B R.symm hU hW hc hcoup
    (Submodule.mem_span_singleton.2 ⟨a, rfl⟩) hw

end Ladder

section ToeplitzRankTwo

variable {E : Type*} [AddCommGroup E] [Module ℝ E]

/-- **Toeplitz + rank-two structure:** `Q = T + 2 ℓ_c² − 2 ℓ_s²` (for ζ: `T` = archimedean − prime
truncated convolution, `ℓ_c = ⟨cosh(x/2), ·⟩`, `ℓ_s = ⟨sinh(x/2), ·⟩`). If `Q ≥ 0`, then on the
hyperplane `ℓ_c = 0`, `T ≥ 2 ℓ_s² ≥ 0`: `T` has at most one negative square, and `Q ≤ T` there
(so `λ_min(Q)` is bounded by the Toeplitz part's spectrum on that hyperplane). PROVED. -/
theorem toeplitz_on_ker_cosh (T : E → ℝ) (lc ls : E →ₗ[ℝ] ℝ)
    (hQ : ∀ v, 0 ≤ T v + 2 * (lc v) ^ 2 - 2 * (ls v) ^ 2) {v : E} (hv : lc v = 0) :
    2 * (ls v) ^ 2 ≤ T v ∧ 0 ≤ T v ∧ T v + 2 * (lc v) ^ 2 - 2 * (ls v) ^ 2 ≤ T v := by
  have h := hQ v
  rw [hv] at h ⊢
  refine ⟨by nlinarith, by nlinarith [sq_nonneg (ls v)], by nlinarith [sq_nonneg (ls v)]⟩

end ToeplitzRankTwo

section QuantitativeMargin

open MeasureTheory GppArithmeticReflection
open scoped ContDiff

/-- **Quantitative rung margin** (named Prop, CONJECTURE): on every window `L ≥ log 2` the ζ Weil
form is `≥ exp(−C·L·e^L)·‖f‖²`. This is strictly stronger than ARP. It is what an induction would
actually have to deliver, and a uniform coupling margin cannot deliver it
(`geometric_of_uniformCouplingMargin`). NUMERICAL (high precision, Rayleigh–Ritz upper bounds,
`num/rung3/`): `−ln λ(L)/(L e^L) = 5.07, 5.03, 4.87, 4.59` at `L = log 3, log 4, log 5, log 7`, so
the data so far are consistent with `C = 6`. -/
def QuantitativeRungMargin (C : ℝ) : Prop :=
  ∀ L : ℝ, Real.log 2 ≤ L → ∀ f : ℝ → ℝ, ContDiff ℝ ∞ f → HasCompactSupport f →
    (∀ y : ℝ, L / 2 < |y| → f y = 0) →
      Real.exp (-(C * L * Real.exp L)) * ∫ y, f y ^ 2 ≤ weilForm (autocorr f)

/-- **The quantitative margin implies ARP.** PROVED (conditional on the named Prop). -/
theorem arp_of_quantitativeRungMargin {C : ℝ} (h : QuantitativeRungMargin C) :
    ArithmeticReflectionPositivity := by
  have hr : ∀ L : ℝ, ReflectionPositiveOnRung L := by
    intro L f hf hc hs
    have := h (max L (Real.log 2)) (le_max_right _ _) f hf hc fun y hy =>
      hs y (lt_of_le_of_lt (by have := le_max_left L (Real.log 2); linarith) hy)
    have hint : 0 ≤ ∫ y, f y ^ 2 := integral_nonneg fun y => sq_nonneg _
    exact le_trans (mul_nonneg (Real.exp_pos _).le hint) this
  intro f hf hc
  obtain ⟨r, hr'⟩ := hc.isCompact.isBounded.subset_closedBall 0
  refine hr (2 * r) f hf hc fun y hy => ?_
  apply image_eq_zero_of_notMem_tsupport
  intro hy'
  have := hr' hy'
  rw [Metric.mem_closedBall, dist_zero_right, Real.norm_eq_abs] at this
  linarith

end QuantitativeMargin

end GppWeilRungLadder

#print axioms GppWeilRungLadder.sq_le_mul_of_forall_nonneg
#print axioms GppWeilRungLadder.add_nonneg_of_sq_le
#print axioms GppWeilRungLadder.add_ge_of_sq_le
#print axioms GppWeilRungLadder.quad_expand
#print axioms GppWeilRungLadder.add_expand
#print axioms GppWeilRungLadder.nonneg_sup_iff_coupling
#print axioms GppWeilRungLadder.add_ge_of_coupling
#print axioms GppWeilRungLadder.margin_sup_of_coupling
#print axioms GppWeilRungLadder.ladder_iff
#print axioms GppWeilRungLadder.geometric_of_uniformCouplingMargin
#print axioms GppWeilRungLadder.slice_margin
#print axioms GppWeilRungLadder.toeplitz_on_ker_cosh
#print axioms GppWeilRungLadder.arp_of_quantitativeRungMargin
