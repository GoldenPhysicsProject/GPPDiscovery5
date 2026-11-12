/-
Copyright (c) 2026 Daniel Toupin. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Daniel Toupin (with AI assistance: Heart of the Cosmos / Grok)
-/
import GppVerify.RiemannHypothesis.ArithmeticReflectionPositivity
import GppVerify.RiemannHypothesis.ArchimedeanFloor
import GppVerify.RiemannHypothesis.BTSigmaLocalEuler
import GppVerify.RiemannHypothesis.CasimirCriticalLine
import GppVerify.RiemannHypothesis.GravitonDeltaIdentification
import GppVerify.RiemannHypothesis.XiReflection
import Mathlib.NumberTheory.SumPrimeReciprocals

/-!
# Building the CFT from arithmetic alone: three routes, each killed or advanced

Companion to `SEqDelta1d.lean`, whose remaining GAP is `ShadowIsAdjointOnSpectrum`
(`∀ ρ ∈ nontrivialZeros, 1 - ρ = conj ρ`), equivalent to ARP / `UnitaryRealization`. The task is to
build the Hilbert space and the states `O_f` (with `gram f = ‖O_f‖²`) with no zero as input.

## Route 1: adelic/Tate local-factor product of positive forms
In the multiplier picture the `p`-local piece of the Weil form multiplies `|f̂(t)|²` by
`-σ_p(t)`, where `σ_p(t) = 2 Σ_{k≥1} log p · p^{-k/2} cos(k t log p) = 2 Re(-L_p'/L_p)(1/2 + it)`.
* `primeSymbol_zero`: `σ_p(0) = 2κ_p`, with `κ_p = log p/(√p - 1)` the P-LOG coupling
  (`kappaP_eq_tateCoupling`: `κ_p = -(L_p'/L_p)(1/2)`, in-tree #579 lane).
* `primeSymbol_le`, `primeSymbol_ge`: `|σ_p(t)| ≤ 2κ_p`; the maximum `2κ_p` is attained at `t = 0`.
* `primeSymbol_at_pi`: `σ_p(π/log p) = -2 log p/(√p + 1) < 0`.
* `localPrimeMultiplier_indefinite`: the local `p`-adic form `-σ_p` takes both signs. **KILLED:**
  no place carries a positive local form, so ARP cannot be a restricted product of local positive
  forms. What survives is the local floor `-σ_p ≥ -2κ_p`.
* `kappaP_not_summable`: `Σ_p κ_p = ∞`. The local floors cannot be summed: the restricted
  product of the floors diverges. Positivity, if true, is a global cancellation across primes.

## Route 2: semi-local truncation over a finite set `S` of primes
`weilS S g = pole + A∞ − Σ_{p∈S} prime_p(g)`, with `prime_p(g) = Σ_k log p · p^{-k/2}(g(k log p) +
g(−k log p))`.
* `primeLocal_le`, `primeLocal_ge` (Cauchy–Schwarz bound `|g| ≤ g(0)` as hypothesis): `|prime_p(g)|
  ≤ 2κ_p g(0)`. `primeLocal_const`: the bound is the value on the infinitely broad state.
* `weilS_floor`: for normalised `g` (in-tree `archimedean_floor` hypotheses),
  `weilS S g ≥ pole(g) + C∞ − 2 Σ_{p∈S} κ_p`, where `C∞ ≈ −5.372`. PROVED.
* `semiLocal_floor_unbounded`: this floor is not uniform in `S` (it tends to −∞). PROVED.
* NUMERICAL (`num/wS_indefinite.py`): for fixed finite `S`, broad tests with zero pole term give
  `weilS S g / g(0) → C∞ − 2 Σ_{p∈S} κ_p < 0`. The floor is sharp: **each fixed finite-`S` form is
  indefinite** on large support. So "positivity for each finite `S`" can only hold on the `S`-rung
  window, where `weilS` is the full Weil form (in-tree `weilForm_eq_truncation`). On the same
  tests the full form is `≥ 0` but only `~10⁻⁶ g(0)` (razor-thin global cancellation).
* `arp_iff_all_rungs`: ARP ⇔ positivity on every rung window. PROVED. So the `S → all` step is
  exactly `∀ L, ReflectionPositiveOnRung L`, i.e. ARP itself. Rung 0 is the Connes–Consani
  archimedean theorem (cited, in-tree as a named hypothesis); rungs ≥ 1 are OPEN.

## Route 3: principal-series model
* `bilinear_shadow_pairing_invariant`: the bilinear pairing `χ_s ⊗ χ_{1-s}` is invariant for every
  `s`. This is the pairing the ζ correlator uses, and it carries no positivity.
* `hermitian_invariant_iff`: a nonzero invariant Hermitian form on `χ_s` exists iff `Re s = 1/2`.
  Turning the bilinear shadow pairing into a Hermitian one is the GAP again.
* `casimir_nonneg_real_iff`: a self-adjoint Casimir eigenvalue `s(1−s) ≥ 0` forces `Re s = 1/2`
  (principal) or real `s ∈ [0,1]` (complementary series). `complementary_casimir`: the
  complementary series is really there.
* `gap_of_nonnegCasimir_noReal`, `rh_strip_of_nonnegCasimir_noReal`: a **nonnegative** symmetric
  Casimir realization (the automorphic `Δ ≥ 0` shape, weaker than Hilbert–Pólya `≥ 1/4`), plus
  `NoRealStripZeros`, gives the GAP. `nonnegCasimir_of_principalCasimir` shows the new hypothesis is
  weaker. `NoRealStripZeros` is STANDARD (ζ(σ) < 0 on (0,1)) but not in Mathlib; it is named here.
* Paper (STANDARD, Selberg / Maass–Selberg): in the automorphic model `PGL₂(ℤ)\ℍ` the Eisenstein
  scattering matrix `φ(s) = ξ(2s−1)/ξ(2s)` is unimodular on `Re s = 1/2` regardless of zeros, and the
  ζ-zeros appear as its poles at `s = ρ/2` (resonances, `Re s = 1/4` under RH), not as unitary
  dimensions. So in that model the zeros are not matrix-coefficient dimensions of a unitary
  representation. **KILLED** as a zero-location route. The arithmetic input a working model would
  need is exactly a self-adjoint `A ≥ 0` with spectrum `{ρ(1−ρ)}` (NonnegCasimirRealization), which
  is Hilbert–Pólya again.

Labels: every theorem PROVED (standard axioms). Named hypotheses `NonnegCasimirRealization`,
`NoRealStripZeros` (STANDARD, not in Mathlib), and rungs ≥ 1 are OPEN. ARP OPEN. RH OPEN.

Daniel Toupin / Golden Physics Project … This seat: Heart of the Cosmos / Grok (GPPDiscovery5).
Audit: Turbo.
-/

open Complex MeasureTheory
open scoped InnerProductSpace ComplexConjugate ContDiff
open GppWeilCriterion GppArithmeticReflection GppCasimirCriticalLine

noncomputable section

namespace GppArithmeticCFTRoutes

/-! ## Route 1: local multipliers -/

/-- `q_p = p^{-1/2}`. -/
def qp (p : ℕ) : ℝ := (Real.sqrt p)⁻¹

/-- The P-LOG coupling `κ_p = log p / (√p − 1)`. -/
def kappaP (p : ℕ) : ℝ := Real.log p / (Real.sqrt p - 1)

theorem one_lt_sqrt' {p : ℕ} (hp : 1 < p) : 1 < Real.sqrt p := by
  have : (1 : ℝ) < p := by exact_mod_cast hp
  rw [show (1 : ℝ) = Real.sqrt 1 by simp]
  exact Real.sqrt_lt_sqrt zero_le_one this

theorem qp_pos {p : ℕ} (hp : 1 < p) : 0 < qp p := by
  unfold qp; have := one_lt_sqrt' hp; positivity

theorem qp_lt_one {p : ℕ} (hp : 1 < p) : qp p < 1 := by
  unfold qp; exact inv_lt_one_of_one_lt₀ (one_lt_sqrt' hp)

theorem log_pos' {p : ℕ} (hp : 1 < p) : 0 < Real.log p :=
  Real.log_pos (by exact_mod_cast hp)

theorem kappaP_pos {p : ℕ} (hp : 1 < p) : 0 < kappaP p := by
  unfold kappaP; have := one_lt_sqrt' hp; have := log_pos' hp
  exact div_pos this (by linarith)

/-- `κ_p` is the P-LOG/Tate coupling `−(L_p'/L_p)(1/2)` (in-tree). PROVED. -/
theorem kappaP_eq_tateCoupling {p : ℕ} (hp : 1 < p) :
    kappaP p = GppBTSigmaLocalEuler.tateCoupling p (1 / 2) :=
  GppBTSigmaLocalEuler.kappa_eq_tateCoupling hp

/-- `log p · Σ_{k≥0} q^{k+1} = κ_p`. -/
theorem geom_qp {p : ℕ} (hp : 1 < p) :
    HasSum (fun k : ℕ => qp p ^ (k + 1)) (qp p / (1 - qp p)) := by
  have h := hasSum_geometric_of_lt_one (qp_pos hp).le (qp_lt_one hp)
  have := h.mul_left (qp p)
  simpa [pow_succ, mul_comm, div_eq_mul_inv] using this

theorem qp_div {p : ℕ} (hp : 1 < p) : qp p / (1 - qp p) = (Real.sqrt p - 1)⁻¹ := by
  have hs := one_lt_sqrt' hp
  unfold qp
  field_simp

theorem geom_neg_qp {p : ℕ} (hp : 1 < p) :
    HasSum (fun k : ℕ => (-qp p) ^ (k + 1)) (-qp p / (1 + qp p)) := by
  have hq := qp_pos hp; have hq1 := qp_lt_one hp
  have h := hasSum_geometric_of_abs_lt_one (r := -qp p) (by rw [abs_neg, abs_of_pos hq]; exact hq1)
  have := h.mul_left (-qp p)
  have e : (1 - -qp p) = 1 + qp p := by ring
  simpa [pow_succ, mul_comm, div_eq_mul_inv, e] using this

/-- Term of the local prime symbol. -/
def symTerm (p : ℕ) (t : ℝ) (k : ℕ) : ℝ :=
  Real.log p * qp p ^ (k + 1) * Real.cos ((k + 1) * t * Real.log p)

/-- The local prime symbol `σ_p(t) = 2 Σ_{k≥1} log p · p^{-k/2} cos(k t log p)`. -/
def primeSymbol (p : ℕ) (t : ℝ) : ℝ := 2 * ∑' k, symTerm p t k

theorem summable_symTerm {p : ℕ} (hp : 1 < p) (t : ℝ) : Summable (symTerm p t) := by
  refine Summable.of_norm_bounded ((geom_qp hp).summable.mul_left (Real.log p)) fun k => ?_
  have hl := log_pos' hp; have hq := qp_pos hp
  rw [Real.norm_eq_abs, symTerm, abs_mul, abs_mul, abs_of_pos hl, abs_of_pos (pow_pos hq _)]
  calc Real.log p * qp p ^ (k + 1) * |Real.cos _| ≤ Real.log p * qp p ^ (k + 1) * 1 := by
        gcongr; exact Real.abs_cos_le_one _
    _ = _ := by ring

theorem hasSum_log_geom {p : ℕ} (hp : 1 < p) :
    HasSum (fun k : ℕ => Real.log p * qp p ^ (k + 1)) (kappaP p) := by
  have := (geom_qp hp).mul_left (Real.log p)
  rw [qp_div hp] at this
  simpa [kappaP, div_eq_mul_inv] using this

/-- **`σ_p(0) = 2κ_p`.** PROVED. -/
theorem primeSymbol_zero {p : ℕ} (hp : 1 < p) : primeSymbol p 0 = 2 * kappaP p := by
  unfold primeSymbol
  congr 1
  have : symTerm p 0 = fun k : ℕ => Real.log p * qp p ^ (k + 1) := by
    funext k; simp [symTerm]
  rw [this, (hasSum_log_geom hp).tsum_eq]

/-- **Upper floor `σ_p ≤ 2κ_p`** (sharp at `t = 0`). PROVED. -/
theorem primeSymbol_le {p : ℕ} (hp : 1 < p) (t : ℝ) : primeSymbol p t ≤ 2 * kappaP p := by
  unfold primeSymbol
  have hl := log_pos' hp; have hq := qp_pos hp
  gcongr
  refine hasSum_le (fun k => ?_) (summable_symTerm hp t).hasSum (hasSum_log_geom hp)
  unfold symTerm
  have : 0 ≤ Real.log p * qp p ^ (k + 1) := by positivity
  nlinarith [Real.cos_le_one ((k + 1) * t * Real.log p)]

/-- **Lower floor `−2κ_p ≤ σ_p`.** PROVED. -/
theorem primeSymbol_ge {p : ℕ} (hp : 1 < p) (t : ℝ) : -(2 * kappaP p) ≤ primeSymbol p t := by
  unfold primeSymbol
  have hl := log_pos' hp; have hq := qp_pos hp
  have h : -kappaP p ≤ ∑' k, symTerm p t k := by
    refine hasSum_le (fun k => ?_) (hasSum_log_geom hp).neg (summable_symTerm hp t).hasSum
    unfold symTerm
    have : 0 ≤ Real.log p * qp p ^ (k + 1) := by positivity
    nlinarith [Real.neg_one_le_cos ((k + 1) * t * Real.log p)]
  linarith

/-- **`σ_p(π / log p) = −2 log p / (√p + 1)`.** PROVED. -/
theorem primeSymbol_at_pi {p : ℕ} (hp : 1 < p) :
    primeSymbol p (Real.pi / Real.log p) = -(2 * Real.log p / (Real.sqrt p + 1)) := by
  have hl := log_pos' hp
  have hs := one_lt_sqrt' hp
  unfold primeSymbol
  have hterm : symTerm p (Real.pi / Real.log p) = fun k : ℕ => Real.log p * (-qp p) ^ (k + 1) := by
    funext k
    unfold symTerm
    have e : ((k : ℝ) + 1) * (Real.pi / Real.log p) * Real.log p = ((k + 1 : ℕ) : ℝ) * Real.pi := by
      field_simp; push_cast; ring
    rw [e, Real.cos_nat_mul_pi, neg_pow, mul_comm ((-1 : ℝ) ^ (k + 1))]; ring
  rw [hterm, ((geom_neg_qp hp).mul_left (Real.log p)).tsum_eq]
  unfold qp
  field_simp

/-- **Route 1 kill: the local `p`-adic form is indefinite.** The multiplier `−σ_p` is `−2κ_p < 0`
at `t = 0` and `> 0` at `t = π / log p`. PROVED. -/
theorem localPrimeMultiplier_indefinite {p : ℕ} (hp : 1 < p) :
    -primeSymbol p 0 < 0 ∧ 0 < -primeSymbol p (Real.pi / Real.log p) := by
  have hk := kappaP_pos hp; have hl := log_pos' hp; have hs := one_lt_sqrt' hp
  refine ⟨by rw [primeSymbol_zero hp]; linarith, ?_⟩
  rw [primeSymbol_at_pi hp, neg_neg]; positivity

/-- **The local floors cannot be summed:** `Σ_p κ_p = ∞`. PROVED (comparison with `Σ 1/p`). -/
theorem kappaP_not_summable : ¬ Summable (fun p : Nat.Primes => kappaP p) := by
  intro h
  apply Nat.Primes.not_summable_one_div
  refine Summable.of_nonneg_of_le (fun p => by positivity) (fun p => ?_) h
  have hp : 1 < (p : ℕ) := p.2.one_lt
  have hl := log_pos' hp; have hs := one_lt_sqrt' hp
  have hp2 : (2 : ℝ) ≤ (p : ℕ) := by exact_mod_cast p.2.two_le
  have hl2 : (1 / 2 : ℝ) ≤ Real.log (p : ℕ) := by
    have h2 : Real.log 2 ≤ Real.log (p : ℕ) := Real.log_le_log (by norm_num) hp2
    have := Real.log_two_gt_d9; linarith
  have hsq : Real.sqrt ((p : ℕ) : ℝ) ^ 2 = (p : ℕ) := Real.sq_sqrt (by positivity)
  unfold kappaP
  rw [div_le_div_iff₀ (by positivity) (by linarith)]
  nlinarith [sq_nonneg (Real.sqrt ((p : ℕ) : ℝ) - 1)]

/-! ## Route 2: semi-local truncation -/

/-- The `p`-local prime functional `Σ_{k≥1} log p · p^{-k/2} (g(k log p) + g(−k log p))`. -/
def primeLocal (p : ℕ) (g : ℝ → ℝ) : ℝ :=
  ∑' k : ℕ, Real.log p * qp p ^ (k + 1) *
    (g ((k + 1) * Real.log p) + g (-((k + 1) * Real.log p)))

theorem summable_primeLocal {p : ℕ} (hp : 1 < p) {g : ℝ → ℝ} (hb : ∀ x, |g x| ≤ g 0) :
    Summable fun k : ℕ => Real.log p * qp p ^ (k + 1) *
      (g ((k + 1) * Real.log p) + g (-((k + 1) * Real.log p))) := by
  refine Summable.of_norm_bounded ((hasSum_log_geom hp).summable.mul_right (2 * g 0)) fun k => ?_
  have hl := log_pos' hp; have hq := qp_pos hp
  rw [Real.norm_eq_abs, abs_mul, abs_of_pos (by positivity : 0 < Real.log p * qp p ^ (k + 1))]
  have := abs_add_le (g ((k + 1) * Real.log p)) (g (-((k + 1) * Real.log p)))
  have h1 := hb ((k + 1) * Real.log p); have h2 := hb (-((k + 1) * Real.log p))
  have : 0 ≤ Real.log p * qp p ^ (k + 1) := by positivity
  nlinarith

/-- **`prime_p(g) ≤ 2κ_p g(0)`** under the Cauchy–Schwarz bound `|g| ≤ g(0)`. PROVED. -/
theorem primeLocal_le {p : ℕ} (hp : 1 < p) {g : ℝ → ℝ} (hb : ∀ x, |g x| ≤ g 0) :
    primeLocal p g ≤ 2 * kappaP p * g 0 := by
  have hs := (hasSum_log_geom hp).mul_right (2 * g 0)
  rw [show 2 * kappaP p * g 0 = kappaP p * (2 * g 0) by ring]
  refine hasSum_le (fun k => ?_) (summable_primeLocal hp hb).hasSum hs
  have hl := log_pos' hp; have hq := qp_pos hp
  have h1 := (abs_le.1 (hb ((k + 1) * Real.log p))).2
  have h2 := (abs_le.1 (hb (-((k + 1) * Real.log p)))).2
  have : 0 ≤ Real.log p * qp p ^ (k + 1) := by positivity
  nlinarith

/-- **`−2κ_p g(0) ≤ prime_p(g)`.** PROVED. -/
theorem primeLocal_ge {p : ℕ} (hp : 1 < p) {g : ℝ → ℝ} (hb : ∀ x, |g x| ≤ g 0) :
    -(2 * kappaP p * g 0) ≤ primeLocal p g := by
  have hs := ((hasSum_log_geom hp).mul_right (2 * g 0)).neg
  rw [show -(2 * kappaP p * g 0) = -(kappaP p * (2 * g 0)) by ring]
  refine hasSum_le (fun k => ?_) hs (summable_primeLocal hp hb).hasSum
  have hl := log_pos' hp; have hq := qp_pos hp
  have h1 := (abs_le.1 (hb ((k + 1) * Real.log p))).1
  have h2 := (abs_le.1 (hb (-((k + 1) * Real.log p)))).1
  have : 0 ≤ Real.log p * qp p ^ (k + 1) := by positivity
  nlinarith

/-- The envelope is attained by the infinitely broad state `g ≡ c`. PROVED. -/
theorem primeLocal_const {p : ℕ} (hp : 1 < p) (c : ℝ) :
    primeLocal p (fun _ => c) = 2 * kappaP p * c := by
  unfold primeLocal
  have := (hasSum_log_geom hp).mul_right (2 * c)
  rw [show 2 * kappaP p * c = kappaP p * (2 * c) by ring, ← this.tsum_eq]
  congr 1; funext k; ring

/-- The semi-local truncated Weil form over a finite set `S` of primes. -/
def weilS (S : Finset Nat.Primes) (g : ℝ → ℝ) : ℝ :=
  poleTerm g + GppArchimedeanFloor.archimedeanTerm g - ∑ p ∈ S, primeLocal p g

/-- **Semi-local floor.** For a normalised `g` (`g 0 = 1`, `|g| ≤ 1`, archimedean integrand
integrable): `weilS S g ≥ pole(g) + C∞ − 2 Σ_{p∈S} κ_p`. PROVED. -/
theorem weilS_floor (S : Finset Nat.Primes) {g : ℝ → ℝ} (h0 : g 0 = 1)
    (hb : ∀ x, |g x| ≤ g 0)
    (hint : IntegrableOn (GppArchimedeanFloor.archimedeanIntegrand g) (Set.Ioi 0)) :
    poleTerm g + GppArchimedeanFloor.archimedeanFloor - 2 * ∑ p ∈ S, kappaP p ≤ weilS S g := by
  have harch := GppArchimedeanFloor.archimedean_floor g h0
    (fun u _ => by have := (abs_le.1 (hb u)).2; linarith) hint
  have hpr : ∑ p ∈ S, primeLocal p g ≤ 2 * ∑ p ∈ S, kappaP p := by
    rw [Finset.mul_sum]
    refine Finset.sum_le_sum fun p _ => ?_
    have := primeLocal_le p.2.one_lt hb
    rw [h0, mul_one] at this; exact this
  unfold weilS; linarith

/-- **The floor is not uniform in `S`:** `Σ_{p∈S} κ_p` is unbounded over finite `S`. PROVED. -/
theorem semiLocal_floor_unbounded (M : ℝ) : ∃ S : Finset Nat.Primes, M < ∑ p ∈ S, kappaP p := by
  by_contra h
  simp only [not_exists, not_lt] at h
  apply kappaP_not_summable
  exact summable_of_sum_le (fun p => (kappaP_pos p.2.one_lt).le) h

/-- **ARP ⇔ positivity on every rung window.** The `S → all` limit is exactly ARP. PROVED. -/
theorem arp_iff_all_rungs :
    ArithmeticReflectionPositivity ↔ ∀ L : ℝ, ReflectionPositiveOnRung L := by
  refine ⟨fun h L => rung_of_reflectionPositivity h L, fun h f hf hc => ?_⟩
  obtain ⟨r, hr⟩ := hc.isCompact.isBounded.subset_closedBall 0
  refine h (2 * r) f hf hc fun y hy => ?_
  apply image_eq_zero_of_notMem_tsupport
  intro hy'
  have := hr hy'
  rw [Metric.mem_closedBall, dist_zero_right, Real.norm_eq_abs] at this
  linarith

/-! ## Route 3: principal-series model -/

/-- The log-scale character `χ_s(a) = e^{(s − 1/2) a}`. -/
def chi (s : ℂ) (a : ℝ) : ℂ := Complex.exp ((s - 1 / 2) * a)

/-- **The bilinear shadow pairing is invariant for every `s`.** PROVED. -/
theorem bilinear_shadow_pairing_invariant (s : ℂ) (a : ℝ) : chi s a * chi (1 - s) a = 1 := by
  unfold chi; rw [← Complex.exp_add]; ring_nf; exact Complex.exp_zero

/-- **A nonzero invariant Hermitian form on `χ_s` exists iff `Re s = 1/2`.** PROVED. -/
theorem hermitian_invariant_iff (s : ℂ) {b : ℂ} (hb : b ≠ 0) :
    (∀ a : ℝ, conj (chi s a) * chi s a * b = b) ↔ s.re = 1 / 2 := by
  have key : ∀ a : ℝ, conj (chi s a) * chi s a = ((Real.exp ((2 * s.re - 1) * a) : ℝ) : ℂ) := by
    intro a
    unfold chi
    rw [← Complex.exp_conj, ← Complex.exp_add, Complex.ofReal_exp]
    congr 1
    apply Complex.ext <;> simp; ring
  constructor
  · intro h
    have h1 := h 1
    rw [key 1] at h1
    have h2 : ((Real.exp ((2 * s.re - 1) * 1) : ℝ) : ℂ) = 1 := by
      have := mul_right_cancel₀ hb (h1.trans (one_mul b).symm)
      exact this
    have h3 : Real.exp ((2 * s.re - 1) * 1) = 1 := by exact_mod_cast h2
    rw [Real.exp_eq_one_iff] at h3
    linarith
  · intro hs a
    rw [key a, hs]; norm_num

/-- **Nonnegative real Casimir ⇔ principal or complementary.** `s(1−s)` is real and `≥ 0` iff
`Re s = 1/2`, or `s` is real in `[0, 1]`. PROVED. -/
theorem casimir_nonneg_real_iff (ρ : ℂ) :
    ((casimir ρ).im = 0 ∧ 0 ≤ (casimir ρ).re) ↔
      ρ.re = 1 / 2 ∨ (ρ.im = 0 ∧ 0 ≤ ρ.re ∧ ρ.re ≤ 1) := by
  rw [casimir_real_iff, casimir_re]
  constructor
  · rintro ⟨h | h, hre⟩
    · exact Or.inl h
    · rw [h] at hre
      refine Or.inr ⟨h, ?_, ?_⟩ <;> nlinarith
  · rintro (h | ⟨h, h0, h1⟩)
    · refine ⟨Or.inl h, ?_⟩; rw [h]; nlinarith [sq_nonneg ρ.im]
    · refine ⟨Or.inr h, ?_⟩; rw [h]; nlinarith

/-- **The complementary series is really there:** for real `s ∈ (0,1)`, `s(1−s) ∈ (0, 1/4]` is
real. So a nonnegative Casimir realization does not by itself exclude real zeros. PROVED. -/
theorem complementary_casimir (σ : ℝ) (h0 : 0 < σ) (h1 : σ < 1) :
    (casimir (σ : ℂ)).im = 0 ∧ 0 < (casimir (σ : ℂ)).re ∧ (casimir (σ : ℂ)).re ≤ 1 / 4 := by
  rw [casimir_im, casimir_re]
  simp only [Complex.ofReal_re, Complex.ofReal_im]
  refine ⟨by ring, by nlinarith, by nlinarith [sq_nonneg (σ - 1 / 2)]⟩

/-- The GAP of `SEqDelta1d.lean`, restated (definitionally the same Prop). -/
def ShadowIsAdjointOnSpectrum : Prop := ∀ ρ ∈ nontrivialZeros, 1 - ρ = conj ρ

/-- **Named (STANDARD, not in Mathlib):** no nontrivial zero is real. -/
def NoRealStripZeros : Prop := ∀ ρ ∈ nontrivialZeros, ρ.im ≠ 0

/-- **Automorphic shape of the principal-series side:** one symmetric operator `A ≥ 0` (not
`≥ 1/4`) with every `casimir ρ` an eigenvalue. Named hypothesis: OPEN. -/
def NonnegCasimirRealization : Prop :=
  ∃ (E : Type) (_ : NormedAddCommGroup E) (_ : InnerProductSpace ℂ E) (A : E →ₗ[ℂ] E),
    (∀ x y, ⟪A x, y⟫_ℂ = ⟪x, A y⟫_ℂ) ∧ (∀ x, 0 ≤ (⟪A x, x⟫_ℂ).re) ∧
    ∀ ρ ∈ nontrivialZeros, ∃ v : E, v ≠ 0 ∧ A v = casimir ρ • v

/-- Hilbert–Pólya shape `A ≥ 1/4` (as in `SEqDelta1d.lean`). -/
def PrincipalCasimirRealization : Prop :=
  ∃ (E : Type) (_ : NormedAddCommGroup E) (_ : InnerProductSpace ℂ E) (A : E →ₗ[ℂ] E),
    (∀ x y, ⟪A x, y⟫_ℂ = ⟪x, A y⟫_ℂ) ∧ (∀ x, (1 / 4 : ℝ) * ‖x‖ ^ 2 ≤ (⟪A x, x⟫_ℂ).re) ∧
    ∀ ρ ∈ nontrivialZeros, ∃ v : E, v ≠ 0 ∧ A v = casimir ρ • v

/-- The automorphic hypothesis is weaker than Hilbert–Pólya. PROVED. -/
theorem nonnegCasimir_of_principalCasimir (h : PrincipalCasimirRealization) :
    NonnegCasimirRealization := by
  obtain ⟨E, i1, i2, A, hs, hl, he⟩ := h
  exact ⟨E, i1, i2, A, hs, fun x => le_trans (by positivity) (hl x), he⟩

/-- Eigenvalues of a symmetric `A ≥ 0` are real and `≥ 0`. PROVED. -/
theorem eigen_real_nonneg {E : Type} [NormedAddCommGroup E] [InnerProductSpace ℂ E]
    (A : E →ₗ[ℂ] E) (hsym : ∀ x y, ⟪A x, y⟫_ℂ = ⟪x, A y⟫_ℂ) (hpos : ∀ x, 0 ≤ (⟪A x, x⟫_ℂ).re)
    {c : ℂ} {v : E} (hv0 : v ≠ 0) (hAv : A v = c • v) : c.im = 0 ∧ 0 ≤ c.re := by
  have hn : (0 : ℝ) < ‖v‖ ^ 2 := by positivity
  have hin : ⟪A v, v⟫_ℂ = conj c * ((‖v‖ : ℂ) ^ 2) := by
    rw [hAv, inner_smul_left, inner_self_eq_norm_sq_to_K]; rfl
  have hin2 : ⟪v, A v⟫_ℂ = c * ((‖v‖ : ℂ) ^ 2) := by
    rw [hAv, inner_smul_right, inner_self_eq_norm_sq_to_K]; rfl
  have hcc : conj c = c := by
    have := hsym v v
    rw [hin, hin2] at this
    have hne : ((‖v‖ : ℂ) ^ 2) ≠ 0 := by exact_mod_cast hn.ne'
    exact mul_right_cancel₀ hne this
  have him : c.im = 0 := by
    have := congrArg Complex.im hcc
    simp at this; linarith
  refine ⟨him, ?_⟩
  have hl := hpos v
  rw [hin, hcc] at hl
  have hre : (c * ((‖v‖ : ℂ) ^ 2)).re = c.re * ‖v‖ ^ 2 := by
    rw [← Complex.ofReal_pow, Complex.mul_re, Complex.ofReal_re, Complex.ofReal_im]; ring
  rw [hre] at hl
  by_contra hneg; rw [not_le] at hneg
  nlinarith

/-- **Nonnegative Casimir realization + no real zeros ⇒ the GAP.** PROVED (conditional). -/
theorem gap_of_nonnegCasimir_noReal (h : NonnegCasimirRealization) (hr : NoRealStripZeros) :
    ShadowIsAdjointOnSpectrum := by
  obtain ⟨E, _, _, A, hsym, hpos, heig⟩ := h
  intro ρ hρ
  obtain ⟨v, hv0, hAv⟩ := heig ρ hρ
  have hc := eigen_real_nonneg A hsym hpos hv0 hAv
  rcases (casimir_nonneg_real_iff ρ).1 hc with hre | ⟨him, -, -⟩
  · apply Complex.ext <;> simp [hre]; linarith
  · exact absurd him (hr ρ hρ)

/-- **… ⇒ RH (strip form) ⇒ Mathlib `RiemannHypothesis`.** PROVED (conditional). -/
theorem rh_of_nonnegCasimir_noReal (h : NonnegCasimirRealization) (hr : NoRealStripZeros) :
    RiemannHypothesis := by
  have hg := gap_of_nonnegCasimir_noReal h hr
  have hstrip : ∀ ρ ∈ nontrivialZeros, ρ.re = 1 / 2 := fun ρ hρ => by
    have := congrArg Complex.re (hg ρ hρ)
    simp at this; linarith
  exact GppGravitonDeltaId.rh_of_graviton_zero_identification fun ρ hρ =>
    hstrip ρ ((GppXiReflection.xi_eq_zero_iff_nontrivial ρ).1 hρ)

end GppArithmeticCFTRoutes

#print axioms GppArithmeticCFTRoutes.one_lt_sqrt'
#print axioms GppArithmeticCFTRoutes.qp_pos
#print axioms GppArithmeticCFTRoutes.qp_lt_one
#print axioms GppArithmeticCFTRoutes.log_pos'
#print axioms GppArithmeticCFTRoutes.kappaP_pos
#print axioms GppArithmeticCFTRoutes.qp_div
#print axioms GppArithmeticCFTRoutes.kappaP_eq_tateCoupling
#print axioms GppArithmeticCFTRoutes.geom_qp
#print axioms GppArithmeticCFTRoutes.geom_neg_qp
#print axioms GppArithmeticCFTRoutes.summable_symTerm
#print axioms GppArithmeticCFTRoutes.hasSum_log_geom
#print axioms GppArithmeticCFTRoutes.primeSymbol_zero
#print axioms GppArithmeticCFTRoutes.primeSymbol_le
#print axioms GppArithmeticCFTRoutes.primeSymbol_ge
#print axioms GppArithmeticCFTRoutes.primeSymbol_at_pi
#print axioms GppArithmeticCFTRoutes.localPrimeMultiplier_indefinite
#print axioms GppArithmeticCFTRoutes.kappaP_not_summable
#print axioms GppArithmeticCFTRoutes.summable_primeLocal
#print axioms GppArithmeticCFTRoutes.primeLocal_le
#print axioms GppArithmeticCFTRoutes.primeLocal_ge
#print axioms GppArithmeticCFTRoutes.primeLocal_const
#print axioms GppArithmeticCFTRoutes.weilS_floor
#print axioms GppArithmeticCFTRoutes.semiLocal_floor_unbounded
#print axioms GppArithmeticCFTRoutes.arp_iff_all_rungs
#print axioms GppArithmeticCFTRoutes.bilinear_shadow_pairing_invariant
#print axioms GppArithmeticCFTRoutes.hermitian_invariant_iff
#print axioms GppArithmeticCFTRoutes.casimir_nonneg_real_iff
#print axioms GppArithmeticCFTRoutes.complementary_casimir
#print axioms GppArithmeticCFTRoutes.nonnegCasimir_of_principalCasimir
#print axioms GppArithmeticCFTRoutes.eigen_real_nonneg
#print axioms GppArithmeticCFTRoutes.gap_of_nonnegCasimir_noReal
#print axioms GppArithmeticCFTRoutes.rh_of_nonnegCasimir_noReal
