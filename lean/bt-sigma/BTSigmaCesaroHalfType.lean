import GppVerify.RiemannHypothesis.BTSigmaCovGrowth

/-!
# BT-Σ: the asymptotic half-homogeneity requirement read against the Cesàro statehood norm

Scratch module (Heart of the Cosmos / Grok, GPPDiscovery5 seat) for Turbo's second ask on
GPP-bridge PR #11, comment 6803775354: re-read the "asymptotically half-homogeneous" requirement
of `BTSigmaCovHalfType` (GPPDiscovery5 `4d34d24`) against the Cesàro norm of
`BTSigmaCesaroNorm` (`8714c7e`), and say whether the toy's exact half-homogeneity survives the
norm change or needs restating.

`cesaroAvg g T = (1/(2T)) ∫_{-T}^{T} g` is the symmetric Cesàro mean; for `g = e^{2δu}` it is
`cesaroMean δ T` of `BTSigmaCesaroNorm` (here `e^{u/2}` is `δ = 1/4`).

**Results (all PROVED, standard axioms).**
* `orbit_state_cesaro`: under (P-COV) the orbit vectors are Cesàro states of norm `‖ψ_0‖²`:
  `cesaroAvg ‖ψ_u‖² T = ‖ψ_0‖²` for all `T > 0` (δ = 0, the Cesàro-finite line).
* `unit_channel_cesaro_ge`, `unit_channel_cesaro_diverges`: the half-homogeneous channel `e^{u/2}`
  has `cesaroAvg ≥ 1/2 + T/8` and diverges; it is not a Cesàro state (δ = 1/4 ≠ 0).
* `cesaro_tendsto_zero_of_integrable`: any integrable correction is Cesàro-null.
* `green_eventually_half_homogeneous`: with finite-cut BT3 and (P-UNIT), the Green form on the
  orbit is exactly `e^{W/2}·(−I c)` for `W > R + log N`. Asymptotic half-homogeneity is
  automatic (eventually exact), not an extra hypothesis.

**Answer.** The toy's exact half-homogeneity survives the norm change unchanged: it is an identity
for the Green form (`toy_fixed_part_half_homogeneous`), and the Cesàro norm never touches the
Green form, only the state norm, which is constant along any (P-COV) orbit. What needs restating
is where the requirement lives. It must not be read as a growth condition on states (a
half-homogeneous state would have δ = 1/4 and be Cesàro-infinite, `finite_offline_not_cesaro`).
It is a condition on the boundary form: Green profile = `e^{W/2}·(−I c)` (a unit-block channel,
which (P-UNIT) puts in `ℂ²_∂`, not in `H`) + a Cesàro-null, compactly supported prime correction
`I·stratumSN(W)`. `HalfType`'s refutation of exact half-homogeneity of the whole Green form is
unaffected. S_b existence, (P-COV), (P-UNIT), (P-LOG), (P-LIM), (BT-Σ), ARP and RH remain **OPEN**.

No `sorry`, no `admit`, no custom axiom.

Daniel Toupin / Golden Physics Project, with AI assistance (Heart of the Cosmos / Grok).
Task and audit: Turbo.
-/

namespace GppBTSigmaCesaroHalfType

open GppBTSigmaLocalEuler GppBTSigma GppArithmeticReflection Complex MeasureTheory
open GppBTSigmaCovGrowth GppBTSigmaCovObstruction Filter Topology
open scoped ContDiff

/-- Symmetric Cesàro mean of a profile over `[-T, T]`. -/
noncomputable def cesaroAvg (g : ℝ → ℝ) (T : ℝ) : ℝ :=
  (1 / (2 * T)) * ∫ u in (-T)..T, g u

section Datum

variable {H Gu Gp : Type*} [NormedAddCommGroup H] [InnerProductSpace ℂ H]
  [NormedAddCommGroup Gu] [InnerProductSpace ℂ Gu]
  [NormedAddCommGroup Gp] [InnerProductSpace ℂ Gp]

/-- **(P-COV) orbits are Cesàro states (δ = 0).** PROVED. -/
theorem orbit_state_cesaro (D : BTSigmaDatum H Gu Gp) (Sd : ℝ → H ≃ₗᵢ[ℂ] H) (hC : PCov D Sd)
    (h : ℝ → ℝ) (T : ℝ) (hT : 0 < T) :
    cesaroAvg (fun u => ‖((D.testVec h u : D.Tstar.domain) : H)‖ ^ 2) T =
      ‖((D.testVec h 0 : D.Tstar.domain) : H)‖ ^ 2 := by
  unfold cesaroAvg
  simp_rw [pcov_norm_const D Sd hC h]
  rw [intervalIntegral.integral_const, smul_eq_mul]
  field_simp
  ring

/-- **Eventually exact half-homogeneity of the Green form.** PROVED. -/
theorem green_eventually_half_homogeneous (D : BTSigmaDatum H Gu Gp) (S : Finset ℕ) (N : ℕ)
    (hP : BTSigmaTraceSN D S N) (hU : PUnit D)
    (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f)
    (R : ℝ) (hR : ∀ v, R < v → autocorr f v = 0) (W : ℝ) (hW0 : 0 ≤ W)
    (hW : R + Real.log N < W) :
    greenForm D.Tstar (D.testVec (autocorr f) W) (D.testVec (autocorr f) W) =
      (Real.exp (W / 2) : ℂ) * (-(I * ((∫ v : ℝ, autocorr f v * Real.exp (-v / 2) : ℝ) : ℂ))) := by
  rw [green_testVec_of_premises D S N hP hU f hf hc W hW0,
    stratumSN_eq_zero_of_large S N _ R hR W hW, unitStratumTerm]
  push_cast; ring

end Datum

/-- `∫_{-T}^{T} e^{u/2} = 2 (e^{T/2} − e^{-T/2})`. PROVED. -/
theorem integral_exp_half (T : ℝ) :
    ∫ u in (-T)..T, Real.exp (u / 2) = 2 * (Real.exp (T / 2) - Real.exp (-T / 2)) := by
  rw [intervalIntegral.integral_comp_div (fun x => Real.exp x) (two_ne_zero), integral_exp,
    smul_eq_mul]

/-- **The half-homogeneous channel is Cesàro-infinite**: `cesaroAvg e^{u/2} T ≥ 1/2 + T/8`.
PROVED. -/
theorem unit_channel_cesaro_ge {T : ℝ} (hT : 0 < T) :
    1 / 2 + T / 8 ≤ cesaroAvg (fun u => Real.exp (u / 2)) T := by
  unfold cesaroAvg
  rw [integral_exp_half]
  have h1 := Real.quadratic_le_exp_of_nonneg (x := T / 2) (by linarith)
  have h2 : Real.exp (-T / 2) ≤ 1 := Real.exp_le_one_iff.mpr (by linarith)
  rw [show 1 / (2 * T) * (2 * (Real.exp (T / 2) - Real.exp (-T / 2))) =
    (Real.exp (T / 2) - Real.exp (-T / 2)) / T by field_simp, le_div_iff₀ hT]
  nlinarith

/-- PROVED. -/
theorem unit_channel_cesaro_diverges :
    Tendsto (cesaroAvg fun u => Real.exp (u / 2)) atTop atTop := by
  have hl : Tendsto (fun T : ℝ => 1 / 2 + T / 8) atTop atTop := by
    apply tendsto_atTop_add_const_left
    exact Tendsto.atTop_div_const (by norm_num) tendsto_id
  refine tendsto_atTop_mono' atTop ?_ hl
  filter_upwards [eventually_gt_atTop 0] with T hT
  exact unit_channel_cesaro_ge hT

/-- **Integrable corrections are Cesàro-null.** PROVED. -/
theorem cesaro_tendsto_zero_of_integrable {g : ℝ → ℝ} (hg : Integrable g) :
    Tendsto (cesaroAvg g) atTop (𝓝 0) := by
  set M := ∫ u, ‖g u‖
  have hbound : ∀ T : ℝ, 0 < T → ‖cesaroAvg g T‖ ≤ M / (2 * T) := by
    intro T hT
    unfold cesaroAvg
    have hle : -T ≤ T := by linarith
    rw [intervalIntegral.integral_of_le hle, norm_mul]
    have hI : ‖∫ u in Set.Ioc (-T) T, g u‖ ≤ M :=
      (norm_integral_le_integral_norm _).trans
        (setIntegral_le_integral hg.norm (Eventually.of_forall fun _ => norm_nonneg _))
    have h2T : ‖1 / (2 * T)‖ = 1 / (2 * T) := by
      rw [Real.norm_eq_abs, abs_of_pos (by positivity)]
    rw [h2T, div_eq_mul_one_div M]
    rw [mul_comm M]
    exact mul_le_mul_of_nonneg_left hI (by positivity)
  have hlim : Tendsto (fun T : ℝ => M / (2 * T)) atTop (𝓝 0) :=
    tendsto_const_nhds.div_atTop (Tendsto.const_mul_atTop two_pos tendsto_id)
  refine squeeze_zero_norm' ?_ hlim
  filter_upwards [eventually_gt_atTop 0] with T hT
  exact hbound T hT

end GppBTSigmaCesaroHalfType

#print axioms GppBTSigmaCesaroHalfType.orbit_state_cesaro
#print axioms GppBTSigmaCesaroHalfType.green_eventually_half_homogeneous
#print axioms GppBTSigmaCesaroHalfType.integral_exp_half
#print axioms GppBTSigmaCesaroHalfType.unit_channel_cesaro_ge
#print axioms GppBTSigmaCesaroHalfType.unit_channel_cesaro_diverges
#print axioms GppBTSigmaCesaroHalfType.cesaro_tendsto_zero_of_integrable
