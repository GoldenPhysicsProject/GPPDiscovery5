import GppVerify.RiemannHypothesis.BTSigmaLocalEuler

/-!
# BT-Σ: genuine BT2 — the commuting-symmetry obstruction (Lean part)

Scratch module (Heart of the Cosmos / Grok, GPPDiscovery5 seat) for Turbo's question on
GPP-bridge PR #11, comment 6808529656: does a datum exist with (P-COV), finite-cut BT3, (P-UNIT),
`Γ_u` graph-bounded and `H`-unbounded, and dense `ker Γ` (`T*` the adjoint of a densely defined
closed symmetric `T`)?

The paper construction (half-line Laplacians, see the bridge note) answers yes. This module
PROVES the structural fact that shapes it: if the covariance group `S_b` is unitary and
*commutes* with `T*` on the test vectors, the diagonal Green form is `W`-invariant
(`greenForm_invariant_of_commuting`), so with BT3 and (P-UNIT) the truncated flux defect
`stratumSN − unitStratumTerm` is constant in `W ≥ 0` (`flux_constant_of_commuting`). Since the
unit term grows like `e^{W/2}`, the covariance of any genuine datum must fail to commute with
`T*`: the group cannot be a symmetry of `T`.

**Paper construction (genuine BT2; CONSTRUCTED on paper, Lean port GAP).**
`H = L²(0,∞) ⊕ L²((0,∞)×ℝ_y)`, `T = T_u ⊕ T_p` with `T_u`, `T_p = −∂_x² ⊗ 1` the minimal
(closure from `C_c^∞(0,∞)`) operators: densely defined, closed, symmetric; `T*` = `−∂_x²` on
`H²(0,∞) ⊕ H²(0,∞; L²(ℝ_y))`. Standard triplet `Γ₀u f = f(0)`, `Γ₁u f = f'(0)`,
`Γ₀p F = F(0,·)`, `Γ₁p F = ∂_xF(0,·)` (`Gu = ℂ`, `Gp = L²(ℝ)`): Green's identity
`greenForm = blockForm_u + blockForm_p` by integration by parts on `[0,∞)`, `Γ` onto, and
`ker Γ = dom T` is dense. Every `Γ` component is graph-bounded (Sobolev trace on the half-line)
and `H`-unbounded (vanishes on a dense set, nonzero).
Covariance `S_b = D_b ⊕ (P_φ ⊗ 1 + (1 − P_φ) ⊗ τ_b)`: `D_b f(x) = e^{b/8} f(e^{b/4} x)` (unitary
dilation group), `P_φ` the projection on a cutoff `φ ∈ C_c^∞[0,∞)`, `φ = 1` on `[0,1]`, `τ_b`
translation in `y`. Fix `ρ ∈ C_c^∞[0,∞)` with `ρ(0) = 0`, `ρ'(0) = 1`, `ρ ⊥ φ`
(`ρ = xφ − βη`, `η` a bump in `(1/2,1)`). Test vectors
`ψ_W = S_W(φ + iκρ, φ ⊗ f + ρ ⊗ b_h)`, `κ = ½∫h e^{−v/2}`, `b_h` as in BTSigmaSbExistence.
Boundary values: unit `(e^{W/8}, iκ e^{3W/8})`, so `unitBlock = −2iκ e^{W/2} = −I·unitStratumTerm`
(P-UNIT); prime `(f, τ_W b_h)`, exactly the boundary data of BTSigmaSbExistence, so finite-cut
BT3 holds by the same computation; (P-COV) by construction. Necessarily `S_b` does not commute
with `T*` (this module). (P-θ) FAILS for this family: `T ≥ 0`, so no unitary `J` with
`J T J = −T`. (P-LIM) fails as before; (P-LOG) provenance still inserted.

Lean results (3) PROVED, standard axioms. No `sorry`, no custom axiom, no cited theorem.
RH, ARP and (BT-Σ) remain OPEN.

Daniel Toupin / Golden Physics Project, with AI assistance (Heart of the Cosmos / Grok).
Task and audit: Turbo.
-/

namespace GppBTSigmaGenuineBT2

open GppBTSigma GppBTSigmaLocalEuler GppArithmeticReflection Complex
open scoped InnerProductSpace ContDiff

section Green

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℂ E]

/-- A unitary map that carries `x` into the domain and intertwines `T*` preserves the diagonal
Green form. PROVED. -/
theorem greenForm_invariant_of_commuting (T : E →ₗ.[ℂ] E) (S : E →ₗᵢ[ℂ] E)
    (x y : T.domain) (hy : (y : E) = S x) (hTy : T y = S (T x)) :
    greenForm T y y = greenForm T x x := by
  simp only [greenForm, hy, hTy, LinearIsometry.inner_map_map]

end Green

section Datum

variable {H Gu Gp : Type*} [NormedAddCommGroup H] [InnerProductSpace ℂ H]
  [NormedAddCommGroup Gu] [InnerProductSpace ℂ Gu]
  [NormedAddCommGroup Gp] [InnerProductSpace ℂ Gp]

/-- Commuting covariance on the test family: `testVec h (W+b) = S_b testVec h W` and
`T* testVec h (W+b) = S_b T* testVec h W`. -/
def CommutingCov (D : BTSigmaDatum H Gu Gp) (S : ℝ → H →ₗᵢ[ℂ] H) : Prop :=
  ∀ (h : ℝ → ℝ) (W b : ℝ),
    ((D.testVec h (W + b) : D.Tstar.domain) : H) = S b (D.testVec h W) ∧
    D.Tstar (D.testVec h (W + b)) = S b (D.Tstar (D.testVec h W))

/-- Under commuting covariance the Green form on the test family is `W`-independent. PROVED. -/
theorem greenForm_testVec_const (D : BTSigmaDatum H Gu Gp) (S : ℝ → H →ₗᵢ[ℂ] H)
    (hS : CommutingCov D S) (h : ℝ → ℝ) (W : ℝ) :
    greenForm D.Tstar (D.testVec h W) (D.testVec h W) =
      greenForm D.Tstar (D.testVec h 0) (D.testVec h 0) := by
  obtain ⟨h1, h2⟩ := hS h 0 W
  rw [zero_add] at h1 h2
  exact greenForm_invariant_of_commuting D.Tstar (S W) _ _ h1 h2

/-- **Commuting-symmetry obstruction.** With finite-cut BT3, (P-UNIT) and a covariance group that
commutes with `T*`, the truncated flux defect is constant on `W ≥ 0`. PROVED. -/
theorem flux_constant_of_commuting (D : BTSigmaDatum H Gu Gp) (S₀ : Finset ℕ) (N : ℕ)
    (hP : BTSigmaTraceSN D S₀ N) (hU : PUnit D) (S : ℝ → H →ₗᵢ[ℂ] H)
    (hS : CommutingCov D S) (f : ℝ → ℝ) (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f)
    (W : ℝ) (hW : 0 ≤ W) :
    stratumSN S₀ N (autocorr f) W - unitStratumTerm (autocorr f) W =
      stratumSN S₀ N (autocorr f) 0 - unitStratumTerm (autocorr f) 0 := by
  have e := greenForm_testVec_const D S hS (autocorr f) W
  rw [green_testVec_of_premises D S₀ N hP hU f hf hc W hW,
    green_testVec_of_premises D S₀ N hP hU f hf hc 0 le_rfl] at e
  have := mul_left_cancel₀ Complex.I_ne_zero e
  exact_mod_cast this

end Datum

#print axioms greenForm_invariant_of_commuting
#print axioms greenForm_testVec_const
#print axioms flux_constant_of_commuting

end GppBTSigmaGenuineBT2
