import GppVerify.RiemannHypothesis.BTSigmaLocalEuler
import GppVerify.RiemannHypothesis.BTSigmaEulerSkeleton

/-!
# BT-Σ provenance step: what the vertex law produces, and where it stops

Scratch module (Heart of the Cosmos / Grok, GPPDiscovery5 seat) for Turbo's question on
GPP-bridge PR #11, comment 6800991395 §4(b): build `ψ` from `f` and the vertex law so that
the trace jumps are `cutWeight` by the vertex law, or name the exact step that fails.

Setting. The vertex law at spectral parameter `s` is `U_p(s) = p^{1/2 − s} · shift` on the
valuation model, so `⟨ξ_p, U_p(s)^j ξ_p⟩ = (p^{-s})^j` and at the unitary point `s = 1/2`
this is the (E1) value `p^{-j/2}`. Write `C_m` for the cross term on the piece after `m`
cuts of the `p`-channel; the jump at the `(m+1)`-st cut is `C_m − C_{m+1}`, and (J) needs it
to be `Λ(p^{m+1})(p^{m+1})^{-1/2}`.

* **(V1) PROVED. The log-determinant of the vertex law produces the cut weights, with no
  coupling inserted.** With `logDetTail p m N s = Σ_{m<j≤N} (p^{-s})^j / j` (the tail of
  `log L_p(s) = Σ_j ⟨ξ, U_p(s)^j ξ⟩ / j`), `−∂_s logDetTail` at `s = 1/2` is the tail charge
  `Σ_{m<j≤N} Λ(p^j)(p^j)^{-1/2}` (`hasDerivAt_logDetTail_half`), and its drop across one cut is
  exactly the next weight (`logDet_cut_jump`). `κ_p` appears only as the `N → ∞`, `m = 0` value
  (E4); it is not put into any amplitude.
* **(V2) PROVED. The `1/j` is load-bearing.** Without it, `−∂_s ⟨ξ, U_p(s)^j ξ⟩` at `1/2` is
  `j · Λ(p^j)(p^j)^{-1/2}` (`hasDerivAt_pairing_half`): the derivative of `L_p`, not of
  `log L_p`, and wrong by the repetition number `j` for every `j ≥ 2`.
* **(V3) PROVED. Obstruction for sesquilinear vertex laws of degree ≤ 1.** A `ψ` built
  linearly from `f`, `ξ_p`, `U_p(s)` and one `s`-derivative has cross terms
  `C_m = α p^{-m/2} + β · m log p · p^{-m/2}` (the `ξ` and `∂_s ξ` channels). If one
  prime-independent pair `(α, β)` gave the weights at every cut of `p = 2` and `p = 3`, then
  (`vertexJump_forces_kappa`, from cuts 1 and 2 alone) `α = κ_2` and `α = κ_3`, contradicting
  `κ_2 > 8/5 > κ_3` (`not_PLog_vertexJump`).

**The step that fails** (stated in the bridge note, not typed here): a Green-form trace jump
is sesquilinear in `ψ`, so it yields matrix elements of words in `U_p(s)` and `∂_s`. The cut
weights are matrix elements of `−∂_s log(1 − U_p(s))^{-1}`, whose `1/j` weighting
(primitive-orbit normalisation, `Λ(p^j) = log p`) is not a finite sesquilinear word; at the
unitary point `1 − U_p` is not trace class, so the log-determinant has no operator meaning
there. A `ψ_{f,W}` with these jumps by the vertex law therefore needs a boundary space where
`log det` is itself a quadratic form (for example exponential vectors in a Fock space over
`V_p`, where `log ⟨e(x), e(y)⟩ = ⟨x, y⟩`). That is a CONJECTURE-level direction, not a result.
(P-LOG), (BT-Σ), ARP and RH remain **OPEN**.

No `sorry`, no `admit`, no custom axiom.

Daniel Toupin / Golden Physics Project, with AI assistance (Heart of the Cosmos / Grok).
Task and audit: Turbo.
-/

namespace GppBTSigmaVertexProvenance

open GppBTSigmaLocalEuler GppBTSigmaC1 ArithmeticFunction

/-! ## (V1), (V2): derivatives of the vertex-law pairings -/

/-- `d/ds (p^{-s})^j = −j log p (p^{-s})^j`. -/
theorem hasDerivAt_pow_rpow {p : ℕ} (hp : 1 < p) (j : ℕ) (s : ℝ) :
    HasDerivAt (fun t : ℝ => ((p : ℝ) ^ (-t)) ^ j)
      (-((j : ℝ) * Real.log p * ((p : ℝ) ^ (-s)) ^ j)) s := by
  have h1 : HasDerivAt (fun t : ℝ => (p : ℝ) ^ (-t))
      ((p : ℝ) ^ (-s) * Real.log p * (-1)) s :=
    ((Real.hasStrictDerivAt_const_rpow (cast_pos_of_one_lt hp) (-s)).hasDerivAt).comp s
      (hasDerivAt_neg s)
  have h2 : HasDerivAt (fun t : ℝ => ((p : ℝ) ^ (-t)) ^ j)
      ((j : ℝ) * ((p : ℝ) ^ (-s)) ^ (j - 1) * ((p : ℝ) ^ (-s) * Real.log p * (-1))) s :=
    h1.pow j
  convert h2 using 1
  rcases j with _ | j
  · simp
  · rw [Nat.add_sub_cancel, pow_succ]
    ring

theorem rpow_neg_half_eq_qp (p : ℕ) : (p : ℝ) ^ (-(1 / 2 : ℝ)) = qp p := (qp_eq_rpow p).symm

/-- **(V2).** `−∂_s ⟨ξ_p, U_p(s)^j ξ_p⟩` at `s = 1/2` is `j · Λ(p^j)(p^j)^{-1/2}`. PROVED. -/
theorem hasDerivAt_pairing_half {p : ℕ} (hp : p.Prime) {j : ℕ} (hj : j ≠ 0) :
    HasDerivAt (fun t : ℝ => ((p : ℝ) ^ (-t)) ^ j) (-((j : ℝ) * halfWeight (p ^ j))) (1 / 2) := by
  have h := hasDerivAt_pow_rpow hp.one_lt j (1 / 2)
  rw [rpow_neg_half_eq_qp] at h
  rw [halfWeight_primePow hp hj]
  convert h using 1
  ring

/-- Tail of `log L_p(s) = Σ_j (p^{-s})^j / j` beyond `m` cuts, truncated at `N`. -/
noncomputable def logDetTail (p m N : ℕ) (s : ℝ) : ℝ :=
  ∑ j ∈ Finset.Ico (m + 1) (N + 1), ((p : ℝ) ^ (-s)) ^ j / j

/-- Tail charge `Σ_{m<j≤N} Λ(p^j)(p^j)^{-1/2}`. -/
noncomputable def tailW (p m N : ℕ) : ℝ := ∑ j ∈ Finset.Ico (m + 1) (N + 1), halfWeight (p ^ j)

/-- **(V1).** `−∂_s logDetTail` at the centre is the tail charge. PROVED. -/
theorem hasDerivAt_logDetTail_half {p : ℕ} (hp : p.Prime) (m N : ℕ) :
    HasDerivAt (logDetTail p m N) (-tailW p m N) (1 / 2) := by
  have h : ∀ j ∈ Finset.Ico (m + 1) (N + 1),
      HasDerivAt (fun t : ℝ => ((p : ℝ) ^ (-t)) ^ j / j) (-halfWeight (p ^ j)) (1 / 2) := by
    intro j hj
    have hj0 : j ≠ 0 := by simp only [Finset.mem_Ico] at hj; omega
    have hjr : (j : ℝ) ≠ 0 := by exact_mod_cast hj0
    have h' := (hasDerivAt_pairing_half hp hj0).div_const (j : ℝ)
    rwa [neg_div, mul_div_cancel_left₀ _ hjr] at h'
  have hs := HasDerivAt.sum h
  unfold logDetTail tailW
  rw [Finset.sum_neg_distrib] at hs
  have e : (fun t : ℝ => ∑ j ∈ Finset.Ico (m + 1) (N + 1), ((p : ℝ) ^ (-t)) ^ j / j) =
      ∑ j ∈ Finset.Ico (m + 1) (N + 1), fun t : ℝ => ((p : ℝ) ^ (-t)) ^ j / j := by
    funext t
    simp [Finset.sum_apply]
  rw [e]
  exact hs

/-- **(V1), cut form.** The log-determinant tail drops by exactly the next weight. PROVED. -/
theorem logDet_cut_jump (p : ℕ) {m N : ℕ} (hmN : m < N) :
    tailW p m N - tailW p (m + 1) N = halfWeight (p ^ (m + 1)) := by
  unfold tailW
  rw [Finset.sum_eq_sum_Ico_succ_bot (by omega : m + 1 < N + 1)]
  ring

/-! ## (V3): degree-≤1 sesquilinear vertex laws cannot carry the coupling -/

/-- Cross term on the piece after `m` cuts for a `ψ` built from `ξ_p` and `∂_s ξ_p` with
amplitude products `(α, β)`: `α p^{-m/2} + β · m log p · p^{-m/2}`. -/
noncomputable def vertexCross (p : ℕ) (ab : ℝ × ℝ) (m : ℕ) : ℝ :=
  ab.1 * qp p ^ m + ab.2 * ((m : ℝ) * Real.log p) * qp p ^ m

/-- Jump at the `(m+1)`-st cut. -/
noncomputable def vertexJump (p m : ℕ) (ab : ℝ × ℝ) : ℝ :=
  vertexCross p ab m - vertexCross p ab (m + 1)

/-- Matching the first two cuts forces `α = κ_p`. PROVED. -/
theorem vertexJump_forces_kappa {p : ℕ} (hp : p.Prime) (ab : ℝ × ℝ)
    (h0 : vertexJump p 0 ab = halfWeight (p ^ 1)) (h1 : vertexJump p 1 ab = halfWeight (p ^ 2)) :
    ab.1 = kappa p := by
  have hq0 := qp_pos hp.one_lt
  have hq1 : qp p < 1 := by
    unfold qp
    exact inv_lt_one_of_one_lt₀ (one_lt_sqrt hp.one_lt)
  have hk : kappa p * (Real.sqrt p - 1) = Real.log p := kappa_mul hp.one_lt
  have hsq : Real.sqrt p * qp p = 1 := mul_inv_cancel₀ (sqrt_ne_zero' hp.one_lt)
  rw [halfWeight_primePow hp (by norm_num)] at h0 h1
  unfold vertexJump vertexCross at h0 h1
  push_cast at h0 h1
  set q := qp p
  set L := Real.log p
  set x := ab.1 - kappa p
  -- `log p · q = κ (1 − q)` from `κ(√p − 1) = log p` and `√p q = 1`.
  have hL : L * q = kappa p * (1 - q) := by linear_combination q * hk.symm + kappa p * hsq
  -- the two cut equations in the variables `x = α − κ`, `y = β log p`
  have e1 : x - x * q - ab.2 * L * q = 0 := by linear_combination h0 + hL
  have e2 : x * q - x * q ^ 2 - 2 * (ab.2 * L) * q ^ 2 + ab.2 * L * q = 0 := by
    linear_combination h1 + q * hL
  -- eliminate `y`: `x (1 − q)² = 0`
  have e3 : x * (1 - q) ^ 2 = 0 := by linear_combination e2 + (1 - 2 * q) * e1
  have hne : (1 - q) ^ 2 ≠ 0 := pow_ne_zero _ (by linarith)
  have hx : x = 0 := (mul_eq_zero.mp e3).resolve_right hne
  simp only [x] at hx
  linarith

/-- **(P-LOG) for the degree-≤1 vertex law, coupling-free form: refuted.** No single
prime-independent `(α, β)` reproduces every cut weight at `p = 2` and `p = 3`. PROVED. -/
theorem not_PLog_vertexJump :
    ¬ ∃ ab : ℝ × ℝ, ∀ p : ℕ, p.Prime → ∀ m : ℕ, vertexJump p m ab = halfWeight (p ^ (m + 1)) := by
  rintro ⟨ab, h⟩
  have h2 := vertexJump_forces_kappa Nat.prime_two ab (h 2 Nat.prime_two 0) (h 2 Nat.prime_two 1)
  have h3 := vertexJump_forces_kappa Nat.prime_three ab (h 3 Nat.prime_three 0)
    (h 3 Nat.prime_three 1)
  have k2 : (8 / 5 : ℝ) < kappa 2 := GppBTSigmaEuler.kappa_two_gt
  have k3 : kappa 3 < 8 / 5 := GppBTSigmaEuler.kappa_three_lt
  linarith

end GppBTSigmaVertexProvenance
