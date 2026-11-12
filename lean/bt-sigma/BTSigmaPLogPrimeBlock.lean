import GppVerify.RiemannHypothesis.LocalTateEvaluations
import GppVerify.RiemannHypothesis.BTSigmaHjumpWitness

/-!
# BT-Σ (P-LOG): the prime-block coupling from the local Tate scaling — D5 companion module

Heart of the Cosmos / Grok, GPPDiscovery5 seat. Task: Turbo, GPP-bridge PR #11, comments
6810320113 (audited spec) and 6811509097. Companion to GPPVerify #579 (`BTSigmaPrimeBlock.lean`,
head 98e6312); this file adds what #579 lists as GAP or does not cover. Compiled against
GPPVerify main 9242027 (`lake env lean`), Lean 4.33.1.

Unramified sector, multiplicative Tate measure `d^×x`, `vol(ℤ_p^×) = 1` (correction 1): a
`ℤ_p^×`-invariant function supported in `ℤ_p` is a shell profile `a k` on `p^k ℤ_p^×`, and
`Z(f, s) = Σ_k a_k p^{-ks}` (`padicTate`; equal to the in-tree `LocalTate.Finite.localZetaPrime`
for `1_{ℤ_p}`). The shell-series form of the integral is CHECKED on paper; a Lean Haar
measure on `ℚ_p^×` is not built (GAP).

PROVED (32 theorems, axioms propext / Classical.choice / Quot.sound only, 0 sorry,
0 native stubs):
* `local_scaling_trace_critical`: `Z(log|·|_p 1_{ℤ_p}, 1/2)/Z(1_{ℤ_p}, 1/2) = −κ_p`, shell series.
* `padicTate_charZp_half`: `Z(1_{ℤ_p}, 1/2) = √p/(√p − 1) > 0` (correction 4);
  `primeBlockJump_unramified_recovers_kappa`.
* `tate_scaling_is_deriv_shell`: `Z(log|·|_p f, s) = ∂_s Z(f, s)`, `Re s > 0`, for **every**
  unramified Schwartz–Bruhat profile supported in `ℤ_p` (eventually constant shell values),
  not only `1_{ℤ_p}` (correction 3). Ramified / non-`ℤ_p`-supported `f`: GAP (obligation).
* `localScaling_charZp_not_schwartz`: `log|·|_p 1_{ℤ_p}` is not Schwartz–Bruhat (correction 2).
* `additive_ratio_ne_kappa`: under additive Haar the ratio is `−log p/(p^{3/2} − 1)` and is
  **not** `−κ_p` (correction 1 guard).
* Rewiring: `derivedCoupling := −Re(primeBlockJump/Z)` (no `κ_p` in any definition; an
  environment reachability gate confirms `derivedCutWeight`, `derivedCoupling`,
  `primeBlockJump` do not reach `kappa`, while the tree `cutWeight` does);
  `derivedCoupling_eq_kappa` is a conclusion; `finiteCut_BT3_derived`;
  `finiteCut_BT3_derived_witness`: an explicit piecewise `ψ` whose trace jumps ARE the derived
  jumps, with Green form `= I · stratumSN` (end to end, no `hjump` hypothesis left).

Not claimed: the witness `ψ` is weight-built (`psiW`), not the spec's `ψ_{f,W}`; no
`BTSigmaDatum` is constructed; the centre `1/2` is the evaluation point (weakened (P-θ) OPEN).
Whether this meets the (P-LOG) closure condition is the auditors' ruling. ARP, CGC, RH OPEN.

Daniel Toupin / Golden Physics Project, with AI assistance (Heart of the Cosmos / Grok).
Task and audit: Turbo; cross-audit: Flint.
-/

namespace GppBTSigmaPLog

open GppBTSigmaLocalEuler GppBTSigmaHjumpWitness GppBTSigmaC1 GppBTSigma Complex

/-! ## The unramified local sector, multiplicative Tate measure -/

section Local

variable (p : ℕ)

/-- `r(s) = p^{-s} = |x|_p^s` on the first shell `p ℤ_p^×`. -/
noncomputable def rp (s : ℂ) : ℂ := (p : ℂ) ^ (-s)

/-- A `ℤ_p^×`-invariant function supported in `ℤ_p` is a shell profile `a : ℕ → ℂ`
(`a k` = value on `p^k ℤ_p^×`). It is Schwartz–Bruhat iff it is locally constant at `0`,
i.e. eventually constant in `k`. -/
def IsShellSchwartz (a : ℕ → ℂ) : Prop := ∃ K : ℕ, ∃ c : ℂ, ∀ k, K ≤ k → a k = c

/-- `1_{ℤ_p}` as a shell profile. -/
def charZp : ℕ → ℂ := fun _ => 1

/-- Multiplication by `log|x|_p`; on `p^k ℤ_p^×` it is `−k log p`. -/
noncomputable def localScaling (a : ℕ → ℂ) : ℕ → ℂ :=
  fun k => -((k : ℂ) * (Real.log p : ℂ)) * a k

/-- The local Tate integral `∫ f(x) |x|_p^s d^×x` with `vol(ℤ_p^×) = 1`, as the shell series
`Σ_k a_k p^{-ks}` (each shell `p^k ℤ_p^×` has multiplicative volume `1`). -/
noncomputable def padicTate (a : ℕ → ℂ) (s : ℂ) : ℂ := ∑' k : ℕ, a k * rp p s ^ k

/-- The prime block jump: the Tate trace of the local scaling at the centre. No `κ_p` input. -/
noncomputable def primeBlockJump (a : ℕ → ℂ) : ℂ := padicTate p (localScaling p a) (1 / 2)

end Local

theorem norm_rp_lt_one {p : ℕ} (hp : 1 < p) {s : ℂ} (hs : 0 < s.re) : ‖rp p s‖ < 1 := by
  unfold rp
  rw [Complex.norm_natCast_cpow_of_pos (by omega)]
  simp only [Complex.neg_re]
  exact Real.rpow_lt_one_of_one_lt_of_neg (by exact_mod_cast hp) (by linarith)

theorem one_sub_rp_ne_zero {p : ℕ} (hp : 1 < p) {s : ℂ} (hs : 0 < s.re) : 1 - rp p s ≠ 0 := by
  intro h
  have h1 : rp p s = 1 := (sub_eq_zero.mp h).symm
  have := norm_rp_lt_one hp hs
  rw [h1, norm_one] at this
  exact lt_irrefl _ this

/-- `Z(1_{ℤ_p}, s) = (1 − p^{-s})^{-1}`, the in-tree `LocalTate.Finite.localZetaPrime`. -/
theorem padicTate_charZp {p : ℕ} (hp : 1 < p) {s : ℂ} (hs : 0 < s.re) :
    padicTate p (charZp) s = (1 - rp p s)⁻¹ := by
  unfold padicTate charZp
  simp only [one_mul]
  exact tsum_geometric_of_norm_lt_one (norm_rp_lt_one hp hs)

theorem padicTate_charZp_eq_localZetaPrime (p : ℕ) (s : ℂ) :
    padicTate p (charZp) s = LocalTate.Finite.localZetaPrime p s := by
  unfold padicTate charZp LocalTate.Finite.localZetaPrime rp
  simp only [one_mul]

/-- `Z(log|·|_p · 1_{ℤ_p}, s) = −log p · r/(1 − r)²`, termwise from the shell series. -/
theorem padicTate_localScaling_charZp {p : ℕ} (hp : 1 < p) {s : ℂ} (hs : 0 < s.re) :
    padicTate p (localScaling p (charZp)) s =
      -(Real.log p : ℂ) * (rp p s / (1 - rp p s) ^ 2) := by
  unfold padicTate localScaling charZp
  rw [← tsum_coe_mul_geometric_of_norm_lt_one (norm_rp_lt_one hp hs), ← tsum_mul_left]
  congr 1
  funext k
  ring

theorem hasDerivAt_rp {p : ℕ} (hp : 1 < p) (s : ℂ) :
    HasDerivAt (rp p) (-(Real.log p : ℂ) * rp p s) s := by
  have hp0 : (p : ℂ) ≠ 0 := by exact_mod_cast (by omega : p ≠ 0)
  have h := (hasDerivAt_neg s).const_cpow (c := (p : ℂ)) (Or.inl hp0)
  unfold rp
  convert h using 1
  rw [Complex.ofReal_log (Nat.cast_nonneg p), Complex.ofReal_natCast]
  ring

/-- **`tate_scaling_is_deriv`, unramified sector.** `Z(log|·|_p f, s) = ∂_s Z(f, s)` for
`f = 1_{ℤ_p}`, `Re s > 0`. PROVED. -/
theorem tate_scaling_hasDerivAt_charZp {p : ℕ} (hp : 1 < p) {s : ℂ} (hs : 0 < s.re) :
    HasDerivAt (padicTate p (charZp)) (padicTate p (localScaling p (charZp)) s) s := by
  have hd := ((hasDerivAt_rp hp s).const_sub 1).inv (one_sub_rp_ne_zero hp hs)
  have hopen : IsOpen {z : ℂ | 0 < z.re} := isOpen_lt continuous_const Complex.continuous_re
  have hev : padicTate p (charZp) =ᶠ[nhds s] fun z => (1 - rp p z)⁻¹ := by
    filter_upwards [hopen.mem_nhds hs] with z hz
    exact padicTate_charZp hp hz
  rw [padicTate_localScaling_charZp hp hs]
  have heq : -(Real.log p : ℂ) * (rp p s / (1 - rp p s) ^ 2) =
      - -(-(Real.log p : ℂ) * rp p s) / (1 - rp p s) ^ 2 := by ring
  rw [heq]
  exact hd.congr_of_eventuallyEq hev

theorem tate_scaling_is_deriv_charZp {p : ℕ} (hp : 1 < p) {s : ℂ} (hs : 0 < s.re) :
    padicTate p (localScaling p (charZp)) s = deriv (padicTate p (charZp)) s :=
  (tate_scaling_hasDerivAt_charZp hp hs).deriv.symm

/-- **Correction 2, PROVED.** `log|·|_p · 1_{ℤ_p}` is not Schwartz–Bruhat: it is not locally
constant at `0` (its shell values `−k log p` never stabilise). -/
theorem localScaling_charZp_not_schwartz {p : ℕ} (hp : 1 < p) :
    ¬ IsShellSchwartz (localScaling p (charZp)) := by
  rintro ⟨K, c, hc⟩
  have h1 := hc K le_rfl
  have h2 := hc (K + 1) (by omega)
  have hlog : (0 : ℝ) < Real.log p := Real.log_pos (by exact_mod_cast hp)
  have : localScaling p (charZp) K = localScaling p (charZp) (K + 1) := h1.trans h2.symm
  unfold localScaling charZp at this
  simp only [mul_one, Nat.cast_add, Nat.cast_one] at this
  have hl : (Real.log p : ℂ) = 0 := by linear_combination this
  exact hlog.ne' (by exact_mod_cast hl)

theorem charZp_schwartz : IsShellSchwartz (charZp) := ⟨0, 1, fun _ _ => rfl⟩

/-- `r(1/2) = p^{-1/2} = qp p`. -/
theorem rp_half {p : ℕ} : rp p (1 / 2) = ((qp p : ℝ) : ℂ) := by
  unfold rp
  rw [qp_eq_rpow, Complex.ofReal_cpow (Nat.cast_nonneg p)]
  push_cast
  ring_nf

/-- **Correction 4.** `Z(1_{ℤ_p}, 1/2) = √p/(√p − 1) > 0`, in particular `≠ 0`. PROVED. -/
theorem padicTate_charZp_half {p : ℕ} (hp : 1 < p) :
    padicTate p (charZp) (1 / 2) = ((Real.sqrt p / (Real.sqrt p - 1) : ℝ) : ℂ) := by
  rw [padicTate_charZp hp (by norm_num), rp_half]
  have hs := one_lt_sqrt hp
  have h1 : Real.sqrt p - 1 ≠ 0 := by linarith
  have h0 : Real.sqrt p ≠ 0 := by linarith
  unfold qp
  rw [← Complex.ofReal_one, ← Complex.ofReal_sub, ← Complex.ofReal_inv]
  congr 1
  field_simp

theorem padicTate_charZp_half_pos {p : ℕ} (hp : 1 < p) :
    0 < Real.sqrt p / (Real.sqrt p - 1) := by
  have hs := one_lt_sqrt hp
  apply div_pos <;> linarith

theorem padicTate_charZp_half_ne_zero {p : ℕ} (hp : 1 < p) :
    padicTate p (charZp) (1 / 2) ≠ 0 := by
  rw [padicTate_charZp_half hp]
  exact_mod_cast (padicTate_charZp_half_pos hp).ne'

/-- **`local_scaling_trace_critical`.** `Z(log|·|_p 1_{ℤ_p}, 1/2) / Z(1_{ℤ_p}, 1/2) = −κ_p`,
from the shell series alone. PROVED. -/
theorem local_scaling_trace_critical {p : ℕ} (hp : 1 < p) :
    primeBlockJump p (charZp) / padicTate p (charZp) (1 / 2) = -((kappa p : ℝ) : ℂ) := by
  unfold primeBlockJump
  rw [padicTate_localScaling_charZp hp (by norm_num), padicTate_charZp hp (by norm_num), rp_half]
  have hs := one_lt_sqrt hp
  have h1 : Real.sqrt p - 1 ≠ 0 := by linarith
  have h0 : Real.sqrt p ≠ 0 := by linarith
  have hq1 : (1 : ℝ) - qp p ≠ 0 := by
    unfold qp; intro h; field_simp at h; exact h1 (by linarith)
  unfold kappa qp
  have hq1' : (1 : ℂ) - ((Real.sqrt p)⁻¹ : ℝ) ≠ 0 := by exact_mod_cast hq1
  push_cast
  have hc0 : ((Real.sqrt p : ℝ) : ℂ) ≠ 0 := by exact_mod_cast h0
  have hc1 : ((Real.sqrt p : ℝ) : ℂ) - 1 ≠ 0 := by exact_mod_cast h1
  push_cast at hq1'
  field_simp

/-- **`primeBlockJump_unramified_recovers_kappa`.** PROVED. -/
theorem primeBlockJump_unramified_recovers_kappa {p : ℕ} (hp : 1 < p) :
    primeBlockJump p (charZp) = -((kappa p : ℝ) : ℂ) * padicTate p (charZp) (1 / 2) := by
  rw [← local_scaling_trace_critical hp, div_mul_cancel₀ _ (padicTate_charZp_half_ne_zero hp)]

/-! ## General unramified Schwartz–Bruhat profiles (supported in `ℤ_p`) -/

theorem summable_rp_pow {p : ℕ} (hp : 1 < p) {s : ℂ} (hs : 0 < s.re) :
    Summable (fun k : ℕ => rp p s ^ k) :=
  summable_geometric_of_norm_lt_one (norm_rp_lt_one hp hs)

/-- Shell decomposition: `a = (a − c) + c·1_{ℤ_p}` with `a − c` finitely supported. -/
theorem padicTate_shell_decomp {p : ℕ} (hp : 1 < p) {a : ℕ → ℂ} {K : ℕ} {c : ℂ}
    (ha : ∀ k, K ≤ k → a k = c) {s : ℂ} (hs : 0 < s.re) :
    padicTate p a s = ∑ k ∈ Finset.range K, (a k - c) * rp p s ^ k + c * (1 - rp p s)⁻¹ := by
  unfold padicTate
  have hg : Summable (fun k : ℕ => c * rp p s ^ k) := (summable_rp_pow hp hs).mul_left c
  have hz : ∀ k ∉ Finset.range K, (a k - c) * rp p s ^ k = 0 := by
    intro k hk
    rw [Finset.mem_range, not_lt] at hk
    rw [ha k hk, sub_self, zero_mul]
  have hb : Summable (fun k : ℕ => (a k - c) * rp p s ^ k) := summable_of_ne_finset_zero hz
  have he : (fun k : ℕ => a k * rp p s ^ k) =
      fun k => (a k - c) * rp p s ^ k + c * rp p s ^ k := by
    funext k; ring
  rw [he, hb.tsum_add hg, tsum_eq_sum hz, tsum_mul_left,
    tsum_geometric_of_norm_lt_one (norm_rp_lt_one hp hs)]

theorem padicTate_localScaling_shell_decomp {p : ℕ} (hp : 1 < p) {a : ℕ → ℂ} {K : ℕ} {c : ℂ}
    (ha : ∀ k, K ≤ k → a k = c) {s : ℂ} (hs : 0 < s.re) :
    padicTate p (localScaling p a) s =
      ∑ k ∈ Finset.range K, -((k : ℂ) * (Real.log p : ℂ)) * (a k - c) * rp p s ^ k +
        c * (-(Real.log p : ℂ) * (rp p s / (1 - rp p s) ^ 2)) := by
  rw [← padicTate_localScaling_charZp hp hs]
  unfold padicTate localScaling charZp
  have hg : Summable (fun k : ℕ => c * (-((k : ℂ) * (Real.log p : ℂ)) * 1 * rp p s ^ k)) := by
    have h1 := (hasSum_coe_mul_geometric_of_norm_lt_one (norm_rp_lt_one hp hs)).summable
    exact (h1.mul_left (-(c * (Real.log p : ℂ)))).congr (fun k => by ring)
  have hz : ∀ k : ℕ, k ∉ Finset.range K →
      -((k : ℂ) * (Real.log p : ℂ)) * (a k - c) * rp p s ^ k = 0 := by
    intro k hk
    rw [Finset.mem_range, not_lt] at hk
    rw [ha k hk, sub_self, mul_zero, zero_mul]
  have hb : Summable (fun k : ℕ => -((k : ℂ) * (Real.log p : ℂ)) * (a k - c) * rp p s ^ k) :=
    summable_of_ne_finset_zero hz
  have he : (fun k : ℕ => -((k : ℂ) * (Real.log p : ℂ)) * a k * rp p s ^ k) =
      fun k : ℕ => -((k : ℂ) * (Real.log p : ℂ)) * (a k - c) * rp p s ^ k +
        c * (-((k : ℂ) * (Real.log p : ℂ)) * 1 * rp p s ^ k) := by
    funext k; ring
  rw [he, hb.tsum_add hg, tsum_eq_sum hz, tsum_mul_left]

theorem hasDerivAt_rp_pow {p : ℕ} (hp : 1 < p) (s : ℂ) (k : ℕ) :
    HasDerivAt (fun z => rp p z ^ k) (-((k : ℂ) * (Real.log p : ℂ)) * rp p s ^ k) s := by
  have h := (hasDerivAt_rp hp s).pow k
  have heq : -((k : ℂ) * (Real.log p : ℂ)) * rp p s ^ k =
      (k : ℂ) * rp p s ^ (k - 1) * (-(Real.log p : ℂ) * rp p s) := by
    cases k with
    | zero => simp
    | succ m =>
      simp only [Nat.add_sub_cancel]
      push_cast
      ring
  rw [heq]
  exact h

/-- **`tate_scaling_is_deriv` on every unramified Schwartz–Bruhat profile supported in `ℤ_p`.**
`Z(log|·|_p f, s) = ∂_s Z(f, s)` for `Re s > 0`. PROVED (termwise on the finite part, closed
form on the tail; no differentiation under the integral sign). -/
theorem tate_scaling_hasDerivAt_shell {p : ℕ} (hp : 1 < p) {a : ℕ → ℂ} (ha : IsShellSchwartz a)
    {s : ℂ} (hs : 0 < s.re) :
    HasDerivAt (padicTate p a) (padicTate p (localScaling p a) s) s := by
  obtain ⟨K, c, hK⟩ := ha
  have hopen : IsOpen {z : ℂ | 0 < z.re} := isOpen_lt continuous_const Complex.continuous_re
  have hev : padicTate p a =ᶠ[nhds s] fun z =>
      ∑ k ∈ Finset.range K, (a k - c) * rp p z ^ k + c * (1 - rp p z)⁻¹ := by
    filter_upwards [hopen.mem_nhds hs] with z hz
    exact padicTate_shell_decomp hp hK hz
  have hfin : HasDerivAt (fun z => ∑ k ∈ Finset.range K, (a k - c) * rp p z ^ k)
      (∑ k ∈ Finset.range K, (a k - c) * (-((k : ℂ) * (Real.log p : ℂ)) * rp p s ^ k)) s :=
    HasDerivAt.fun_sum (fun k _ => (hasDerivAt_rp_pow hp s k).const_mul (a k - c))
  have htail := (((hasDerivAt_rp hp s).const_sub 1).inv (one_sub_rp_ne_zero hp hs)).const_mul c
  have hsum := hfin.add htail
  rw [padicTate_localScaling_shell_decomp hp hK hs]
  have heq : ∑ k ∈ Finset.range K, -((k : ℂ) * (Real.log p : ℂ)) * (a k - c) * rp p s ^ k +
        c * (-(Real.log p : ℂ) * (rp p s / (1 - rp p s) ^ 2)) =
      ∑ k ∈ Finset.range K, (a k - c) * (-((k : ℂ) * (Real.log p : ℂ)) * rp p s ^ k) +
        c * (- -(-(Real.log p : ℂ) * rp p s) / (1 - rp p s) ^ 2) := by
    congr 1
    · exact Finset.sum_congr rfl (fun k _ => by ring)
    · ring
  rw [heq]
  exact hsum.congr_of_eventuallyEq hev

theorem tate_scaling_is_deriv_shell {p : ℕ} (hp : 1 < p) {a : ℕ → ℂ} (ha : IsShellSchwartz a)
    {s : ℂ} (hs : 0 < s.re) :
    padicTate p (localScaling p a) s = deriv (padicTate p a) s :=
  (tate_scaling_hasDerivAt_shell hp ha hs).deriv.symm

/-! ## Correction 1 guard: the additive measure gives the wrong number -/

/-- The local Tate integral against **additive** Haar measure (`vol(ℤ_p) = 1`):
`dx = (1 − p^{-1}) |x|_p d^×x`, so `Z_add(f, s) = (1 − p^{-1}) Z(f, s + 1)`. -/
noncomputable def padicTateAdd (p : ℕ) (a : ℕ → ℂ) (s : ℂ) : ℂ :=
  (1 - (p : ℂ)⁻¹) * padicTate p a (s + 1)

/-- Normalised scaling trace for `1_{ℤ_p}`: `−log p · r/(1 − r)`. -/
theorem scaling_ratio_charZp {p : ℕ} (hp : 1 < p) {s : ℂ} (hs : 0 < s.re) :
    padicTate p (localScaling p (charZp)) s / padicTate p (charZp) s =
      -(Real.log p : ℂ) * (rp p s / (1 - rp p s)) := by
  rw [padicTate_localScaling_charZp hp hs, padicTate_charZp hp hs]
  have h := one_sub_rp_ne_zero hp hs
  field_simp

theorem rp_three_halves {p : ℕ} : rp p (3 / 2) = ((qp p ^ 3 : ℝ) : ℂ) := by
  have h : rp p (3 / 2) = rp p (1 / 2) ^ 3 := by
    unfold rp
    rw [← Complex.cpow_nat_mul]
    congr 1
    push_cast
    ring
  rw [h, rp_half]
  push_cast
  ring

/-- **Correction 1, PROVED.** Under additive Haar measure the normalised scaling trace at `1/2`
is `−log p · p^{-3/2}/(1 − p^{-3/2}) = −log p/(p^{3/2} − 1)`, and this is **not** `−κ_p`. -/
theorem additive_ratio_ne_kappa {p : ℕ} (hp : 1 < p) :
    padicTateAdd p (localScaling p (charZp)) (1 / 2) / padicTateAdd p (charZp) (1 / 2) =
      -(((Real.log p * (qp p ^ 3 / (1 - qp p ^ 3))) : ℝ) : ℂ) ∧
    padicTateAdd p (localScaling p (charZp)) (1 / 2) / padicTateAdd p (charZp) (1 / 2) ≠
      -((kappa p : ℝ) : ℂ) := by
  have hp1 : (1 : ℂ) - (p : ℂ)⁻¹ ≠ 0 := by
    intro h
    have : (p : ℂ) = 1 := by
      have h' : (p : ℂ)⁻¹ = 1 := (sub_eq_zero.mp h).symm
      simpa using congrArg (·⁻¹) h'
    have : p = 1 := by exact_mod_cast this
    omega
  have hval : padicTateAdd p (localScaling p (charZp)) (1 / 2) / padicTateAdd p (charZp) (1 / 2) =
      -(((Real.log p * (qp p ^ 3 / (1 - qp p ^ 3))) : ℝ) : ℂ) := by
    unfold padicTateAdd
    rw [mul_div_mul_left _ _ hp1, show (1 / 2 : ℂ) + 1 = 3 / 2 by norm_num,
      scaling_ratio_charZp hp (by norm_num), rp_three_halves]
    push_cast
    ring
  refine ⟨hval, ?_⟩
  rw [hval]
  intro h
  have hr : Real.log p * (qp p ^ 3 / (1 - qp p ^ 3)) = kappa p := by
    exact_mod_cast neg_injective h
  have hlog : 0 < Real.log p := Real.log_pos (by exact_mod_cast hp)
  have hq0 := qp_pos hp
  have hsq := one_lt_sqrt hp
  have hq1 : qp p < 1 := by unfold qp; exact inv_lt_one_of_one_lt₀ hsq
  have hk : kappa p = Real.log p * (qp p / (1 - qp p)) := by
    unfold kappa qp
    have h1 : Real.sqrt p - 1 ≠ 0 := by linarith
    have h0 : Real.sqrt p ≠ 0 := by linarith
    have h2 : 1 - (Real.sqrt p)⁻¹ ≠ 0 := by
      rw [sub_ne_zero]; intro h; rw [eq_comm, inv_eq_one] at h; linarith
    field_simp
  rw [hk] at hr
  have hr' : qp p ^ 3 / (1 - qp p ^ 3) = qp p / (1 - qp p) := mul_left_cancel₀ hlog.ne' hr
  have hq3 : qp p ^ 3 < 1 := pow_lt_one₀ hq0.le hq1 (by norm_num)
  rw [div_eq_div_iff (by linarith) (by linarith)] at hr'
  have : qp p ^ 3 = qp p := by linarith
  have : qp p ^ 2 = 1 := by
    have h := this
    have : qp p * (qp p ^ 2 - 1) = 0 := by ring_nf; linarith
    rcases mul_eq_zero.mp this with h0 | h0
    · linarith
    · linarith
  nlinarith

/-! ## Rewiring the (BT-Σ) prime block to `primeBlockJump` -/

/-- The prime-block coupling **derived** from the local Tate action: minus the normalised
`primeBlockJump`. Its definition contains no `κ_p`. -/
noncomputable def derivedCoupling (p : ℕ) : ℝ :=
  -(primeBlockJump p (charZp) / padicTate p (charZp) (1 / 2)).re

/-- Bridge to the banked value: `derivedCoupling p = κ_p`. PROVED. -/
theorem derivedCoupling_eq_kappa {p : ℕ} (hp : 1 < p) : derivedCoupling p = kappa p := by
  unfold derivedCoupling
  rw [local_scaling_trace_critical hp]
  simp

/-- The rewired cut jump: the cross-term drop at the derived coupling. No `κ_p` in the
definition. -/
noncomputable def derivedCutWeight (c : ℕ × ℕ) (g : ℝ → ℝ) (W : ℝ) : ℝ :=
  crossTerm (derivedCoupling c.1) c.1 (g (W - c.2 * Real.log c.1)) (c.2 - 1) -
    crossTerm (derivedCoupling c.1) c.1 (g (W - c.2 * Real.log c.1)) c.2

theorem prime_of_mem_cutIndex {S : Finset ℕ} {N : ℕ} {c : ℕ × ℕ} (hc : c ∈ cutIndex S N) :
    c.1.Prime := by
  simp only [cutIndex, Finset.mem_filter, Finset.mem_product, Finset.mem_Icc] at hc
  exact hc.1.1.2

theorem derivedCutWeight_eq_cutWeight {S : Finset ℕ} {N : ℕ} {c : ℕ × ℕ}
    (hc : c ∈ cutIndex S N) (g : ℝ → ℝ) (W : ℝ) :
    derivedCutWeight c g W = cutWeight c g W := by
  unfold derivedCutWeight cutWeight
  rw [derivedCoupling_eq_kappa (prime_of_mem_cutIndex hc).one_lt]

/-- The derived jump at each cut is the BT3 weight `Λ(p^k)(p^k)^{-1/2} g(W − log p^k)`. -/
theorem derivedCutWeight_eq {S : Finset ℕ} {N : ℕ} {c : ℕ × ℕ} (hc : c ∈ cutIndex S N)
    (g : ℝ → ℝ) (W : ℝ) :
    derivedCutWeight c g W = halfWeight (c.1 ^ c.2) * g (W - Real.log ((c.1 ^ c.2 : ℕ) : ℝ)) :=
  (derivedCutWeight_eq_cutWeight hc g W).trans (cutWeight_eq hc g W)

theorem sum_derivedCutWeight (S : Finset ℕ) (N : ℕ) (g : ℝ → ℝ) (W : ℝ) :
    ∑ c ∈ cutIndex S N, derivedCutWeight c g W = stratumSN S N g W :=
  Finset.sum_congr rfl fun _ hc => derivedCutWeight_eq hc g W

/-- **(J) finite-cut BT3 with the rewired prime block.** If each trace jump is the derived
jump, the Green form is `I ·` the `(S, N)` stratum. PROVED. -/
theorem finiteCut_BT3_derived (S : Finset ℕ) (N : ℕ) (g : ℝ → ℝ) (W : ℝ) (n : ℕ)
    (σ : Fin n ≃ (cutIndex S N)) (ψ : ℕ → ℝ → ℂ) (x : ℕ → ℝ)
    (hψ : ∀ j, ContDiff ℝ 1 (ψ j)) (h0 : ψ 0 (x 0) = 0) (hn : ψ n (x (n + 1)) = 0)
    (hjump : ∀ j : Fin n,
      inner ℂ (ψ j (x (j + 1))) (ψ j (x (j + 1))) -
        inner ℂ (ψ (j + 1) (x (j + 1))) (ψ (j + 1) (x (j + 1))) =
      ((derivedCutWeight (σ j) g W : ℝ) : ℂ)) :
    ∑ j ∈ Finset.range (n + 1), (∫ t in x j..x (j + 1), greenIntegrand (ψ j) t) =
      I * (stratumSN S N g W : ℂ) :=
  finiteCut_BT3 S N g W n σ ψ x hψ h0 hn
    (fun j => by rw [hjump j, derivedCutWeight_eq_cutWeight (σ j).2])

/-- **End to end.** An explicit piecewise `ψ` whose trace jumps are exactly the derived jumps
(`derivedCoupling`, from `primeBlockJump`), and finite-cut BT3 holds for it. PROVED. -/
theorem finiteCut_BT3_derived_witness (S : Finset ℕ) (N : ℕ) (g : ℝ → ℝ) (W : ℝ) (n : ℕ)
    (σ : Fin n ≃ (cutIndex S N)) (x : ℕ → ℝ) (h01 : x 0 ≠ x 1) (hnn : x n ≠ x (n + 1)) :
    ∃ ψ : ℕ → ℝ → ℂ, (∀ j, ContDiff ℝ 1 (ψ j)) ∧ ψ 0 (x 0) = 0 ∧ ψ n (x (n + 1)) = 0 ∧
      (∀ j : Fin n,
        inner ℂ (ψ j (x (j + 1))) (ψ j (x (j + 1))) -
          inner ℂ (ψ (j + 1) (x (j + 1))) (ψ (j + 1) (x (j + 1))) =
        ((derivedCutWeight (σ j) g W : ℝ) : ℂ)) ∧
      ∑ j ∈ Finset.range (n + 1), (∫ t in x j..x (j + 1), greenIntegrand (ψ j) t) =
        I * (stratumSN S N g W : ℂ) := by
  let w : ℕ → ℝ := fun j => if h : j < n then derivedCutWeight (σ ⟨j, h⟩) g W else 0
  obtain ⟨hψ, h0, hn, hj⟩ := hjump_witness n w x h01 hnn
  have hJ : ∀ j : Fin n,
      inner ℂ (psiW n w x j (x (j + 1))) (psiW n w x j (x (j + 1))) -
        inner ℂ (psiW n w x (j + 1) (x (j + 1))) (psiW n w x (j + 1) (x (j + 1))) =
      ((derivedCutWeight (σ j) g W : ℝ) : ℂ) := by
    intro j
    rw [hj j j.isLt]
    simp [w, j.isLt]
  exact ⟨psiW n w x, hψ, h0, hn, hJ, finiteCut_BT3_derived S N g W n σ _ x hψ h0 hn hJ⟩

end GppBTSigmaPLog

#print axioms GppBTSigmaPLog.padicTate_charZp_eq_localZetaPrime
#print axioms GppBTSigmaPLog.tate_scaling_is_deriv_charZp
#print axioms GppBTSigmaPLog.tate_scaling_is_deriv_shell
#print axioms GppBTSigmaPLog.localScaling_charZp_not_schwartz
#print axioms GppBTSigmaPLog.padicTate_charZp_half
#print axioms GppBTSigmaPLog.padicTate_charZp_half_ne_zero
#print axioms GppBTSigmaPLog.local_scaling_trace_critical
#print axioms GppBTSigmaPLog.primeBlockJump_unramified_recovers_kappa
#print axioms GppBTSigmaPLog.additive_ratio_ne_kappa
#print axioms GppBTSigmaPLog.derivedCoupling_eq_kappa
#print axioms GppBTSigmaPLog.derivedCutWeight_eq
#print axioms GppBTSigmaPLog.sum_derivedCutWeight
#print axioms GppBTSigmaPLog.finiteCut_BT3_derived
#print axioms GppBTSigmaPLog.finiteCut_BT3_derived_witness
