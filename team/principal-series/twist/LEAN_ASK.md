# Lean ask: `twist_cos_sum` and `twist_cosh_krein` (Principal series seat; the lead will route it to Claude)

Daniel Toupin / Golden Physics Project, with AI assistance from Grok. **RH is OPEN.** These are elementary kernel identities. They are not RH-strength, and nothing here implies RH.

**Target.** GPPVerify, new file `GppVerify/RiemannHypothesis/WeilTwist.lean`, namespace `GppWeilTwist`. Checked against GPPVerify main e1c1aa9.

## Why
Session 1 of `work/grok/team/principal-series.md` (bridge) found the following:
- Λ_G = Λ(s+iθ)Λ(s−iθ) has window Weil kernel 2cos(θr)·W_ζ(r).
- Λ_F = Λ(s−θ)Λ(s+θ) has kernel 2cosh(θr)·W_ζ(r).

Every downstream claim then rests on two algebraic identities:
- G_θ inherits ζ's certificates.
- F_θ has a Kreĭn J-form, and its onset law is λ_min = θ²κ.

The identities hold at the level of the **autocorrelation** g = f ⋆ f̃, *before* any Weil functional is applied. So the analytic content is only "a linear functional applied to equal functions gives equal values". The explicit-formula step (that the Weil functional of Λ_G is W_ζ applied to 2cos(θu)·g) stays labelled PROVED-in-literature (Weil 1952; Burnol math/9809119). It is **not** part of this ask and should enter as a hypothesis or definition.

## A. Finite-atom version (pure algebra; reuse `GppNoGhost.weilPairing` / `weilForm`, FiniteNoGhost.lean)
`weilForm P x c = Σ_i Σ_j conj(c i) · c j · P(x i − x j)`, with `G = ℝ`.

```lean
namespace GppWeilTwist
open GppNoGhost Real

/-- Real multiplier applied to coefficients: (m • c) i = m (x i) * c i. -/
noncomputable def twistCoef {n : ℕ} (m : ℝ → ℝ) (x : Fin n → ℝ) (c : Fin n → ℂ) : Fin n → ℂ :=
  fun i => (m (x i) : ℂ) * c i

theorem twist_cos_sum (P : ℝ → ℝ) (θ : ℝ) {n : ℕ} (x : Fin n → ℝ) (c : Fin n → ℂ) :
    weilForm (fun r => 2 * cos (θ * r) * P r) x c =
      2 * weilForm P x (twistCoef (fun t => cos (θ * t)) x c) +
      2 * weilForm P x (twistCoef (fun t => sin (θ * t)) x c)

theorem twist_cosh_krein (P : ℝ → ℝ) (θ : ℝ) {n : ℕ} (x : Fin n → ℝ) (c : Fin n → ℂ) :
    weilForm (fun r => 2 * cosh (θ * r) * P r) x c =
      2 * weilForm P x (twistCoef (fun t => cosh (θ * t)) x c) -
      2 * weilForm P x (twistCoef (fun t => sinh (θ * t)) x c)
```

**Hypotheses.** None: P is an arbitrary real function and positivity is not used.

**Proof sketch.**
- Unfold `weilForm` and `weilPairing`. Push the sums together with `Finset.sum_add_distrib` / `Finset.mul_sum`.
- Termwise, use `cos_sub`: cos(θ(a−b)) = cos θa cos θb + sin θa sin θb (here r = x i − x j).
- For the second identity, use `cosh_sub`: cosh(θ(a−b)) = cosh θa cosh θb − sinh θa sinh θb.
- The multipliers are real, so `conj ((m a : ℂ) * c i) = (m a : ℂ) * conj (c i)` (`Complex.conj_ofReal`, `map_mul`).
- Close with `push_cast; ring` per term (`Finset.sum_congr rfl`).

**Corollaries (nice to have, one line each from `weilForm_nonneg`).**
- `twist_cos_nonneg`: if `PositiveTypeOn P`, then `0 ≤ weilForm (fun r => 2 * cos (θ * r) * P r) x c`. This is the finite shadow of "G_θ' inherits".
- `twist_cosh_krein_indefinite` (optional): P r = 1 + cos r (positive type), θ = 10, x = ![0, 1/2], c = ![1, -1]. The cosh-twisted form is 8 − 4cosh(5)(1 + cos(1/2)) < 0, since cosh 5 > 74 and 1 + cos(1/2) > 1.8. This is checked numerically: the 2×2 eigenvalues are ≈ −275, +283. It shows that the cosh twist destroys positive type, i.e. the J-form is genuinely indefinite.

## B. Continuous autocorrelation version (reuse the pattern of OddWindowAutocorrelation.lean / BpyDoubledKmsWedge.`autocorr_even`)
```lean
/-- Autocorrelation A_f(u) = ∫ f(y) · conj (f (y - u)) dy. -/
noncomputable def acorr (f : ℝ → ℂ) (u : ℝ) : ℂ := ∫ y, f y * starRingEnd ℂ (f (y - u))

theorem acorr_twist_cos (f : ℝ → ℂ) (hf : Integrable (fun y => ‖f y‖ ^ 2)) -- or: f bounded measurable, support ⊆ [-L, L]
    (θ u : ℝ) :
    (cos (θ * u) : ℂ) * acorr f u =
      acorr (fun y => (cos (θ * y) : ℂ) * f y) u + acorr (fun y => (sin (θ * y) : ℂ) * f y) u

theorem acorr_twist_cosh (f : ℝ → ℂ) (L : ℝ) (hsupp : ∀ y, L < |y| → f y = 0)
    (hf : Measurable f) (hb : ∃ M, ∀ y, ‖f y‖ ≤ M) (θ u : ℝ) :
    (cosh (θ * u) : ℂ) * acorr f u =
      acorr (fun y => (cosh (θ * y) : ℂ) * f y) u - acorr (fun y => (sinh (θ * y) : ℂ) * f y) u
```

Then, for **any** ℂ-linear functional `W : (ℝ → ℂ) →ₗ[ℂ] ℂ` (the Weil distribution on the window), with
`QW f := W (acorr f)` and `QW_twist k f := W (fun u => 2 * k u * acorr f u)`:
- `twist_cos_sum'`: `QW_twist (cos θ·) f = 2·QW (f cos θ·) + 2·QW (f sin θ·)`.
- `twist_cosh_krein'`: `QW_twist (cosh θ·) f = 2·QW (f cosh θ·) − 2·QW (f sinh θ·)`.

Both follow from `funext` with the pointwise lemmas and `map_add` / `map_sub` / `map_smul`.

**Hypotheses.**
- For the cos case, an L² integrability condition is enough: each integrand is ‖·‖-dominated by |f(y)||f(y−u)|.
- For the cosh case, use compact support plus boundedness, because cosh and sinh are unbounded. That matches the window setting.

**Proof sketch (pointwise).**
- Expand cos(θu) = cos(θy)cos(θ(y−u)) + sin(θy)sin(θ(y−u)), and the cosh analogue.
- Rewrite the integrand of the RHS sum as the LHS integrand (`integral_add` / `integral_sub`, using integrability of each piece), then pull out the constant (`integral_mul_left`).
- The pattern is already in `OddWindowAutocorrelation.autocorr` (sinh product-to-sum, then integrate).

**Parity lemma (for the class-transfer statement of notes §2; optional C).**
- If W is even (`W (g ∘ neg) = W g`) and real on real-even functions, then QW(f) = QW(f_e) + QW(f_o), where f_e and f_o are the even and odd parts.
- Consequently, the bounds `QW ≥ c_e‖·‖²` on even and `QW ≥ c_o‖·‖²` on odd functions give `QW_twist (cos θ·) f ≥ 2·min(c_e, c_o)·‖f‖²`.
- Hypothesis: closure of the certified class under multiplication by cos θ· and sin θ· (smooth, support-preserving). In Lean, state that closure as a hypothesis on the class predicate.

## Names to reuse / not to duplicate
- `GppNoGhost.weilPairing`, `weilForm`, `weilForm_nonneg`, `PositiveTypeOn` (FiniteNoGhost.lean, TruncatedTransport.lean).
- `GppOddWindowAutocorrelation.sinh_mul_sinh` and the `autocorr` proof pattern.
- BpyDoubledKmsWedge `autocorr_even` (translation invariance).
- WeilSupportLadder `convolution_hasSupportIn` (support of f ⋆ f̃ ⊆ [−2L, 2L]).

**Not in scope:** the explicit formula and the identification of Λ_G's Weil distribution (`open_weil_explicit_formula` remains a `True := trivial` placeholder in SpectralWeil.lean). The ask gives the algebra only.
