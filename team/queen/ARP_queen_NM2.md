# ARP_queen_NM2 — Part VI (2026-10-28): the zero-density route to (NM₂), made explicit; is β < 2 reachable unconditionally?

Daniel Toupin / Golden Physics Project, with AI assistance from Grok. **RH is OPEN; (T3) is OPEN.** Analytic only, no numerics. **No RH is used anywhere.** Following the lead, opinions are labelled REMARK. This is the continuation of ARP_queen_Hkappa.md Part V (H.17), which is audited SOUND (main dd3ca88, fixes 7447660). Part V's text is untouched; the lead's fixes (i)–(vi) are kept. (MB)(c) and (WK) belong to Spectrum and are not done here; only their loss *form* is used.

**Headline.**
- **Verdict: β < 2 is NOT reachable unconditionally by the zero-density route with currently citable inputs (PROVED for the route as structured, H.19.5).**
  - For every window length a > 0, the exponent Φ(a) of the (IM₂) bound satisfies Φ(a) ≥ (52/21)(1 − θ_d/(8 + θ_d)).
  - With Conrey's θ_d = 8/7 (any value below it) this is ≥ **13/6 > 2**. With Selberg's θ_d = 1/4 it is ≥ 2.40.
  - The Kovrijkine loss 2κ_K a only adds to Φ(a), so it is **not** the bottleneck.
  - The bottleneck is the worst-case floor on the *exceptional* windows (those near off-line zeros). In the current argument it inherits the subconvexity exponent 13/84 through log⁺|ζ|.
- **What would suffice (named, with status; H.19.6).** One explicit inequality decides the route:
  2κ_K a + 8c₂(a)·4A(a)/(4A(a) + θ_d/2) < 2.
  - Here θ_d is the near-line density exponent, c₂ the floor exponent on exceptional windows, A the window-penalty constant, and κ_K the Kovrijkine exponent (all defined in H.19.1–H.19.4).
  - Since A ≥ 1, a *necessary* condition is c₂ < 1/4 + θ_d/32 ≤ 2/7 (with Conrey).
  - The improvements needed:
    - (I1) a graded treatment of log⁺|ζ| on exceptional windows, removing the 13/84 (unconditional; plausibly doable; next step);
    - (I2) an explicit, small constant C₁ in Titchmarsh 9.6(B) / H.2′, with C₁ log(e/a) < 1/7 − α;
    - (I3) a small Kovrijkine exponent;
    - optionally (I4) a larger θ_d. Density-hypothesis strength near ½ (θ_d = 2) is a CONJECTURE; it is implied by Lindelöf, so it is not circular, but it is not unconditional.
  - **The route is OPEN pending (I1)–(I3).** With present constants it fails.
- **Where the Euler product enters (flagged):**
  - [EP1] the near-line zero-density (mollifier mean values);
  - [EP2] the Dirichlet-polynomial expansion of log ζ off exceptional windows (Selberg's lemma, Λ(n));
  - [EP3] (optional) Ingham's fourth moment for (I1).
  Everything else is archimedean or functional-equation-only, and therefore fake-blind.
- **Falsifiers (H.19.7).** DH and Epstein violate [EP1]: they have ≫ T zeros with β ≥ σ₀ > ½. So the route correctly fails for them; the input is arithmetic but RH-free. RH would make the exceptional set empty and give (NM₂) with exponent o(1) (HEURISTIC sketch). That route is **BARRED**.

## H.19.0 Inputs (all RH-free; cited exactly; bibliographic data web-checked where marked)
- **[Z1] Riemann–von Mangoldt, explicit.** |N(T) − (T/2π)log(T/2πe) − 7/8| ≤ 0.11 log T + 0.29 log log T + 2.29 + 1/(5T) for T ≥ e. Platt–Trudgian 2015, Cor. 1; Trudgian 2014. Quoted via Simonič; web-checked.
- **[Z2] Local zero count.** N(t+1) − N(t−1) ≪ log t (Titchmarsh Thm 9.2). By [Z1]: ≤ (1/π)log(t/2π) + 0.44 log t + O(log log t). Uses only FE + Stirling, so it is fake-blind.
- **[Z3] Local factorization.** For −1 ≤ σ ≤ 2: log ζ(s) = Σ_{|t−γ|<1} log(s − ρ) + O(log t) (Titchmarsh Thm 9.6(B)). Call the constant in O(log t) **C₁**; it is untracked in the literature I have. H.2′ is derived from this, so it is FE + growth only and fake-blind.
- **[Z4] Near-line zero density.** N(σ, T) ≪ T^{1−θ_d(σ−½)} log T, uniformly for ½ ≤ σ ≤ 1:
  - θ_d = 1/4: Selberg 1946; Titchmarsh Thm 9.19(C). Explicit version: Simonič, J. Math. Anal. Appl. (2020), arXiv 1910.08274 (expdb records the constant 5874.051 for N(σ, 2T) − N(σ, T) with ½ ≤ σ ≤ 0.778).
  - any θ_d < 1: Jutila 1983.
  - any θ_d < 8/7: Conrey 1989 (claimed), full proof by Baluyot.
  - θ_d = 2 near ½ is density-hypothesis strength.
  - Normalisation (web-checked): expdb and Baluyot write T^{1+θ(1−2σ)} with θ < 4/7, so θ_d = 2θ.
  - **[EP1]:** every proof bounds ∫|ζM|² with a mollifier M(s) = Σ_{n≤T^{θ}} μ(n)P(log n)n^{−s}, i.e. the Euler product.
- **[Z5] Positive proportion on the line.** Levinson 1/3; Conrey > 2/5; later improvements. **REMARK:** this does not enter the exponents. Floors are hurt by *all* nearby zeros, on-line or not, and the only gain comes from *excluding* off-line zeros, which is [Z4]. The lead's fix (ii) applies only to the H.17.2 lower-tail sketch.
- **[Z6] Subconvexity.** ζ(½+it) ≪ t^{13/84+ε} (Bourgain 2017). FE + exponential sums; fake-blind in form.
- **[K] Kovrijkine.** For f ∈ L²(ℝ) with Fourier support in an interval of length ℓ, and E (γ, a)-thick: ‖f‖ ≤ (C_K/γ)^{C_K(aℓ+1)}‖f‖_{L²(E)} (Proc. AMS 2001). An explicit C_K is given by Egidi–Veselić (2018) (REMARK: constant not re-verified).

## H.19.1 Precise statement of (NM₂) and the exponent bookkeeping (PROVED, bookkeeping)

**Notation.**
- Heights s ∈ [−2πΛ², 2πΛ²], T := 2πΛ². Windows W of length a tile this range.
- Window floor: η_W := sup{η : |{t ∈ W : |ζ(½+it)| ≥ η}| ≥ |W|/2}, and ω := η_W² on W.

**(NM₂)(ν):** (1/T)∫_{|s|≤T} η_{W(s)}^{−4} ds ≤ Λ^{ν}·(log Λ)^{O(1)}.

**(WK₂) loss form (Spectrum).** ‖Gh‖ ≥ c_W(a)‖P_B M_ω h‖ with c_W(a)^{−2} ≤ Λ^{2κ_K a}.
- Reason: Kovrijkine with γ = ½ and ℓ = 2 log Λ + O(1) (H.1) gives the factor K_W = (2C_K)^{C_K(2a log Λ + O(a) + 1)} = Λ^{2C_K log(2C_K)·a}·O_a(1).
- Set κ_K := 4C_K log(2C_K), so K_W² = Λ^{κ_K a}. (WK₂) loses K_W² once per side, hence c_W^{−2} ≤ K_W⁴ = Λ^{2κ_K a}.
- If (WK₂) is proved with a different power of K_W, κ_K rescales; nothing below changes qualitatively.

**Bookkeeping.**
- By H.17.3, (IM₂) ⇐ (WK₂) + (NM₂). Its loss is c_W^{−2}·Λ^{ν}·r_K·(log Λ)^{O(1)}, which must be ≤ Λ^{2−δ}.
- So the closing condition is **Φ(a) := 2κ_K a + ν(a) < 2**.
- (A_β) then holds with β = Φ/2 (H.17.1, Cauchy–Schwarz).
- **Unconditional worst case (Part V).** η_W^{−4} ≤ T^{4c₂} = Λ^{8c₂+o(1)}, with c₂ = 2(θ_sub + ε + C₁ log(e/a)) and θ_sub = 1/6, so ν ≥ 8/3. With [Z6], θ_sub = 13/84 and ν_worst ≥ 52/21 ≈ 2.476.

## H.19.2 Off-line-free windows: the floor at level h (PROVED (sketch); constants partly explicit)

Fix h ∈ (0, 1] and σ₁ := ½ + h. Call W **h-clean** if no zero ρ = β + iγ with β − ½ ≥ h/2 has |γ − t| ≤ 2 for some t ∈ W.

**Lemma H.19.2.** On an h-clean window, η_W ≥ T^{−A h}·e^{−O_h(1)}, with A := A_Z + A_E and:
- **A_Z = (1 + 2/a)/(2f)** (explicit). Here f ∈ (0, ½) is the fraction of W conceded to the zero penalty; the remaining ½ − f is conceded to the [EP2] and log⁺ terms.
- **A_E**: the constant in the O((σ₁ − ½)log t) error of [Z3] differentiated (Titchmarsh 9.6(A)), plus Selberg's-lemma error. Untracked.

*Proof sketch.* Integrate Re ζ′/ζ along [½, σ₁] using 9.6(A):

  log|ζ(½+it)| = log|ζ(σ₁+it)| − Z_h(t) + O(A_E h log t),  where Z_h(t) := Σ_{|t−γ|<1} log(|σ₁+it−ρ|/|½+it−ρ|).

1. **Zero penalty.** For zeros with |β − ½| < h/2, the summand is ≥ −log 2 on the clean window, and ≤ ½log(1 + 4h²/(t−γ)²) + O(h). So
   ∫_W Z_h ≤ n_W ∫_ℝ ½log(1 + h²/x²)dx + O(n_W h) = n_W·πh + O(n_W h),
   with n_W ≤ (a + 2)(log T)/(2π)(1 + o(1)) by [Z1]–[Z2].
   Chebyshev: |{t ∈ W : Z_h > V}| ≤ fa once V = (a + 2)h log T/(2fa)·(1 + o(1)). This is A_Z.
2. **The value at σ₁ [EP2].** There is no zero with β > ½ + h/2 within distance 2. Selberg's unconditional lemma (Selberg 1946, with his σ_{x,t}; Tsang's thesis), with x = e^{c/h}, gives log|ζ(σ₁+it)| = Re P_x(t) + O(A_E h log T) + O(|Q_x(t)|/log x). Here P_x and Q_x are Dirichlet polynomials of length x³ = e^{3c/h}, which is bounded for fixed h. So |P_x|, |Q_x| = O_h(1) pointwise.
   - **This step uses Λ(n), i.e. the Euler product [EP2].**
3. **Combine.** On ≥ ½ of W, log|ζ(½+it)| ≥ −(A_Z + A_E)h log T − O_h(1). ∎

**REMARK (why A_Z ≥ 1 is not an artifact).** For h ≫ 1/log T the penalty Z_h is essentially deterministic: Z_h ≈ (log T/2π)·πh = h log T/2. In the exact Hadamard form this mean is cancelled by +(h/2)log(t/2π), and what remains is the S(t)-fluctuation. But the local form [Z3] does not see the cancellation, and the worst-case S(t) = O(log t) (Titchmarsh 9.4) restores a loss of the same order. So with the local form, **A ≥ (1+2/a)/(2f) > 1** is genuine. Exploiting the Hadamard cancellation in the worst case is item (I2′) in H.19.6.

## H.19.3 Measure of non-clean windows [EP1] (PROVED (sketch))

A window fails to be h-clean only if a zero with β − ½ ≥ h/2 has γ within distance 2 + a. By [Z4], the fraction of such windows among those in [T/2, T] is

  ≤ (a + 4)/(T/2)·N(½ + h/2, T) ≪ (a + 4) T^{−θ_d h/2} log T.

Lower heights contribute at most the same after dyadic summation. **[EP1].**

## H.19.4 The window-floor negative moment: two strategies (PROVED (sketch))

Write T^{x} for exponents; ν = 2·(T-exponent), since T ≍ Λ².

1. **Graded.** Give each window its level h(W) := min(1, 2·max{β − ½ : ρ within distance 2 of W}), or h₀ → 0⁺ if there is none. Then η_W^{−4} ≤ T^{4A h(W)}, and the fraction with h(W) ≥ h is ≪ T^{−θ_d h/2}log T. So
   mean η^{−4} ≪ log²T·Σ_{h dyadic ≤ 1} T^{(4A − θ_d/2)h} ⇒ ν_grad = 2·max(0, 4A − θ_d/2) = max(0, 8A − θ_d).
   Since h ≤ 1 (β ≤ 1), no worst-case floor is needed. But **A ≥ 1 forces ν_grad ≥ 8 − 8/7 = 48/7** even with Conrey.
2. **Balanced.** Use level h* on clean windows and the worst case c₂ on the rest. Then
   mean η^{−4} ≪ T^{4Ah*} + T^{4c₂ − θ_d h*/2}log T.
   Optimizing at h* = 4c₂/(4A + θ_d/2) gives
   **ν_bal = 8c₂·4A/(4A + θ_d/2) = 8c₂(1 − θ_d/(8A + θ_d)).**
   - This is a strict improvement on the worst case 8c₂, by the factor 1 − θ_d/(8A + θ_d), but only that factor.

The best of the two, ν(a) = min(ν_grad, ν_bal), is ν_bal in every admissible regime (A ≥ 1, θ_d ≤ 8/7).

## H.19.5 Joint optimisation in a, and the verdict (PROVED for this route)

The parameter a enters in three places:
- A(a) = (1 + 2/a)/(2f) + A_E, which decreases in a and tends to 1/(2f) + A_E > 1;
- c₂(a) = 2(θ_sub + ε + C₁ log(e/a)) for a ≤ 1, and ≥ 2(θ_sub + ε) for every a (the log⁺ term in Part V's Chebyshev step does not shrink with a);
- the Kovrijkine term 2κ_K a, which increases in a.

So:

  **Φ(a) = 2κ_K a + 8c₂(a)·(1 − θ_d/(8A(a) + θ_d)) ≥ 16θ_sub·(1 − θ_d/(8 + θ_d)).**

The inequality uses only κ_K ≥ 0, C₁ ≥ 0, A ≥ 1, ε ≥ 0. With θ_sub = 13/84 [Z6]:

| density input | θ_d | inf_a Φ(a) ≥ | β = Φ/2 ≥ |
|---|---|---|---|
| Selberg 1946 (Titchmarsh 9.19(C)) | 1/4 | (52/21)(32/33) ≈ 2.40 | 1.20 |
| Jutila 1983 | 1⁻ | (52/21)(8/9) ≈ 2.20 | 1.10 |
| Conrey 1989 / Baluyot | 8/7⁻ | (52/21)(7/8) = **13/6 ≈ 2.167** | 1.083 |
| density hypothesis near ½ (CONJECTURE) | 2 | (52/21)(4/5) ≈ 1.98 | 0.99 |

**Verdict (PROVED for the route as structured).** With every currently proved density exponent, inf_a Φ(a) > 2, **whatever the untracked constants C₁, A_E, κ_K are** (they only increase Φ). So **β < 2 is not reachable unconditionally by this route as it stands.**

The obstruction is the worst-case floor on exceptional windows. It inherits 2θ_sub from the log⁺ term, and the density gain multiplies it only by 1 − θ_d/(8A + θ_d) ≥ 7/8. Kovrijkine is not the bottleneck: the arithmetic part alone exceeds 2 at every a.

Even the density-hypothesis row (≈ 1.98) clears 2 only in the idealized limit C₁ = A_E = κ_K = 0 with A → 1 (which needs a → ∞, contradicting κ_K a → 0). REMARK: with honest constants, DH-strength density alone would not close it either.

## H.19.6 What would suffice (named precisely)

The exact closing condition is

  (★) 2κ_K a + 8c₂(a)·4A(a)/(4A(a) + θ_d/2) < 2 for some a > 0.

Since A ≥ 1, a necessary condition is **c₂ < 1/4 + θ_d/32**, i.e. c₂ < 0.258 (Selberg) or c₂ < 2/7 ≈ 0.286 (Conrey). Equivalently, on the exceptional windows the half-window floor must satisfy η_W ≥ T^{−2/7}, and that floor exponent must be the *only* loss.

The improvements, by item:
- **(I1) Graded log⁺ on exceptional windows (unconditional; expected doable; next step).**
  - Split exceptional windows further by mean_W|ζ|² ≤ T^{2α}. Ingham's fourth moment ∫_0^T|ζ|⁴ ≪ T log⁴T [EP3] gives a fraction ≪ T^{−4α}log⁴T for the complement, and the log⁺ term becomes α log T.
  - This removes θ_sub from c₂ at the cost of one more balanced class.
  - Then c₂ → 2(α + C₁ log(e/a)), and the necessary condition becomes C₁ log(e/a) + α < 1/7 (Conrey). **Label: CONJECTURE that it works; not done.**
- **(I2) Explicit C₁ (Titchmarsh 9.6(B)) and A_E (9.6(A) and Selberg's lemma).** Needed: C₁ log(e/a) < 1/7 − α at the chosen a. Explicit versions in the Trudgian/Simonič line are the natural source. **Status: OPEN** (not in the literature I have; REMARK: plausible only for a ≳ 1, where log(e/a) ≈ 1).
- **(I2′) Hadamard cancellation in the worst case.** Replace the local factorization by the exact Hadamard form and bound the residual through windowed S(t). The explicit |S(t)| ≤ 0.11 log t + … [Z1] would make the worst-case constant ≈ 0.22·(kernel mass) instead of C₁. **Status: OPEN; this is the most promising lever** (REMARK).
- **(I3) Kovrijkine exponent.** κ_K = 4C_K log(2C_K) explicit (Egidi–Veselić), or a weighted (WK₂) with better a-dependence. Spectrum's item; the condition only needs 2κ_K a < 2 − ν(a).
- **(I4) Density.** Larger θ_d helps only through the factor 1 − θ_d/(8A + θ_d):
  - θ_d = 2 (density hypothesis near ½): CONJECTURE, implied by Lindelöf (REMARK: Ingham-type LH ⇒ DH). Not circular for RH, but not unconditional.
  - θ_d → ∞ is RH: BARRED.

**Overall status of "β < 2 unconditionally": OPEN.** It fails with present constants (PROVED, H.19.5) and is decided by (★) after (I1)–(I3). REMARK: under Lindelöf (θ_sub = ε), (★) reduces to C₁ log(e/a) + A-terms, which is a pure explicit-constants question. LH is weaker than RH and not circular, but it is a CONJECTURE.

## H.19.7 Euler-product flags and fake behaviour (R2, R5)

| Step | Uses | DH | Epstein (class number > 1) |
|---|---|---|---|
| [Z1]–[Z3], H.2′, [Z6] | FE + Stirling + growth | holds (same Γ-factor and N(T) shape) | holds |
| [K], (WK₂) | archimedean | holds | holds |
| [EP1] near-line density [Z4] | mollifier = Euler product | **fails** (see below) | **fails** |
| [EP2] Selberg's lemma | Λ(n) | fails (no Euler product; log f_DH has no prime expansion) | fails |
| [EP3] Ingham's 4th moment (for I1) | approximate FE + divisor sum | analogue holds up to constants (REMARK) | analogue holds |

- **DH fails [Z4].** f_DH has infinitely many zeros with σ > 1 (Davenport–Heilbronn 1936), and ≫ T zeros in ½ < σ₁ < σ < σ₂ < 1 (Voronin; Bombieri–Ghosh, Russian Math. Surveys 66 (2011); REMARK: not re-verified). At fixed σ > ½ this contradicts T^{1−θ_d(σ−½)}.
  - Near the line, zeros of linear combinations still *cluster* (Selberg's announced results; Bombieri–Hejhal), so the failure is at fixed distance h/2 > 0. That is exactly where H.19.3 uses [Z4].
- **Epstein:** the same (Davenport–Heilbronn 1936 for Epstein; Voronin).
- **Consequence.** The route is **not fake-blind**: it correctly fails for DH and Epstein at [EP1]/[EP2]. The input is genuinely arithmetic.
- **It is not secretly RH.** [Z4] is proved unconditionally, and RH ⇔ θ_d = ∞ is used nowhere. RH would make every window clean (no exceptional set) and give ν = o(1) (HEURISTIC sketch). That route is **BARRED**.
- Lindelöf (I4) is not RH-circular, but it is a CONJECTURE.
- R-line H.19: R1 Λ-uniform | R2 Euler product at [EP1], [EP2], ([EP3]) | R3 n/a | R4 n/a | R5 the DH and Epstein failures recorded; no false positive possible, because the route fails for them | R6 the bookkeeping identities are exact.

## H.20 Summary table

| Claim | Label |
|---|---|
| (NM₂)(ν) stated; closing condition Φ(a) = 2κ_K a + ν(a) < 2; β = Φ/2 | PROVED (bookkeeping) |
| Kovrijkine loss c_W^{−2} = Λ^{2κ_K a}, κ_K = 4C_K log(2C_K) | PROVED-conditional on the (WK₂) loss form (Spectrum) |
| Lemma H.19.2: h-clean floor η_W ≥ T^{−Ah}, A = (1+2/a)/(2f) + A_E | PROVED (sketch) [EP2] |
| A ≥ 1 is genuine in the local form | REMARK (Hadamard cancellation lost; worst S(t)) |
| Non-clean fraction ≪ (a+4)T^{−θ_d h/2} log T | PROVED (sketch) [EP1] |
| ν_grad = max(0, 8A − θ_d) ≥ 48/7; ν_bal = 8c₂(1 − θ_d/(8A+θ_d)) | PROVED (sketch) |
| inf_a Φ ≥ (52/21)(1 − θ_d/(8+θ_d)): ≥ 2.40 (Selberg), 2.20 (Jutila), 13/6 (Conrey) | PROVED (for this route, constants-independent lower bound) |
| **β < 2 unconditionally via zero density: not reachable as it stands** | PROVED (route); overall OPEN pending (I1)–(I3) |
| Necessary: c₂ < 1/4 + θ_d/32 (≤ 2/7) | PROVED |
| (I1) graded log⁺ via Ingham's 4th moment | CONJECTURE (expected doable; next) |
| (I2) explicit C₁, A_E; (I2′) Hadamard + explicit S(t) | OPEN |
| (I4) θ_d = 2 near ½ (⇐ Lindelöf) | CONJECTURE; not circular |
| RH ⇒ (NM₂) with ν = o(1) | HEURISTIC (sketch); **BARRED** |
| DH and Epstein fail [Z4], so the route is arithmetic, RH-free, and not fake-blind | PROVED (sketch; DH/Epstein facts cited, not re-verified) |
| (T3) | **OPEN** |

Not cited: ArithmeticTemperedPrincipalSeries, rh_iff_weil_pairedForm_nonneg, or any True-placeholder.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
