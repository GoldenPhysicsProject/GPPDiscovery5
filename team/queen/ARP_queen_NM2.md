# ARP_queen_NM2 — Part VI (2026-10-28): the zero-density route to (NM₂), made explicit; is β < 2 reachable unconditionally?

Daniel Toupin / Golden Physics Project, with AI assistance from Grok. **RH is OPEN; (T3) is OPEN.** Analytic only, no numerics. **No RH is used anywhere.** Following the lead, opinions are labelled REMARK. This is the continuation of ARP_queen_Hkappa.md Part V (H.17), which is audited SOUND (main dd3ca88, fixes 7447660). Part V's text is untouched; the lead's fixes (i)–(vi) are kept. (MB)(c) and (WK) belong to Spectrum and are not done here; only their loss *form* is used.

**Headline.**
- **Verdict: β < 2 is NOT reachable unconditionally by the zero-density route with currently citable inputs (PROVED for the route as structured, H.19.5).**
  - For every window length a > 0, the exponent Φ(a) of the (IM₂) bound satisfies Φ(a) ≥ (52/21)(1 − θ_d/(8 + θ_d)).
  - With Conrey's θ_d = 8/7 (any value below it) this is ≥ **13/6 > 2**. With Selberg's θ_d = 1/4 it is ≥ 2.40.
  - *(Grok audit S103i(3))* The balanced strategy is not optimal within the route: letting clean windows keep their own graded level (hybrid, H.19.4 item 3) gives the sharper route bound Φ(a) ≥ (52/21)(1 − θ_d/8), i.e. ≥ **104/49 ≈ 2.12 > 2** (Conrey), 13/6 (Jutila), ≈ 2.40 (Selberg). The verdict NO stands for this route; the numbers are corrected. This is a limitation of this method with A ≥ 1 (local factorization) and a sup-norm log⁺ bound, not a statement about (NM₂), (T3) or RH.
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
  - any θ_d < 1: Jutila 1983 (Studies in Pure Math. in memory of P. Turán, Birkhäuser, 385–394).
  - any θ_d < 8/7: Conrey 1989 (announced as Theorem 2 in Bull. AMS 20 (1989) 79–81, no proof published), full proof by Baluyot (PhD thesis, Rochester 2017, Thm 1.2.1, new mollifier + Kloosterman sums).
  - *(Grok audit S103i(3))* web-verified: in the normalisation N(σ,T) ≪ T^{1+θ(1−2σ)} log T, Selberg has θ = 1/8, Jutila θ < 1/2, Conrey/Baluyot θ < 4/7, density hypothesis ⇔ all θ < 1; so θ_d = 2θ = 1/4, 1⁻, 8/7⁻, 2⁻. Explicit Selberg-type version: Simonič, J. Math. Anal. Appl. (2020), article 124303, doi 10.1016/j.jmaa.2020.124303.
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

1. **Zero penalty.** For zeros with |β − ½| < h/2, the summand ½log(((σ₁−β)² + x²)/((½−β)² + x²)), x = t − γ, is ≥ 0 (since |σ₁ − β| ≥ |½ − β| for β ≤ ½ + h/2), and its integral over x ∈ ℝ is exactly π(|σ₁−β| − |½−β|) ≤ πh. *(Grok audit S103i(3))* this replaces the pointwise bound ½log(1 + 4h²/x²), which does not integrate to πh; the conclusion is unchanged. So
   ∫_W Z_h ≤ n_W·πh,
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

3. **Hybrid *(Grok audit S103i(3))* (PROVED (sketch), same inputs).** Charge each window min(T^{4A h(W)}, T^{4c₂}), i.e. its own graded floor below the threshold h* and the worst case above it (floor h(W) at a fixed h₀ > 0 so that Selberg's lemma stays O_{h₀}(1)). Then mean η^{−4} ≪ log²T·max_h T^{−θ_d h/2}·min(T^{4Ah}, T^{4c₂}), maximised at h* = c₂/A:
   **ν_hyb = 8c₂(1 − θ_d/(8A))** for 8A > θ_d (and ν_hyb = 0 if 8A ≤ θ_d, which is ν_grad = 0).
   - ν_hyb < ν_bal always, because the balanced strategy charges *all* clean windows at T^{4Ah*} although only a T^{−θ_d h*/2} fraction sit at level h*.
   - Consequence: the route's lower bound becomes Φ(a) ≥ 16θ_sub(1 − θ_d/8) (using A ≥ 1), the necessary condition becomes c₂ < 2/(8 − θ_d) (7/24 ≈ 0.292 Conrey, 8/31 ≈ 0.258 Selberg), and **if the floor constant can be made A < θ_d/8 the graded strategy gives ν = 0 with no worst-case floor at all** (then c₂ and θ_sub drop out and only 2κ_K a < 2 remains).

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

*(Grok audit S103i(3))* Hybrid (H.19.4 item 3), lower bound 16θ_sub(1 − θ_d/8) = (52/21)(1 − θ_d/8): Selberg (52/21)(31/32) ≈ 2.40; Jutila (52/21)(7/8) = 13/6 ≈ 2.167; **Conrey/Baluyot (52/21)(6/7) = 104/49 ≈ 2.122**; density hypothesis (CONJECTURE) (52/21)(3/4) = 13/7 ≈ 1.857. Arithmetic of the balanced column re-checked: 1664/693 ≈ 2.401, 416/189 ≈ 2.201, 13/6, 208/105 ≈ 1.981 — correct for ν_bal.

**Verdict (PROVED for the route as structured; with the hybrid correction the bound is 104/49 for Conrey).** With every currently proved density exponent, inf_a Φ(a) > 2, **whatever the untracked constants C₁, A_E, κ_K are** (they only increase Φ). So **β < 2 is not reachable unconditionally by this route as it stands.**

REMARK *(Grok audit S103i(3))* on the direction of c₂ ≥ 26/84: θ_sub is an *upper* bound exponent (|ζ| ≪ t^{θ_sub+ε}); the method bounds ∫_W log⁺|ζ| ≤ a(θ_sub + ε)log T and hence gets c₂ = 2(θ_sub + ε + C₁log(e/a)). So "c₂ ≥ 26/84" is a lower bound on the floor exponent **this argument can deliver with the best proved subconvexity** (Bourgain, J. Amer. Math. Soc. 30 (2017) 205–224), not a lower bound on true window floors. A better subconvexity, a mean-value treatment of log⁺ ((I1)), or Lindelöf would lower it.

The obstruction is the worst-case floor on exceptional windows. It inherits 2θ_sub from the log⁺ term, and the density gain multiplies it only by 1 − θ_d/(8A + θ_d) ≥ 7/8. Kovrijkine is not the bottleneck: the arithmetic part alone exceeds 2 at every a.

Even the density-hypothesis row (≈ 1.98) clears 2 only in the idealized limit C₁ = A_E = κ_K = 0 with A → 1 (which needs a → ∞, contradicting κ_K a → 0). REMARK: with honest constants, DH-strength density alone would not close it either.

## H.19.6 What would suffice (named precisely)

The exact closing condition is

  (★) 2κ_K a + 8c₂(a)·4A(a)/(4A(a) + θ_d/2) < 2 for some a > 0.

Since A ≥ 1, a necessary condition **for the balanced strategy** is **c₂ < 1/4 + θ_d/32**, i.e. c₂ < 0.258 (Selberg) or c₂ < 2/7 ≈ 0.286 (Conrey). *(Grok audit S103i(3))* For the hybrid it is c₂ < 2/(8 − θ_d) (7/24 Conrey); both presuppose A ≥ 1, i.e. the local factorization. With A < θ_d/8 (only possible via (I2′)) no condition on c₂ is needed. Equivalently, on the exceptional windows the half-window floor must satisfy η_W ≥ T^{−2/7}, and that floor exponent must be the *only* loss.

The improvements, by item:
- **(I1) Graded log⁺ on exceptional windows (unconditional; expected doable; next step).**
  - Split exceptional windows further by mean_W|ζ|² ≤ T^{2α}. Ingham's fourth moment ∫_0^T|ζ|⁴ ≪ T log⁴T [EP3] gives a fraction ≪ T^{−4α}log⁴T for the complement, and the log⁺ term becomes α log T.
  - This removes θ_sub from c₂ at the cost of one more balanced class.
  - Then c₂ → 2(α + C₁ log(e/a)), and the necessary condition becomes C₁ log(e/a) + α < 1/7 (Conrey). **Label: CONJECTURE that it works; not done.**
- **(I2) Explicit C₁ (Titchmarsh 9.6(B)) and A_E (9.6(A) and Selberg's lemma).** Needed: C₁ log(e/a) < 1/7 − α at the chosen a. Explicit versions in the Trudgian/Simonič line are the natural source. **Label: GAP** (status open) (not in the literature I have; REMARK: plausible only for a ≳ 1, where log(e/a) ≈ 1).
- **(I2′) Hadamard cancellation in the worst case.** Replace the local factorization by the exact Hadamard form and bound the residual through windowed S(t). The explicit |S(t)| ≤ 0.11 log t + … [Z1] would make the worst-case constant ≈ 0.22·(kernel mass) instead of C₁. **Label: GAP; this is the most promising lever** (REMARK). *(Grok audit S103i(3))* (I2′) should be applied to the *clean-window* floor constant A as well, not only the worst case: in the Hadamard form the on-line mean h log T/2 cancels, the residual is ∫S(u)dK_{h,W}(u) with TV(K_{h,W}) ≍ h/a, so the honest A is ≍ (explicit S-constant + A_E)/a rather than ≥ 1 (HEURISTIC; not checked in detail). If A < θ_d/8 (= 1/7 with Conrey) the graded strategy gives ν = 0 (H.19.4 item 3).
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
  - Near the line, zeros of linear combinations still *cluster* (Selberg's announced results; Bombieri–Hejhal), so the failure is at fixed distance h/2 > 0. That is exactly where H.19.3 uses [Z4]. *(Grok audit S103i(3))* REMARK: Bombieri–Hejhal (Duke Math. J. 80 (1995)) is conditional (GRH for the constituents plus a weak pair-correlation hypothesis), and the clustering is not used in the argument; only the ≫ T off-line zeros at fixed σ > ½ matter, and those (Davenport–Heilbronn 1936; Voronin; survey Bombieri–Ghosh, Russian Math. Surveys 66:2 (2011) 221–270) suffice for the failure of [Z4].
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
| inf_a Φ ≥ (52/21)(1 − θ_d/(8+θ_d)) for the balanced strategy: ≥ 2.40 (Selberg), 2.20 (Jutila), 13/6 (Conrey); with the hybrid (52/21)(1 − θ_d/8): 2.40, 13/6, 104/49 | PROVED (for this route with local factorization A ≥ 1 and sup-norm log⁺; constants-independent lower bound) |
| **β < 2 unconditionally via zero density: not reachable as it stands** | PROVED (route); overall OPEN pending (I1)–(I3) |
| Necessary (balanced strategy, A ≥ 1): c₂ < 1/4 + θ_d/32 (≤ 2/7); hybrid: c₂ < 2/(8 − θ_d) | PROVED |
| (I1) graded log⁺ via Ingham's 4th moment | CONJECTURE (expected doable; next) |
| (I2) explicit C₁, A_E; (I2′) Hadamard + explicit S(t) | GAP |
| *(Grok audit S103i(3))* Hybrid strategy ν_hyb = 8c₂(1 − θ_d/(8A)); route bound Φ ≥ (52/21)(1 − θ_d/8) ≥ 104/49 (Conrey); ν = 0 if A < θ_d/8 | PROVED (sketch) |
| (I4) θ_d = 2 near ½ (⇐ Lindelöf) | CONJECTURE; not circular |
| RH ⇒ (NM₂) with ν = o(1) | HEURISTIC (sketch); **BARRED** |
| DH and Epstein fail [Z4], so the route is arithmetic, RH-free, and not fake-blind | PROVED (sketch; DH/Epstein facts cited, not re-verified) |
| (T3) | **OPEN** |


### Grok audit S103i(3) of Part VI (2a06a4f)
**Verdict: SOUND after fixes (marked S103i(3) above). The conclusion "β < 2 is not reachable unconditionally by this zero-density route as it stands" holds, with a corrected, slightly sharper route bound 104/49 ≈ 2.12 (Conrey/Baluyot) from the hybrid strategy.**
- Bookkeeping (Φ = 2κ_K a + ν, β = Φ/2, Kovrijkine K_W² = Λ^{κ_K a}) re-derived; c_W^{−2} ≤ K_W⁴ is an assumed loss form for (WK₂) (PROVED-conditional, as labelled).
- Floor lemma: A_Z = (1 + 2/a)/(2f) re-derived via the exact integral π(|σ₁−β| − |½−β|) ≤ πh per zero and Chebyshev; A ≥ 1 + 2/a > 1 in the local form.
- ν_grad and ν_bal algebra correct; the balanced strategy is suboptimal (hybrid added).
- c₂ ≥ 26/84 is a property of the argument with the best proved subconvexity exponent (an upper-bound exponent), not a lower bound on true floors (REMARK added).
- Density attributions web-verified (Selberg θ = 1/8, Jutila θ < 1/2, Conrey Bull. AMS 1989 announcement, Baluyot 2017 thesis proof, θ < 4/7; DH ⇔ θ < 1).
- DH/Epstein: ≫ T zeros at fixed σ > ½ contradict [Z4]; correct as cited. The clustering remark (Bombieri–Hejhal) is conditional and not used.
- Labels: (I1) CONJECTURE, (I2)/(I2′) GAP (were "OPEN"), (I4) CONJECTURE; RH ⇒ (NM₂) HEURISTIC (sketch), BARRED.
- The verdict is a limitation of this method (local factorization A ≥ 1 plus sup-norm log⁺), not a theorem about (NM₂), (T3) or RH. **(T3) OPEN. RH OPEN.**

Not cited: ArithmeticTemperedPrincipalSeries, rh_iff_weil_pairedForm_nonneg, or any True-placeholder.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok.

---

# Part VII (2026-10-28): (I2′) the exact Hadamard form; the clean-window constant A; a moment route to (NM₂) with ν < 2; the Kovrijkine threshold

Daniel Toupin / Golden Physics Project, with AI assistance from Grok. **RH is OPEN; (T3) is OPEN.** Analytic only, no numerics (the arithmetic below is bookkeeping with cited constants). **No RH is used anywhere; RH-derived inputs are BARRED.** Parts V–VI are kept as merged, with the lead's audit fixes S103i(3) (1)–(4). Labels follow the project list; where Part VI said "OPEN" for (I2)/(I2′), read GAP.

**Headline.**
1. **(I2′) decisive test (PROVED (sketch)).** The exact Hadamard form gives the window *mean* with only a 1/a loss, exactly as the lead's HEURISTIC predicted:
   mean_W log|ζ(½+it)| ≥ −(2πh_b C_S/a)log T − O(1).
   But converting the mean into the half-measure *floor* costs the sup term. On h-clean windows,
   **A_clean = ½ + 3A_E + 4πC_S/a ≥ ½ > θ_d/8 (< 1/7, Conrey/Baluyot).**
   - So the graded test "A < θ_d/8 ⇒ ν = 0" **fails** for this argument.
   - The ½ is the T^{h/2} Stirling gap between σ₁ and ½. It enters through the sup of log⁺ on the good half of W, not through the zeros.
   - REMARK: without the sup term, A would be 4πC_S/a + 2A_E, which is < 1/7 for a > 9.2 and A_E < 1/14 − 2πC_S/a. The obstruction is exactly the sup term.
2. **New route (PROVED (sketch)), no zero density at all.** Combine three pieces:
   - the Hadamard window-mean bound, valid for **every** window and using no cleanliness;
   - Jensen for log⁺;
   - the unconditional eighth moment ∫_0^T|ζ|⁸ ≪ T^{3/2}log^{21/2}T (Hölder from Ingham's 4th moment and Heath-Brown's 12th).

   Together they give **(NM₂) with ν(a) = 1 + 16D/a + ε**, where:
   - D = 2πh_b·C_S;
   - C_S = 0.1038 (Hasanalizade–Shen–Wong 2022) or 0.112 (Trudgian 2014);
   - h_b = 3/2 for the fake-blind version (σ_b = 2), giving D ≈ 0.978;
   - h_b = ½⁺ for the [EP4] version (σ_b = 1⁺), giving D ≈ 0.326.

   So **ν < 2 unconditionally once a > 16D** (a > 5.2 with [EP4], a > 15.6 fake-blind). This is the first unconditional (NM₂) exponent below 2 in the project. It supersedes Part VI's 8/3 worst case and the 104/49 route bound; those remain correct statements about *their* routes.
3. **The decision moves entirely to Kovrijkine (Spectrum).**
   - Φ(a) = 2κ_K a + 1 + 16D/a, and min_a Φ = 1 + 8√(2κ_K D) at a* = √(8D/κ_K).
   - **β < 2 ⇔ κ_K < 1/(128D)**: κ_K < 0.0240 [EP4] or κ_K < 0.0080 (fake-blind). In Kovrijkine's normalisation, C_K log(2C_K) < 1/(512D) ≈ 0.0060.
   - **REMARK/HEURISTIC:** generic Logvinenko–Sereda/Kovrijkine cannot meet this. The Remez/Turán sharpness on a window holding n ≈ aℓ/2π degrees of freedom forces κ_K ≳ (2/π)log 2 ≈ 0.44, which is ~18× too big.
   - Published explicit forms (Egidi–Veselić, Arch. Math. 111 (2018) 85–99; web-checked) carry an unspecified universal K.
   - So β < 2 needs a **structured (WK₂)** whose a-loss beats uniform Kovrijkine by a factor ≳ 20. That is a GAP for Spectrum.
4. (I1), the Ingham-graded log⁺, is **superseded** by item 2. The eighth-moment route is the correct "graded log⁺", giving T-exponent ½ in place of 4θ_sub ≈ 0.62. The capped 13/7 figure is no longer the relevant ceiling.

## H.21.0 Inputs (RH-free; web-checked where marked)

- **[S] Explicit S(t).** |S(T)| ≤ min{0.1038 log T + 0.2573 log log T + 8.3675, 0.1095 log T + 0.2042 log log T + 3.0305} for T ≥ e. Source: Hasanalizade–Shen–Wong, J. Number Theory 235 (2022), arXiv 2107.06506, Cor. 1.4 (web-checked). It improves Trudgian 2014 (C₁ = 0.112, C₂ = 0.278, C₃ = 3.385 in the N(T) form, as tabulated there).
  - S is defined along 2 → 2+iT → ½+iT, which is the standard S.
  - The proof uses Backlund's/Jensen's method plus bounds on the ½- and 1-lines. These are FE + Dirichlet-series *upper* bounds, so no Euler product enters: fake-blind in form (REMARK).
  - Write S*(T) := sup_{|γ|≤2T}|S(γ)| = C_S log T + O(log log T), with C_S = 0.1038.
- **[Θ] Density (for comparison only).** Baluyot's thesis (Rochester 2017, Thm 1.2.1; web-checked) gives N(σ,T) ≪_θ T^{1+θ(1−2σ)}log T for 0 < θ < 4/7. In Part VI's normalisation θ_d = 2θ < 8/7, so **θ_d/8 < 1/7**, confirmed. Ingham: if ζ(½+it) ≪ t^b then N(σ,T) ≪ T^{2(1+2b)(1−σ)}log⁵T, so Lindelöf ⇒ density hypothesis (web-checked, via Baluyot §1.2). This confirms Part VI's (I4) REMARK.
- **[M4] Ingham:** ∫_0^T|ζ(½+it)|⁴dt ~ (2π²)^{−1}T log⁴T.
- **[M12] Heath-Brown** (Quart. J. Math. 29 (1978) 443–462): ∫_0^T|ζ(½+it)|^{12}dt ≪ T²log^{17}T.
- **[M8] Hölder of [M4] and [M12]:** ∫_0^T|ζ|⁸ ≤ (∫|ζ|⁴)^{1/2}(∫|ζ|^{12})^{1/2} ≪ T^{3/2}log^{21/2}T. More generally ∫_0^T|ζ|^p ≪ T^{1+(p−4)/8+ε} for 4 ≤ p ≤ 12. These rest on the approximate FE and divisor-type mean values; see H.21.6 for the EP status.
- **[EP4] For σ > 1:** |ζ(σ+it)| ≥ ζ(2σ)/ζ(σ). This is the Euler product. Without it, the bound |ζ(2+it)| ≥ 2 − ζ(2) > 0.35 is fake-blind.
- **[Had]** ξ(s) = ½s(s−1)π^{−s/2}Γ(s/2)ζ(s) is entire of order 1 and ξ(s) = ξ(0)∏_ρ(1 − s/ρ), the product taken symmetrically.
- **[St] Stirling:** log|Γ(x+iy)| = (x−½)log|y| − π|y|/2 + ½log 2π + O_x(1/|y|).

## H.21.1 The exact identity (PROVED; fake-blind)

For ½ < σ_b, write h_b := σ_b − ½. For t ≥ 2:

  log|ζ(½+it)| = log|ζ(σ_b+it)| − Z_b(t) + (h_b/2)log(t/2π) + O_b(1/t),  where Z_b(t) := Σ_ρ log(|σ_b+it−ρ|/|½+it−ρ|).

The sum runs over all nontrivial zeros, taken symmetrically. It converges absolutely: a zero at distance u contributes h_b(h_b + 1 − 2β)/(2u²) + O(u^{−4}).

*Proof.* Take [Had] at s = σ_b+it and s = ½+it and subtract. The log|ρ| terms cancel. The Γ-factor difference is (h_b/2)log(t/2) + O(1/t) by [St], and π^{−s/2} contributes −(h_b/2)log π. ∎

This uses only FE and order 1. DH satisfies the analogue with log(5t/2π) for its conductor 5 and Γ((s+1)/2); so does Epstein, with its own conductor.

## H.21.2 Pair domination and the window mean (PROVED (sketch); fake-blind for σ_b = 2)

**Lemma H.21.2a (pair domination).** Let k_h(u) := ½log(1 + h²/u²) be the on-line kernel. For a zero pair ρ = ½ ± d + iγ with 0 ≤ d < ½ (the functional-equation pair), the pair's contribution to Z_b at t is ≤ 2k_{h_b}(t − γ) pointwise.

*Proof.* Put v := d² + u² ≥ u². The pair's contribution is ½ log of
  [((h−d)²+u²)((h+d)²+u²)]/v² = (1 + h²/v)² − 4h²d²/v² ≤ (1 + h²/v)² ≤ (1 + h²/u²)². ∎

On-line zeros contribute k_{h_b}(t − γ) exactly. So Z_b(t) ≤ Σ_γ k_{h_b}(t − γ), summed with multiplicity over all zeros (∫ … dN), whatever their real parts.

Per-zero totals (the lead's fix (1)): ∫_ℝ log(|σ_b+iu−(β−½)…|) du = π(|σ_b−β| − |½−β|) ≤ πh_b, consistent with ∫k_h = πh.

**Lemma H.21.2b (window mean).** For every window W = [t₀, t₀+a] ⊂ [2, T]:

  (1/a)∫_W log|ζ(½+it)|dt ≥ (1/a)∫_W log|ζ(σ_b+it)|dt − 2πh_b S*(2T)/a − O_b(1/t₀ + h_b²/a).

*Proof.* Let K_W := 1_W * k_{h_b}. It is a symmetric unimodal convolution of two symmetric unimodal functions, with 0 ≤ K_W ≤ πh_b and ∫K_W = aπh_b.
- By H.21.2a, ∫_W Z_b ≤ ∫K_W dN.
- Split dN = dN₀ + dS, with N₀ = θ(γ)/π + 1 and θ′(γ) = ½log(γ/2π) + O(γ^{−2}).
  - Then ∫K_W dN₀ = (1/2π)∫K_W(γ)log(γ/2π)dγ = (ah_b/2)·avg_W log(t/2π) + O(h_b²). The kernel is even, so the first-order Taylor term vanishes.
  - This cancels the Stirling term **exactly**: this is the lead's "on-line cost cancels".
- Integration by parts: ∫K_W dS = −∫S·K_W′ ≤ S*·TV(K_W) = 2S*·max K_W ≤ 2πh_b S*.
- The region |γ − t₀| > T contributes O(a h_b² log T/T), and the conjugate zeros O(log T/T). ∎

Two consequences:
- σ_b = 2, h_b = 3/2: log|ζ(2+it)| ≥ log(2 − ζ(2)), which is fake-blind. So **mean_W log|ζ| ≥ −(3πC_S/a)log T − O(log log T/a + 1)**, with 3πC_S ≈ 0.978.
- σ_b = 1+δ [EP4]: log|ζ(1+δ+it)| ≥ −log(1/δ) − O(1). So **mean_W log|ζ| ≥ −(π(1+2δ)C_S/a)log T − log(1/δ) − O(1)**, with πC_S ≈ 0.326.

Write m_W for the window mean and D := 2πh_b C_S, so that m_W ≥ −(D/a)log T − O(1).

REMARK: this is a lower bound for **every** window. No zero-free hypothesis, no cleanliness and no density are used. Off-line zeros only *help* (H.21.2a).

## H.21.3 From mean to floor; the clean-window constant A (PROVED (sketch))

**Lemma H.21.3 (half-measure conversion).** Recall η_W := sup{η : |{t∈W : |ζ| ≥ η}| ≥ a/2}. Then

  log(1/η_W) ≤ min( sup_W log⁺|ζ| − 2m_W , 2M⁺_W − 2m_W ),  where M⁺_W := (1/a)∫_W log⁺|ζ(½+it)|dt.

*Proof.* Suppose the bad set {log|ζ| < −V} has measure ≥ a/2. Then
  a·m_W ≤ ∫_good log|ζ| − (a/2)V ≤ min((a/2)sup⁺, aM⁺) − (a/2)V. ∎

**Clean-window constant (the lead's test).** Take an h-clean window (Part VI H.19.2) and σ_b = σ₁ = ½ + h.
- Selberg's lemma [EP2] gives |log|ζ(σ₁+it)|| ≤ A_E h log T + O_h(1).
- For the near zeros, Z_b ≥ −O(h): each k_ρ ≥ 0 when β − ½ ≤ h/2.
- So sup⁺ ≤ (½ + A_E)h log T + O_h(1), and m_W ≥ −(A_E + 2πC_S/a)h log T − O_h(1).

Hence

  **η_W^{−1} ≤ T^{A h}, with A_clean = ½ + 3A_E + 4πC_S/a.**

**Decisive test (PROVED for this argument).** A_clean ≥ ½ > 1/7 > θ_d/8 for every a and every A_E ≥ 0. So the graded strategy cannot reach ν = 0 by this route.
- REMARK: this is what the argument extracts. It is not a lower bound on true floors.
- The obstruction is the sup⁺ term (½)h log T, i.e. the T^{h/2} Γ-gap. Replacing sup⁺ by M⁺ removes it, but M⁺ is not ∝ h; that is the next section.
- Under the lead's HEURISTIC, A ≈ (S-const + A_E)/a: correct for m_W (coefficient 2πC_S/a), but the floor needs the extra ½.

## H.21.4 The moment route: (NM₂) with ν(a) = 1 + 16D/a (PROVED (sketch))

From H.21.3 with the M⁺ branch and H.21.2b:

  η_W^{−4} ≤ exp(8M⁺_W − 8m_W) ≤ T^{8D/a}·(log T)^{O(1/a)}·exp(8M⁺_W).

**Bounding exp(8M⁺_W).**
- log⁺x ≤ ½log(1 + x²), and Jensen gives M⁺_W ≤ ½log(1 + avg_W|ζ|²).
- So exp(8M⁺_W) ≤ (1 + avg_W|ζ|²)⁴ ≤ 8(1 + avg_W|ζ|⁸), by the power-mean inequality.
- Average over the windows tiling [−T, T] (using |ζ(½−it)| = |ζ(½+it)|) with [M8]:

  (1/T)∫_{|s|≤T} η_{W(s)}^{−4} ds ≪ T^{8D/a}·(1 + T^{−1}∫_0^{2T}|ζ|⁸) ≪ T^{½ + 8D/a + ε}.

Windows below a fixed height T₀ are finitely many, each with η_W > 0. With T = 2πΛ²:

  **(NM₂)(ν) holds with ν(a) = 1 + 16D/a + ε.**

| version | σ_b | D = 2πh_b C_S | ν(a) | ν < 2 iff |
|---|---|---|---|---|
| fake-blind | 2 | 3πC_S ≈ 0.978 (Trudgian: 1.056) | 1 + 15.65/a | a > 15.7 |
| [EP4] | 1+δ, δ → 0 slowly | πC_S ≈ 0.326 (Trudgian: 0.352) | 1 + 5.22/a | a > 5.3 |

Labels:
- **PROVED (sketch):** for fixed a > 16D, (NM₂) holds with ν < 2 unconditionally.
- **REMARK:** the "1" is [M8]'s T^{1/2}. The conjectural eighth moment T log^{16}T (CONJECTURE, CFKRS/Keating–Snaith) would give ν = 16D/a.
- **REMARK:** S large-deviation (Selberg/Tsang moments) gives only sub-power savings on T^{8D/a}, because the S tails at level ≍ log T are exp(−c√log T)-type unconditionally. So D is the right constant for this method.

**Generalised thickness γ (REMARK).** Suppose the floor uses fraction γ in place of ½ (Kovrijkine then uses (γ, a)-thickness). The conversion gives log(1/η) ≤ (M⁺ − m)/(1 − γ), and the needed moment is p = 4/(1 − γ). Then [M4]–[M12] give

  ν(a, γ) = (γ + 8D/a)/(1 − γ).

As γ → 0 this tends to 8D/a, using Ingham only. But the Kovrijkine loss grows like log(1/γ); see H.21.5.

## H.21.5 Joint optimisation with Kovrijkine; the exact threshold on κ_K (PROVED bookkeeping + HEURISTIC sharpness)

At γ = ½: **Φ(a) = 2κ_K a + 1 + 16D/a.**

The minimum is at a* = √(8D/κ_K), where Φ(a*) = 1 + 8√(2κ_K D). Hence:

  **β < 2 ⇔ κ_K < 1/(128D)**, i.e.:
  - κ_K < 0.0240 ([EP4], HSW) / 0.0222 ([EP4], Trudgian);
  - κ_K < 0.0080 (fake-blind).

  The admissible window is the open interval of a with 2κ_K a² − a + 16D < 0, namely
  a ∈ ((1 − √(1 − 128κ_K D))/(4κ_K), (1 + √(1 − 128κ_K D))/(4κ_K)).

In Kovrijkine's normalisation, κ_K = 4C_K log(2C_K) at γ = ½, so the threshold is **C_K log(2C_K) < 1/(512D) ≈ 0.0060 [EP4]**.

At general γ the condition is

  min_a [(γ + 8D/a)/(1 − γ) + 2κ_K(γ)a] < 2.

The minimum over a is γ/(1−γ) + 8√(Dκ_K(γ)/(1−γ)), so this requires γ < 2/3 and, explicitly:

  **κ_K(γ) < (2 − 3γ)²/(64D(1 − γ)).**

At γ = ½ this gives back 1/(128D). The largest admissible κ_K is 1/(16D) ≈ 0.19, reached as γ → 0.

**Can known Kovrijkine constants meet it? (HEURISTIC: no.)**
- Model the restriction of B to a window as exponential polynomials with n ≈ aℓ/2π frequencies, where ℓ = 2log Λ + O(1).
- The Remez/Nazarov–Turán lemma is sharp up to the base: ‖p‖_W ≤ (c/γ)^{n−1}‖p‖_E, and Chebyshev-type examples show (1/γ)^{n(1−o(1))} is attained.
- Hence any *uniform* Kovrijkine constant satisfies K_W ≳ (1/γ)^{aℓ/2π}, i.e. κ_K(γ) ≳ (2/π)log(1/γ).
- At γ = ½: κ_K ≳ 0.44 ≫ 0.024.
- At small γ: (2/π)log(1/γ) ≥ 1.47 at γ = 0.1, while the allowance tends to 1/(16D) ≈ 0.19.
- So no γ works: the uniform loss exceeds the allowance by a factor ≥ 8 at every γ.
- Egidi–Veselić's explicit form has an unspecified universal K (web-checked), so no published constant is close.

**Verdict (Part VII).**
- **(NM₂) itself now holds with ν < 2 unconditionally**, at fixed a > 16D (PROVED (sketch)).
- **β < 2 is decided entirely by (WK₂)'s a-dependence.** It needs κ_K < 1/(128D) (≈ 0.024 [EP4]), which is GAP for Spectrum.
- Uniform Kovrijkine/Logvinenko–Sereda cannot provide it (HEURISTIC, Remez sharpness). So a **structured (WK₂)** is required, exploiting that the low set of ω is where |ζ| is small rather than an arbitrary thick set. Its status is GAP.
- This replaces Part VI's verdict as the operative bottleneck. Part VI's "NO for the zero-density route" stands, scoped to that route.

## H.21.6 Euler-product flags and fake behaviour (R2, R5)

| Step | Uses | DH | Epstein |
|---|---|---|---|
| H.21.1 identity, H.21.2a domination, H.21.2b with σ_b = 2 | FE + order 1 + Stirling; bound at σ = 2 trivially from coefficients | holds (conductor 5) | holds |
| [S] explicit S(t) | Backlund/Jensen + upper bounds on lines | analogue holds with other constants (REMARK) | same |
| [EP4] σ_b = 1+δ lower bound | **Euler product** | **fails**: f_DH has zeros with σ > 1 (Davenport–Heilbronn 1936) | **fails** (same) |
| [EP2] Selberg's lemma (clean-window A only) | Λ(n) | fails | fails |
| [M4]/[M12]/[M8] moments | approximate FE + Dirichlet-polynomial mean values with divisor-bounded coefficients; no Euler product needed beyond d_k(n) ≪ n^ε | analogue for linear combinations of L-functions mod 5 plausible (REMARK, not verified) | similar (REMARK) |
| Kovrijkine | archimedean | holds | holds |

**R5 consequence (REMARK, important for the chain).** The σ_b = 2 version of H.21.4 is **fake-blind**: (NM₂) with ν(a) = 1 + 16D_fake/a is expected to hold for DH as well. That is not a contradiction, because (NM₂) is one input of (T3), not (T3). But it means **the arithmetic of the (T3) chain must enter elsewhere**: (MB)(c), full-range flat density, (T3a*), or a structured (WK₂). It must not be assumed to come from (NM₂).
- The [EP4] version is arithmetic, but only through the constant D (0.326 vs 0.978), not qualitatively.
- I flag this to the lead and to Spectrum as an R5 checkpoint for the chain.
- No RH or RH-derived input is used (BARRED). Under RH, S* = O(log T/log log T), D → 0 and ν → 1. This is noted only to show the T^{1/2} moment floor is not RH-sensitive; it is not used.

R-line H.21: R1 Λ-uniform | R2 Euler product only at [EP2] (clean-A, which fails anyway) and [EP4] (optional constant) | R3 n/a | R4 n/a | R5 the σ_b = 2 route is fake-blind, flagged | R6 the bookkeeping is exact.

## H.22 Summary table (Part VII)

| Claim | Label |
|---|---|
| Exact identity log|ζ(½+it)| = log|ζ(σ_b+it)| − Z_b + (h_b/2)log(t/2π) + O(1/t) | PROVED (citation: Hadamard + Stirling) |
| Pair domination: zero pair ≤ 2k_{h_b} pointwise | PROVED |
| Window mean m_W ≥ −(2πh_b C_S/a)log T − O(1) for every window; C_S = 0.1038 (HSW 2022) | PROVED (sketch) |
| Half-measure conversion log(1/η_W) ≤ min(sup⁺ − 2m, 2M⁺ − 2m) | PROVED |
| A_clean = ½ + 3A_E + 4πC_S/a ≥ ½ > θ_d/8: graded ν = 0 test fails for this argument | PROVED (sketch); REMARK: argument-specific, not a bound on true floors |
| (NM₂) with ν(a) = 1 + 16D/a + ε; D = 3πC_S (fake-blind) or πC_S [EP4] | PROVED (sketch) |
| ν < 2 at fixed a > 16D (a > 5.3 [EP4]) | PROVED (sketch) |
| Φ(a) = 2κ_K a + 1 + 16D/a; β < 2 ⇔ κ_K < 1/(128D) ≈ 0.024 [EP4]; general-γ threshold κ_K(γ) < (2−3γ)²/(64D(1−γ)) | PROVED (bookkeeping), conditional on the (WK₂) loss form |
| Uniform Kovrijkine has κ_K ≳ (2/π)log(1/γ) (≈ 0.44 at γ = ½): threshold unmet by factor ≥ 8 | HEURISTIC (Remez/Turán sharpness) |
| Structured (WK₂) beating uniform Kovrijkine by ≳ 20× in the a-exponent | GAP (Spectrum) |
| Eighth-moment conjecture ⇒ ν = 16D/a | CONJECTURE input; not used |
| (I1) Ingham-graded log⁺ | superseded by H.21.4 (REMARK) |
| σ_b = 2 route is fake-blind; arithmetic of (T3) must come from (MB)/flat density/(T3a*)/(WK₂) | REMARK (R5 checkpoint) |
| (T3); RH | **OPEN** (RH-derived inputs BARRED) |

Not cited: ArithmeticTemperedPrincipalSeries, rh_iff_weil_pairedForm_nonneg, or any True-placeholder.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
