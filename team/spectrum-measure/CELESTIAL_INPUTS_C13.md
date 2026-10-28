# CELESTIAL_INPUTS_C13 — Spectrum & Measure inputs for Celestial §C.13, Prop C.13.A (IF)

Answers team/celestial-cft 407715d, §C.13.3 "named inputs … From §S / prolate data at Λ = 1, support L: a lower bound for σ_S(t)
and the spectral margin 1 − κ_L of K_L for L slightly above ½ log 2." **RH is OPEN. Nothing here proves anything about RH.**
Analytic only; the only computations are evaluations of closed-form elementary/special functions (no exploratory numerics).
Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok. Date: 2026-10-28.
Track: **[L]** (Λ = 1 fixed, 2L ∈ (log 2, log 3]), Celestial's. Labels PROVED / via restatement / NUMERICAL (CC's, not replayed) /
CONJECTURE / GAP.
Note (Celestial §C.13.0, agreed): "The Λ↔L dictionary (L = log Λ) is CONJECTURE, a convention fixed by Q.8.E. It is not forced by CC."
This matches our withdrawal (LAMBDA_L_dictionary.md, 60e900f/618c63b): Λ and L are independent; here Λ = 1 throughout.

Sources (verified verbatim, CC = Connes–Consani arXiv:2006.13771, Selecta 27 (2021) 77):
- CC (153)/(154): "h∞(τ) = −log π + ℜ(ψ(1/4 + iτ/2)) … It is the derivative of 2θ(τ), where θ is the Riemann–Siegel angular function".
- CC (102): "E∞(Q∞f) = −2ε′(1⁺)f(0) + ∫_0^∞ (f(x) + f(−x))(Qε)(exp(x))dx" (for all f ∈ C_c^∞(ℝ); (101) E∞(f) := ∫ f(x)ε(exp|x|)dx).
- CC Lemma 5.4 / text: ε′(1⁺) = Σ t(n), t(n) = λ(n)²/(1−λ(n)²)·ξ_n(1)², "t(0) = 11.9719, t(1) = 8.77574, t(2) = 2.20528,
  t(3) = 0.0433983, t(4) = 0.000125459 … and the total value is of the order of 22.9965" (NUMERICAL, CC).
- CC Prop 5.5: N_I = −2ε′(1⁺)(Id − K_I), "⟨η | K_I(ξ)⟩ = (1/(2ε′(1⁺))) ∫_{−log 2}^{log 2}∫ η(x)ξ(x+v)(Qε)(exp(|v|))dxdv" (104), K_I Hilbert–Schmidt;
  "Since the length of I is ≤ log 2, the function f(v) vanishes outside [−log 2, log 2] and thus we obtain (104)."
- CC Remark 5.6: "The function Qε(ρ) is 0 for ρ = 1."
- CC §6: "while for intervals of length a bit smaller than log 2 the largest eigenvalue λmax of K_I is less than 1 … this no longer holds
  true when the interval length gets closer to log 2 … the next to largest eigenvalue stays always smaller than 1!"; "λ ∼ 1.05177,
  slightly exceeds 1 when one considers the full interval [0, log 2]" (NUMERICAL / computer-assisted, discretization ω = 10⁻³).
- CC Thm 6.11 (141): W∞(g∗g*) ≥ Tr(ϑ(g)Sϑ(g)*) − c|ĝ(0)|², "c = 4γ/log 2", with γ the constant of CC Lemma 6.10 (not Euler's);
  CC: "the best constant c fulfilling (141) is such that 13 < c < 17". Proof uses "⟨ξ0 | ξ⟩ = (log 2)^{−1/2} ξ̂(0)".
Notation (Celestial's): q(v) := (Qε)(e^{|v|}), κ(v) := q(v)/(2ε′(1⁺)), ξ ∈ C_c^∞(−L, L) real, h = ξ⋆ξ*, ξ̂(t) = ∫ξ(x)e^{−itx}dx,
c₂ = (log 2)/√2, 2c₂ = √2 log 2 = 0.980258.

## 0. A structural identity first (PROVED; changes how input (1) should be used)
**(0.1) Symbol of K_L and of ε.** (102) holds for every f ∈ C_c^∞(ℝ); the length restriction in Prop 5.5 only fixes the integration
range in (104). Hence for every L, with ⟨ξ, K_Lξ⟩ := (1/(2ε′(1⁺)))∫∫ξ(x)ξ(x+v)q(v)dxdv (Celestial's K_L, a compression of convolution
by κ to L²(−L, L)), and for every ξ ∈ C_c^∞(−L, L):
  ⟨ξ, K_L ξ⟩ = (1/2π)∫ κ̂(t)|ξ̂(t)|² dt,   and, as tempered distributions (q tempered), **(¼ + t²) ε̂(t) = −2ε′(1⁺) + q̂(t)**.
Proof: apply (102) to f = h and pass to Fourier (ĥ = |ξ̂|², (Qh)^ = (¼ + t²)ĥ). ∎ Hence
  **σ_S(t) = 2θ′(t) − 2ε′(1⁺)(1 − κ̂(t))/(¼ + t²)**,   and   **κ̂(t) > 1 ⇔ ε̂(t) > 0 ⇔ σ_S(t) > 2θ′(t)**.
**(0.2) (IF) is W∞(Qh) ≥ W₂(Qh), one symbol.** Substituting (0.1) into (IF-F), the σ_S- and ε-terms recombine:
  (IF-F) ⇔ (1/2π)∫ (¼ + t²)[2θ′(t) − 2c₂ cos(t log 2)] |ξ̂(t)|² dt ≥ 0  for all ξ ∈ C_c^∞(−L, L).
(Directly: 𝒮_L − ℰ(Qh) = W∞(Qh) by T3.) So the Sonin split σ_S + (−ε̂) is a *choice of partition of one symbol*; the "split (i)/(ii)"
of C.13.3 is sufficient but lossy. Consequences used below:
- High band: the unsplit symbol 2θ′ is **better** than σ_S, because there κ̂(t) → 0 (Riemann–Lebesgue) and ε̂ < 0, i.e. σ_S < 2θ′; the
  ε-term contributes +2ε′(1⁺)(1 − κ̂(t)) ≥ 0 where κ̂ ≤ 1. Use σ_S only where it beats 2θ′, i.e. on {κ̂ > 1}.
- Low band: 2θ′(t) < 0 exactly for |t| < t_θ = 6.28984 (zero of Re ψ(¼ + it/2) − log π; Re ψ(a + iy) is increasing in y > 0 by the
  series Re ψ(a+iy) = −γ_E + Σ_n [1/(n+1) − (n+a)/((n+a)² + y²)]). Since σ_S ≥ 0 (T3, Bochner, PROVED), (0.1) gives
  **κ̂(t) > 1 for all |t| < t_θ** (ε̂(t) ≥ −2θ′(t) > 0 there). PROVED, given CC (102) and q tempered.
R-line §0: R1 every L | R2 FE via W∞, p = 2 line explicit | R3 used with sign | R4 pole spent at T4 | R5 fake-blind | R6 identity.

## 1. Input (1): explicit lower bounds for σ_S (and for 2θ′)
**(1.1) Stirling with explicit remainder. PROVED.** Binet: ψ(z) = log z − 1/(2z) − 2∫_0^∞ s ds/((s² + z²)(e^{2πs} − 1)), Re z > 0.
For z = ¼ + it/2: |s² + z²| ≥ |Im z²| = |t|/4 and ∫_0^∞ s/(e^{2πs} − 1)ds = 1/24, so the integral term is ≤ 1/(3|t|) in modulus;
Re log z ≥ log(|t|/2), |Re 1/(2z)| ≤ 1/(2t²)·… ≤ 1/(2t²). Hence for t ≠ 0:
  **log(|t|/2π) − 1/(3|t|) − 1/(2t²) ≤ 2θ′(t) ≤ ½ log(1/16 + t²/4) − log π + 1/(3|t|).**
(Check of the inequality, not exploration: 2θ′ − lower bound ≥ 1.6·10⁻³ at t = 1, 2, 5, 10, 17, 50, 200.)
**(1.2) The unsplit high-band symbol. PROVED.** 2θ′(t) ≥ 2c₂ for |t| ≥ 17.11 (from (1.1); exact crossing of 2θ′ = 2c₂ is t = 16.748).
2θ′(17) = 0.99519.
**(1.3) σ_S. PROVED (conditional constant).** From (0.1) and (1.1), if q ∈ L¹(ℝ):
  **σ_S(t) ≥ log(|t|/2π) − 1/(3|t|) − 1/(2t²) − B/(¼ + t²),  B := 2ε′(1⁺) + ‖q‖_{L¹}** (and σ_S ≥ 0 always).
Threshold t₀(B) with bound ≥ 2c₂ for all |t| ≥ t₀, and t₁(B) with bound ≥ 2c₂|cos(t log 2)| for all |t| ≥ t₁:
| B | t₀ (target 0.9803) | t₁ (target 0.9803·|cos(t log 2)|) |
|---|---|---|
| 46 (‖q‖₁ = 0, lower limit since 2ε′(1⁺) ≈ 45.993) | 19.30 | 18.58 |
| 56 | 19.70 | 18.68 |
| 66 | 20.08 | 18.76 |
| 96 | 21.12 | 18.97 |
**(1.4) Answer on the 0.98 target at |t| ≈ 17. Does NOT follow; what holds is (1.2)/(1.3).** Exactly, by (0.1):
σ_S(17) ≥ 2c₂ ⇔ κ̂(17) ≥ 1 − (2θ′(17) − 2c₂)(¼ + 17²)/(2ε′(1⁺)) ≈ 1 − 4.32/45.99 ≈ **0.906**, i.e. q̂(17) ≳ 41.7.
- With q̂(17) ≈ 0 one gets σ_S(17) ≈ 0.995 − 0.159 = 0.836 < 0.98 (and σ_S(17) ≥ 2c₂|cos(17 log 2)| = 0.695 would still need
  κ̂(17) ≥ 1 − (0.995 − 0.695)·289.25/46 < 0, i.e. it holds there iff q̂(17) ≥ −40.8 — fine for that t, but the target at the cos
  peaks needs κ̂ ≈ 0.9).
- CONJECTURE (heuristic): κ̂(17) is small (κ is a modest smooth function on short scale, CC Fig. 10), so **0.98 at |t| ≈ 17 fails for
  σ_S**; the heuristic "≈ 17" of C.13.3 is the 2θ′ threshold (16.75), not the σ_S threshold.
- **What holds (PROVED):** σ_S ≥ 2c₂|cos(t log 2)| for |t| ≥ t₁(B) ∈ [18.58, 18.97] for B ≤ 96, and σ_S ≥ 2c₂ for |t| ≥ t₀(B)
  (19.3–21.1); unconditionally only in terms of B. **GAP:** an analytic bound on ‖q‖_{L¹} (and q ∈ L¹ itself; CC Lemma F.1 is a
  computer approximation of Qε to 10⁻¹¹ by 11 terms, not an analytic L¹ bound). Recommended instead (by §0.2): use 2θ′ ≥ 2c₂ for
  |t| ≥ 17.11 (PROVED, no GAP) and put the ε-margin 2ε′(1⁺)(1 − κ̂) on the high band, where it is ≥ 0 iff κ̂ ≤ 1.
- Constant 2ε′(1⁺): its formula (Lemma 5.4) is PROVED (CC); the value 45.993 is NUMERICAL (CC). Analytic lower bound: ε′(1⁺) ≥ t(0)
  (terms positive). Upper bound needs the tail: with CC (100) |ξ_n(1)| ≤ √(2n + ½) and PROLATE_BOUNDS.md §5 at c = 2π
  (λ̃_k ≤ 10 e^{−(k−10)/1.0975}, k ≥ 4), the terms with λ(n)² ≤ ½ sum to ≤ 2Σ λ(n)²(2n + ½), explicit; the first few terms need
  rigorous evaluation (GAP, computational).
R-line §1: R1 all |t| ≥ t₀, one rule | R2 archimedean only | R3 — | R4 — | R5 fake-blind (any fake with ζ's W∞) | R6 explicit constants.

## 2. Input (2): the margin 1 − κ_L, interlacing, and the eigenvalue > 1
Define μ₁(K_L) ≥ μ₂(K_L) ≥ … (positive eigenvalues, with multiplicity) and κ_L := μ₂(K_L) (the margin below 1 after the one CC
direction). K_L is Hilbert–Schmidt for each L since q is smooth on v ≠ 0 and bounded near 0 (Celestial T9 from CC's formulas; q(0) = 0
by Remark 5.6).
**(2.1) Interlacing / monotonicity. PROVED (verifies Celestial T5).** For L < L′, K_L = J*K_{L′}J with J: L²(−L, L) → L²(−L′, L′) the
isometric extension by 0. Courant–Fischer (max over j-dimensional subspaces) gives μ_j(K_L) ≤ μ_j(K_{L′}) for every j. Hence each
μ_j(K_L), in particular κ_L, is nondecreasing in L, and #{j : μ_j(K_L) > 1} is nondecreasing. (Cauchy interlacing proper,
μ_j(K_{L′}) ≥ μ_j(K_L) ≥ μ_{j+r}(K_{L′}), needs finite codimension r and does not apply to an infinite-dimensional enlargement; the
monotonicity is what holds.) So the margin 1 − κ_L **decreases** past ½ log 2.
**(2.2) Explicit continuity modulus. PROVED.** K_{L′} − (K_L ⊕ 0) has kernel κ(y − x) on [−L′, L′]² ∖ [−L, L]² (area 4(L′² − L²)),
so by Weyl's inequality
  |μ_j(K_{L′}) − μ_j(K_L)| ≤ ‖κ‖_{L^∞(|v| ≤ 2L′)} · 2√(L′² − L²),
hence **1 − κ_{L′} ≥ 1 − κ_{½log 2} − 2‖κ‖_{L^∞(|v|≤2L′)} √(L′² − ¼log²2)**, L′ = ½ log 2 + η: √(·) = √(η log 2 + η²).
The margin persists for η < (1 − κ_{½log2})²/(4‖κ‖_∞² log 2) (to first order).
- Inputs: 1 − κ_{½log 2} > 0 is CC NUMERICAL ("next to largest … always smaller than 1"; value not stated in CC); ‖κ‖_∞ on
  [0, log 2] is ≈ 3 read off CC Fig. 10 (NUMERICAL, not a bound; on (log 2, log 3] not plotted). **GAP:** rigorous values of both.
- Weyl/HS is crude (O(√η)); a sharper rate would use smoothness of κ (not done).
**(2.3) The eigenvalue > 1 and the symbol. PROVED (given §0).**
- μ_j(K_L) ≤ sup_t κ̂(t) for all L, j (K_L is a compression of convolution by κ, norm sup|κ̂|); μ₁(K_L) → sup κ̂ as L → ∞ (Rayleigh
  quotients of long modulated bumps). By §0.1, {κ̂ > 1} ⊃ (−t_θ, t_θ), t_θ = 6.28984: the bad directions are exactly low-frequency, where
  W∞'s own symbol 2θ′ is negative and only σ_S ≥ 0 pays.
- CONJECTURE (Szegő/Landau-type count, not verified here): #{j : μ_j(K_L) > 1} ≈ (2L/2π)·|{κ̂ > 1}| ≥ (2L/2π)·2t_θ ≈ 4.0 L for large L.
  At 2L = log 2 this gives ≈ 1.39, consistent with CC's "exactly one" (NUMERICAL). Past ½ log 2 the second eigenvalue rises toward 1
  (2.1); the count suggests a second bad direction only near 2L ≈ 2 log 2 (heuristic), i.e. not on (log 2, log 3] — CONJECTURE.
- **Paying the ĝ(0) direction from 𝒮_L — exact requirement (PROVED reformulation).** Let ξ₁ be a top eigenvector, μ₁ > 1. On ξ₁, the
  ε-term is −2ε′(1⁺)(μ₁ − 1)‖ξ₁‖² < 0, so (IF) on ξ₁ requires
    𝒮_L(ξ₁) ≥ 2ε′(1⁺)(μ₁ − 1)‖ξ₁‖² + W₂(Qh₁)
  and, for mixed ξ = αξ₁ + η (η ⊥ ξ₁), the 2×2 form with the cross terms of 𝒮_L and W₂ must be ≥ 0. At L = ½ log 2 (W₂ = 0) the CC
  numbers give the cost 2ε′(1⁺)(μ₁ − 1) ≈ 45.993 × 0.05177 ≈ 2.381 (NUMERICAL), and 2.381/⟨ζ|ξ₀⟩² ≈ 2.381/0.94865² ≈ 2.65 lies in CC's
  range γ = c·(log 2)/4 ∈ (2.25, 2.95) for 13 < c < 17 (consistency check of CC's own numbers, no new numerics).
  Sufficient analytic route (PROVED as an implication; its input is a GAP): if σ_S(t) ≥ s₀ > 0 on |t| ≤ τ₀, then for ξ ∈ C_c^∞(−L, L),
  |ξ̂(t)| ≥ |ξ̂(0)| − |t|·L√(2L)‖ξ‖, so 𝒮_L(ξ) ≥ (s₀/8π)∫_{|t|≤τ₀}(|ξ̂(0)| − |t|L√(2L)‖ξ‖)₊² dt, which controls ⟨ξ₀|ξ⟩ =
  (2L)^{−1/2}ξ̂(0). **GAP:** an explicit s₀ (σ_S near 0 is not accessible analytically: 2θ′(0) = −5.372, so σ_S(0) ≥ 0 forces
  ε̂(0) ≥ 5.372, but no explicit lower bound for σ_S near 0 beyond 0 is known to us).
**(2.4) Prolate inputs (from PROLATE_BOUNDS.md, c = 2π, Λ = 1).** 2c/π = 4; β₁ = (2/π²)log 225 = 1.0975; even trace 2.23748;
Σλ(1−λ) = D(1) = 0.48633 (exact up to 3.5·10⁻⁴); KRD: λ̃_k ≤ 10e^{−(k−10)/1.0975} (k ≥ 4), Σ_{k≥K}λ̃_k ≤ 10.975e^{−(K−11)/1.0975}.
These control the convergence of CC's series for ε, ε′(1⁺) (= Σ t(n)) and Qε (CC (97)/(99)) and hence are the analytic route to the
GAPs in (1.4) and (2.2) (‖q‖₁, ‖κ‖_∞, ε′(1⁺)); the finitely many head terms still need certified evaluation.
R-line §2: R1 one rule for all L in (½log 2, ½log 3] | R2 — (archimedean K_L) | R3 operator ordering (compression), not weights |
R4 the ĝ(0) direction is the pole-adjacent CC direction, paid once | R5 fake-blind | R6 identity + sign.

## 3. Summary for Celestial
| item | label |
|---|---|
| (0.1) (¼+t²)ε̂ = −2ε′(1⁺) + q̂; σ_S = 2θ′ − 2ε′(1⁺)(1−κ̂)/(¼+t²) | PROVED (from CC (102)) |
| (0.2) (IF-F) ⇔ ∫(¼+t²)[2θ′ − 2c₂cos(t log 2)]|ξ̂|² ≥ 0 on C_c^∞(−L,L) | PROVED (identity) |
| κ̂ > 1 on |t| < 6.28984 | PROVED (σ_S ≥ 0) |
| Stirling bound for 2θ′; 2θ′ ≥ 2c₂ for |t| ≥ 17.11 | PROVED |
| σ_S ≥ 2c₂|cos| for |t| ≥ t₁(B) ≈ 18.6–19.0; ≥ 2c₂ for |t| ≥ 19.3–21.1 | PROVED given B (‖q‖₁ GAP) |
| σ_S ≥ 0.98|cos(t log 2)| from |t| ≈ 17 | not obtainable; needs κ̂ ≈ 0.9 at the cos peaks; CONJECTURE: false |
| μ_j(K_L) nondecreasing in L (T5) | PROVED (Courant–Fischer) |
| 1 − κ_{L′} ≥ 1 − κ_{½log2} − 2‖κ‖_∞√(L′² − ¼log²2) | PROVED; inputs NUMERICAL/GAP |
| ĝ(0) direction: exact payment condition; low-frequency nature | PROVED (reformulation); explicit s₀ GAP |
| count of K_L eigenvalues > 1 ≈ 4L | CONJECTURE |
| analytic L₀ > ½ log 2 | GAP (Celestial's; unchanged) |
| RH | OPEN |
Fakes: everything in §§0–2 is archimedean at Λ = 1 and holds for any L-function with ζ's γ-factor; the only arithmetic input is the
p = 2 line (Celestial T7/T9). Not stronger than RH, not equivalent to it.
