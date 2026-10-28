# SECTOR_COERCIVITY_p2 — sector-resolved archimedean coercivity for the p = 2 step (main G.4g / Prime N.12, N.G2)

Seat: Spectrum & Measure, jointly with Prime (§N N.12). **RH is OPEN. Nothing here proves or implies RH** (fixed window
2L ∈ (log 2, log 3], where ζ ≡ P_X for X ≥ 2). First principles; no numerics (the only arithmetic is evaluating explicit closed-form
thresholds, marked "eval."). Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok. Date: 2026-10-28.
Read verbatim: main 25f8f83 derivation/ARP.md G.4g; team/prime-number-theory N.12.0–N.12.7 (as merged on main 25f8f83);
Celestial §C.13 (407715d); our CELESTIAL_INPUTS_C13.md (d20c838). Track **[L]** (Λ = 1, L past ½ log 2).
Labels: PROVED / PROVED-conditional (hypotheses named) / NUMERICAL (CC's, not replayed) / CONJECTURE / GAP.

## 0. Notation (exactly main G.4g / N.12.0)
a := log 2; c₂ := (log 2)/√2 (2c₂ = √2 log 2); Q := −∂² + ¼; ξ ∈ C_c^∞(−L, L) real, h = ξ⋆ξ*; ℓ := 2L − a ∈ (0, log(3/2)];
S₋ := [−L, L − a], S₊ := S₋ + a, E := S₋ ∪ S₊ (|E| = 2ℓ); J: x ↦ −a − x on S₋, P_± := ½(1 ± J); ‖u‖²_{H(I)} := ∫_I (u′² + ¼u²);
**R_L(ξ) := 𝒮_L(ξ) + 2ε′(1⁺)(‖ξ‖² − ⟨ξ, K_Lξ⟩) = W_∞(Qh)**; B(ξ) := ∫_{S₋}[ξ′ξ′(·+a) + ¼ξξ(·+a)] = Qh(a); (IF): 2c₂B(ξ) ≤ R_L(ξ).
(K2) odd: R_L(ξ_o) + 2c₂‖P₊ξ_o‖²_{H(S₋)} ≥ 2c₂‖P₋ξ_o‖²_{H(S₋)}. (K3) even: R_L(ξ_e) + 2c₂‖P₋ξ_e‖²_{H(S₋)} ≥ 2c₂‖P₊ξ_e‖²_{H(S₋)}.
Additional (ours): ξ̂(t) = ∫ξe^{−itx}dx; for T > 0, ξ_< := F⁻¹(1_{|t|<T}ξ̂), ξ_> := ξ − ξ_< (functions on ℝ, not supported in [−L, L]).
h_∞(t) := 2θ′(t) = Re ψ(¼ + it/2) − log π (CC (153)); σ_S = 2θ′ + ε̂ (Celestial T3); κ := q/(2ε′(1⁺)), q(v) := (Qε)(e^{|v|}) (CC (102)/(104)).
Sector margins: μ₁ ≥ μ₂ ≥ … eigenvalues of K_L; K_L commutes with x ↦ −x (κ even), so they split by parity.
  δ_o(L) := 1 − sup spec(K_L|odd),  μ₁(L) := sup spec(K_L) (even, CC's ĝ(0) direction),  δ_e(L) := 1 − sup spec(K_L|even ⊖ ζ₁),
with ζ₁ a normalized top eigenvector.

## 1. Two elementary lemmas (PROVED)
**Lemma 1.1 (orthogonality and parity bookkeeping).** J is an isometric involution of H(S₋) (it preserves ∫u′² and ∫u²), so P_± are
complementary orthogonal projections in H(S₋). For odd ξ, x ↦ −x maps S₋ onto S₊, so ‖ξ‖²_{H(E)} = 2‖ξ‖²_{H(S₋)}; same for even ξ.
Hence on either parity, for σ ∈ {e, o},
  c₂‖ξ_σ‖²_{H(E)} = 2c₂(‖P₊ξ_σ‖²_{H(S₋)} + ‖P₋ξ_σ‖²_{H(S₋)}),
so **R_L(ξ_σ) ≥ c₂‖ξ_σ‖²_{H(E)} implies (K2) (σ = o) and (K3) (σ = e)**, with slack 4c₂‖P_{helped}ξ_σ‖²_{H(S₋)} (Prime's N.12.4(ii)
overpayment, accepted here). Edge order (N.12.3): c₂‖ξ‖²_{H(E)} = 2c₂s²ℓ + O(ℓ²), i.e. this is exactly the **margin 2c₂s²ℓ** asked for.
**Lemma 1.2 (low band on the strips).** For every ξ ∈ L²(ℝ) and T > 0:
  ‖ξ_<‖²_{H(E)} ≤ (2ℓ T(¼ + T²)/π) ‖ξ‖²_{L²}.
Proof: |ξ_<′(x)|² + ¼|ξ_<(x)|² ≤ (2π)⁻²(∫_{|t|<T}(t² + ¼)^{1/2}|ξ̂|)² (Minkowski in ℝ²) ≤ (2π)⁻²·2T·(¼ + T²)·2π‖ξ‖² (Cauchy–Schwarz,
Plancherel); integrate over E, |E| = 2ℓ. ∎ Also ‖ξ_>‖²_{H(E)} ≤ ‖ξ_>‖²_{H(ℝ)} = (1/2π)∫_{|t|≥T}(¼ + t²)|ξ̂|².

## 2. The coercivity theorem (PROVED-conditional)
**Theorem 2.1.** Fix a parity σ. Suppose, for some θ ∈ (0, 1], δ > 0, η > 0, T > 0:
- (H_σ) R_L(ξ) ≥ θ 𝒮_L(ξ) + 2ε′(1⁺) δ ‖ξ‖² for all ξ of parity σ in C_c^∞(−L, L);
- (HB) σ_S(t) ≥ c₂(1 + η)/θ for |t| ≥ T.
If **ℓ ≤ ℓ*(θ, δ, η, T) := π ε′(1⁺) δ η / (c₂ (1 + η) T (¼ + T²))**, then R_L(ξ_σ) ≥ c₂‖ξ_σ‖²_{H(E)}; hence (K2) (σ = o) or (K3) (σ = e).
*Proof.* σ_S ≥ 0 (T3) and (HB): θ𝒮_L(ξ) ≥ θ(1/2π)∫_{|t|≥T}(¼+t²)σ_S|ξ̂|² ≥ c₂(1+η)‖ξ_>‖²_{H(ℝ)}. By ‖x+y‖² ≤ (1+η)‖x‖² + (1+η⁻¹)‖y‖²
and Lemma 1.2, c₂‖ξ‖²_{H(E)} ≤ c₂(1+η)‖ξ_>‖²_{H(ℝ)} + c₂(1+η⁻¹)(2ℓT(¼+T²)/π)‖ξ‖² ≤ θ𝒮_L(ξ) + 2ε′(1⁺)δ‖ξ‖² ≤ R_L(ξ) under the
ℓ-condition. Lemma 1.1. ∎
R-line: R1 one rule for all ℓ ≤ ℓ* | R2 FE only via W∞ (in R_L) | R3 uses the J-sign of B through Lemma 1.1, not weight signs |
R4 pole pair spent at T4; the ĝ(0) direction is handled in (H_e) | R5 fake-blind: any datum with ζ's W∞ and the same p = 2 line.

## 3. (a) ODD sector
**(3.1) (H_o) holds with θ = 1, δ = δ_o(L). PROVED** (spectral theorem: K_L|odd ≤ (1 − δ_o)·Id, so 2ε′(1⁺)⟨ξ,(1 − K_L)ξ⟩ ≥ 2ε′(1⁺)δ_o‖ξ‖²).
Note ξ̂(0) = 0 for odd ξ, so CC's eigenvalue > 1 (whose eigenvector has ⟨ζ|ξ₀⟩ ≈ 0.94865 ≠ 0, CC §6.7, NUMERICAL) is even: the odd sector
carries no ĝ(0) direction (PROVED that the ĝ(0) functional vanishes on odd ξ; parity of ζ₁ is CC NUMERICAL).
**(3.2) Margin near ½ log 2. PROVED (CELESTIAL_INPUTS_C13 §2.1–2.2, sector-wise).** δ_o(L) is nonincreasing in L (Courant–Fischer on the
odd subspace) and δ_o(L) ≥ δ_o⁰ − 2‖κ‖_{L^∞(|v|≤2L)} √((ℓ/2)(a + ℓ/2)), δ_o⁰ := δ_o(½ log 2).
**(3.3) Odd result. PROVED-conditional.** Let δ_o⁰ > 0 and k_∞ := ‖κ‖_{L^∞(|v| ≤ log 3)}. If
  ℓ ≤ ℓ₀^odd := min{ ℓ*(1, δ_o⁰/2, η, T), (δ_o⁰)²/(12 a k_∞²) },  with T = T(c₂(1+η), B) from §5,
then (K2) holds, with margin R_L(ξ_o) ≥ c₂‖ξ_o‖²_{H(E)} (= 2c₂s²ℓ + O(ℓ²)). (Second term: forces δ_o(L) ≥ δ_o⁰/2, using a + ℓ/2 ≤ 3a/2.)
Eval. (formula only; ε′(1⁺) = 22.9965 CC NUMERICAL): best η ≈ 0.6 gives ℓ*(1, δ, η, T) ≈ 0.0120·δ (B = 46, T = 16.6) and
≈ 0.0104·δ (B = 66, T = 17.5). So ℓ₀^odd ≈ min{0.006 δ_o⁰, (δ_o⁰)²/(8.32 k_∞²)}.
Inputs and their status: δ_o⁰ > 0 — CC NUMERICAL ("next to largest eigenvalue stays always smaller than 1"; value not stated) → **GAP
(certified value)**; k_∞ — **GAP** (≈ 3 read off CC Fig. 10 on [0, log 2], not a bound); B — **GAP** (‖q‖_{L¹}). Everything else PROVED.
**Verdict (a): N.G2-odd is reduced to three explicit numbers (δ_o⁰, k_∞, B); given them, an analytic ℓ₀^odd > 0 follows (PROVED-conditional).**
Size: ℓ₀^odd is small (≲ 10⁻² δ_o⁰), far below CC's numerical μ_a (Prime N.12.4(iii)); the route is not sharp (S2-type overpayment).

## 4. (b) EVEN sector, including the ĝ(0) eigenvalue
**(4.1) Spectral split. PROVED.** On even ξ, with Π₁ = |ζ₁⟩⟨ζ₁|: 1 − K_L ≥ δ_e(1 − Π₁) − (μ₁ − 1)Π₁, so
  R_L(ξ) ≥ 𝒮_L(ξ) + 2ε′(1⁺)[δ_e‖ξ‖² − (μ₁ − 1 + δ_e)|⟨ζ₁, ξ⟩|²].
**(4.2) Payment inequality (P) ⇒ (H_e) with θ = ½, δ = δ_e/2. PROVED.** (P): ½𝒮_L(ξ) + ε′(1⁺)δ_e‖ξ‖² ≥ 2ε′(1⁺)(μ₁ − 1 + δ_e)|⟨ζ₁,ξ⟩|².
(Subtract (P) from (4.1).) This is exactly "the K_L eigenvalue above 1 is paid from 𝒮_L".
**(4.3) A sufficient condition for (P). PROVED-conditional.** If σ_S(t) ≥ σ_min > 0 for all t, then for ξ ∈ C_c^∞(−L, L)
𝒮_L(ξ) ≥ σ_min‖ξ‖²_{H(ℝ)} ≥ σ_min((π/2L)² + ¼)‖ξ‖² (Dirichlet Poincaré on (−L, L)), and |⟨ζ₁,ξ⟩| ≤ ‖ξ‖, so (P) holds if
  **σ_min ≥ 2ε′(1⁺)(2(μ₁ − 1) + δ_e) / ((π/2L)² + ¼).**
Eval. at 2L = log 2 with CC's μ₁ ≈ 1.05177 (NUMERICAL): (π/log 2)² + ¼ ≈ 20.79, so σ_min ≳ 0.229 + 2.21·δ_e.
Caveat (PROVED): ζ₁ has nonzero boundary values, so 𝒮_L(ζ₁) = +∞ (σ_S ~ log|t|, ζ̂₁ ~ 1/t); ζ₁ is not in the form domain and (P) is a
statement about C_c^∞ approximants. The crude |⟨ζ₁,ξ⟩| ≤ ‖ξ‖ step ignores this and is lossy.
**GAP (b-1): a positive lower bound σ_min of σ_S on the low band** (|t| ≲ 6.3, where 2θ′ < 0; CC Figs 5–6 say δ̂ "compensates almost
exactly", so σ_S may be small there). σ_S ≥ 0 (Bochner) is all that is proved. Equivalent sharp form of (P): the rank-one Schur condition
2ε′(1⁺)(μ₁ − 1 + δ_e)⟨ζ₁, (½𝒮̃_L + ε′(1⁺)δ_e)⁻¹ζ₁⟩ ≤ 1, 𝒮̃_L the Friedrichs operator of 𝒮_L on H¹₀(−L, L) (PROVED equivalence).
**(4.4) Even result. PROVED-conditional.** Under (P) at L and (HB) with θ = ½ (σ_S ≥ 2c₂(1 + η) for |t| ≥ T), (K3) holds for
ℓ ≤ ℓ*(½, δ_e/2, η, T). Eval.: best η ≈ 0.3, T ≈ 24.6 (B = 46): ℓ* ≈ 0.0023·(δ_e/2). Continuity in ℓ: μ₁(L), 1 − δ_e(L) increase by at
most 2k_∞√((ℓ/2)(a + ℓ/2)) (Weyl), which enters (4.3) linearly.
Remark (PROVED, N.12.3): on the even nodeless profile the p = 2 term **helps** (B(ξ_e) = −s²ℓ + O(ℓ²)); Theorem 2.1 discards this help.
A sharper even argument would credit 2c₂‖P₋ξ_e‖²_{H(S₋)} against the rank-one defect; it is O(ℓ) while the defect 2ε′(μ₁ − 1) ≈ 2.38 is
O(1), so for small ℓ the payment must come from 𝒮_L, not from the p = 2 term (PROVED as an order comparison).
**Verdict (b): N.G2-even is reduced to (P) (exactly: the rank-one Schur condition), plus (δ_e, μ₁, k_∞, B). Sufficient: σ_min ≥ 0.229 +
2.21δ_e near 2L = log 2 (GAP b-1).**

## 5. (c) High band: σ_S ≥ c₂ (and Celestial's 0.98|cos(t log 2)|)
**(5.1) Stirling with remainder. PROVED** (CELESTIAL_INPUTS_C13 §1.1, Binet): for t ≠ 0,
  log(|t|/2π) − 1/(3|t|) − 1/(2t²) ≤ 2θ′(t) ≤ ½ log(1/16 + t²/4) − log π + 1/(3|t|).
Also 2θ′ is increasing in |t| and 2θ′(0) = ψ(¼) − log π = −γ_E − π/2 − 3 log 2 − log π (= −5.3722).
**(5.2) Exact link. PROVED** (CELESTIAL_INPUTS_C13 §0.1, from CC (102)): σ_S(t) = 2θ′(t) − 2ε′(1⁺)(1 − κ̂(t))/(¼ + t²), hence with
B := 2ε′(1⁺) + ‖q‖_{L¹}: **σ_S(t) ≥ log(|t|/2π) − 1/(3|t|) − 1/(2t²) − B/(¼ + t²).**
**(5.3) Thresholds** (eval. of (5.1)/(5.2); T such that the bound ≥ target for all |t| ≥ T):
| target | 2θ′ (no B) | σ_S, B = 46 | B = 56 | B = 66 | B = 96 |
|---|---|---|---|---|---|
| c₂ = 0.4901 (N.12 (K1)) | 10.63 | 13.54 | 14.01 | 14.44 | 15.58 |
| 2c₂ = 0.9803 (Celestial) | 17.10 | 19.30 | 19.70 | 20.08 | 21.12 |
| 2c₂|cos(t log 2)| | — | 18.58 | 18.68 | 18.76 | 18.97 |
B ≥ 2ε′(1⁺) ≈ 46 always (CC NUMERICAL for ε′(1⁺)). The (K1) heuristic 2πe^{c₂} ≈ 10.26 (N.12.4(iv)) is the 2θ′ threshold up to the
Stirling remainder (10.63 provable); for σ_S the provable threshold is 13.5–15.6 given B. Celestial's "0.98 from |t| ≈ 17" holds for 2θ′
(17.10 provable) but not for σ_S (needs κ̂ ≈ 0.9 at |t| ≈ 17; CONJECTURE: false). **GAP (c-1): analytic ‖q‖_{L¹}** (CC Lemma F.1 is a
computer approximation of Qε, not an L¹ bound). Labels: (5.1)–(5.2) PROVED; thresholds PROVED for each stated B.

## 6. Summary
| item | label |
|---|---|
| Lemma 1.1 (parity/J bookkeeping; R_L ≥ c₂‖·‖²_{H(E)} ⇒ (K2)/(K3), margin 2c₂s²ℓ) | PROVED |
| Lemma 1.2 (low-band strip energy ≤ 2ℓT(¼+T²)/π · ‖ξ‖²) | PROVED |
| Thm 2.1 (sector coercivity for ℓ ≤ ℓ*) | PROVED-conditional on (H_σ), (HB) |
| (a) odd: (H_o) with δ_o; ℓ₀^odd ≈ min{0.006δ_o⁰, (δ_o⁰)²/(8.32k_∞²)} | PROVED-conditional; inputs δ_o⁰ (CC NUMERICAL), k_∞, B: GAP |
| (b) even: payment (P) ⇔ rank-one Schur condition; sufficient σ_min ≥ 0.229 + 2.21δ_e | PROVED-conditional; σ_min: GAP (b-1) |
| (c) σ_S ≥ c₂ for |t| ≥ 13.5–15.6 (B = 46–96); 2θ′ ≥ c₂ for |t| ≥ 10.63 | PROVED (σ_S given B; B: GAP c-1) |
| margin 1 − κ_L near ½ log 2, sector-wise, Weyl modulus 2k_∞√((ℓ/2)(a+ℓ/2)) | PROVED; base values CC NUMERICAL |
| N.G2 | still GAP; now reduced to the numbers δ_o⁰, δ_e⁰, μ₁⁰, k_∞, B and the even payment (P) |
| RH | OPEN |
Joint with Prime: Lemma 1.1 is N.12.3 in the form we use; Theorem 2.1 is S2 restricted to a sector and to ℓ ≤ ℓ*, which is consistent
with N.12.4 (S2 cannot reach log 3, and is wasteful: here ℓ* is O(10⁻²)).
