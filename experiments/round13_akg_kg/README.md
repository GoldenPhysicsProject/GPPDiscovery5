# Round 13 — AKG through the Pasterski–Shao Klein–Gordon construction: results

Heart of the Cosmos (Grok / GPPDiscovery5), 2026-11-18 (Toronto). Plan written first: `PLAN.md` (commit `ed87c7a`).
Labels: PROVED / CHECKED / NUMERICAL / CONJECTURE / GAP / OPEN. AKG (GPPVerify `GppAKG.AKG`, #610, `akg_iff_rh`), H_ψ, ARP, RH: **OPEN**.
Zeros: Platt table `zeros_14.dat` (GPP Drive "Zeros"; first 4520 zeros, matches mpmath.zetazero at n = 1, 2, 100) — **comparison only**. Every construction below uses primes, the pole and the archimedean kernel only.

## 0. Normalisation (CHECKED, NUMERICAL)
W(g) = C_∞A(0) + 2∫₀^∞ w(t)[A(0) − A(t)]dt + 2∫A(t)cosh(t/2)dt − 2Σ_n Λ(n)n^{−1/2}A(log n), A = g ⋆ g̃, w(t) = e^{−t/2}/(1 − e^{−2t}), C_∞ = ψ(1/4) − log π.
Against the zero side Σ_ρ|ǧ(γ_ρ)|² (ǧ(r) = ∫g e^{iry}): relative agreement 9e-15 (L = 0.5), 5e-13 (L = 1), 3e-10 (L = 2; zero-tail limited), for an even and a non-even polynomial bump at each width.

## 1. The construction (Pasterski–Shao transplanted)
PS (1705.01027 §3.2 eq. 3.6–3.8, §4): the KG product is positive on positive-energy solutions because those already form a Hilbert space; Mellin (conformal-primary) modes are δ-normalizable iff Re Δ = d/2. Arithmetic transplant, with u = log x as bulk time:
- **Bulk field (prime/pole/arch only).** For real g ∈ C_c² supported in [−L/2, L/2] and u > L/2:
  Φ_g(u) = e^{−u/2}ǧ(i/2) + e^{u/2}ǧ(−i/2) − Σ_n Λ(n)n^{−1/2} g(log n − u) − ∫ w(x) g(x − u) dx.
  By the explicit formula applied to g(· − u), Φ_g(u) = Σ_ρ ǧ(γ_ρ)e^{iγ_ρ u} (CHECKED; NUMERICAL pointwise at u = 3, 7, 11: prime side = zero side to all 7 printed digits, all six tests). It is e^{−u/2}(ψ(e^u) − e^u + trivial-zero terms), smoothed by g.
- **Hilbert space.** H = B², the Besicovitch space (mean ⟨F, G⟩ = lim (1/U)∫^U F Ḡ, quotient by null functions, completed). A mode e^{(ρ−1/2)u} has B² norm 1 if Re ρ = 1/2 and ∞ if Re ρ > 1/2 — Daniel's "Cesàro: 1 on the line, diverges everywhere else", which is the PS δ-normalizability statement.
- **Map.** ι(g) = Φ_g, linear in g, exhibited from prime/arch data — not an abstract GNS witness.

## 2. Candidates, tested at L = 0.5, 1, 2 (U = 13, primes to e^{14.05}; `akg_kg.py` → `akg_kg_results.json`)
### C-B2 (Cesàro/Besicovitch KG norm) — RIGHT as a conditional exhibited witness; does NOT prove AKG
| L | test g | W(g) (primes+arch) | Hann–Cesàro mean of Φ_g² on [2, 13] | rel. err |
|---|---|---|---|---|
| 0.5 | even | 3.11325e-2 | 3.11324e-2 | −5.5e-6 |
| 0.5 | asym | 3.42617e-2 | 3.42610e-2 | −1.9e-5 |
| 1 | even | 7.70264e-4 | 7.70265e-4 | +1.3e-6 |
| 1 | asym | 1.74540e-3 | 1.74541e-3 | +3.7e-6 |
| 2 | even | 2.26737e-6 | 2.26737e-6 | +1.6e-6 |
| 2 | asym | 3.18168e-6 | 3.18169e-6 | +4.0e-6 |
NUMERICAL: the B² norm of the prime-built field equals the Weil LHS at all three widths, from primes alone, to ≤ 2e-5.
Exact status:
- **(R13.1, CHECKED, classical almost-periodic Parseval):** under RH, Φ_g is uniformly almost periodic with absolutely summable spectrum (ǧ(γ) = O(γ^{−2}), N(T) ~ T log T), so ‖Φ_g‖²_{B²} = Σ_{distinct γ} m_γ²|ǧ(γ)|². Hence ‖ι g‖² = W(g) = Σ m_γ|ǧ(γ)|² for all L **iff, in addition, every zero with ǧ(γ) ≠ 0 is simple.** So RH ∧ (simple zeros) ⇒ AKG with an *exhibited* witness (B², Φ), which is the witness-level content `AkgStatement.lean` says is owed.
- **(R13.2, CHECKED, Landau-type):** Φ_g ∈ B² (finite Cesàro mean square) for all g ⇔ RH. ⇒: if ∫^U Φ_g² ≤ CU, then ∫e^{−2εu}Φ_g² < ∞ for every ε > 0, so the Laplace transform of Φ_g (a multiple of −ζ′/ζ(s + 1/2) times ǧ, plus pole/arch terms, with residue m_ρ ǧ(γ_ρ) at s = ρ − 1/2) is analytic in Re s > 0; choosing g with ǧ(γ_ρ) ≠ 0 excludes every ρ with Re ρ > 1/2. ⇐: R13.1.
- **Plain verdict:** positivity is free (a mean of squares) and O8 is automatic (Φ_g has no spectrum below 14.13 — the prime signal's power spectrum *is* the zero set). But *membership* ι(g) ∈ H is exactly the RH-equivalent content: a mean-square bound ψ(x) − x = O(√x) in the Cesàro sense (Cramér/Landau). The KG route relocates AKG to normalizability of the arithmetic bulk field; it does not prove it. **O9 (new):** in the PS transplant, Re Δ = 1/2 comes out of normalizability only once Φ_g ∈ H is known a priori; PS get that from the KG Cauchy problem, and the arithmetic analogue of that a-priori bound is RH-equivalent. Multiplicity gap: **C-Mult OPEN** (the B² witness carries m², Weil carries m; numerics cannot see it since known zeros are simple).
- Falsifier check against the D2 Mellin twist (CHECKED): for the shifted Ξ_θ the field has modes e^{(ρ ± θ − 1/2)u}, which are off the line, so ‖Φ^θ‖_{B²} = ∞ and C-B2 does **not** certify the negative F_{1/4} box. It is not twist-invariant, so it passes the D2 falsifier.
- All-L detection (NUMERICAL, hypothetical toy; zeros as background only): one off-line quartet β = 0.1, γ₀ = 28 inserted. The Cesàro ratio ‖Φ‖²/W grows monotonically with U at L = 0.5 (1.028 → 1.051 for U = 5 → 13) and at L = 1 (1.0028 → 1.0069). At L = 2 the toy amplitude |ǧ(28 ∓ 0.1i)|² is below the window leakage for U ≤ 13 (ratio 1.0009 → 1.00014, not monotone). Asymptotically it grows like e^{2βU}/U for every L with ǧ(γ₀ ∓ iβ) ≠ 0. Contrast s2 (Round 11): fixed-width forms are blind to the pair below L*.
### C-Abel (one-sided weight e^{−2εu}) — RIGHT but weaker; same limit, O(ε) bias
Weighted means at ε = 0.5 / 0.25 / 0.1 vs W: L = 0.5 even 0.844 / 0.920 / 0.967 × W; L = 1 even 1.044 / 1.023 / 1.010 × W; L = 2 even 0.976 / 0.987 / 0.995 × W. They converge to W as ε → 0, with Lorentzian leakage ε²/(ε² + Δγ²). This is the `CesaroAbel.abelForm` shape of GPPVerify's `WeilIsGram`, with F = Φ_g, which is a *prime-built infinite* family, not a finite character family. That evades Flint's three falsifications of the finite and ε-independent shapes. Finiteness for small ε is again R13.2 (RH-equivalent). NUMERICAL.
### C-Window (finite-U KG norm, no limit) — KILLED as a certificate
Q_U − W relative to W, over windows of length 1 / 2 / 4 placed across [2, 13]: L = 0.5: [−0.35, +0.41] / [−0.043, +0.042] / [−0.0033, +0.0051]; L = 1: [−0.12, +0.13] / [−0.0083, +0.0080] / [−0.0011, +0.0009]; L = 2: [−0.072, +0.077] / [−0.0046, +0.0036] / [−0.0006, +0.0006]. Exact obstruction: the leakage Σ_{γ≠γ′} ǧ(γ)ǧ(γ′)*σ̂_U(γ − γ′) is indefinite, so Q_U ≥ 0 (always true) says nothing about the sign of W at any finite U. NUMERICAL.
### C-Pair (prime-pair form of the identity) — CHECKED reformulation, numbers below
Mean Φ_g² = diagonal Σ_n Λ(n)²n^{−1}g(log n − u)² (which grows ~ u∫g² by PNT) + off-diagonal prime pairs |log(m/n)| ≤ L (with the main-term cross). Hann means: diagonal 1.111 / 2.222 / 4.444 (L = 0.5 / 1 / 2, even g), against W = 3.1e-2 / 7.7e-4 / 2.3e-6. So the off-diagonal prime pairs cancel the diagonal to 1 part in 36 / 2900 / 2,000,000. The B² identity is an averaged Hardy–Littlewood pair statement over multiplicative windows of width L; for all g it is RH-equivalent by R13.2. NUMERICAL (cancellation) / CHECKED (reformulation).

## 3. Best lead: bulk time = width, so a *defect-tolerant* AKG suffices
For real g, h supported in [−L/2, L/2], Φ_{g⋆h̃}(u) = W(g, h(· − u)), the bilinear Weil form of two tests whose joint support has width u + L (explicit formula, CHECKED). So **evolution of the bulk field in time u is the Weil form across widths.** Cauchy–Schwarz on the shifted form W + δ(L′)‖·‖²_{L²} gives
  **(R13.3, CHECKED, Heart paper argument, audit owed).** If for some κ ≥ 0 and every ε > 0 the width-L′ Weil form satisfies W_{L′}(f) ≥ −C_ε e^{(κ+ε)L′}‖f‖²_{L²} for all L′ (on C_c² tests of width L′), then |Φ_{g⋆h̃}(u)| = O_ε(e^{(κ+ε)u}). By the Laplace argument of R13.2, every zero then has Re ρ − 1/2 ≤ κ.
- κ = 0 ("subexponential-defect AKG", **SD-AKG**) ⇒ RH. AKG itself is the case δ ≡ 0. SD-AKG ⇔ RH (⇐ is trivial), so it is RH-equivalent. But it tolerates any e^{o(L)} negative defect, which no exact-positivity certificate does.
- Calibration — **CORRECTED in `ROUND13B.md` §2**: the two-space bound gives κ = 1/2, not 1 (ϑ(N) cancels), and κ = 1/2 corresponds to the trivial half-plane Re ρ ≤ 1, not PNT. Any κ < 1/2 from primes and arch alone would be a new zero-free strip.
- **Relabelled (ROUND13B §3): R13.3 is a classical restatement** (Guinand–Weil explicit formula + Cauchy–Schwarz + the Ingham/Landau abscissa theorem; shape as in Bombieri–Lagarias 1999), not a new theorem.

## 4. Exact next step
1. Turbo/Claude audit of R13.1–R13.3. R13.3 is the one new implication; it is short and Lean-shaped (Cauchy–Schwarz plus Laplace/Landau).
2. Compute the actual defect curve δ(L′) := −min(0, λ_min(W_{L′}; L²)) and the *certified* lower bounds the two-space/Schur machinery gives, for L′ = 0.5 … 8. Find the exponent κ that the current fixed-width method yields when run uniformly in L′ (expected κ = 1 from the crude t_N, which is the first thing to beat), and which single term (prime Gram growth ϑ(N) vs arch coercivity log n) sets it.
3. Try to reach κ = 1/2 from the two-space form alone, which would reproduce PNT inside the AKG framework. That is the first non-trivial calibration of whether the framework can ever see a zero-free strip.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5). Audit: Turbo.
