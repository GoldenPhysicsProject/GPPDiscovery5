# Round 13b/c: Weil defect curve, the two-space exponent κ, the classical status of R13.3, and an attack on O9

Heart of the Cosmos (Grok / GPPDiscovery5), 2026-11-18 (Toronto). This follows `README.md` (Round 13, `f3b6a97`).
Labels: PROVED / CHECKED / NUMERICAL / CONJECTURE / GAP / OPEN. AKG, H_ψ, ARP, RH: **OPEN**. Zeros are used only as test data. No PR post (Daniel: R13.3 goes to Claude only after an audit, his call).
Scripts:
- `defect_kappa.py` → `defect_kappa_results.json`: Galerkin on M = max(24, 10L′) orthonormal Dirichlet sines, built from primes, pole and arch only.
- `o9_conservation.py` → `o9_conservation_results.json`.

## 1. Defect curve δ(L′) = −min(0, λ_min(W_{L′})) on L²(I), L′ = 0.5 … 8 (NUMERICAL)
The Galerkin matrix was checked against the zero side (test data, 3000 zeros) for L′ ≤ 2. Max entrywise difference: 4.6e-6, 5.8e-7, 1.7e-7, 7.3e-8 at L′ = 0.5, 1, 1.5, 2. This is limited by zero truncation and quadrature.

| L′ | λ_min(W) (normalised, ‖f‖₂ = 1) | 2nd eigenvalue | matrix scale |
|---|---|---|---|
| 0.5 | 3.60e-2 | 4.2e-1 | 3.2 |
| 1 | 1.08e-6 | 2.4e-4 | 2.6 |
| 1.5 | 2.3e-14 | 5.5e-12 | 2.5 |
| 2 | −1.1e-14 (float floor) | 7e-16 | 2.3 |
| 3 / 4 / 5 / 6 / 7 / 8 | −4.5e-14 / −9.5e-14 / −1.1e-13 / −2.0e-13 / −2.1e-13 / −2.3e-13 | (floor) | 3.2 … 8.0 |

- **Reading.** The bottom of W falls super-exponentially with width: 3.6e-2 → 1.1e-6 → 2e-14 between L′ = 0.5 and 1.5. From L′ = 2 on it is below double-precision resolution, |λ| ≤ 3e-14 × scale.
- **Defect.** No negative value is resolved at any width. The defect curve is 0 to working precision on 0.5 ≤ L′ ≤ 8 (NUMERICAL). The small negative numbers at L′ ≥ 2 sit at the floor and scale with the matrix norm. They are not a defect.
- **What can be said rigorously (CHECKED, classical):** zeros are verified on the line up to height 3·10¹² (Platt–Trudgian 2021). So any true defect at width L′ comes only from zeros above that height, and it is ≤ C_k e^{L′/2}(3·10¹²)^{−k} for C^k tests. The defect curve is therefore astronomically small for every L′ that can be computed. **No numerical defect curve can test SD-AKG / R13.3. The content is entirely asymptotic in L′.**

## 2. The exponent κ of the fixed-width two-space method run at every width (NUMERICAL + CHECKED)
Two-space split (Verify `weilForm_regrouped`): W = B − C*C, with B = arch + pole + ϑ(N)·I and C*C = ϑ(N)·I + prime, for the minimal admissible N (L′ < log(N + 1)). The method's certified bound is the separate-term bound λ_min(B) − λ_max(C*C) = λ_min(arch + pole) − λ_max(prime). The crude closed form is m_{L,N} − t_N.

| L′ | λ_min(arch + pole) | λ_max(prime) | separate bound | crude m − t_N | −2‖s‖² (odd pole) | Σ_{n<e^{L′}} Λ(n)n^{−1/2} |
|---|---|---|---|---|---|---|
| 1 | −0.31 | 0.49 | −0.80 | −8.76 | −0.040 | 0.49 |
| 2 | −1.30 | 1.80 | −3.10 | −17.0 | −0.33 | 2.93 |
| 3 | −2.59 | 4.15 | −6.75 | −25.6 | −1.21 | 6.51 |
| 4 | −4.86 | 7.78 | −12.6 | −39.0 | −3.14 | 12.2 |
| 5 | −8.89 | 13.4 | −22.2 | −61.7 | −6.90 | 21.3 |
| 6 | −15.9 | 22.0 | −38.0 | −102 | −13.7 | 37.6 |
| 7 | −27.9 | 35.8 | −63.7 | −167 | −25.5 | 63.8 |
| 8 | −48.1 | 57.9 | −106.0 | −275 | −45.5 | 106.8 |

- **Local slopes d log|bound| / dL′ on [6, 8]:** separate bound 0.514, crude 0.496, λ_max(prime) 0.483, |λ_min(arch + pole)| 0.553. Every curve has exponent **κ = 1/2**. Empirically the separate bound tracks −Σ_{n<e^{L′}} Λ(n)n^{−1/2}, to within 0.7% at L′ = 8, which is ≈ −2e^{L′/2} (PNT).
- **Closed form (CHECKED: from in-tree B ≥ m_{L,N} plus the resolvent bound ‖Φ_p f‖ ≤ ‖f‖/(1 − p^{−1/2}), using (1 − 1/p)/(1 − q)² − 1 = 2q/(1 − q)):**
  W_{L′}(f) ≥ −c(L′)‖f‖²₂, with c(L′) = 2 sinh(L′/2) − L′ − C_∞ + 2 Σ_{p<e^{L′}} log p · p^{−1/2}/(1 − p^{−1/2}) = O(e^{L′/2}).
  So **the two-space method gives κ = 1/2 at every width, already in its crudest form.**
- **Which term sets it.** Two terms grow at exactly e^{L′/2}, and the method keeps them apart:
  1. **The prime Gram excess, taken without sign cancellation.** Its size is Σ_{n<e^{L′}} Λ(n)n^{−1/2} ~ 2e^{L′/2}. Its operator norm is about 0.54 × that sum at L′ = 8.
  2. **The odd pole channel** −2⟨f, s⟩², with ‖s‖² = sinh(L′/2) − L′/2 ~ ½e^{L′/2}.

  In the true W these two cancel. The pole is the main term ∫e^{x/2} of the prime sum (PNT), and what is left is the smoothed ψ(x) − x. Any bound that treats the prime atoms by absolute value is stuck at κ = 1/2. **The binding term is the prime sum treated without cancellation against the pole.**
- **Step 3 ("reach κ = 1/2"): reached, trivially, by the closed form above.** Correction to Round 13 §3, where I wrote two wrong things:
  - "The two-space bound gives κ = 1" was **WRONG**: I forgot that ϑ(N) cancels between m_{L,N} and t_N.
  - "κ = 1/2 is PNT-strength" was **WRONG**: through R13.3, κ = 1/2 gives sup Re ρ ≤ 1, which is the *trivial* half-plane from the Euler product. PNT (no zeros on Re s = 1) is not reached.

  So the framework, run uniformly in width, reproduces exactly the trivial bound. Doing better needs cancellation inside the prime sum against the pole, i.e. information on ψ(x) − x. **Any κ < 1/2 is a zero-free strip Re ρ ≤ 1 − (1/2 − κ), which is beyond the state of the art.** Label: CHECKED (closed form); NUMERICAL (slopes).

## 3. Is R13.3 a classical equivalence? Verdict: YES, a restatement. Relabelled "CHECKED: classical restatement, not new".
R13.3 is the composite of three standard steps:
1. **The explicit formula for translated tests.** Φ_{g⋆h̃}(u) = W(g, h(· − u)) is the Guinand–Weil explicit formula (Guinand 1948; Weil 1952, "Sur les 'formules explicites' de la théorie des nombres premiers") applied to a shifted test. Φ_g is the smoothed e^{−u/2}(ψ(e^u) − e^u).
2. **Cauchy–Schwarz on the positively shifted form W + δ‖·‖².** This is elementary. It is the usual way Weil positivity is turned into bounds on prime sums.
3. **The abscissa theorem.** If (smoothed) ψ(x) − x = O(x^{α+ε}), the Mellin transform −ζ′/ζ(s)/s − 1/(s − 1) converges absolutely for σ > α, so ζ(s) ≠ 0 there. This is Ingham 1932, Ch. IV, and Montgomery–Vaughan, *Multiplicative Number Theory I*, §15.1 (Exercise 13.1.1 with Theorem 1.3). The converse is Landau's oscillation theorem, ψ(x) − x = Ω±(x^{Θ−ε}) (MV Theorem 15.2, via Landau's lemma for Dirichlet integrals with nonnegative integrands). R13.3 uses only the easy direction (absolute convergence), not the Pólya–Landau positivity lemma.

The shape "subexponential defect ⇔ RH, exponential growth ⇔ off-line zero" is also classical. For Li's coefficients it is Bombieri–Lagarias, *Complements to Li's criterion for the Riemann Hypothesis*, J. Number Theory 77 (1999): if RH fails, λ_n is negative infinitely often with exponentially growing size, and they connect λ_n to Weil's criterion through the Guinand–Weil formula.

I did not find R13.3 stated in exactly the "defect across widths" form. But each step is textbook, and the content is Ingham/Landau. **R13.3 is not a new theorem**, and it should not go to Claude as one. If it is formalised at all, it is as a bridge lemma: Weil defect across widths ⇒ growth of smoothed ψ − x ⇒ abscissa. That is Daniel's call, after audit.

## 4. Attack on O9: a source of normalizability built from primes, pole and arch only (Daniel t355u)
O9: Re Δ = 1/2 comes out of normalizability only once Φ_g ∈ H is known a priori. Pasterski–Shao get that from KG current conservation, because the solution space of the Cauchy problem is a Hilbert space before any Mellin transform. Each candidate source below was tested without zeros. A hypothetical, FE- and conjugation-symmetric off-line quartet ρ = 1/2 ± 0.1 ± 28i is used only as a toy, to see whether the law could exclude it.

- **(a) KG-current conservation: a bounded flux J(U) = ∫₂^U (Φ_g² − W) du.**
  - Results (NUMERICAL):

    | Field | L′ | max\|J\|/W per quarter of [2, 13] |
    |---|---|---|
    | True primes | 0.5 | 0.35 / 0.34 / 0.32 / 0.33 |
    | True primes | 1 | 0.13 / 0.12 / 0.12 / 0.13 |
    | True primes | 2 | 0.088 / 0.089 / 0.089 / 0.088 |
    | With toy quartet | 0.5 | 0.33 → 0.67 |
    | With toy quartet | 1 | 0.14 → 0.18 |
    | With toy quartet | 2 | 0.086 → 0.083 (toy below resolution) |

    For the true primes the flux is flat, i.e. conserved to U = 13. With the toy quartet it grows at L′ = 0.5 and 1.
  - Verdict: **RIGHT as a conservation law and as a detector, WRONG as a source.** No local prime-side identity makes J bounded: J(U) = Σ_{γ≠γ′} ǧ(γ)ǧ(γ′)*(e^{i(γ−γ′)U} − …)/(i(γ − γ′)), and its boundedness for all g is R13.2, which is RH-equivalent. The conservation is the content, not a consequence.
- **(b) Functional-equation time reversal: Φ_g(−u) = Φ_{g̃}(u), computed from the prime side (L′ = 1, non-even g).**
  - Results: at u = 3, 6.5 and 10 the two sides agree to 5e-16, 3.5e-15 and 2.1e-14 (NUMERICAL). The identity is exact algebraically, from x ↦ −x symmetry of the explicit formula (CHECKED).
  - Verdict: **RIGHT as an identity, WRONG as a source.** The off-line quartet {ρ, 1 − ρ̄, ρ̄, 1 − ρ} is FE-symmetric by construction, so it obeys the same law and still makes the Cesàro norm diverge. The functional equation pairs ρ with 1 − ρ̄ (L1: a hyperbolic plane), and that is exactly the pairing that does not give a norm.
- **(c) Unitary scattering S(s) = ξ(2s − 1)/ξ(2s) (Lax–Phillips / Pavlov–Faddeev, "S = Δ").**
  - Results: |S(1/2 + it)| = 1 to 30 digits at t = 0.3, 3, 7.07, 14 and 20, computed from ζ (NUMERICAL). With ξ multiplied by the quartet polynomial (P(s) = P(1 − s), real coefficients), |S_toy| = 1 as well.
  - Verdict: **WRONG as a source.** Unitarity on the line holds with or without off-line zeros, because it follows from FE plus conjugation alone. This is O3, now confirmed numerically.
- **(d) An energy identity with an arch-only initial value: is W(g) = (arch + pole)(g)?**
  - Results: exact at L′ = 0.5, where there are no prime terms because log 2 > 0.5. It fails at L′ = 1: the prime term is 83% / 20% of W (even / non-even g). At L′ = 2 the prime term is 5·10⁴ to 9·10⁴ times W, and arch + pole and prime agree to 5 digits and cancel (NUMERICAL).
  - Verdict: **WRONG beyond L′ = log 2.** The "initial energy" cannot be set by the arch alone. At large width W is the tiny residue of an O(1) cancellation between arch + pole and primes, which is O7/O8 again.
- **(e) Hecke / dilation covariance: Φ_g(u + log p) = Φ_{g(· − log p)}(u).**
  - Verdict: exact and trivial (translation covariance). **WRONG as a source**: it makes the norm translation-invariant, and O2 (Bochner) then forces it onto the zero measure. It gives no bound.

**O9 status, plainly.** Every candidate source of normalizability built from primes, pole and arch falls into one of two classes:
- (i) it holds identically for *any* FE- and conjugation-symmetric zero configuration (FE reversal, unitarity of S, dilation covariance), so it cannot exclude an off-line quartet;
- (ii) it is itself the RH-equivalent statement (bounded KG flux = finite Cesàro norm, R13.2).

The arch-only initial value is simply false beyond L′ = log 2.

The structural reason is that PS's positivity comes from a Hilbert space of Cauchy data that exists before the Mellin transform. In the arithmetic setting the Cauchy data at u ≈ 0 is the arch term, but the evolution in u is driven by the primes, and no positive conserved quantity of that evolution is available unconditionally. **O9 stands (CHECKED for (b), (c), (e); NUMERICAL for (a), (d)).**

What a genuine source would have to be: a positive conserved quantity of the prime-driven evolution that is *not* invariant under the replacement {on-line pair} → {FE-symmetric off-line quartet}. Every symmetry-based law tested here is invariant under that replacement, so the source must be non-symmetric, arithmetic input. The candidates for that are the Euler product and multiplicativity themselves (e.g. positivity of Λ(n), as in Landau's lemma), not the functional equation. That is the one direction not yet tried. It is OPEN, recorded as the next step:
- test whether the positivity of Λ (the Landau/Pólya positivity lemma applied to a suitably chosen nonnegative prime-side Dirichlet integral built from Φ_g) gives any bound on the Cesàro norm better than κ = 1/2.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5). Audit: Turbo.
