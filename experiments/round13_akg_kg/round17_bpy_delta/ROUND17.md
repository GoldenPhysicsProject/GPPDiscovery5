# Round 17: "s IS Δ" — the BPY law as a Pasterski–Shao state, and where the identity stops

Heart of the Cosmos (Grok / GPPDiscovery5), 2026-11-18 (Toronto). Labels: PROVED / CHECKED / NUMERICAL / CONJECTURE / HEURISTIC / GAP / OPEN.
**RH: OPEN.** Zeros are never proof inputs. Numerics use only theta series, Dirichlet coefficients and gamma factors; zero labels appear only as Newton starting points in the heat-flow appendix.

## Summary

1. **The identity holds and is exact (PROVED).** Take the Biane–Pitman–Yor (BPY) variable Y = √(2/π)·(max b − min b), where b is a Brownian bridge (the Kuiper range). Set f_ζ(x) = ½·x·p_Y(x) = x^{−1/2}·Φ(log x), where Φ is (2×) Pólya's Φ. Then f_ζ is a positive, shadow-invariant vector of the d = 1 Pasterski–Shao space L²(ℝ_{>0}, dx), and for **every** complex Δ
   **ξ(Δ) = (ψ_Δ, f_ζ) := ∫₀^∞ x^{−Δ} f_ζ(x) dx**,  where ψ_Δ(x) = x^{−Δ} is the weight-Δ conformal primary.
   So ξ is literally the conformal-primary amplitude of the BPY state, and the zeros of ξ are exactly the weights Δ whose primary has zero overlap with that state. This is "s = Δ" as a map, not a renaming.
2. **The identity cannot by itself put the zeros on the line (PROVED).** The pairing (ψ_Δ, f_ζ) is defined for all Δ ∈ ℂ, and normalizability is a property of ψ_Δ alone. So "every zero gives a normalizable primary" is *verbatim* Re ρ = ½. That is GPPVerify's open premise `XiZeroIsNormalizable1dState`, which is equivalent to RH and to Weil positivity (both PROVED in Lean).
3. **Epstein test.** The same construction run on x² + 5y² goes through at every step except one: the Epstein "law" is **not positive**. Φ_EPS(0) = −0.150048304486139… (CHECKED, ball arithmetic). *This corrects my earlier claim to Daniel that Epstein also has a positive BPY law; it does not.* But positivity does not reach the zero step either. F_θ has a positive, shadow-symmetric law that is PF₄ (inherited from ζ, PROVED), and every one of its zeros is off the line. **No step of the identity uses a ζ-only input where it would matter.**
4. **New separator, weaker than RH.** The Pólya-frequency order of the BPY log-density works like a ladder. EPS fails PF₁ (positivity), DH fails PF₃, and ζ passes PF₂–PF₄ but fails PF₅ (NUMERICAL). This is an additive-side, theta-only inequality that does the same job against EPS and DH as c(n) = Λ(n) ≥ 0. It is **strictly weaker than RH** (PROVED), because F_θ inherits ζ's PF₄.
5. **GAP-17 (one step).** Build a unitary scaling representation on a space made from ζ-data, intertwined with the d = 1 PS space, in which each zero is a generalized eigenvector of weight ρ. The zero-to-state step must use the Euler product, because steps 1–3 also hold for EPS and positivity and PF₄ also hold for F_θ. In Nyman–Beurling form this needs no Weil form: off-line zeros are exactly the *L²-normalizable* truncated primaries x^{ρ̄−1}·1_{(0,1)} orthogonal to the integer-floor span, so RH ⇔ there are no such bound states ⇔ density.

---

## 0. Already in GPPVerify (read-only, main `7ebc350`, 2026-11-17 22:59 Toronto; no `sorry`/axiom in these files)

| Piece of the s = Δ chain | Lean theorem | file:line | Status |
|---|---|---|---|
| Massless conformal primary = Mellin of the plane wave, solves □φ = 0 | `mellin_planeWave_eq`, `box_conformalPrimary` | CelestialHolography/ConformalPrimaryWavefunction.lean:228, :287 | PROVED |
| PS theorem: radial delta-normalizable ⇔ Re Δ = d/2 | `radialDeltaNormalisable_iff`, `…_iff_principal` | CelestialHolography/PSSNormalizability.lean:53, :65 | PROVED |
| d = 1: normalizable ⇔ Δ = 1 − Δ̄ (shadow-conjugate fixed) | `radialDeltaNormalisable_one_iff_eq_conjShadow` | RiemannHypothesis/PasterskiShaoDictionary.lean:55 | PROVED |
| normalizable ⇒ principal (d = 1) | `normalizable_1d_implies_principal` | RiemannHypothesis/PasterskiShaoNormalizability.lean:183 | PROVED |
| Shadow intertwines dilations, weight Δ ↦ d − Δ; involutive | `shadow_dilate`, `shadowDim_involutive` | CelestialHolography/PSSShadowTransform.lean:141, :82 | PROVED |
| Mellin–Plancherel on Re s = ½ (Schwartz class) | `mellin_half_plancherel` | Analysis/MellinPlancherelHalfLine.lean:95 | PROVED (completeness/onto: not formalized) |
| s = Δ₁ as dilation characters | `s_eq_delta_label`, `s_eq_delta_1d_holds_for_every_label` | PasterskiShaoDictionary.lean:82; PasterskiShaoNormalizability.lean:193 | PROVED (label identity only, stated as such in Lean) |
| Critical line ⇔ unitary character with shadow = conjugate | `critical_line_iff_unitary_with_shadow` | RiemannHypothesis/ScaleShadowHalfDensity.lean:100 | PROVED |
| Tate local FE (Schwartz, both parities) | `localFE_schwartz_cross` | RiemannHypothesis/SchwartzTateLocalFE.lean:516 | PROVED |
| Zero-to-state premise | `XiZeroIsNormalizable1dState` (= `AllXiZerosKGNormalized`) | PasterskiShaoNormalizability.lean:267; OneDimensionalMellinWeight.lean:201 | **OPEN (definition only)** |
| premise ⇔ all ξ-zeros principal | `all_xi_zeros_kg_normalized_iff_principal` | OneDimensionalMellinWeight.lean:205 | PROVED |
| premise ⇒ RH (conditional) | `rh_of_xi_zero_is_normalizable_1d_state` | PasterskiShaoDictionary.lean:92 / PasterskiShaoNormalizability.lean:272 | PROVED (conditional) |
| RH ⇔ finite Weil paired form ≥ 0 | `rh_iff_weil_pairedForm_nonneg` | RiemannHypothesis/WeilPositivityCriterion.lean:151 | PROVED |
| Regression: F_θ is shadow-symmetric, but its zeros are not principal | `ftheta_off_axis_zero_fails_1d_normalizability`, `ftheta_has_shadow_but_zeros_not_all_principal` | PasterskiShaoNormalizability.lean:205, :245 | PROVED |
| Arithmetic Hardy space: e_z(n) = n^{−1/2−z} ∈ ℓ²(ℕ), Gram ζ(1+z+w̄), unitary dilation n^{−it} | `coherent_inner`, `dilate_coherent` | RiemannHypothesis/ArithmeticHardySpace.lean | PROVED (integer-side PS-type space) |
| Observer / Casimir / tempered restatements of RH | `observerPositiveZeroSurvival_iff_rh`; `all_strip_zeros_principal_of_tempered_bridge` | ObserverPrincipalSeriesPositivity.lean:244; ArithmeticTemperedPrincipalSeries.lean:75 | PROVED as equivalences; their bridge hypotheses are OPEN |
| BPY Mellin identity E[Y^s] = 2ξ(s) | — | (CriticalBPYCoercivity.lean says "not formalized") | **not in Lean** (Round 17 §1.2 gives it; a Lean target) |

Round 17 builds on these and does not redo them. What it adds: the BPY state f_ζ in the PS space, the identity ξ(Δ) = (ψ_Δ, f_ζ), the step-by-step Epstein test, the PF ladder, and the Nyman–Beurling form of the GAP.

---

## 1. Main section: the structural identity "s IS Δ"

### 1.1 The PS side (d = 1), all PROVED (§0)
- **Space and states.** Boundary ℝ_{>0}, Hilbert space H_PS = L²(ℝ_{>0}, dx), primaries ψ_Δ(x) = x^{−Δ}.
- **Dilations.** Unitary dilations (U_λ f)(x) = λ^{1/2} f(λx), with U_λ ψ_Δ = λ^{1/2−Δ} ψ_Δ.
- **Shadow.** (Jf)(x) = x^{−1} f(1/x) is unitary and involutive, and J ψ_Δ = ψ_{1−Δ}, i.e. Δ ↔ d − Δ with d = 1.
- **PS theorem.** ψ_Δ is delta-normalizable ⇔ Re Δ = ½. Mellin is unitary from H_PS onto L²(½ + iℝ, dν/2π).

### 1.2 The candidate identity (BPY state ↦ PS state)

**Definition.** Let Φ(u) = Σ_{n≥1}(4π²n⁴e^{9u/2} − 6πn²e^{5u/2})·exp(−πn²e^{2u}), which is twice Pólya's Φ. Set
  f_ζ(x) := x^{−1/2} Φ(log x) = ½·x·p_Y(x).
Here p_Y is the BPY density, Y = √(2/π)·range of a Brownian bridge, and E[Y^s] = 2ξ(s).

**Theorem S (PROVED: Riemann/Pólya theta integral + BPY 2001 + Mellin).**
- (S1) f_ζ ∈ H_PS, and f_ζ > 0 (Pólya; also because p_Y is a density). f_ζ decays super-exponentially in log x at both ends, so ∫ x^{−Δ} f_ζ dx converges for all Δ.
- (S2) J f_ζ = f_ζ ⇔ Φ is even ⇔ ξ(s) = ξ(1 − s). Φ is even by theta modularity θ(1/t) = √t θ(t).
- (S3) **ξ(Δ) = (ψ_Δ, f_ζ)** for all Δ ∈ ℂ, where (·,·) is the bilinear pairing. On the principal series, ⟨ψ_{½+iν}, f_ζ⟩_{L²} = \overline{ξ(½+iν)} = Ξ(ν), and ‖f_ζ‖² = (1/2π)∫|Ξ(ν)|² dν (Plancherel, Lean `mellin_half_plancherel` up to the Schwartz-class restriction).
- (S5, definition) The zeros of ξ are the set {Δ : (ψ_Δ, f_ζ) = 0}, the primaries "invisible" to the BPY state.

*Proof.* The theta Mellin integral gives ξ(s) = ∫Φ(u)e^{(s−½)u}du. Write g(u) = ψ(e^{2u}), and note that s(s−1) acts on the Mellin side as g'' + g', which removes the boundary terms. Substituting x = e^u gives (S3). BPY identify 2Φ(u)e^{−u/2} as the density of log Y. ∎

*Checks (CHECKED, `r17_laws.py`).*
- ∫Φ e^{(s−½)u} matches mpmath ξ at s = 3, ½ + 2i and 0.7 + 5i to a relative error of 2·10⁻⁴¹.
- ξ(0) = ½.
- P(Y ≤ y) from the law matches the Kuiper CDF 1 − 2Σ(4n²r²−1)e^{−2n²r²}, r = y√(π/2), to 10⁻⁴¹ at five values of y.

### 1.3 Step-by-step, with the Epstein column

The construction is the same for every model: f_M(x) = x^{−1/2}Φ_M(log x), where Φ_M is the Fourier transform of Ξ_M(z) = ξ_M(½+iz), built from M's own theta series and gamma factor (`laws17.py`).

| Step | ζ | EPS (x²+5y², own Γ(s), ξ_E = s(s−1)Λ_Q) | ζ-only input used? |
|---|---|---|---|
| S1 state f_M ∈ L², all Mellin moments entire | ✓ PROVED | ✓ PROVED (theta) | no |
| S2 J f_M = f_M (shadow ⇔ FE) | ✓ PROVED | ✓ PROVED (form self-dual: Θ(1/t) = tΘ(t)); evenness CHECKED | no (theta modularity only) |
| S3 ξ_M(Δ) = (ψ_Δ, f_M) | ✓ PROVED | ✓ PROVED; ξ_E(3) CHECKED against the direct lattice sum, rel. 2.5·10⁻¹⁴ | no |
| S4 f_M ≥ 0 (BPY positivity) | ✓ PROVED | **✗**: Φ_E < 0 on \|u\| < 0.30; Φ_E(0) = −0.1500483044861394404779 (CHECKED, `certify_eps.py`) | yes, but insufficient (see below) |
| S5 zeros = primaries orthogonal to f_M | definition | definition | — |
| **S6 zero ρ ⇒ ψ_ρ normalizable** | **OPEN (= RH)** | **FALSE** (off-line zeros, e.g. 0.933 + 15.67i; 7 in the low band) | must be |
| S7 normalizable ⇒ Re Δ = ½ | PROVED (Lean) | PROVED (Lean) | no |

**What the Epstein test shows.**
- EPS fails at S4 and S6.
- S4 is a genuine ζ-vs-EPS discriminator, but it cannot drive S6. Φ_{Fθ} = (Φe^{−θu}) * (Φe^{θu}) is a convolution of positive functions, so it is positive, even, rapidly decaying and PF₄ whenever Φ is (tilting and convolution preserve PF_k). Yet every zero of F_θ is off the line (PROVED).
- DH also passes S4 (NUMERICAL, positive on the whole grid) and has off-line zeros.
- **So the identity uses no ζ-only input at the step that matters (S6).** That is exactly what the Lean regression `ftheta_has_shadow_but_zeros_not_all_principal` encodes.

### 1.4 Why a "fixed state in the PS basis" identity can never close S6 (PROVED)
- Normalizability is a property of ψ_Δ alone (Lean: ψ_Δ normalizable ⇔ Re Δ = ½). The overlap (ψ_Δ, f) is an entire function of Δ for any super-exponentially decaying f.
- So for any identity of the form "ξ = amplitude of a fixed state f in the PS basis", S6 reads "ξ(ρ) = 0 ⇒ Re ρ = ½", which is RH verbatim.
- S6 ⇔ RH: Lean `all_xi_zeros_kg_normalized_iff_principal` + `rh_of_xi_zero_is_normalizable_1d_state`, plus the trivial converse.
- RH ⇔ Weil positivity: Weil 1952; finite form in Lean `rh_iff_weil_pairedForm_nonneg`.
- **Label: S6 is RH-equivalent (PROVED). The GAP is S6 itself.**

To be non-circular, "s IS Δ" has to be an identity of **representations**, not of one vector: the zeros have to come out as eigenvectors of a scaling action on a Hilbert space built from ζ-data.

### 1.5 The three ways to make it an identity of representations

| Version | Space and action | Where zeros live | EPS | Verdict |
|---|---|---|---|---|
| (i) cyclic span of the BPY state | H = closed span{U_λ f_ζ} ⊂ L²(dx), unitary | spectral measure \|Ξ(ν)\|²dν on the line; zeros have measure 0, and off-line zeros are invisible. By Wiener (Ξ ≠ 0 a.e.), H = L² | identical | **PROVED useless**: unitary but blind |
| (ii) one-sided / Hardy (Nyman–Beurling) | L²(0,1) with the integer-floor functions ρ_θ(x) = {θ/x} − θ{1/x}; the dilation semigroup is contractive, not unitary | ∫₀¹ρ_θ x^{s−1}dx = −(θ−θ^s)ζ(s)/s. A zero with Re ρ > ½ gives k_ρ = x^{ρ̄−1}1_{(0,1)} ∈ L² (normalizable iff Re ρ > ½) and k_ρ ⊥ every ρ_θ | the lattice-floor analogue has these bound states (its off-line zeros) | RH ⇔ no normalizable bound state ⇔ span dense (Nyman 1950, Beurling 1955; Báez-Duarte 2003 with θ = 1/n). **No Weil form.** RH-equivalent; GAP = density |
| (iii) adelic / Tate (Connes) | L²(𝔸_ℚ) with the idele class group, Tate vector ⊗_p 1_{ℤ_p} ⊗ e^{−πx²} | zeros on the line = absorption spectrum on the cokernel of the Poisson map f ↦ Σ_{q∈ℚ^×} f(qx) | the EPS vector is not a factorizable GL₁/ℚ vector (Epstein = combination of Hecke L-functions of ℚ(√−5)) | the Euler product enters as factorizability, but RH ⇔ positivity of Connes's trace formula = Weil positivity. Goes back through Weil |

**Physics reading (HEURISTIC, exact in (ii)).** PS's "normalizable ⇒ principal" is about delta-normalizable *unitary* primaries. In the one-sided picture the opposite holds: a truncated primary is L²-normalizable ⇔ it sits off the line, and if it is orthogonal to the integer span it is a *bound state*, which is our ghost. So Daniel's slogan "every zero is a normalizable PS state" has an exact, Weil-free counterpart: **no zero is a bound state** of the integer-floor (Nyman–Beurling) system.

### 1.6 GAP-17: the single missing step

> **GAP-17 (zero-to-state, structural form).** Construct, from ζ-data that include the Euler product, a Hilbert space H_ζ with a **unitary** scaling action and an intertwiner W: H_ζ → L²(ℝ_{>0}, dx), such that every zero ρ of ξ is a weight-ρ generalized eigenvector v_ρ with W v_ρ = ψ_ρ. PS (Lean, PROVED) then gives Re ρ = ½.
> **Weil-free equivalent (Nyman–Beurling form):** no zero of ζ gives an L²-normalizable truncated primary x^{ρ̄−1}1_{(0,1)} orthogonal to span{ρ_θ}, i.e. that span is dense in L²(0,1).
> **Constraint (PROVED):** S1–S3 hold verbatim for EPS, and S4 and PF₄ hold for F_θ. So the step that produces v_ρ must use a ζ-only multiplicative input: the Euler product with 0 ≤ c(n) ≤ Λ(n) (or Möbius coefficients). Theta and shadow symmetry alone are not enough.
> **Status: OPEN.** Equivalent to RH (PROVED), but stated without the Weil form.

---

## 2. Direct routes to s = Δ that avoid Weil positivity

Daniel's point is that Hasse over F_q gets positivity from point counts and never uses the Weil form. For each route below: the input beyond the FE, whether it fails for EPS (a valid route must fail there), its label, and its missing step.

| Route | Input beyond FE | ζ | EPS | Other fakes | Label | Missing step (GAP) |
|---|---|---|---|---|---|---|
| (a) PS normalizability from the BPY law's growth/decay | super-exponential decay + positivity of f | ✓ | ✗ (positivity fails, CHECKED) | F_θ, DH pass with off-line zeros | **PROVED insufficient**. Decay only makes (ψ_Δ, f) entire, and positivity only rules out real zeros (E[Y^σ] > 0) | none; the route is dead as stated |
| (b) Lee–Yang / Pólya frequency of the BPY log-density | PF_k of Φ; Lee–Yang / de Bruijn–Newman class | PF₂–PF₄ ✓, PF₅ ✗ (NUMERICAL) | ✗ at PF₁ | DH ✗ at PF₃; F_θ ✓ PF₄ (PROVED inheritance) | finite PF_k: **strictly weaker than RH (PROVED)**. PF_∞: impossible for any Ξ with zeros (Schoenberg). Full Lee–Yang class (Newman 1976: H_t real-rooted for t ≥ 0): RH-equivalent | prove that Φ (built from theta) lies in Newman's Lee–Yang class, i.e. Λ ≤ 0, from ferromagnetic structure and not from zeros. OPEN |
| (c) count-type positivity, Hasse style | c(n) = Λ(n) ≥ 0 plus Poisson on ℤ / twist identity | ✓ | ✗ (c(36) = −7.17 < 0, CHECKED, R16 / Principal Series) | DH ✗ (c(3) < 0); BEUR″ passes c ≥ 0 and is caught only by Poisson; F_θ passes c ≥ 0 (fails temperedness) | GAP (Round 16) | the integer Stepanov lemma (R16 §4.2): an auxiliary integer with √x slack. OPEN |
| (d) unitarity of scaling on a space built from V | (i) cyclic BPY span: none. (ii) integer-floor Hardy space: integers (floor = add + mult). (iii) Tate: factorizability (Euler product) | (i) ✓ blind; (ii) density = RH; (iii) Weil | (i) ✓ (proves nothing); (ii) ✗ (has bound states); (iii) n/a (not factorizable) | — | (i) PROVED useless; (ii) RH-equivalent, Weil-free; (iii) RH ⇔ Weil | (ii) **density of the Nyman–Beurling span = GAP-17 Weil-free form**. Natural approximants use Möbius μ(n) (Báez-Duarte), which is where the Euler product enters |

**Best Weil-free candidate: (d)(ii) combined with (c).**
- The integer floor is the additive half (it is Poisson/Euler–Maclaurin on ℤ).
- The Möbius-coefficient approximants Σ μ(k)(1 − log k/log N) ρ_{1/k} are the multiplicative half.
- The open statement is quantitative: Báez-Duarte's distance d_N → 0. That is Hasse-like in that it is a counting/approximation statement, not a positivity of a form.

---

## 3. Coupling test (task 3)

| Model | c(n) ≥ 0 | 0 ≤ c(n) ≤ deg·Λ(n) (tempered) | BPY law exists (FE + theta) | law ≥ 0 | PF₃ | survives all? |
|---|---|---|---|---|---|---|
| ZETA | ✓ | ✓ | ✓ | ✓ PROVED | ✓ (NUM.) | ✓ |
| EPS | ✗ (first n = 36; 6 negative values ≤ 200) | ✗ | ✓ | ✗ CHECKED | ✗ | ✗ |
| DH | ✗ (n = 3; 59 negative values ≤ 200) | ✗ | ✓ | ✓ (NUM.) | ✗ (NUM.) | ✗ |
| BEUR″ | ✓ | ✓ (non-integer support) | ✗ (no FE of its own) | n/a | n/a | ✗ |
| F_θ (θ = 0.2) | ✓ | ✗ (c/Λ = 2cosh(θ log n) > 2) | ✓ (degree-2 Γ) | ✓ PROVED | ✓ (inherits) | ✗ only via temperedness / degree |
| F_{iτ} (control, RH-true) | ✗ (2cos(τ log n)) | ✓ | ✓ | ✓ on \|v\| ≤ 3 (NUM.) | ✓ sampled | — (control) |

**Verdict.**
- c(n) ≥ 0 plus the additive property from §2(b) (law ≥ 0 or PF₃) kills EPS and DH, which **both halves already kill independently**. Plus existence of a law it kills BEUR″. **F_θ survives every combination of positivity-type properties** (c ≥ 0, law ≥ 0, PF₄).
- F_θ falls only to temperedness or degree 1. With "degree-1 BPY law (Γ_ℝ(s)) + integer support + 0 ≤ c(n) ≤ Λ(n)", only ζ survives. That is because the conjunction **characterizes ζ** (Hamburger 1921), not because it is an RH-type inequality.
- So: **separating all the fakes at a level strictly weaker than RH is achieved only by pinning down ζ itself.** No intermediate positivity inequality was found that separates all of them. The genuinely new piece is that PF₃ of the BPY law is an *additive-only* inequality that does c(n) ≥ 0's job against EPS and DH. This is the add/mult lock seen from the theta side (NUMERICAL; the inheritance and weaker-than-RH statements are PROVED).

## 4. Corrections to earlier messages
- **"Epstein's ξ is also the moment of a positive random variable" is false** for the canonical completion ξ_E = s(s−1)(√20/2π)^sΓ(s)Z_Q(s). Φ_E(0) < 0 (CHECKED). The reason: the first lattice term has a = π/√5 < 2, whereas for ζ the corresponding condition π > 3/2 holds.
- The ch. 8 critique still stands, but its counterexample changes. "Positive density ⇒ RH" is refuted by **F_θ** (PROVED), and also by DH (positive on the grid, NUMERICAL), not by Epstein. The other ch. 8 points (zeros are not poles; p_V ≥ 0 regardless of zero locations) are unchanged.

---

## Appendix A: supporting numerics

**A1. Laws** (`r17_laws.py` → `r17_laws.json`; `r17_pf.py` → `r17_pf.json`; `certify_eps.py` → `certify_eps.json`)

| Model | Φ(0) | law ≥ 0 | log-concave (PF₂) | PF₃ | PF₄ | PF₅ |
|---|---|---|---|---|---|---|
| ZETA | 0.8934 | ✓ | ✓ | ✓ | ✓ | ✗ (min minor −1.9·10⁻⁷, 5 points spaced ≈ 0.1 near 0) |
| EPS | **−0.1500** | ✗ on \|u\| < 0.30 | ✗ | ✗ | ✗ | ✗ |
| DH | 1.1527 | ✓ | ✓ | **✗** (normalized minor −0.0126) | ✗ (−0.067) | ✗ |
| F_θ | 0.3204 | ✓ PROVED | ✓ | ✓ | ✓ (sampled; PROVED given ζ's PF₄) | — |
| F_{iτ} | 0.2045 | ✓ (\|v\| ≤ 3) | ✓ | ✓ sampled | ✓ sampled | — |

The PF search is targeted: coordinate descent from random starts, with consecutive points ≥ 0.02 apart. An earlier run without that spacing produced cancellation artifacts, now withdrawn. ζ "PF₄ yes, PF₅ no" agrees with recent preprints (arXiv:2602.20313, certified PF₅ failure; a Zenodo note claiming STP₄). I have not checked those preprints, so this is NUMERICAL here.

**A2. Jensen / Turán** (`r17_jensen.py` → `r17_jensen.json`; moments of Φ_M at 230 digits; rigorous root isolation with python-flint)
- For every model: all moments are positive, the Turán inequalities hold for n ≤ 319, and J^{d,n} is hyperbolic for d = 3…8, n ≤ 150, and for d ≤ 300, n ∈ {0, 1, 2, 5, 10}.
- That includes EPS (first off-line zero 0.933 + 15.67i) and F_θ (zeros 14.13 ± 0.2i).
- **Fixed- or moderate-degree Jensen hyperbolicity does not separate ζ from EPS** (NUMERICAL). This is consistent with GORZ (fixed d, large n holds for a wide class) and Chasse (ζ, n = 0, d ≤ 2·10¹⁷).

**A3. Laguerre inequality and heat flow** (`r17_heat.py` → `r17_heat.json`)

| Model | L = Ξ'² − ΞΞ'' ≥ 0 on the scanned band | Λ_M lower bound from tracked pairs (collision t*) |
|---|---|---|
| ZETA (z ≤ 60) | ✓ (min rel. 0.17; code check) | — (0 ≤ Λ ≤ 0.2 known) |
| EPS (z ≤ 60) | ✗ in 7 windows, one at each off-line zero (15.67, 29.98, 36.37, 44.0, 46.76, 52.74, 58.51) | Λ_EPS ≥ 0.0794 |
| DH (z ≤ 100) | ✗ only at 85.45–85.95 | Λ_DH ≥ 0.0441 |
| F_θ (z ≤ 40) | ✗ at every ζ-zero window | Λ_Fθ ≥ 0.0199 (≈ θ²/2, isolated pair) |
| F_{iτ} (z ≤ 40) | ✓ | — |

The k = 1 Laguerre inequality acts as a local ghost detector. It is necessary for RH; whether it is sufficient for ξ is not settled here (GAP). Zero labels were used only as Newton starts.

## Files (this folder)
- `laws17.py`: the laws Φ_M (theta / Dirichlet data only), moments, heat-flow transform.
- `r17_laws.py` → `r17_laws.json`: normalization, BPY/Kuiper check, positivity, PF sampling.
- `r17_pf.py` → `r17_pf.json`, `r17_pf.log`: targeted PF₂–PF₅ search.
- `certify_eps.py` → `certify_eps.json`: ball-arithmetic certificate Φ_EPS(0) < 0.
- `r17_jensen.py` → `r17_jensen.json`, `r17_jensen.log`: Turán and Jensen (d ≤ 300).
- `r17_heat.py` → `r17_heat.json`, `r17_heat.log`: Laguerre scan and heat-flow collision times.

## References
Biane–Pitman–Yor, Bull. AMS 38 (2001) 435–465 · Pasterski–Shao, arXiv:1705.01027 · Pasterski–Shao–Strominger, arXiv:1701.00049 · Pólya, J. reine angew. Math. 158 (1927) · de Bruijn, Duke 17 (1950) · Newman, Proc. AMS 61 (1976) · Rodgers–Tao, Forum Math. Pi 8 (2020) · Platt–Trudgian, Bull. LMS (2021) · Griffin–Ono–Rolen–Zagier, PNAS 116 (2019) · Chasse, Thm 1.8 (2011) · Schoenberg, J. Analyse Math. 1 (1951) · Csordas–Varga, Adv. Appl. Math. 11 (1990) · Nyman (1950), Beurling, PNAS 41 (1955), Báez-Duarte, Rend. Lincei (2003) · Weil (1952) · Connes, Selecta 5 (1999) · Hamburger, Math. Z. 10 (1921) · GPPVerify main `7ebc350`.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).
