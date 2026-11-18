# Round 14B — the Euler route. Part 1: Daniel's No-Ghosts theorem (physics first)

Idea and name: **No-Ghosts theorem — Daniel Toupin.** Daniel's clue: "same time". The integer state space (primon Fock space, i.e. the Euler product) together with the time reflection x → 1/x (Poisson / theta) should be ghost-free, and then the energies are real and equal to ξ's zeros.

Code (this directory): `weil_galerkin.py` (closed-form Galerkin), `models14b.py` (the models), `check_models.py` (validation), `noghosts.py` (the pairing, per model and width). Data: `check_models.json`, `noghosts.json` (Weil-reflected pairing), `noghosts_B.json` (Euclidean OS pairing). Zeros are never inputs. They appear below only as labels for where a ghost sits, taken from `../epstein_zeros*.json` and `../principal_series/ROUND14A.md`. Candidates 1–4 of Round 14B (Pólya Φ, Nyman–Beurling Gram, Selberg squares, Jensen polynomials) come in Part 2.

Labels: PROVED / CHECKED / NUMERICAL / CONJECTURE / GAP / OPEN.

## 1. The reflected pairing, built from integer data only

- **State space.** The space is the test functions f ∈ C_c^∞(ℝ_{>0}) on the scale line, with Euclidean time u = log x.
  - Mellin transform: f̂(s) = ∫ f(x) x^{s−1} dx.
  - In u-coordinates the state is g(u) = f(e^u) e^{u/2}, so ĝ(γ) = f̂(½ + iγ).
- **Time reflection.** (Θf)(x) := x^{−1} · conj f(1/x). This is x ↦ 1/x with the Poisson weight, and it is unitary on L²(dx). In u-coordinates it is g(u) ↦ conj g(−u). Its Mellin transform is (Θf)^(s) = conj f̂(1 − s̄).
- **Pairing.** ⟨f, f⟩_Θ := 𝒲(f ∗ Θf), where ∗ is multiplicative convolution and 𝒲 is the explicit-formula functional evaluated on its **prime side only**:
  𝒲(h) = [pole terms] + [archimedean term] − [prime atoms c(n) at x = n^{±1}]. The exact normalisation is the u-coordinate formula below.
  - c(n) = Λ(n) for ζ. These are the integers' multiplicative data, i.e. the Euler product.
  - The archimedean term (C, w) comes from Γ_ℝ(s). That is Poisson summation over ℤ ⊂ ℝ, i.e. the theta function.
  - The pole terms come from s = 0, 1.
  - In u-coordinates this is exactly the Round-13 form W(g) = C·A(0) + 2∫w(t)[A(0) − A(t)]dt + 2∫A(t)cosh(t/2)dt − 2Σ c(n) n^{−1/2} A(log n), with A = g ⋆ g.
- **Exact identity (PROVED; Weil 1952, Bombieri 2000).**
  ⟨f, f⟩_Θ = Σ_ρ f̂(ρ) · conj f̂(1 − ρ̄)
  - This is the explicit formula applied to h = f ∗ Θf, whose Mellin transform is f̂(s)·conj f̂(1 − s̄).
  - So the pairing **equals** the Weil form. There is no defect term.
  - The same holds with each model's own data for every model that has a functional equation (EPS, DH, F_θ).
- **Galerkin check.** The closed-form Galerkin (Dirichlet sines, A_ij(t) in closed form) agrees with the Round-13 harness to 4.0e-14 (arch), 6.2e-15 (pole), and 3.3e-15 (prime) (CHECKED, `check_models.json`).

## 2. What "no ghosts" gives, and what a ghost is

- **(i) No ghosts ⇒ real energies equal to the zeros (PROVED, standard).** Suppose ⟨f, f⟩_Θ ≥ 0 for all f.
  - Then the kernel k(u) of the pairing is a positive-definite distribution on ℝ.
  - By Bochner–Schwartz, k̂ = Σ_ρ m_ρ δ_{γ_ρ} is a positive measure, so every γ_ρ is real.
  - GNS construction: quotient by the null vectors and complete. This gives H_Θ ≅ L²(Σ m_ρ δ_{γ_ρ}).
  - Dilation acts unitarily on H_Θ. Its generator has **pure point spectrum {γ_ρ}**: the zeros are honest L² eigenvalues of i d/du on H_Θ.
  - Conversely, RH ⇒ no ghosts. So **No-Ghosts ⇔ RH**: this is Weil's criterion read as a no-ghost statement.
- **(ii) Anatomy of a ghost (PROVED).** Take real test functions and set a = f̂(ρ), b = f̂(1 − ρ̄).
  - An **on-line** zero pair (ρ = 1 − ρ̄) contributes 2|a|². That is positive.
  - An **off-line** quartet {ρ, ρ̄, 1 − ρ, 1 − ρ̄} contributes 4 Re(a · conj b).
    - The "diagonal" norms vanish: ρ pairs only with its reflection 1 − ρ̄, never with itself.
    - The whole contribution is the **cross term** between a state and its reflected partner.
    - In the 4 real coordinates (Re a, Im a, Re b, Im b) the form 4(Re a·Re b + Im a·Im b) is a pair of hyperbolic planes, with **signature (2, 2)**.
  - Hence the ghost index is at most 2 × (number of off-line quartets) (PROVED). Equality holds once the test space resolves each quartet (NUMERICAL, §3).
- **What remains.** "ζ has no ghosts" is RH itself: **OPEN**. The theorem is the equivalence plus the exact identity. The content is still the positivity.

## 3. Numbers: ghost count and ghost location (NUMERICAL, `noghosts.json`)

- Basis: M = ⌈F·L/π⌉ + 8 sines on [−L/2, L/2]. F = 45 for most models; F = 125 for DH and for the ζ re-run.
- "Ghost" means an eigenvalue below −1e-7 × max|W|.
- ghost f = peak frequency of the negative eigenvector.
- Model data:
  - ζ uses Λ(n).
  - EPS (x² + 5y²) uses its own r_Q, Γ(s), and poles at 0, 1.
  - BEUR″ (spliced Beurling, x₀ = 55) has no arch of its own and borrows ζ's.
  - DH uses its own coefficients, Γ_ℝ(s+1), conductor 5, and no pole.
  - F_0.2 = Λ(s−0.2)Λ(s+0.2) uses weights 2cosh(0.2 log n)Λ(n), Γ_ℝ(s±0.2), and poles at ±0.2, 1±0.2. It was checked against the polarisation 2B_ζ(e^{−θu}g, e^{θu}g) to 4e-10 relative.

| L | ζ (F = 45) | ζ (F = 125) | EPS | BEUR″ | DH (F = 125) | F_0.2 |
|---|---|---|---|---|---|---|
| 1 | +1.1e-6, 0 | – | +1.92, 0 | = ζ | +0.31, 0 | −2.5e-5, 1 (f ≈ 0, marginal) |
| 2 | −1.9e-14, 0 | −1.9e-14, 0 | +1.09, 0 | = ζ | +3.9e-5, 0 | −0.091, 7 |
| 3 | −4.8e-14, 0 | – | −1.7e-5, 1 (onset) | = ζ | −2.1e-13, 0 | −0.36, 16 |
| 4 | −1.4e-13, 0 | −2.7e-13, 0 | −1.62, **4** (14.5, 14.5, 31.2, 31.1) | = ζ (0.0 difference, log 55 > 4) | −0.76, **2** (84.5, 84.5) | −0.75, 20 |
| 4.5 | −2.4e-13, 0 | – | −2.96, 4 | **−0.59, 8** (31.6, 31.7, 19.8, 42.2, …) | −1.34, **4** (84.7, 84.7, 115.3, 115.3) | −1.08, 20 |
| 5 | −2.4e-13, 0 | – | −4.43, 8 | −1.13, 16 | −1.86, 4 | −1.39, 20 |
| 6 | −3.4e-13, 0 | −1.0e-12, 0 | −8.64, 8 | −2.37, 30 | −3.68, 4 | −2.53, 20 |
| 7 | −5.3e-13, 0 | – | −16.2, **10** (15.1, 16.3, 30.6, 30.6, 43.3, 43.3, 37.0, 37.0, 45.9, 45.9) | −4.06, 42 | −6.31, 4 | −3.26, 20 |
| 8 | −3.5e-13, 0 | −1.1e-9 (rel. 1.4e-10), 0 | −29.0, 10 | −5.59, 52 | −10.6, **4** (85.2, 86.2, 114.7, 114.7) | −5.55, 20 |

Each cell gives λ_min and the ghost count, with ghost frequencies in brackets where listed.

- **ζ: no ghosts.** This holds at every width up to L = 8 and every resolved frequency up to 125. λ_min stays at the float floor, which at M = 327 is 1e-10 relative (NUMERICAL).
- **Epstein: exactly 2 ghosts per off-line quartet, sitting at the quartet.**
  - The quartets inside the basis band are at 15.67, 29.98, 36.37, 44.00, 46.76 (test data). The ghost pairs are at 14.5–16.5, 30.5–31.2, 35.5–37.0, 43.0–43.5, 45.9–46.1.
  - The count 10 matches 5 quartets × 2.
  - Rayleigh parts of the ghosts (L = 8, first ghost): arch +4.9, pole +0.6, prime +34.5. **The prime term creates the ghost.**
- **DH: exactly 2 ghosts per quartet.**
  - The quartets below 125 are at 85.70 and 114.16 (test data). Ghosts appear at 84.5–86.6 from L = 4 and at 113.4–115.3 from L = 4.5.
  - The count 4 is stable from L = 4.5 to 8.
  - The arch is ≈ +4.2, the pole is 0 (DH has none), and the prime term is larger.
- **F_θ (control): 2 ghosts per ζ zero.** Every zero is off the line by θ.
  - The count saturates at 20 = 2 × 10 ζ ordinates inside the basis band (14.1 … 49.8).
  - This happens despite an Euler product with **positive** weights and a functional equation.
- **Spliced Beurling: ghost proliferation, not a ghost count.**
  - It is identical to ζ up to L = 4, with 0.0 difference (PROVED, splice at 55 > e⁴).
  - From L = 4.5 the ghost count is 8, 16, 30, 42, 52, growing with L. The frequencies are spread out: the quartet at 28 shows up only as one pair (27.4) among many.
  - Reason: with no functional equation, the borrowed arch is not BEUR″'s own explicit formula. The pairing is **not** Σ_ρ over BEUR″'s zeros, and the defect is the whole arch mismatch.

## 4. The literal Osterwalder–Schrader version is false for ζ (PROVED + NUMERICAL, `noghosts_B.json`)

- **The Euclidean form.** OS reflection positivity uses the bilinear form ⟨Θf, f⟩_OS = ∫∫ f(u) f(v) k(u + v) du dv. The test functions are supported at positive time, u, v ∈ (0, T). Here Θ is the reflection u ↦ −u without complex conjugation, and k is the same integer-built kernel.
- **PROVED: for 2T < log 2, ζ has at most 2 positive directions and everything else is a ghost.** No prime atom enters at that size, so the kernel is −w(u+v) + 2cosh((u+v)/2).
  - w(s) = Σ_{k≥0} e^{−(½+2k)s} is completely monotone. So ∫∫ f f w(u+v) = Σ_k (∫ f e^{−(½+2k)u})² ≥ 0, and −w is negative semidefinite.
  - 2cosh((u+v)/2) = e^{u/2}e^{v/2} + e^{−u/2}e^{−v/2} has rank 2 and is positive.
  - Numerics at T = 0.3: 7 negative, 1 positive.
- **NUMERICAL, ζ:** at T = 0.5 / 1 / 2 / 3 / 4 the counts (negative, positive) are 10/6, 12/11, 13/16, 12/15, 13/15. EPS and DH look the same (EPS 32/33 at T = 4, DH 24/27).
- **Reason.** On-line zeros are *real frequencies*, e^{iγu}. In a bilinear Hankel form, cos γ(u+v) = cos γu cos γv − sin γu sin γv has signature (1, 1). OS positivity needs real decay e^{−E(u+v)}, and that describes the archimedean and pole terms, not the zeros.
- **Physics reading (Wick rotation).** The zeros live on the real-time (Minkowski) side. The correct reflection is the Weil Θ = x ↦ 1/x **together with complex conjugation**, i.e. the antiunitary time reversal. Its positivity is Bochner / Wightman positivity, and that is §1–§3.
- **Verdict.** "No ghosts" in the OS-Euclidean sense is **WRONG** as a statement about ζ. In the Weil-reflected sense it is **exactly** W ≥ 0.

## 5. Which physical axiom each fake breaks ("same time")

- **The single state with both properties at once (PROVED, Tate 1950 / Weil representation).** It is f₀ = e^{−πx²} ⊗ ∏_p 1_{ℤ_p}.
  - At ∞, the Gaussian is the ground state of the harmonic oscillator. The Fourier transform is the oscillator's **quarter-period time evolution** e^{−iπN/2}, so the ground state is Fourier-fixed.
  - At each p, 1_{ℤ_p} is the unramified ("ground") vector, and it is fixed by the p-adic Fourier transform.
  - **Euler product** means the state is a tensor product over places.
  - **Poisson** means the state is invariant under the global Fourier transform, i.e. the same quarter-period "time" applied at every place simultaneously.
  - Its Mellin transform is Λ(s). This is Daniel's "same time": the two properties are one property, "ground state of the oscillator at every place at once".
- **Fakes, by the axiom they break, with the ghost signature observed:**

| model | product over primes (Euler) | Θ-symmetric (Poisson / FE) | constituents unitary | pure (one product state) | ghosts (§3) |
|---|---|---|---|---|---|
| ζ | yes | yes | yes | yes | none (L ≤ 8, f ≤ 125) |
| Beurling | yes | **no** (no additive group, no Fourier, no own arch) | yes | yes | proliferate after the splice; the pairing is not even its own Weil form |
| Epstein x²+5y² | no | yes (ℤ² Poisson) | yes | **no**: Z_Q = ζL(χ₋₂₀) + L(χ₋₄)L(χ₅) is a superposition of two product states (the genus average, Siegel–Weil, is the product) | 2 per quartet, prime-driven |
| DH | no | yes (exact Tate integral, 14A) | yes | **no**: a superposition of the χ and χ̄ product states, mixed at p = 5 | 2 per quartet (85.7, 114.2) |
| F_θ | **yes** (Λ-weights ≥ 0) | **yes** (degree 2) | **no**: e^{±θu} is a Euclidean-time shift, i.e. the quasi-characters \|·\|^{±θ} | yes | 2 per ζ zero |

- **F_θ is the decisive control.** It has an Euler product with positive weights and a functional equation, and it still has ghosts. So "Euler product + Poisson" alone, as ingredients, does not exclude ghosts. A source has to use **pure + unitary + degree-one**. That is the Selberg-class degree-1 situation: Kaczorowski–Perelli (Acta Math. 182, 1999) show that the degree-1 functions are exactly ζ(s) and L(s + iθ, χ). This **selects** the function. It is not a positivity mechanism (CHECKED, literature; consistent with 14A §7).
- **GAP.** None of the axioms in the table is a positivity statement about the pairing. Every one of them is satisfied by ζ, and for each axiom there is a fake whose violation produces ghosts. That shows the axioms are *necessary* for no ghosts on our test cases. It does not show they are *sufficient*. Sufficiency is RH.

## 6. Relation to the 2025 draft

The 2025 draft (ch. 8) has two versions of a no-ghost argument. Each has one gap, and the construction here is designed to avoid both.

- **(i) The adelic version.**
  - The gap: it takes Im ρ to be an L² eigenvalue of i d/d log|a| on L²(𝔸^×/ℚ^×). On that space the operator has purely continuous spectrum. Also, the (non-)normalisability of |a|^{σ−1/2} is kinematic: it comes from Haar measure and holds equally for Beurling and Epstein (13i, 14A).
  - Here the Hilbert space is not L²(C_ℚ). It is the GNS space H_Θ of the Weil-reflected pairing. In H_Θ the zeros *are* eigenvalues (§2(i)), but only once its positivity is established. So positivity is the stated hypothesis, not something inherited from Haar measure.
- **(ii) The CFT version.**
  - The gap: Q[f] = ‖Tf‖² is positive by construction (a Gram form), but it is never shown to equal the Weil form. Step 5 also sets ⟨ψ₊|ψ₋⟩ = 0.
  - Here the pairing is shown to equal W **exactly**, through the explicit formula with no defect term (§1). It is not positive by construction.
  - The computation in §2(ii) shows that the off-line ghost is precisely the cross term between a state and its reflection: ⟨ψ_ρ|ψ_ρ⟩ = 0 and ⟨ψ_ρ|ψ_{1−ρ̄}⟩ ≠ 0, forming a hyperbolic plane. So the orthogonality assumption is exactly where an off-line zero would hide. The numerics confirm it: EPS and DH ghosts come in pairs at their quartets.

## 7. Labels (Part 1)

- Pairing = W, an exact identity from integer data (prime side + own arch + poles): PROVED (explicit formula; Weil 1952, Bombieri 2000). Galerkin vs the 13e harness: CHECKED (≤ 4e-14).
- No ghosts ⇒ real spectrum = zeros, pure point in H_Θ; No-Ghosts ⇔ RH: PROVED (standard).
- Off-line quartet = signature (2, 2) cross term; ghost index ≤ 2 × #quartets: PROVED. Equality at the quartet heights for EPS (10 = 5 × 2) and DH (4 = 2 × 2): NUMERICAL.
- ζ ghost-free for L ≤ 8, f ≤ 125: NUMERICAL.
- BEUR″ ghost proliferation after the splice: NUMERICAL (identity to L = 4 PROVED).
- Euclidean OS form indefinite for ζ: PROVED (2T < log 2) and NUMERICAL.
- Tate vector = oscillator ground state at every place: PROVED (standard). Axiom table: CHECKED (structural) + NUMERICAL ghosts. Sufficiency of the axioms: GAP.
- **No-Ghosts for ζ: OPEN** (it is RH). Nothing here closes RH, ARP or H_ψ.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5). Audit: Turbo.
