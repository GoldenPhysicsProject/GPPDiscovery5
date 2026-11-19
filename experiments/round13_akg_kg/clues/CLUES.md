# Clue hunt: hints that might show the way (Round 14 → clues)

Heart of the Cosmos (Grok / GPPDiscovery5), 2026-11-18 (Toronto). Daniel asked for clues "in anything": hints that point somewhere, not more no-go theorems.

**Labels:** PROVED / CHECKED / NUMERICAL / CONJECTURE / HEURISTIC.
**Status:** RH, ARP and H_ψ are all **OPEN**. Nothing below proves RH. Zeros enter only as test data, and only in probes A and C.

**Markers:**
- ★ = separates ζ from **every** fake we use: Epstein x²+5y² (EPS), Davenport–Heilbronn (DH), spliced Beurling (BEUR″), and the non-unitary control F_θ. Any proof has to use a ★ clue somewhere.
- ◆ = echoes one of Daniel's ideas: integers fundamental, add/mult coupling, "same time", x→1/x, no actual infinity, s = Δ.

**Probes run for this file** (all in this folder):

| probe | files | what it tests |
|---|---|---|
| A | `hko_interlace.py`, `hko_runs.py` → `hko_interlace.json`, `hko_runs.json` | interlacing of the two pure pieces of EPS and DH |
| B | `ff_epstein.py` → `ff_epstein.json` | an exactly solvable, fully finite Epstein in a function field |
| C | `lagarias_pos.py` → `lagarias_pos.json` / `.log` | Re Λ′/Λ > 0 to the right of the line, for ζ, EPS, DH |

**Our own rounds, cited below:**

| round | where | commit |
|---|---|---|
| 13 | `../README.md` | main |
| 14B | `../round14b_euler/ROUND14B.md` | `ea7490f` |
| 14C | `../round14c_prolate/ROUND14C.md` | `b581f8c` |
| 14D | `../round14d_global/ROUND14D.md` | `7e89769` |
| 14E | `../round14e_gap/ROUND14E.md` | `cbcf08a` |
| 14E2 | `../round14e_gap/ROUND14E2.md` | `24e17c3` |
| 14F | team/principal-series `principal_series/ROUND14F.md` | `623cbc9` |
| Golden Eye audit | team/golden-eye `golden_eye/AUDIT14.md` | `eae661d` |

---

## Ranked clues

### 1. ★◆ The one proven RH gets its positivity from a norm form in two integers: deg(m − nF) ≥ 0. Our finite toy reproduces the Epstein failure exactly.

**Observation (classical).** For an elliptic curve E/F_q with Frobenius F, the degree map is a positive-definite quadratic form on ℤ[F]:

  deg(m − nF) = m² − a·mn + q·n² ≥ 0 for all integers m, n, where a = q + 1 − #E(F_q).

A positive binary form has discriminant ≤ 0, so a² ≤ 4q. That is Hasse's RH. Weil's general proof is the same thing on C×C: the Hodge index theorem, i.e. positivity of the Rosati involution, which encodes F·F† = q. Bombieri–Stepanov give a proof of the curve RH that uses only finite counting with auxiliary polynomials.

**New probe B** (`ff_epstein.py`; PROVED by Riemann–Roch, CHECKED by brute force):
- Setting: A = F_q[E − O] has class group Cl(A) = E(F_q), with h classes. The ideal-class partial zetas, the exact analogue of Epstein, come out in closed form:
  - Z₀(T) = (1 − qT + qT²)/(1 − qT);
  - Z_c(T) = T/(1 − qT) for c ≠ 0;
  - Σ_c Z_c = (1 − aT + qT²)/(1 − qT) = ζ_A.
- Brute-force degree-1 and degree-2 counts agree with these formulas on 6 curves, p = 5, 7, 11, 13, 17, 19.
- Z₀ satisfies the functional equation (the numerator is palindromic under T → 1/(qT)) but has no Euler product. For q ≥ 4 its zeros are real, with moduli e.g. 0.276 and 0.724 against the RH circle 0.447 at q = 5. So they are **off the circle**.
- The class sum restores the Euler product, and its zeros are on the circle exactly by Hasse.
- Pencil 1 − qT + qT² + λ(h − 1)T: the zeros stay on the circle iff λ ∈ [(q − 2√q)/(h − 1), (q + 2√q)/(h − 1)], which contains λ = 1. At the window's edge the two circle zeros collide at T = ±q^{−1/2} and leave the circle. This is the exact, solvable version of 14F's collision λ* ∈ (0.1400, 0.1425) at height 15.24.

**Why it might point the way.**
- In the one case where RH is proved, positivity is not attached to the zeta function as a functional. It is a positive form on an **integer lattice of operators**, ℤ·1 + ℤ·F: additive integer combinations of the identity and the multiplicative (Frobenius) operator, measured at the same time by a multiplicative and positive function, the degree.
- Mapping the pieces:
  - the FE T → 1/(qT) is F ↦ q/F, i.e. Rosati duality, i.e. Daniel's x → 1/x;
  - the Euler product is the class average;
  - unitarity is |eigenvalues of F/√q| = 1, which is exactly what the non-unitary control F_θ breaks.
- The missing object for ℤ is the "ℤ[F]" on which the Weil form is a degree.

**Echoes:** integers fundamental; add/mult at the same time; x → 1/x; no infinity (everything here is finite).

**Cheap next test.**
1. Form the "Hasse Gram" for ℤ: the Weil form W on span{g, τ_p g, τ_p² g, …}, where τ_p is the shift by log p with weight p^{−1/2} and g is a bump at L = 1–3. Compare ζ / EPS / DH / BEUR″ / F_θ with 14B's `weil_galerkin.py`.
2. Ask whether the 2×2 and 3×3 determinants already separate the fakes, and whether ζ's matrix looks like m² − a mn + q n² with |a| ≤ 2√q, a "local Hasse bound at p".

---

### 2. ★◆ The axioms that separate ζ from every fake: Euler product AND Poisson on ℤ AND unitary AND pure. In degree 1 these force the function to be ζ.

**Observation.**
- 14B §5 (`ROUND14B.md` lines 120–138, `ea7490f`): each fake breaks exactly one leg.

| model | breaks | lacks |
|---|---|---|
| EPS | purity | Euler product: a(6) = 2 ≠ a(2)a(3) = 0 (14F) |
| DH | purity | Euler product: 906 non-prime-power atoms of −F′/F up to n ≤ 3000 (14F) |
| BEUR″ | Poisson / FE | its "integers" are not an additive group, so there is no Fourier transform and no arch of its own |
| F_θ | unitarity | quasi-characters \|·\|^{±θ}; it keeps the Euler product with positive weights and the FE, and still has ghosts |

- Rigidity:
  - Hamburger (Math. Z. 10, 1921): a Dirichlet series over ℕ with ζ's functional equation (plus mild growth) is c·ζ.
  - Kaczorowski–Perelli (Acta Math. 182, 1999): the Selberg class in degree 1 is exactly ζ(s) and L(s + iθ, χ).
- PROVED (literature, plus 14B/14F checks).

**Why it might point the way.**
- Only ζ (and its pure, unitary relatives) satisfies all four legs. So a positivity proof must use them **inside one inequality**, not one after another.
- Hamburger says the additive face (theta and Poisson on ℤ) alone already pins ζ down. The Euler product is not extra input: it follows from ℤ being ℤ.
- This is Daniel's "integers, not primes, are fundamental" in theorem form. Every fake with an FE (EPS, DH) is a superposition of pure Euler states (14F). The one fake with an Euler product but no Poisson (BEUR″) has fake integers.
- The single state that has all four legs is Tate's vector f₀ = e^{−πx²} ⊗ ∏_p 1_{ℤ_p}, the oscillator ground state at every place at once (14B §5).

**Echoes:** integers fundamental; same time; x → 1/x.

**Cheap next test.**
- Run three one-leg homotopies from ζ, one per leg:
  - Euler leg: ζ_K + λL_ψ, from 14F;
  - unitarity leg: ζ(s + θ)ζ(s − θ) with θ going from imaginary to real;
  - Poisson leg: the BEUR″ splice point x₀ moving down from ∞.
- For each, record the λ (or θ, or x₀) at which the first ghost is born, as a function of height. A strictly positive onset threshold at every height would say positivity is *robust* along that leg, so a perturbative/stability proof might work there.
- 14F already shows a finite threshold (0.14) on the Euler leg at height 15.

---

### 3. ◆ NUMERICAL, new: ghosts are born exactly where the two pure pieces' zeros get out of step (runs of ≥ 3 same-piece zeros)

**Observation, probe A** (`hko_interlace.py`, `hko_runs.py`):
- The two pencil fakes split into real, FE-symmetric pure pieces on the line:
  - EPS = A + B, with A = Λ(ζ·L(χ₋₂₀)) and B = Λ(L(χ₋₄)L(χ₅));
  - DH = A + B, with A = Z_χ and B = Z_χ̄. The pencil ratio comes out exactly r = 1.000 (imaginary part 2e-31).
- Hermite–Kakeya–Obreschkoff: A + λB is real-rooted for every real λ iff the zeros of A and B interlace.
- Single interlacing defects are common (EPS: 41 below t = 60; DH: 38 below t = 130) and create no ghost.
- **Every off-line zero sits at a run of ≥ 3 consecutive zeros of one piece with none of the other.**
  - EPS: 7/7 off-line zeros (15.67, 29.98, 36.37, 44.00, 46.76, 52.74, 58.51) lie within 0.45 of such a run.
  - DH: 2/2 (85.70, 114.16) lie within 0.30.
  - Conversely, all 10 EPS runs and both DH runs lie within 1.5 of an off-line zero.
- The fakes' own zeros are confirmed: |EPS(ρ)| = 8e-11 and |DH(ρ)| = 3e-11 at the 14F zeros.

**Why it might point the way.**
- A ghost is a local **phase slip** between two pure (Euler) pieces. It is a Krein collision of two modes of opposite signature, the same mechanism as PT-symmetry breaking (Bender–Brody–Müller, PRL 118, 2017). It matches 14B's "(2,2) cross term" and 14F's collision homotopy.
- ζ is a single pure piece, so it has no partner to slip against. The de Branges / Hermite–Biehler form of No-Ghosts for ζ is that the phase of ζ's own E-function is monotone.
- It is equivalent to RH, but it is a **local, sign-counting** statement, which is the kind finite calculations can attack.

**Label:** NUMERICAL (HKO is PROVED; the run ↔ ghost correspondence is CONJECTURE).

**Echoes:** No-Ghosts as unbroken PT symmetry; "ghosts are interference" (14F).

**Cheap next test.**
1. Extend to t ≤ 300 (DH) and t ≤ 150 (EPS) and see whether 1-to-1 survives.
2. Prove the easy half: a run of ≥ 3 forces a sign-change deficit for A + B in that window, by a Sturm/argument-principle count.
3. Look for the ζ analogue: split e^{iθ(t)}ζ(σ + it) into Re/Im at σ = 0.5 + ε and test interlacing (Lagarias's framework, see clue 9).

---

### 4. ★◆ ζ sits exactly at a critical point: de Bruijn–Newman 0 ≤ Λ ≤ 0.2, with RH ⇔ Λ = 0. Lee–Yang is where this question came from.

**Observation.**
- Rodgers–Tao (Forum Math. Pi 2020, arXiv:1801.05914) prove Λ ≥ 0. Polymath15 (arXiv:1904.12438) gives Λ ≤ 0.22, and Platt–Trudgian (Bull. LMS 2021) give Λ ≤ 0.2.
- Newman (Proc. AMS 61, 1976) came to this from the Lee–Yang circle theorem (Phys. Rev. 87, 1952). In Lee–Yang, ferromagnetic (positive) couplings put the partition-function zeros on the unit circle.
- "RH, if true, is only barely so."
- The fakes have off-line zeros, so Λ_F > 0 strictly. This separation is tautological.

**Why it might point the way.**
- Physicist reading: ζ is the critical member of the heat-flow family H_t. A critical point is scale/conformally invariant, which is Daniel's s = Δ (conformal dimension), and that is where Lee–Yang zeros pinch the real axis.
- The live question is the Lee–Yang one: which positive ("ferromagnetic") structure of Φ(u) = Σ(2π²n⁴e^{9u/2} − 3πn²e^{5u/2})e^{−πn²e^{2u}} pins the zeros? Note that Φ is built from the integers n, via theta.
- Pólya's Φ > 0 holds for EPS too (its theta is positive), so plain positivity of Φ is fake-blind (HEURISTIC).

**Label:** PROVED anchors; HEURISTIC reading.

**Echoes:** s = Δ; physicist-first.

**Cheap next test.** Run the backward-heat zero dynamics dz_j/dt = −2Σ_k 1/(z_j − z_k) on EPS's and DH's first quartets (test data). Get the time t*(quartet) at which each pair reaches the real line, and compare t* with the 14F mixing thresholds λ*. If t* ∝ (1 − λ*)², the "mixing" and "heat" pictures are one picture.

---

### 5. ◆ s = Δ holds literally in Bost–Connes and modular theory: the modular operator is n^{−β}, the partition function is ζ(β), and Tomita–Takesaki supplies a built-in positivity with a "½"

**Observation.**
- Bost–Connes (Selecta Math. 1, 1995): H ε_n = (log n) ε_n, and the Gibbs/modular operator is e^{−βH} = n^{−β}, so Tr = ζ(β).
- There is a phase transition at β = 1: a unique KMS state for β ≤ 1, and symmetry breaking by Ẑ^× = Gal(ℚ^ab/ℚ) for β > 1.
- The primon gas (Julia 1990; Spector, CMP 1990): the Euler product is a free Bose gas of primes, which is Daniel's "Euler product as particles".
- In Tomita–Takesaki theory, the modular operator Δ is **positive**.
  - The positive cone is built with Δ^{1/4}, and the KMS/reflection positivity reads ⟨aξ, J aξ⟩ ≥ 0 (self-polar cone), equivalently ⟨aξ, Δ^{1/2} aξ⟩ = ‖Δ^{1/4}aξ‖² ≥ 0.
  - By Bisognano–Wichmann, J is the wedge reflection, which is the origin of OS reflection positivity.
- Daniel's draft line 749 makes the same connection on the physics side: the Plancherel weight πλ/sinh(πλ) is the boost oscillator at the Unruh temperature (Bisognano–Wichmann).

**Why it might point the way.**
- The "½" in Re s = ½ may be the Δ^{1/2} of modular theory.
- Weil positivity for ζ might be a modular-positivity statement ⟨ξ_f, Jξ_f⟩ ≥ 0 for some state whose J is Daniel's Θ f(x) = x^{−1} f̄(1/x).
- Positivity is free in modular theory once the state exists (it is automatic for any cyclic separating vector). So the whole problem would move to constructing the right von Neumann algebra and state.

**Label:** HEURISTIC; the BC and Tomita facts are PROVED.

**Separation:** partial. BEUR″ also has a free gas, so a Beurling BC system exists; the clue must be combined with clue 2's Poisson leg.

**Echoes:** s = Δ, x → 1/x, Euler as particles, no ghosts.

**Cheap next test.** On BC truncated to n ≤ N, take the KMS_β state at β = 1 and β = ½ and compute ⟨aξ, J aξ⟩ on the 14B Galerkin span (finite matrices). Compare with the 14B Θ-pairing matrix. A match of the archimedean term, or a failure, says whether Θ is a modular conjugation.

---

### 6. ◆ Positivity is proved until the first integer > 1 enters, and then the margin decays like exp(−13·e^L), i.e. by the count of integers, not primes

**Observation.**
- Connes–Consani (Selecta 27, 2021, arXiv:2006.13771) and Yoshida (1992) prove Weil positivity for support width L ≤ log 2, where no prime atom enters.
- 14C confirms the identity is exact until the n = 2 atom at log 2 (`ROUND14C.md` lines 45, 79).
- 14E2 and Golden Eye's 80-digit data: λ_min(W_L) = exp(a + bL − B e^L) with B = 13.13 ≈ 4π (`ROUND14E2.md` line 13, `24e17c3`). Plain exponential and Gaussian fits are rejected.

**Why it might point the way.**
- e^L is the number of integers in the window, not π(e^L). The positivity margin behaves like exp(−4π·#integers), which is the Slepian/prolate law: eigenvalues die beyond the time–bandwidth (Shannon) number.
- So the margin's decay is kinematic, from the degrees of freedom of the integer lattice. It is not prime-driven.
- If that holds, the prime term only has to fit inside a kinematically shrinking but **never negative** window. That is a statement about integers, which is Daniel's view.

**Label:** NUMERICAL fit (5 converged points, L ≤ 2); HEURISTIC reading.

**Echoes:** integers fundamental; no infinity (finite windows).

**Cheap next test.**
1. Refit the existing 80-digit λ_min against ⌊e^L⌋, π(e^L) and ψ(e^L), adding L points just below and above log 3, log 4, log 5.
2. Check whether λ_min jumps at every integer, only at prime powers, or at neither.
3. Do the same for the fakes' positive parts.

---

### 7. ★(tautological)◆ Finite average energy: RH ⇔ the prime-built field has finite Cesàro energy

**Observation.** Round 13 R13.2 (`../README.md` line 32, CHECKED): Φ_g, built from primes, the pole and the arch only, lies in B² (finite Cesàro mean square) for all g ⇔ RH. Its B² norm equals the Weil form to ≤ 2e-5 at L = 0.5, 1, 2.

**Why it might point the way.**
- RH becomes "the arithmetic field has finite energy". That is a physicist's normalizability statement and Daniel's "no actual infinity": an off-line zero is an exponentially growing mode, a runaway.
- In Pasterski–Shao, finiteness comes from the KG Cauchy problem (energy conservation). So the real clue is to look for the arithmetic **equation of motion and conserved energy** of Φ_g. `../o9_conservation.py` started this.

**Label:** CHECKED. It is RH-equivalent, so the separation is tautological.

**Cheap next test.**
1. Measure the growth exponent of the fakes' fields at u ≤ 13 from their own prime sides. EPS should show e^{0.87u}.
2. Search for a quadratic functional E[Φ_g](u) that is non-increasing in u for ζ at the numerical level (Lyapunov/energy). If one exists, test whether it fails exactly for the fakes.

---

### 8. ◆ Nyman–Beurling / Báez-Duarte: RH is a statement about the floor function, where addition and multiplication meet

**Observation.**
- RH ⇔ 1 is in the closed span of ρ_θ(x) = {θ/x} (Nyman 1950, Beurling 1955).
- Báez-Duarte (Rend. Lincei 14, 2003) restricts to integer dilations θ = 1/n.
- Burnol (Adv. Math. 170, 2002) gives an unconditional lower bound d_N² ≥ (Σ_{ρ on line} m_ρ²/|ρ|² + o(1))/log N. The conjecture is equality.
- Repo: `../nyman_beurling.py` / `.json` (Round 13).

**Why it might point the way.**
- {x} = x − ⌊x⌋ is the additive structure of ℤ (the floor) acting under multiplicative dilations. The criterion exists only because ℤ is an additive group.
- BEUR″ has no floor function, so no criterion exists for it. EPS and DH need character twists.
- This is the most literal "add/mult coupling" form of RH, and it is finite N at every step.

**Label:** PROVED (literature); HEURISTIC that this is the right lever.

**Echoes:** integers fundamental; add/mult; no infinity.

**Cheap next test.** Compute Báez-Duarte d_N² log N for N ≤ 2·10⁴ (Möbius-weighted Gram, cheap) and compare with Σ1/|ρ|² = 2 + γ − log 4π = 0.0461. Then compute the DH analogue with χ-twisted fractional parts, to see whether the fake's ghost shows up as a deficit in the constant.

---

### 9. Lagarias positivity Re Λ′/Λ > 0 to the right of the line is a sharp, local ghost detector (probe C)

**Observation.**
- Lagarias (Acta Arith. 89, 1999): RH ⇔ Re ξ′/ξ(s) > 0 for Re s > ½.
- Probe C (`lagarias_pos.py`, t ∈ [2, 120], step 0.25):
  - ζ: positive everywhere sampled, with min 0.0024 at σ = 0.55.
  - EPS: negative only in thin bands at its off-line zeros, for σ ≤ 0.9.
  - DH: negative only at t ≈ 85.5–85.75 and 114.25, for σ ≤ 0.8.
  - For σ ≥ 1 all three are positive. DH does have zeros in σ > 1 (Davenport–Heilbronn 1936), but none show up in this sample (t ≤ 120).

**Why it might point the way.** It is a pointwise positivity condition: no Gram matrix and no test functions. The failure is localized, which matches clue 3. For σ > 1 the Euler product controls the prime part of Re ζ′/ζ, and the FE mirrors it, so the RH strip is exactly where the two legs must act at the same time.

**Label:** NUMERICAL; the equivalence is PROVED.

**Echoes:** same time.

**Cheap next test.** For ζ, bound Re ξ′/ξ(σ + it) from below at σ = 1 + δ using only the Euler product and Stirling. Track how far left an Euler-only bound reaches before it must borrow the FE.

---

### 10. Scattering: zeros as resonances, and RH as "no resonance on the wrong sheet"

**Observation.**
- Pavlov–Faddeev (1972) and Lax–Phillips (Ann. Math. Stud. 87, 1976): for the automorphic wave equation on SL₂(ℤ)\ℍ, the scattering matrix is φ(s) = ξ(2s − 1)/ξ(2s). The zeta zeros are its resonances.
- Connes–Consani ζ-cycles (arXiv:2106.01715) and CCM prolate (arXiv:2310.18423) realize the **critical** zeros at every circle length L. Off-line zeros could only be resonances, as GPT's note and 14F §2c both state.

**Why it might point the way.** Physicists exclude wrong-sheet resonances through unitarity plus causality (analyticity of S), not through a positive norm on a state space. The SL₂(ℤ) case is attached to ζ itself, the pure object. The congruence/class-number versions bring in Epstein-type mixtures. This is Daniel's celestial/KG physics, applied to the right object.

**Label:** PROVED (the scattering identity); HEURISTIC (the route).

**Echoes:** physicist-first; unitarity / no ghosts.

**Cheap next test.** Compare the Lax–Phillips cut-off semigroup's generator spectrum at small cut-off with our 14C/14E finite-L λ_min ladders. Is exp(−13e^L) the Lax–Phillips decay rate?

---

### 11. Where not to look: statistical and asymptotic signatures cannot see ghosts

**Observation.**
- Bombieri–Hejhal (Duke 80, 1995): linear combinations of Euler products have almost all zeros on the line, under GRH plus a spacing hypothesis. Yet DH-type combinations have infinitely many off-line zeros, about T of them against T log T in total.
- Montgomery/GUE statistics, smooth counting (Berry–Keating xp, GPT's Lambert-W carrier), and Griffin–Ono–Rolen–Zagier asymptotic Jensen hyperbolicity (PNAS 116, 2019, arXiv:1902.07321) are all asymptotic or density statements. HEURISTIC: they cannot distinguish ζ from the fakes.
- 14B ghost counts are small and local (EPS 10, DH 4 at L ≤ 8).

**Why it might point the way.** It prunes the search: ghosts form a density-zero, local phenomenon. Any proof must be exact and local in height (clues 3, 9), not statistical. This matches Daniel's "finite calculations".

**Label:** PROVED (BH, DH 1936) + HEURISTIC.

**Echoes:** no infinity.

**Cheap next test.** Compute GORZ Jensen polynomials J^{d,n} for EPS at d ≤ 6, n ≤ 60, and find the first non-hyperbolic (d, n). This confirms the fakes fail only at low n, i.e. locally.

---

### 12. ◆ Add/mult uncertainty: ℤ cannot be structured both ways at once, except through the self-dual Gaussian

**Observation.**
- Results limiting simultaneous structure:
  - sum-product: Erdős–Szemerédi (1983), and Bourgain–Katz–Tao (GAFA 2004) over F_p;
  - ×2×3 rigidity: Furstenberg (1967), Rudolph (ETDS 1990);
  - Bondarenko–Radchenko–Seip (arXiv:2005.02996), whose Fourier interpolation formulas have nodes tied to log n and to zeta zeros.
- All of these say additive and multiplicative structure cannot both be concentrated at once.
- The Tate vector (clue 2) is the minimum-uncertainty (Gaussian) state between them.

**Why it might point the way.** Uncertainty principles are positivity statements: ⟨ψ, A†Aψ⟩ ≥ 0 gives ΔxΔp ≥ ½. If the Weil form is a Robertson–Schrödinger expression in an "additive" operator (Fourier on ℤ) and a "multiplicative" one (dilation / log n), then its positivity would be an uncertainty principle for ℤ. That is the strongest physical form of Daniel's coupling idea.

**Label:** HEURISTIC; the cited theorems are PROVED.

**Echoes:** add/mult coupling, "same time", integers fundamental.

**Cheap next test.** On the 14B Galerkin basis, compute the commutator [X, D], with X = Poisson/theta projection and D = −i d/du on log-scale. Test whether W − c·(Var X · Var D − ¼|⟨[X, D]⟩|²) has a sign for some c, for ζ vs EPS.

---

## Which clues separate ζ from all the fakes

- **Non-tautological ★ (the valuable ones):**
  - **Clue 2:** the conjunction Euler ∧ Poisson-on-ℤ ∧ unitary ∧ pure, plus Hamburger / Kaczorowski–Perelli rigidity.
  - **Clue 1:** the norm-form mechanism. It uses all of clue 2's legs together and is the only proven positivity source we know. Our finite toy shows it fails, in exactly the Epstein way, when the Euler leg is removed.
- **Tautological ★ (RH-equivalent detectors, useful for locating ghosts, not for proving):** clues 4 (Λ), 7 (finite energy), 9 (Lagarias).
- **Partial:** clue 3 (EPS, DH; not BEUR″, which is not a pencil), clue 5 (not BEUR″), clue 8 (BEUR″ has no floor function; EPS/DH untested).

## What I'd do next (one round, all cheap)

1. **Hasse Gram for ℤ** (clue 1): W on span{τ_{p^k} g} for ζ and all fakes.
2. **One-leg onset thresholds** (clue 2): ghost-birth λ* against height along the three homotopies.
3. **Run ↔ ghost** to t = 300 (clue 3), plus the Sturm-count proof of the easy half.

Item 1 is the most direct test of whether the one proven mechanism has a shadow over ℤ.

RH / ARP / H_ψ: **OPEN**.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).
