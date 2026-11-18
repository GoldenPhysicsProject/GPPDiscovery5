# Round 14F — representation theory of the Connes quotient: isobaric sums vs linear combinations

Lane: `experiments/round13_akg_kg/principal_series/`, branch `team/principal-series`.

Inputs:
- `../round14d_global/ROUND14D.md` (main 7e89769);
- `../round14b_euler/ROUND14B.md`;
- Golden Eye `../golden_eye/AUDIT14.md`: pre-audit at team/golden-eye 062e4e2, audit §5 at 74740e6.

14E is not on main at 7e89769.

Code (each script is self-contained, so the outputs do not depend on run order):
- `r14f_quotient.py` → `r14f_quotient.json`: integer identities, DH, interference, homotopy, collision.
- `r14f_logderiv.py` → `r14f_logderiv.json`: atoms of −F′/F.
- `r14f_twist.py` → `r14f_twist.json`: 14B Galerkin for F_θ vs F_{iτ}.

Zeros are test data only. Labels: PROVED / CHECKED / NUMERICAL / CONJECTURE / GAP / OPEN.

> **Revision (Golden Eye audit §5, 74740e6).** Fix commit on top of 104719a.
> 1. §2a: W_{FG} = W_F + W_G is relabelled as an *identity* (log(FG) = log F + log G). Its float check is not content.
> 2. Item 5 is restated as the Tate/Hecke containment fact. The Connes-quotient wording is dropped as vacuous. "Equal moduli" at the off-line zeros is marked as forced by |β/α| = 1. GL_n with n ≥ 2 now cites Weil's criterion.
> 3. Added the elementary Epstein non-Euler proof a(6) = 2 ≠ a(2)a(3) = 0. The ζ / ζ_K ghost entries are now code checks, and Beurling is N.A. / defect. Temperedness is separated from unitarity.
> 4. Added a seventh row, F_{iτ} (τ = 3), the unitary-twist control.
> 5. The collision bracket is tightened to λ* ∈ (0.1400, 0.1425) at height ≈ 15.24. The script bug that tested r = 2(b₁ + b₂) is fixed, `r14f_patch.py` and `r14f_collision.py` are folded into `r14f_quotient.py` and removed, and all JSONs are regenerated.

**Bottom line.**
- "One irreducible vs a sum of two" is the wrong dichotomy.
- The working dividing line has two parts:
  - **products** (isobaric sums) vs **linear combinations** of Dirichlet series, which is the same as Hecke-equivariant vs non-equivariant test vectors;
  - **unitarity** of the constituents.
- The candidate criterion: F = L(s, Π) with Π = π₁ ⊞ … ⊞ π_r an isobaric sum of unitary cuspidal automorphic representations of GL_{n_i}(𝔸_ℚ). This is the classical domain in which GRH is conjectured (Langlands; Selberg class), not a new criterion.
- It sorts all seven test cases. It **selects** objects. It does not prove positivity: W ≥ 0 on that class is GRH, which is CONJECTURE / OPEN.
- The new content here is expository and numerical: the interference picture, the collision homotopy, and the unitary/non-unitary twist pair.

---

## 0. Which group, which "irreducible"

- **The scaling group cannot carry the distinction.** ℝ_{>0} and C_ℚ are abelian, so their irreducibles are characters, and ζ's quotient is already a sum of characters |·|^{iγ}. "Single irreducible" is vacuous there. (PROVED, trivial.)
- **The non-abelian setting.** Automorphic representations π of GL_n(𝔸_ℚ), with the Hecke algebra ℋ = ⊗′_p ℋ_p acting (spherical at unramified p, and GL_n(ℤ_p) or units at ramified p).
  - n = 1 gives Dirichlet / Hecke characters.
  - **Epstein** lives on GL₂: θ_{Q₁} via the theta lift / automorphic induction from GL₁(𝔸_K), K = ℚ(√−5).
  - **DH** lives on GL₁(𝔸_ℚ), mixed in the ramified component at 5.
- **"Irreducible / eigen" is defined algebraically**, since for the fakes W_F is indefinite and the quotient is at best a Krein space.
  - The test vector is ℋ-eigen, and its cyclic ℋ-module is a finite sum of eigen-modules λ_{π_i}.
  - E_F is ℋ-equivariant when its range is ℋ-stable.
  - Strong multiplicity one makes λ_π determine π. (Definition.)

## 1. Exact decompositions

### Epstein Q₁ = x² + 5y²
- **Setting.** K = ℚ(√−5), d = −20, h = 2, w = 2, Q₂ = 2x² + 2xy + 3y², and ψ is the genus character.
- **Decomposition.**
  - Z_{Q₁} = ζ_K + L_K(ψ) = ζ(s)L(s, χ₋₂₀) + L(s, χ₋₄)L(s, χ₅).
  - Z_{Q₂} = ζ(s)L(s, χ₋₂₀) − L(s, χ₋₄)L(s, χ₅).
  - The constants are exactly ±1: Z_{Q_i} = w·ζ_K(s, C_i) = (w/h)(ζ_K ± L_ψ). (PROVED, genus theory.)
- **Checks.**
  - At s = 3, lattice sum (R = 400) vs formula: 2.0783968688555 vs 2.0783968688659, and 0.4218385421188 vs 0.4218385421236. The gap is lattice-truncation level. (NUMERICAL)
  - r_{Q₁} = b₁ + b₂ and r_{Q₂} = b₁ − b₂, with b₁ = 1 ∗ χ₋₂₀ and b₂ = χ₋₄ ∗ χ₅: 0 failures for n ≤ 200. GE independently found 0 failures for n ≤ 300. (CHECKED)
- **No Euler product (elementary proof, from GE).** Set a = r_{Q₁}/2. Then a(2) = a(3) = 0, since 2 and 3 are not of the form x² + 5y², but a(6) = 2, from (±1, ±1). So a(6) ≠ a(2)a(3), a is not multiplicative, and Z_{Q₁} is not L(s, Π) for any Π. (PROVED; CHECKED in `epstein_non_euler`.)
- **Automorphic pieces.** Π₁ = 1 ⊞ χ₋₂₀ and Π₂ = χ₋₄ ⊞ χ₅, both weight-1 Eisenstein of level 20 with nebentypus χ₋₂₀. θ_{Q₁} lies in Π₁ ⊕ Π₂ and in neither summand alone. (PROVED, standard.)

### Davenport–Heilbronn
- **Decomposition.** DH = α L(s, χ) + β L(s, χ̄), where χ mod 5 has χ(2) = i, α, β = (1 ∓ iκ)/2, and κ = (√(10−2√5) − 2)/(√5 − 1) = 0.28408.
- **Root numbers.**
  - β/α = ε(χ) = 1/ε(χ̄) = 0.850651 + 0.525731i, with |β/α| = 1 (CHECKED).
  - The FE of L(χ) holds to 4e-21 (CHECKED).
  - The FE-symmetric combinations form the single complex ray β = α/ε(χ̄). (PROVED, two lines of algebra.)
- **No Euler product.** There are 906 non-prime-power atoms of −F′/F up to n ≤ 3000 (CHECKED).

## 2. W is not bilinear in F: what is actually true

**2a. Additivity over products (an identity, not content).**
- W_F is linear in −F′/F, log Γ_F and the pole divisor, and all three are additive under F ↦ FG. So W_{FG} = W_F + W_G and zeros(FG) = zeros(F) ∪ zeros(G). This is log(FG) = log F + log G. (PROVED, tautological.)
- The value 9.8e-15 for ζ·L(χ₅) is a **float check of the Dirichlet-convolution code**, not evidence.
- Likewise the Epstein defect max |c_EPS − c_{ζ_K} − c_{L_ψ}| = 315 just says log is not additive over sums.
- Consequence: the "cross term" of a sum F₁ + F₂ is **not** a bilinear B(F₁, F₂). (PROVED)
- **Compatibility with 14B (GE §5.1).** 14B's polarisation W_{F_θ}(g) = 2B_ζ(e^{−θu}g, e^{θu}g) is bilinearity in the *test vector*, not in F, so there is no conflict.
  - GE verified the archimedean constants term by term.
  - The θ- and τ-digamma identities ψ(¼ + θ/2) + ψ(¼ − θ/2) − 2ψ(¼) = −4∫w_ζ(cosh θt − 1), and its τ analogue, hold to 2e-15 at θ = 0.1, 0.2, 0.3 and to 2e-11 at τ = 1, 3 (`r14f_twist.json`). (CHECKED)
- **Atoms of −F′/F** (`r14f_logderiv.json`, N = 3000; CHECKED):

| model | non-prime-power atoms | max_p \|c(p)\|/log p |
|---|---|---|
| ζ, L(χ₅) | 0 | 1 |
| ζ_{ℚ(√5)} | 0 | 2 |
| ζ_K(−20), L_ψ | 0 | 2 |
| EPS | **203** (6, 14, 21, 36, …) | 2 |
| DH | **906** (6, 12, 14, 18, …) | 1 |
| F_{0.2} | 0 | **5.16** |

**2b. Off-line zeros of a sum come from interference.**
- **Identity.** For F = αF₁ + βF₂ with αβ ≠ 0: F(ρ) = 0 ⇔ F₁(ρ)/F₂(ρ) = −β/α. A common zero of F and F₁ is a common zero of F₁ and F₂. (PROVED)
- **Equal moduli are forced.** Here |β/α| = 1 (DH: ε(χ); EPS: α = β = 1), so |F₁(ρ)| = |F₂(ρ)| at *every* zero, on or off the line. That is not evidence. (PROVED, tautological.)
- **The content is the size of the pieces** (NUMERICAL, test data):
  - At DH's ρ = 0.808517182457 + 85.6993484854i: |L(χ)| = |L(χ̄)| = 0.70220.
  - At EPS's ρ = 0.932969697485 + 15.6682495313i: |ζ_K| = |L_ψ| = 0.93341, with ratio −1 to 3e-21.
  - So neither piece is near its own zeros there.
  - At ζ's first zero |L_ψ| = 9.41, so it is not a zero of Z_{Q₁}.
- **Collision homotopy Z_λ = ζ_K + λL_ψ** (`r14f_quotient.json`; NUMERICAL). Real λ keeps the FE, gamma factor, conductor and real coefficients.
  - **λ = 0.** On t ∈ [14.5, 16.2] the only on-line zero is L(χ₋₂₀)'s 15.4925. The other zero involved, L(χ₋₂₀)'s 14.338, enters the window as λ grows. (ζ's 14.1347 moves down and out of the window, to 13.40 at λ = 0.25.)
  - **The pair approaches.** 14.79 / 15.45 at λ = 0.06, 14.97 / 15.40 at λ = 0.10, 15.13 / 15.33 at λ = 0.13, 15.22 / 15.26 at λ = 0.140.
  - **No on-line zero in the window** at λ = 0.1425, 0.145, 0.15, 0.16 (step 0.005). So **λ* ∈ (0.1400, 0.1425)**, at height ≈ 15.24. This is consistent with GE's (0.140, 0.145).
  - **After the collision, the off-line track.** Re ρ = 0.566 / 0.595 / 0.633 / 0.718 / 0.890 / 0.933 at λ = 0.145 / 0.15 / 0.16 / 0.2 / 0.5 / 1, with Im ρ going from 15.244 to 15.668.
  - Losing on-line zeros in colliding pairs is classical for real FE-symmetric families. Here it shows the Epstein quartet is born from two zeros of one piece, L(χ₋₂₀), once the second eigen-component is mixed in at about 14%.

**2c. Restated item 5: the Tate/Hecke containment (replaces the vacuous Connes-quotient sentence).**
- **Retracted.** The earlier claim "ghosts are not spectrum of any ℋ-equivariant quotient" was true but vacuous.
  - Connes 1999 Thm 1 works in L²_δ, δ > 1, where the cokernel spectrum consists of **critical** zeros only. In the δ = 0 picture off-line zeros appear only as resonances.
  - So *no* off-line zero is ever spectrum of such a quotient, equivariant or not, and equivariance did no work. Also, Connes 1999 is GL₁ (Grössencharakters) only.
- **Correct statement (PROVED, standard Tate / local theory).**
  - For an ℋ-eigen (isotypic) test vector, the global zeta integral is L(s, π_i) times finitely many local factors, which are polynomials in p^{±s} or archimedean entire factors.
  - So the zeros of any equivariant construction lie in ∪_i zeros(L(π_i)) ∪ (local-factor zeros).
  - Consequently, any off-line zero of an equivariant zeta integral that is not a local-factor zero is an off-line zero of some L(π_i).
- **Application to the fakes (PROVED + NUMERICAL).**
  - The EPS and DH off-line zeros are not zeros of any piece (sizes 0.93 and 0.70, §2b).
  - They are not local-factor zeros either: a finite polynomial in p^{−s} vanishes only on vertical lines, and EPS's and DH's test vectors carry no such factor.
  - So they come from **non-equivariance**, i.e. interference.
- **The non-equivariance itself (CHECKED).**
  - DH evaluates at u = 1 in (ℤ/5)^×: the ratio c(2a)/c(a) = 0.284, −3.520, −3.520, 0.284 is not constant, against i for χ.
  - Epstein evaluates at the class C₁ in Cl_K ≅ ℤ/2. The GL₂ Hecke relation b(p)b(n) = b(pn) + χ₋₂₀(p)b(n/p) gives defect 0 for b₁ and b₂ and −4, −4, −4, −2, −4, −4, −4 for r_{Q₁}/2 at (p, n) = (3,3), (3,7), (7,3), (2,3), (7,7), (23,3), (3,23).
- **The same gcd principle for ζ (PROVED).** The 14A vector has zeta integral ζ(s)(1 − 2^{1−s})·Γ-factor, whose zeros on Re s = 1 are local-factor zeros. Only the gcd over the Hecke orbit is intrinsic.

## 3. The criterion and the seven test cases

**Criterion.** F = L(s, Π), with Π an isobaric sum of unitary cuspidal automorphic π_i on GL_{n_i}(𝔸_ℚ).
- As a *predicate* on the seven cases it is PROVED / CHECKED.
- As a *positivity criterion* it is CONJECTURE (GRH).

| model | Euler product | FE | automorphic (isobaric Π) | product vs linear combination | unitary constituents | tempered at p (Ramanujan) | criterion | ghosts |
|---|---|---|---|---|---|---|---|---|
| ζ | yes | yes | GL₁, trivial | single | yes | yes (1) | passes | none, L ≤ 8: **code check** (zeros in band known on the line) |
| ζ_{ℚ(√5)} = ζ·L(χ₅) | yes | yes | 1 ⊞ χ₅ | product (W additive: identity) | yes | yes (2 = degree) | passes | none (13j): **code check** |
| F_{iτ} = ζ(s−iτ)ζ(s+iτ), τ = 3 | yes (c = 2cos(τ log n)Λ) | yes | \|·\|^{iτ} ⊞ \|·\|^{−iτ} | product | **yes** | yes (\|2cos\| ≤ 2) | passes | **0 at L = 2 (M = 37) and L = 4 (M = 66)**, λ_min = +1.4e-10 / +6.1e-11: **code check** |
| F_θ = ζ(s−θ)ζ(s+θ), θ = 0.2 | yes, Λ ≥ 0 | yes | \|·\|^θ ⊞ \|·\|^{−θ} | product | **no** (non-unitary exponents) | no (5.16 > 2) | fails (unitarity) | 8 at L = 2 (λ_min = −0.091), 20 at L = 4 (λ_min = −0.755): NUMERICAL (built in) |
| Epstein x²+5y² | **no**: a(6) = 2 ≠ a(2)a(3) = 0 | yes (Poisson on ℤ²) | in Π₁ ⊕ Π₂, not eigen | **linear combination** | each piece yes | each piece yes | fails | quartet 0.933+15.67i: NUMERICAL (verified zero) |
| DH | **no** | yes | GL₁ vector mixed at 5 | **linear combination** | each piece yes | each piece yes | fails | 0.8085+85.70i: NUMERICAL (verified zero) |
| spliced Beurling | yes, Λ ≥ 0 | **no** | **no** (no ℚ^× ⊂ 𝔸 realisation) | – | – | – | fails (not automorphic) | **N.A. / defect**: no own explicit formula, so its borrowed-arch form is a defect, not a ghost count (14B erratum c) |

- **The twist pair** (`r14f_twist.json`; own code; pole term by 400-node quadrature, validated against ζ's closed form to 3.3e-14).
  - F_{iτ} and F_θ have the same product structure and differ only in unitarity. The ghost counts are 0 and 0 against 8 and 20.
  - GE reports 7 ghosts for F_θ at L = 2. The λ_min agrees (−0.091). The count differs by one at the 1e-8 relative threshold, which is a threshold artifact.
  - **Statement (PROVED, conditional).** A unitary twist multiplies the test vector by e^{iτu}, which is an isometry. So W_{ζ(·∓iτ)}(g) = W_ζ(e^{±iτu}g), and W_{F_{iτ}}(g) = W_ζ(e^{iτu}g) + W_ζ(e^{−iτu}g) for complex g.
  - It therefore **transports** positivity: W_{F_{iτ}} ≥ 0 if W_ζ ≥ 0 (i.e. under RH). It does not create positivity, and the "0 ghosts" is a code check.
  - A non-unitary twist e^{±θu} is not an isometry. W_{F_θ}(g) = 2B_ζ(e^{−θu}g, e^{θu}g) is a cross term between two different vectors, indefinite by construction (14B erratum b).
  - So unitarity is the axiom that keeps the positive cone invariant under twisting. (PROVED, elementary.)
- **Unitarity vs temperedness.**
  - For GL₁ they coincide (|χ(p)| = 1, or |p^{iτ}| = 1).
  - For GL_n, n ≥ 2, unitary cuspidal π are not known to be tempered (Ramanujan is open, e.g. for Maass forms). So "c(p)/log p > degree" detects **non-temperedness**, not non-unitarity in general.
  - For F_θ both fail, since the constituents |·|^{±θ} are non-unitary characters of GL₁. By Jacquet–Shalika uniqueness of the isobaric decomposition, F_θ has no unitary-constituent representation. (CHECKED, literature.)
- **The pre-audit falsifier is resolved.**
  - ζ_{ℚ(√5)} is two pieces and ghost-free (a product). F_θ is a product and ghosted (non-unitary).
  - "Irreducible ⇔ no ghosts" fails. The isobaric-unitary predicate separates all seven, and the F_{iτ} / F_θ pair isolates unitarity.
  - The agreement with ghost-freeness is conditional on RH/GRH for ζ, ζ_{ℚ(√5)} and F_{iτ} (code checks only). It is NUMERICAL for the ghosted fakes. Beurling is not a ghost case at all.
- **Strong multiplicity one** (Jacquet–Shalika 1981, isobaric classification) makes "F = L(Π)" a property of the coefficient sequence alone, through multiplicativity and the Hecke relations of §2c. Degree ≤ 1 Selberg class: Kaczorowski–Perelli 1999. (CHECKED, literature.)

**What the criterion gives, and what remains.**
- **Selection, not positivity (PROVED as a reduction).** Each axiom is an identity constraint, not an inequality: equivariance, multiplicativity, Hecke relations, isobaric uniqueness, unitarity.
- **The product reduction.** W_Π = Σ_i W_{π_i} (identity) and zeros(L(Π)) = ∪ zeros(L(π_i)). By Weil's criterion for each L(π_i), W_Π ≥ 0 ⇔ W_{π_i} ≥ 0 for all i. (PROVED, given Weil's criterion and its standard extensions to automorphic L-functions; CHECKED, literature.)
- **What §2c excludes.** It removes only *interference* ghosts. A ghost inside a single L(π) is not excluded. (PROVED; the remaining case is OPEN.)
- **The remaining inequality (GAP).** For each unitary cuspidal π (here ζ, χ₅, χ₋₂₀, χ₋₄, χ mod 5, and |·|^{±iτ}-twists of ζ), and all admissible g with h = g ∗ g̃:
  Σ_n Λ_π(n) n^{−1/2} [h(n) + h̄(1/n)] ≤ Pole_π(h) + ∫ |ĝ(½+it)|² w_π(t) dt,
  where w_π(t) = (1/2π) Re Σ_j ψ((½ + it + μ_j)/2) + (1/2π) log(N_π / π^n), the μ_j are the archimedean parameters and N_π is the conductor.
  - This is Weil's criterion, i.e. GRH for π: **CONJECTURE / OPEN**.
  - For GL₁ it is equivalent to the validity of Connes' global trace formula (Connes 1999). For n ≥ 2 the citation is Weil's criterion, not a Connes quotient.

## 4. Pair test (matched local data)
- **ζ_K(−20) vs Z_{Q₁}.** Same gamma factor, conductor and FE. The homotopy produces the ghost at λ* ∈ (0.1400, 0.1425), once about 14% of the second eigen-component is mixed in. (NUMERICAL)
- **L(χ mod 5) vs DH.** Same gamma factor and conductor. L(χ) is ℤ₅^×-eigen (ratio i) and DH is not, and DH's off-line zero sits where both pieces have size 0.70. The difference shows only in the eigen-property at 5. (CHECKED + NUMERICAL)
- **F_{iτ} vs F_θ.** Same product structure, differing only in unitarity: 0 vs 8 / 20 ghosts. (Code check vs NUMERICAL.)

## 5. Labels summary
- **§0:** "irreducible" for C_ℚ is vacuous: PROVED. The ℋ-module definition is a definition, not a claim.
- **§1:**
  - Epstein constants ±1: PROVED; coefficient identity (n ≤ 200): CHECKED.
  - a(6) ≠ a(2)a(3): PROVED / CHECKED.
  - DH root number and FE: PROVED / CHECKED.
- **§2a:**
  - Additivity over products: identity (PROVED, tautological); the 9.8e-15 is a float check.
  - Digamma twist identities: CHECKED.
- **§2b:**
  - Interference identity: PROVED.
  - Equal moduli: forced, so tautological.
  - Piece sizes 0.70 / 0.93: NUMERICAL.
  - λ* ∈ (0.1400, 0.1425): NUMERICAL.
- **§2c:**
  - Tate/Hecke containment: PROVED (standard).
  - Fakes' zeros come from non-equivariance: PROVED + NUMERICAL.
  - Hecke and translation defects: CHECKED.
  - The earlier Connes-quotient sentence: retracted as vacuous.
- **§3:**
  - Predicate on seven cases: PROVED / CHECKED.
  - ζ, ζ_{ℚ(√5)}, F_{iτ} ghost-free: code checks.
  - Fakes' ghosts: NUMERICAL.
  - Beurling: N.A. / defect.
  - Unitary twist transports positivity: PROVED (conditional on W_ζ ≥ 0).
  - Criterion ⇒ positivity: GAP. GRH for each unitary cuspidal π: CONJECTURE / OPEN.
- RH, GRH, ARP and H_ψ: OPEN. Nothing here closes any of them.

— Daniel Toupin / Golden Physics Project. This seat: Grok Bot (principal-series lane). Tasking: Heart of the Cosmos. Audit: Golden Eye (AUDIT14 §5) and Turbo.
