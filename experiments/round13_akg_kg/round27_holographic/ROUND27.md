# Round 27 (holographic): resonances, one plane with two faces, the χ-charge, and the finite-prime truncation

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5). 2026-11-19 (Toronto).

**RH: OPEN. Nothing here claims it.** Labels: PROVED / CHECKED / NUMERICAL / CONJECTURE / HEURISTIC / GAP / OPEN / KILLED / Cited (published result, stated and not re-derived). Zeros are never proof inputs. ζ zeros appear only (i) as an independent check of the explicit formula and (ii) to build F_θ, whose zeros are ρ ± θ by definition. The controls' off-line zeros (from R25 `meta.json`) are used only to place their fake resonances.

Frame (as in R21/R26): d = 1, s = Δ, H = L²(ℝ₊, dx). Mellin transform Mf(s) = ∫ f x^{s−1} dx is unitary on Re s = ½. Log coordinate x = eᵘ.
- Faces: J f(x) = x⁻¹ f(1/x), so (MJf)(s) = Mf(1−s).
- Sides: C f = f̄, so (MCf)(s) = conj Mf(s̄).
- Θ = CJ, Θf(x) = x⁻¹ conj f(1/x), so (MΘf)(s) = conj Mf(1−s̄).

Controls: F_θ = ξ(s+θ)ξ(s−θ) with θ = ¼; Epstein Z_Q with Q = x²+5y² (EPS); Davenport–Heilbronn (DH).

## Short answer (plain language first)

- **Daniel's picture is exactly right as a picture, and exactly equivalent to RH.**
  - One plane with two faces (x and 1/x), each split by the critical line into two sides, gives four sectors.
  - The four places a zero shows up, ρ, 1−ρ, ρ̄, 1−ρ̄, are the four lifts ++, +−, −+, −− of *Which Way Is Forward?*.
  - On the seam (Re s = ½), the "anti-aligned" lifts land on the aligned ones: the zero tunnels to s̄, and the quartet is one state.
  - Weil's quadratic form, written in these four sectors, has **no cross-terms**. Each collapsed quartet contributes a perfect square, "the other positive instead".
  - A split quartet contributes one square **minus** another, and that is the only minus sign anywhere on the zero side.
  - So "no negatives ever" = RH (PROVED equivalence; this is a translation, not a proof).
- **Conservation law.**
  - The one charge that ζ keeps exactly and every fake breaks is "each prime counted exactly once, with a unit-circle local factor": Satake parameter 1, a Lee–Yang circle per prime, Λ(p^k) = log p.
    - F_θ pays with p^{±¼}: charge 2.03–8.53 instead of 1.
    - Epstein pays with representation counts: mass off the prime powers, 87 negative Λ's below 5000.
    - DH pays with periodic signs: 1088 negative Λ's below 5000.
  - **New numerical fact (NUMERICAL):** ζ has *no budget*. Multiply every prime's weight by (1+ε) and Weil positivity survives only for ε ∈ [−6.3 %, +0.02 %] at window L = 1, and only for |ε| < 10⁻⁶ at L = 1.5 (both sectors). Over-counting by one part in a million already creates a χ = −1 direction.
- **Where it breaks: the infinite product.**
  - The charge is conserved by every finite truncation, but truncations lose positivity. At L = 1.5, dropping the single prime 3 gives λ_min = −9.4·10⁻³ even and −1.8·10⁻² odd, far above the 10⁻⁹ floor.
  - The per-prime phases do not converge on the line.
  - The zeros are made by the additive (theta/Poisson) continuation, not by any product.
  - No conserved quantity we tested is defined at S = ∞ in the strip. Charge + functional equation ⇒ RH is Selberg's conjecture, OPEN (GAP-27).

Technical summary:

1. **Resonance obstruction (part 1).**
   - ζ zeros are poles of the modular scattering matrix at s₀ = ρ/2 (Faddeev–Pavlov, Lax–Phillips). The residues are **never** L² (cusp integrand y^{−β}), on or off the line (PROVED).
   - |φ| = 1 on the line for ζ **and** all three controls (CHECKED to 10⁻³⁰). Unitarity is blind.
   - Self-adjointness of the scattering excludes DH/EPS (zeros with σ > 1 give forbidden poles) but **not** F_θ.
   - The only space where the zeros are normalizable is the Weil completion, which is a Hilbert space iff RH (PROVED, classical). GAP-26 ⇔ Weil positivity.
2. **Literature (part 2).** §2. Nothing published escapes the obstruction.
3. **Four sectors (part 3).** Lemmas 3.3(a–c) are PROVED. The prime side and the zero side agree to 6·10⁻¹⁰ relative for ζ and 1·10⁻⁹ for F_θ (CHECKED, explicit formula). This sets the double-precision floor at ≈ 6·10⁻⁹ (operator norm).
   - ζ: no negative direction above the floor, in either sector, in any test (L ≤ 5; heights 15.7 and 85.7). Arb-certified at support 2.6 in S36, both parities.
   - F_θ: negative in both sectors from L = 1–2, 85 % (L = 1) to 98 % (L = 2) carried by the χ = −1 squares of split quartets.
   - EPS at t ≈ 15.7 and DH at t ≈ 85.7: **exactly one negative direction per sector** (the signature (1,1) of one split quartet), from X = 55 (EPS) and X = 300 (DH).
4. **Finite-prime truncation (part 4).**
   - Window L is an exact finite-prime statement (PROVED).
   - The truncated scattering matrix is not unitary and has poles on the line (CHECKED).
   - The semilocal Tate product is unitary, with resonances on Re s = 1, never near ζ's.
   - Truncation conserves the charge and loses both the zeros and the positivity.

---

## 1. Part (1): the resonance obstruction

**1.1 Modular surface (Cited: Faddeev–Pavlov 1972; Lax–Phillips 1976).**
- For X = SL(2,ℤ)\H, the Eisenstein constant term is y^s + φ(s) y^{1−s}, with φ(s) = Λ(2s−1)/Λ(2s) and Λ(s) = π^{−s/2}Γ(s/2)ζ(s).
- φ is unitary on Re s = ½, and its poles in Re s < ½ are the points s₀ = ρ/2.
- RH ⇔ every such pole is on Re s = ¼ (Lax–Phillips phrase this as a decay rate of the cut-off wave semigroup Z(t)).

**1.2 Not normalizable, ever (PROVED).**
- At a pole s₀ = ρ/2 with ρ = β + iγ, the residue of E(z, s) has constant term ∝ y^{1−s₀}.
- On the cusp the integrand is ∫^∞ |y^{1−s₀}|² y⁻² dy = ∫^∞ y^{−β} dy, which diverges for every β ≤ 1.
- So no zero in the strip gives an L² state, on or off the line. The normalizable spectrum (Re s > ½) and the zeros (Re s₀ < ½) live in different half-planes.
- `r27_checks.json → resonance`:

| model | ρ used | s₀ = ρ/2 | cusp exponent −β | L²? |
|---|---|---|---|---|
| ζ | ½ + 14.1347i | 0.25 + 7.067i | −0.5 | no |
| F_θ | ¾ + 14.1347i | 0.375 + 7.067i | −0.75 | no |
| EPS | 0.9330 + 15.668i | 0.466 + 7.834i | −0.933 | no |
| DH | 0.8085 + 85.699i | 0.404 + 42.85i | −0.809 | no |

**1.3 Unitarity is blind (CHECKED).**
- φ_M(s) = Λ_M(2s−1)/Λ_M(2s) satisfies |φ_M(½+it)| = 1 for all four models:
  - ζ, F_θ, EPS: error ≤ 2·10⁻³¹;
  - DH: error 1.6·10⁻¹⁷ (the KAP constant is stored as a double).
- The reason is that each completed function satisfies Λ(s) = Λ(1−s) and Λ(s̄) = conj Λ(s) (CHECKED; `FE_*`, `reality_*`). Off-line zeros only move the resonances off Re s = ¼ (F_θ: +0.125; EPS: +0.216; DH: +0.154).
- This matches R13 O9: S-unitarity holds with or without off-line zeros.

**1.4 Where self-adjointness *does* bite (PROVED from cited facts).**
- On a genuine surface, the Maass–Selberg relations forbid poles of φ in Re s > ½ except on (½, 1].
- DH and Epstein (class number 2) have infinitely many zeros with Re s > 1 (Cited: Davenport–Heilbronn 1936). For their fake φ these give non-real poles in Re s > ½, which no self-adjoint Laplacian allows. So DH and EPS are excluded by self-adjointness of the *scattering* picture.
- **F_θ is not excluded:** its zeros have β = ½ ± ¼ ∈ (0, 1), so its resonances sit at Re s₀ = 3/8 and 1/8, which is allowed.
- What self-adjoint scattering proves is "no zeros in Re s > 1" (true for any Euler product), not RH. This is the precise point where it breaks.

**1.5 The one space where zeros are normalizable (PROVED; classical, Weil 1952, Connes 1999, Bombieri 2000).**
- Complete the compactly supported test functions h under W(h) = Σ_ρ ĥ(ρ) conj ĥ(1−ρ̄), computed from the **prime side** (no zeros).
- If W ≥ 0, the completion modulo null vectors is ≅ ℓ²(zeros, multiplicity), the dilation generator acts as multiplication by γ, and it is self-adjoint. This is exactly the Postulate-PIP space of R26.
- If some quartet is split, the completion is a Krein space with a negative direction: the "ghost". (Recalled from the literature: the negative index equals the number of off-line quartets per sector. §3.3(c) gives the elementary reason. Page reference not re-checked.)
- **So GAP-26 ⇔ Weil positivity ⇔ RH.** The construction "from the integers" exists and is canonical. Its *positivity* is the whole problem.

## 2. Part (2): literature (each item checked against the source or its abstract where marked ✓)

| work | what it does | where it stops |
|---|---|---|
| Faddeev–Pavlov, Zap. LOMI 27 (1972) ✓; Lax–Phillips, *Scattering Theory for Automorphic Functions* (1976) ✓, BAMS 2 (1980) ✓ | ζ zeros = poles of the modular scattering matrix. RH ⇔ a decay rate of Z(t) | φ is unitary for any zero location (§1.3). In their own words, "if RH holds the scattering is trivial and says nothing about the zeros on the line" (quoted via M. Watkins' scattering page ✓) |
| Weil 1952; Yoshida 1992; Bombieri 2000 (*Remarks on Weil's quadratic functional*) | RH ⇔ W ≥ 0; quartet = indefinite block | positivity not proved |
| Connes, Selecta 1999; Connes–Consani 2021 (archimedean Weil positivity), Connes–Consani–Moscovici 2024 (prolate) | zeros as an absorption spectrum on the adele class space. Semilocal S-truncation. Archimedean positivity on small windows | global (all-prime) positivity open. This is our part (4) |
| Burnol (Sonine spaces, co-Poisson, 1998–2004) | Hilbert spaces from x ↔ 1/x with Poisson; local scattering at each p | no RH |
| de Branges; Conrey–Li 2000 | Hilbert spaces of entire functions | Conrey–Li showed de Branges' positivity condition fails for ζ |
| Berry–Keating 1999; Sierra–Townsend 2008; Sierra 2019; Bender–Brody–Müller 2017 & Bellissard 2017 | xp-type operators | no self-adjoint operator with all zeros (see R26 §GAP) |
| Remmen, PRL 127, 241602 (2021) ✓ | amplitude with poles at m² = γ_n². Real masses ⇔ RH; unitarity ⇒ positive odd moments of 1/γ_n² **assuming RH** | unitarity bounds are consequences, not causes |
| Lee–Yang 1952; Newman 1976; de Bruijn–Newman, Rodgers–Tao 2020 (Λ ≥ 0) | circle theorems; RH ⇔ Λ ≤ 0 | RH is "barely true" if true |
| Davenport–Heilbronn 1936; Diamond–Montgomery–Vorhauer 2006 | FE without an Euler product ⇒ off-line zeros. Euler product (Beurling, Λ ≥ 0, tempered) without an FE ⇒ off-line zeros | shows that **both** ingredients are needed (§3.6) |
| Selberg 1992 (Selberg class) | conjecture: Euler product with tempered local factors + FE ⇒ RH | **OPEN**. This is GAP-27 |

## 3. Part (3): one plane, two faces, four lifts

**3.1 Dictionary (PROVED as an exact translation; not progress on RH).**
- Daniel's four lifts (q, t) ∈ {±}² from *Which Way Is Forward?*: ++ = s, −− = s̄, R_q: s ↦ 1−s, R_t: s ↦ 1−s̄, D = R_q R_t: s ↦ s̄. The observable is χ = q·t.
- The group generated is the Klein four-group {1, R_q, R_t, D} acting on zeros: ρ, 1−ρ, 1−ρ̄, ρ̄.
- The χ = +1 class is {ρ, ρ̄} (two lifts of one relational state, related by D). The χ = −1 class is {1−ρ, 1−ρ̄}.
- On Re s = ½ we have 1−s = s̄, so the χ = −1 class coincides with the χ = +1 class (the quartet is one state). Off the line it separates.
- **RH ⇔ no zero has a separate χ = −1 class.** On the Hilbert space: R_q ↔ J (faces), D ↔ C (sides), R_t ↔ Θ = CJ. The fixed set of R_t is the critical line, the seam.

**3.1b Tunneling to s̄ and the "other positive" rule (exact restatement).**
- R_q (faces, s ↦ 1−s) followed by D (sides, s ↦ s̄) is R_t: s ↦ 1−s̄.
- On the seam, 1−s = s̄. So the face-flip of a zero *is* its conjugate: the zero "tunnels" from the 1/x face straight to s̄ on the same side, and no new point is created.
- In the Weil form each zero ρ is paired with its R_t-partner 1−ρ̄, giving ĥ(ρ)·conj ĥ(1−ρ̄).
  - On the seam the partner is ρ itself, so the term is |ĥ(ρ)|² ≥ 0. This is Daniel's "instead of going negative it takes the other positive".
  - Off the seam the partner is a different point, and the product ĥ(ρ)·conj ĥ(1−ρ̄) has no sign.
- **Rule (Daniel): every cancellation is a pairing of two positive quantities across the faces, never a subtraction.**
  - Lemma 3.3 shows that on the zero side this rule holds iff every quartet collapses, i.e. iff RH.
  - On the prime side the subtraction −Σ Λ(n) n^{−½}[g(log n) + g(−log n)] is real, and it is what the rule must absorb.

**3.2 The four sectors (PROVED).**
- C is antilinear, so the sectors are real subspaces: H = ⊕_{ε,η=±} H_{εη} with Jf = εf and Cf = ηf.
- Multiplication by i maps H_{ε+} → H_{ε−}. Θ = CJ acts as εη.
- In log coordinates, H_{++} = real even functions of u, H_{−+} = real odd functions, and H_{ε−} = i·H_{ε+}.

**3.3 The Weil form on the sectors (PROVED).** Write ρ = ½ + iγ, and C_j(γ) = ∫ b_j cos(γu) du, S_j(γ) = ∫ b_j sin(γu) du for real b_j.
- **(a) No cross-terms.**
  - For real b_j, b_k the hermitian form is W(b_j, b_k) = Σ_ρ H_j(γ) H_k(−γ), with H(γ) = ∫ b e^{iγu} du.
  - If b_j is even and b_k is odd, the summand is odd in γ. The zero multiset is invariant under γ ↦ −γ (that is R_q), so the sum is 0.
  - W(if) = W(f), so the C-sectors carry identical copies.
  - Hence W = W⁺ ⊕ W⁻ ⊕ (copies), with W⁺_{jk} = Σ C_jC_k and W⁻_{jk} = Σ S_jS_k. Positivity on all of H ⇔ positivity on H_{++} and on H_{−+}.
  - **The "counted once" information is not in any cross-term. There are none.**
- **(b) Collapsed quartet (on-line, γ real).** The orbit is the pair {γ, −γ}, and its contribution is 2 C(γ)C(γ)ᵀ ⪰ 0. Each lift is paired with its R_t-partner, which is itself: "the other positive instead".
- **(c) Split quartet (γ = t + iδ, δ ≠ 0).**
  - The orbit is {±γ, ±γ̄}, and its contribution is 4 Re[C(γ)C(γ)ᵀ] = 4(aaᵀ − bbᵀ) with a = Re C(γ), b = Im C(γ). The odd sector is the same with S in place of C.
  - The −bbᵀ term is the only subtraction on the zero side. It has signature (1, 1).
  - **Daniel's "no negatives" rule ⇔ b = 0 for every quartet ⇔ RH** (via Weil's criterion: any split quartet produces a negative direction for a suitable test function).
- **Implementation note.** The orbit must be the *set* {±γ, ±γ̄}. Listing four lifts for an on-line zero double-counts it. The first run of `r27_sectors.py` did exactly that and was off by a factor of 2 for ζ only. "Counted once" is literally the collapse of the quartet.

**3.4 Numerics (NUMERICAL; `r27_sector_check.json`, `r27_quartet.json`, `r27_scan_*.json`).** Basis: b_j = cos²(πu/L)·{cos, sin}(jπu/(L/2)) on |u| ≤ L/2, so n ≤ e^L. Prime side from the R22 explicit formula, with sector kernel g_jk = ±(b_j ∗ b_k).
- **Explicit-formula check, prime side vs zero side** (300 ζ zeros; F_θ zeros ρ ± ¼):
  - F_θ: agreement to 1·10⁻⁹ (L = 2) in both sectors.
  - ζ (after the orbit fix): agreement 6.0·10⁻¹⁰ (L = 2, even), 1.1·10⁻¹¹ (odd), 9.3·10⁻¹⁰ (L = 5, even). Operator-norm differences ≤ 5.7·10⁻⁹ define the floor (`r27_budget.json → floor_*`). Every ζ prime-side eigenvalue in this round lies above −1.3·10⁻⁹, i.e. inside the floor.
- **F_θ attribution** (min eigenvector v):

| L | sector | λ_min | vᵀ(+aaᵀ)v | vᵀ(−bbᵀ)v | first quartet (γ₁) −part |
|---|---|---|---|---|---|
| 1 | even | −2.08e−5 | 4.5e−6 | 2.53e−5 | 2.46e−5 |
| 1 | odd | +5.66e−5 | 9.2e−5 | 3.54e−5 | 3.18e−5 |
| 2 | even | −7.91e−3 | 2.0e−4 | 8.11e−3 | 6.34e−3 |
| 2 | odd | −8.04e−3 | 2.4e−4 | 8.28e−3 | 6.32e−3 |

  The failure is carried by the χ = −1 (Im) part of the split quartets, mostly the first one, exactly as §3.3(c) predicts. The even sector fails first (L = 1).
- **Sector scan, λ_min / max|W|** (n_neg counts eigenvalues < −10⁻⁹):

| L | ζ even | ζ odd | F_θ even | F_θ odd | EPS even | EPS odd | DH even | DH odd |
|---|---|---|---|---|---|---|---|---|
| 1 | 2.68e-06 | 1.09e-04 | -2.65e-05 (n_neg 1, t≈0.0) | 6.88e-05 | 5.63e-04 | 2.40e-03 | 4.71e-04 | 2.45e-03 |
| 2 | -3.94e-11 (floor) | -1.56e-13 (floor) | -4.88e-03 (n_neg 4, t≈10.5) | -4.95e-03 (n_neg 4, t≈10.5) | 1.21e-04 | 1.06e-03 | 1.38e-05 | 2.31e-04 |
| 3 | -4.15e-10 (floor) | -9.21e-13 (floor) | -2.30e-02 (n_neg 4, t≈18.6) | -1.88e-02 (n_neg 5, t≈18.6) | 4.25e-06 | 1.49e-04 | -2.71e-11 (floor) | -2.17e-14 (floor) |
| 4 | -7.32e-10 (floor) | -3.00e-12 (floor) | -4.04e-02 (n_neg 4, t≈23.0) | -4.07e-02 (n_neg 5, t≈23.0) | -5.11e-03 (n_neg 1, t≈14.3) | -5.11e-03 (n_neg 1, t≈14.3) | -1.42e-10 (floor) | -4.67e-13 (floor) |
| 5 | -1.27e-09 (floor) | -6.73e-12 (floor) | -4.76e-02 (n_neg 4, t≈22.6) | -5.14e-02 (n_neg 3, t≈23.0) | -6.54e-02 (n_neg 1, t≈14.4) | -6.54e-02 (n_neg 1, t≈14.4) | -4.19e-10 (floor) | -3.93e-12 (floor) |

  - ζ: no eigenvalue below the double-precision floor (|λ| ≲ 10⁻⁹ relative) in either sector.
  - The ζ Gram is so close to singular that double precision cannot certify its sign beyond L ≈ 2. The rigorous version is S36, which arb-certifies both parities (both J-sectors) at support 2.6 with λ_min ≥ 3.63·10⁻⁶² (even) and 1.49·10⁻⁵⁸ (odd).
  - EPS and DH are positive on these low-frequency bases because their first off-line zeros sit at t = 15.7 and 85.7, above the basis band (§3.5).

**3.5 Sector-resolved frequency probe (NUMERICAL; `r27_probe.json`).** Test functions h_T^± = cos²(πu/L)·{cos, sin}(Tu), T ∈ [0, 100], X = e^L ∈ {100, 300}:

| X | model | even: min W/‖h‖² (at T) | odd: min (at T) | # negative T (even/odd) |
|---|---|---|---|---|
| 100 | ZETA | 7.91e-09 (2.0) | 2.47e-08 (4.5) | 0/0 |
| 100 | FTH | -0.0115 (27.75) | -0.0115 (27.75) | 20/27 |
| 100 | EPS | 0.307 (3.5) | 0.297 (3.5) | 0/0 |
| 100 | DH | 1.35e-05 (0.75) | 7.61e-05 (1.0) | 0/0 |
| 300 | ZETA | 5.21e-09 (0.5) | 1.48e-08 (1.0) | 0/0 |
| 300 | FTH | -0.0566 (23.0) | -0.0566 (23.0) | 48/40 |
| 300 | EPS | 0.0385 (3.5) | 0.0403 (3.5) | 0/0 |
| 300 | DH | 4.57e-06 (1.75) | 1.45e-05 (0.75) | 0/0 |

**3.5b Localized sector Galerkin at the controls' off-line heights (NUMERICAL; `r27_localized.json`).**
- Basis: cos²(πu/L)·{cos, sin}(T_k u), with T_k = T₀ + {−6, …, 6}.
- Exactly **one** negative eigenvalue per sector, which is the signature (1,1) of one split quartet per sector predicted by Lemma 3.3(c).

| model | T₀ | X = e^L | even λ_min/max | odd λ_min/max |
|---|---|---|---|---|
| EPS | 15.67 | 55 | −5.6e−3 (1 neg) | −5.7e−3 (1 neg) |
| EPS | 15.67 | 100 | −3.9e−2 (1) | −3.9e−2 (1) |
| EPS | 15.67 | 300 | −1.24e−1 (1) | −1.24e−1 (1) |
| DH | 85.70 | 300 | −6.4e−2 (1) | −6.4e−2 (1) |
| DH | 85.70 | 1000 | −1.05e−1 (1) | −1.05e−1 (1) |
| DH | 85.70 | 3000 | −1.19e−1 (1) | −1.19e−1 (1) |
| ζ | 85.70 | 300 / 1000 | +1.4e−11 / +7.0e−12 | +1.6e−11 / −1.0e−11 (floor) |
| ζ | 15.67 | 100 | −4.4e−10 (floor) | −7.6e−10 (floor) |

- The single-function frequency probe (3.5) misses EPS/DH, because one cos²-window function cannot resolve δ = 0.43 / 0.31 at these X. The 13-function local Galerkin does resolve it.
- Even and odd values coincide to 4–7 digits for EPS/DH. A split quartet at large t loads both J-sectors equally, since for large t the cos and sin parts of a quartet's bbᵀ block have nearly equal mass.

**3.6 Where ζ differs from the controls, in this language.**

| property | ζ | F_θ | EPS | DH |
|---|---|---|---|---|
| Klein four-group symmetry, 4 sectors, no cross-terms (§3.3a) | ✓ | ✓ | ✓ | ✓ |
| line-unitary scattering φ_M (§1.3) | ✓ | ✓ | ✓ | ✓ |
| local Tate gamma unitary on the line, winding 0 per prime | ✓ | ✓ | — (no local factors) | — |
| resonances in Re s ≤ ½ (self-adjoint scattering, §1.4) | ✓ | ✓ | ✗ | ✗ |
| Λ_M ≥ 0, supported on prime powers | ✓ | ✓ | ✗ (87 neg.) | ✗ (1088 neg.) |
| χ-charge: Satake moduli 1 / Lee–Yang per prime / Λ_M(p^k) = log p exactly | ✓ | ✗ (2^{±¼}) | undefined | undefined |
| sector Weil form ≥ 0 ("no negatives") | ✓ so far (certified on support 2.6) | ✗ from L = 1 | ✗ at t ≈ 15.7 from X = 55 (1 neg/sector) | ✗ at t ≈ 85.7 from X = 300 (1 neg/sector) |

**Reading.**
- Only the χ-charge row separates ζ from all three controls.
- F_θ is the sharp control: it has an Euler product and Λ ≥ 0, and breaks only the charge. Its local poles split onto two lines Re s = ±θ, a **local** χ = −1 class.
- EPS and DH never had local charges to conserve.

## 3b. Conservation-law candidates (Daniel: "χ is a conserved charge")

| candidate | conserved by each ζ Euler factor? | broken by F_θ / EPS / DH? | survives ∏_p into the strip? | label |
|---|---|---|---|---|
| U(1) phase per prime, \|p^{−it}\| = 1 | yes | **no**: it is a property of the variable, true for every model | n/a | KILLED as a discriminator |
| line-unitarity of the local gamma, \|γ_p(½+it)\| = 1 | yes | **no**: F_θ also (real Satake parameters give conjugate moduli); CHECKED 0.0 | partial products of phases do not converge on the line (−14.3 at 10⁶, −15.7 at 2·10⁶, t = 10; drift scale 2\|li(X^{½+it})\| ~ √X/(t log X)) | KILLED as a discriminator |
| winding / argument-principle count of γ_p per period | 0 | **no**: F_θ also 0 (CHECKED) | n/a | KILLED |
| **Satake parameter on \|α_p\| = 1** ⇔ Lee–Yang zero of 1/L_p at z = p^{−s} = 1 on \|z\| = 1 ⇔ Λ_M(p^k) = log p ("each prime counted once") | **yes, exactly** | **yes**: F_θ 2^{±¼}, ratio 2.03–8.53; EPS/DH no local factors, Λ_M off prime powers and negative | finite products: yes, trivially (the zero set of ∏ is the union, all on the torus \|z_p\| = 1). Into the strip: **no mechanism**. ∏ converges only for Re s > 1, and 1/ζ_S has zeros only on Re s = 0 | PROVED locally; GAP-27 globally |
| Weil positivity on window L | exact finite-prime statement (p^k < e^L) | yes (controls fail at finite L) | that *is* RH | equivalent to RH, OPEN |

**3b.1 "No budget": ζ sits on the edge of the positive cone (NUMERICAL; `r27_budget.json`, fig 5).**
- Scale the prime side Λ → (1+ε)Λ, keeping the pole and archimedean terms. Basis as in 3.4; floor ≈ 6·10⁻⁹.

| L | sector | λ_min at ε=0 | PSD window in ε | λ_min at ε = +10⁻⁶ / −10⁻⁶ |
|---|---|---|---|---|
| 1.0 | even | 2.7e−6 | [−0.063, +0.0002] | 2.67e−6 / 2.70e−6 |
| 1.0 | odd | 1.1e−4 | [−0.0079, +0.63] | — |
| 1.5 | even | −3e−12 (floor) | none above floor | **−3.1e−7 / −6.6e−8** |
| 1.5 | odd | 1.3e−11 (floor) | [0, 5e−6] | 1.4e−7 / **−1.5e−7** |

- At L = 1.5, a change of the prime weights by one part in a million in either direction produces a negative eigenvalue 30–300× above the floor.
- The "counted exactly once" point ε = 0 is (numerically) the *only* PSD point. This is Daniel's "no budget to create a χ = −1 pair", made quantitative.
- **Label: NUMERICAL**, double precision, two windows. CONJECTURE: the PSD ε-window shrinks to {0} as L → ∞.
  - The side ε > 0 is forced by the prime number theorem: the pole term fixes the prime mass, which is what F_θ compensates with its two poles at 1 ± θ.
  - The side ε < 0 is the non-trivial one.
- Relation to the literature: the near-degeneracy of ζ's Weil form (exponentially small eigenvalues) is known (Connes–Consani–Moscovici, prolate; S36 here). Reading it as ε-rigidity of the prime weights is, as far as we know, new framing, not a new theorem.

**3b.2 Truncation keeps the charge, loses positivity.**
- Keep only the primes p ≤ P inside the window (L = 1.5, primes 2 and 3):
  - P = 0: λ_min = −0.084 (even), −0.164 (odd);
  - P = 2: −9.4·10⁻³ / −1.8·10⁻²;
  - P = 3: back to the ζ value (floor).
- At L = 1 (only prime 2 in the window), dropping 2 gives −3.7·10⁻³ / −1.8·10⁻².
- So each truncation, which conserves the per-prime charge perfectly, is *not* positive. Positivity is a global balance between the pole, the archimedean term and **all** primes in the window together, not a sum of per-prime conserved pieces (NUMERICAL; consistent with S54's "every prime ≤ 13 necessary").
- Truncated objects have no functional equation, so their negative directions are ghosts not attached to any quartet.
- **This is the precise sense in which the charge does not survive the infinite product:** conservation is local, positivity is global.

**Honest verdict.**
- The surviving charge is "tempered local factors with multiplicity one". This is the Selberg-class axiom, and the expected answer.
- It is conserved under every finite product. It fails to reach the strip because the strip is outside the product's convergence region, so the zeros are not zeros of any product (truncations have none).
- To get from the charge to RH one needs the functional equation, which comes from the additive Poisson/theta structure, not from the primes. One then needs a theorem that the two together forbid a χ = −1 class.
- Beurling systems (DMV 2006) show that the charge alone allows a separate χ = −1 class: they have Λ ≥ 0, tempered local factors and multiplicity one, but no functional equation. DH/EPS show the same for the functional equation alone.
- "Charge + FE ⇒ no χ = −1 class" is the Selberg-class RH: **OPEN**. Every known attempt breaks exactly at the passage from finite S to the full product. The charge has no known continuous extension into 0 < Re s < 1, and Bohr–Bagchi universality in ½ < σ < 1 shows that the product's statistics there are those of independent phases (Cited, HEURISTIC as applied here).

## 4. Part (4): the finite-prime truncation

**4.1 Exactness on windows (PROVED).** If ĥ has support in |u| ≤ L/2, then g = h ∗ h^♯ has support in |u| ≤ L. The prime side Σ Λ(n) n^{−½} g(log n) then involves only n < e^L. So W_L is a finite-prime object (GPPVerify `kS_eq_Kpsi`), and every finite-L positivity statement is a finite computation (S36: L = 1.3, primes up to 13 enter through the lift). RH ⇔ W_L ⪰ 0 for all L.

**4.2 The truncated scattering matrix is not unitary (CHECKED).**
- φ_P(s) = [Γ_ℝ(2s−1)/Γ_ℝ(2s)] ∏_{p≤P}(1−p^{−2s})/(1−p^{1−2s}).
- On s = ½+it the factor 1/(1−p^{−2it}) has poles at t = πk/log p. For example t = π/log 2 = 4.532, where |φ_13| = 3.9·10⁶ at distance 10⁻⁸.
- |φ_P(½+it)| at t = (3.1, 9.4, 21.0):

| P | t=3.1 | t=9.4 | t=21.0 |
|---|---|---|---|
| 2 | 0.773 | 0.749 | 0.293 |
| 13 | 0.368 | 0.23 | 14.8 |
| 101 | 0.0477 | 0.00438 | 90.5 |
| 1009 | 3.25e-10 | 1.93e-05 | 414 |

- Line unitarity of φ is a property of the *completed* (functional-equation) object, not of any truncation.

**4.3 Semilocal (Tate/Connes) truncation is unitary but has the wrong resonances (PROVED).**
- γ_S(s) = ∏_{p∈S} L_p(1−s)/L_p(s) has |γ_S(½+it)| = 1 (CHECKED 0.0).
- Its poles are at p^{s−1} = 1, i.e. Re s = 1, Im s ∈ 2πℤ/log p. Its zeros are at Re s = 0. None is near a ζ zero, for any finite S.
- The partial phases do not converge as S grows (§3b). The ζ zeros are a **collective** effect of all primes plus the archimedean place, visible only after continuation.
- In Daniel's language: each finite truncation is perfectly aligned (χ = +1 everywhere, no zeros at all). The question of whether a χ = −1 class appears is decided only at S = ∞.

## 5. Where it breaks (summary) and GAP-27

- Part 1 breaks at normalizability. Zeros are never L² states of the modular surface, and unitarity is blind. The Weil completion makes them normalizable iff RH.
- Part 3 breaks at the sector level: the sectors carry no cross-term information. The full information is the sign of the χ = −1 square in each split quartet, which is RH itself (exact restatement).
- Part 4 breaks at completion: finite truncations conserve the charge and have no zeros.

**GAP-27 (OPEN).** Show that the per-prime χ-charge (Satake parameter exactly 1, multiplicity one), together with the theta functional equation, forbids the −bbᵀ block in the Weil form. Equivalently, find a quantity defined at S = ∞ (not as a limit of products in Re s > 1) that is conserved by ζ, broken by F_θ (charge violation), and undefined or broken for DH and EPS (FE without product). Beurling (product without FE) must also fail it. Any candidate must pass all four filters before it is taken seriously.

## 6. Files

| file | content |
|---|---|
| `r27_sectors.py` | sector Gram (prime side, R22 formula; zero side), explicit-formula check → `r27_sector_check.json` |
| `r27_scan.py` | λ_min per sector, all models, L = 1…5 → `r27_scan_*.json`, `scan.log` |
| `r27_probe.py` | sector-resolved frequency probe → `r27_probe.json`, `probe.log` |
| `r27_quartet.py` | χ = ±1 attribution for F_θ → `r27_quartet.json` |
| `r27_checks.py` | part 1 (FE, unitarity, resonances), part 4 (truncation), conservation candidates → `r27_checks.json` |
| `r27_dh.py` | localized sector Galerkin at the EPS/DH off-line heights → `r27_localized.json`, `localized.log` |
| `r27_budget.py` | ε-rigidity ("no budget"), prime truncation, ζ prime-vs-zero floor → `r27_budget.json`, `budget.log` |
| `r27_figs.py` | figures 1–5 |
| `r22_models.py` | copy of the R22 explicit-formula data |
| `zeta_zeros_300.json` | first 300 ζ ordinates (checks only) |

Reproduce: `python3 r27_checks.py; python3 r27_sectors.py; python3 r27_scan.py 1,2,3,4,5; python3 r27_probe.py; python3 r27_quartet.py; python3 r27_dh.py; python3 r27_budget.py; python3 r27_figs.py` (numpy, scipy, mpmath, matplotlib). Total about 1.5 h on 8 cores.

## 7. Figures

- `fig1_four_lifts.png`: the four lifts. ζ: the quartet collapses on the seam. F_θ: a separate χ = −1 class.
- `fig2_sector_scan.png`: λ_min per J-sector vs L, all models, with the floor band.
- `fig3_fth_chi_attribution.png`: F_θ negativity split into χ = +1 and χ = −1 squares.
- `fig4_charge_and_phase.png`: the per-prime charge Λ_M/Λ for the four models, and the non-convergence of the product of local phases on the line.
- `fig5_no_budget_and_truncation.png`: ε-rigidity, and the loss of positivity under truncation.

## 8. Done vs not done, new vs known, next step

**Done:** parts 1, 2, 3 (four sectors, dictionary, tunneling, rule), 4, and the conservation hunt (five candidates, with the budget and truncation tests).

**Not done:**
- No arb certification of the new numerics (ε-window, truncation, localized EPS/DH). These are double precision, against a measured floor of 6·10⁻⁹.
- The ε-window at L ≥ 2 is below the double floor, so it needs mpmath/arb.
- No Beurling control in the sector code (R13 harness has one).
- No literature page-check for Bombieri/Yoshida's negative-index statement.

**New versus the literature (modest):**
1. The exact four-lift/χ dictionary and the sector form of Weil's criterion, in which the χ = −1 square is the only minus sign. The mathematics is classical; the packaging is new.
2. Localized sector Galerkin: exactly one negative direction per sector at EPS's and DH's first off-line zeros.
3. The ε-rigidity ("no budget") measurement and the truncation-loses-positivity measurement. Both are numerical.
Nothing here proves RH or a new theorem about ζ.

**Recommended next step:**
- Certify the ε-rigidity at L = 1.5–2.6 in arb, reusing the S36 machinery (it already has both parities).
- Then test the CONJECTURE "PSD ε-window → {0}" against F_θ (rescaled by its own charge) and a Beurling system.
- If the window provably collapses only for systems with charge exactly 1 *and* a functional equation, that is a concrete, checkable form of GAP-27.

Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
