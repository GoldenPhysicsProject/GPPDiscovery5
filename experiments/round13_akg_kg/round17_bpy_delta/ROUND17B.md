# Round 17B: "s IS Δ" for the massless spin-2 (graviton) conformal primary

Heart of the Cosmos (Grok / GPPDiscovery5), 2026-11-18 (Toronto). Labels: PROVED / CHECKED / NUMERICAL / CONJECTURE / HEURISTIC / GAP / OPEN.
**RH: OPEN.** Zeros are never proof inputs. Every profile below is built from theta series or Dirichlet coefficients. The one zero in §3 is located by Newton from a Round 17 label, and only to show the failure mode.

Daniel's correction: "s IS Δ" means the **massless (graviton) conformal primary Δ**, not a generic weight-Δ primary x^{−Δ} as in Round 17. This round redoes the identity with the PS graviton wavefunction and tests whether the massless/spin-2 data add a constraint.

## Summary

1. **Graviton facts (PROVED in the literature; exact citations in §1).** In 4D, massless spin-2 conformal primaries are Mellin transforms in the **energy ω** of plane waves (Pasterski–Shao 2017, eq. 6.20). In the inner product (6.16) they are δ-normalizable on Δ ∈ 1 + iℝ and span the nonzero-energy plane waves (PS §6.3). That conclusion is reached "following the same argument used in Section 4.2", i.e. by reduction to the explicit massless scalar computation (4.6)–(4.9). The selection rule is the convergence of ∫dω ω^{2c−d−1+iν} (eq. 4.7), which holds only at c = d/2 = 1.
2. **The ζ map (PROVED, CHECKED to 1e−22).** Let a_ζ(ω) := 2 f_ζ(ω²), where f_ζ is Round 17's BPY state. Take the graviton wavepacket Φ_ζ with momentum-space amplitude a_ζ(ω), any angular profile g(ŵ), and either helicity. Then
   **A_ζ(Δ) := ∫₀^∞ ω^{Δ−1} a_ζ(ω) dω = ξ(Δ/2) for every Δ ∈ ℂ**,
   and its overlap with the graviton primary h^{Δ,±}(·; w) is a nonzero constant × g(w) × ξ(Δ̄/2). So ζ's s is exactly the graviton's Δ/2, and the reason for the 2 is now explicit. The Lorentz-invariant measure d³p/2ω ∝ ω dω d²w is flat Lebesgue measure in y = ω², so PS's line Re Δ = 1 *is* the Mellin–Plancherel line Re s = ½ of L²(ℝ₊, dy). Among the power maps y = ω^k, only k = d = 2 does this (PROVED).
3. **What normalizability gives (PROVED) and what it doesn't (GAP).** (i) The PS theorem says which Δ label normalizable basis elements: Re Δ = 1 ⇔ Re s = ½. (ii) RH needs ξ's zeros to be such Δ. Φ_ζ is normalizable (KG norm = Plancherel norm on Re Δ = 1, CHECKED to 1e−26). That gives only A_ζ ∈ L²(1 + iℝ). It says nothing about where the entire function A_ζ vanishes. The missing link is exactly GPPVerify's `GravitonZeroIdentification` (= `AllXiZerosPrincipalOne`, OPEN, RiemannHypothesis/GravitonDeltaIdentification.lean:154).
4. **The extra graviton structure, tested item by item (§2.5).** Helicity ±2, the factor 2, conformally soft values Δ = 2, 1, 0, −1, BMS Ward identities, memory and shadow: none constrains the zeros. One finding is new and exact. The raw lattice graviton Σ_n e^{−π n² ω⁴} has a soft expansion that **Poisson summation truncates after two terms**, 1/(2ω²) − ½ + O(e^{−π/ω⁴}/ω²). So its only soft poles are at **Δ = 2 (stress tensor) and Δ = 0 (superrotation)**, a shadow pair. It has no Weinberg pole at Δ = 1 (CHECKED).
5. **Epstein verdict.** Every graviton step goes through for x² + 5y²: the identity, normalizability, the shadow, and the same two soft poles. The one exception is positivity of the profile, which is not a physical requirement, and F_θ satisfies it anyway. EPS has a normalizable graviton state that is blind to the **non-normalizable** primary Δ = 1.8659 + 31.3365i (overlap 4e−32 against a scale of 2e−8). The only step that must use the Euler product is the zero-to-normalizable step.
6. **Does "graviton" add anything over generic Δ? Kinematically yes, on zeros no.** It fixes the variable (y = ω²), the coefficient (Δ = 2s), and which ζ-object is the primary's amplitude (ξ(Δ/2), with Γ_R(Δ/2)). It also rules out the naive "Planck gas" reading s = Δ, which would put ζ's critical line at Re Δ = ½, off the normalizable line (CHECKED). But it gives no condition on where the zeros sit beyond what the generic d = 1 statement of Round 17 gave. Under ω ↦ ω² the two are unitarily equivalent (PROVED).

---

## 0. Sources read

- **Pasterski–Shao, "Conformal Primary Basis for Flat Space Amplitudes", arXiv:1705.01027 (PRD 96, 065022)**, PDF supplied by Daniel. Primary reference for everything massless and spinning.
- **Pasterski–Shao–Strominger, "Flat Space Amplitudes and Conformal Symmetry of the Celestial Sphere", arXiv:1701.00049 (PRD 95, 065026)**, PDF supplied by Daniel. §2.2 (massless wavefunction) and App. A (KG inner product).
- Donnay–Puhm–Strominger, arXiv:1810.05219; Kapec–Mitra–Raclariu–Strominger, arXiv:1609.00282; Guevara–Himwich–Pate–Strominger, arXiv:2103.03961 (soft/special values).
- GPPVerify main `7ebc350` (read-only): `RiemannHypothesis/GravitonDeltaIdentification.lean`, `RiemannHypothesis/CelestialRiemannAffineBridge.lean` (`celestialDelta s = 2s`), `QuantumGravity/DiagonalConformalLift.lean`, `NumberTheory/WeightKShadow.lean`, `CelestialHolography/ConformalPrimaryWavefunction.lean`, `PSSNormalizability.lean` (inventory in ROUND17.md §0).
- Daniel's draft: ch. 7 "canonical dictionary" (four proofs of Δ = 2s, lines ≈15380–15530; Ghost Theorem ≈15500) and ch. 5 graviton sections (soft theorem ≈5500, ≈9313; graviton operators ≈7493; stress tensor ≈8762).

## 1. The graviton facts, pinned to the papers

**1.1 What PS 2017 prove, and in which inner product.**

| Case | Wavefunction | Inner product | What is shown | Where |
|---|---|---|---|---|
| Massive scalar | φ^±_{Δ,m} (2.10), bulk-to-boundary propagator on H^{d+1} | KG (3.6); plane waves (3.7) | Inverse transform exists on Δ ∈ d/2 + iℝ (3.1), using the H^{d+1} orthonormality relations (3.2)–(3.4) from their ref. [35]. Explicit δ-normalization (3.8). Basis statement in §3.3. Discrete series not normalizable (3.12). | §3.1–3.3 |
| Massless scalar | φ^±_Δ = ∫₀^∞dω ω^{Δ−1}e^{±iωq·X−εω} = (∓i)^ΔΓ(Δ)/(−q·X∓iε)^Δ | KG (3.6) | Inverse Mellin (4.6) expands every nonzero-energy plane wave on any line c + iℝ, c > 0. The KG product (4.7) ∝ δ^{(d)}(w₁−w₂)∫dω ω^{2c−d+iν₁−iν₂−1} diverges unless **c = d/2**. At c = d/2 it gives (4.8)–(4.9) (± 8π^{d+2} δ(ν₁−ν₂)δ^{(d)}(w₁−w₂)). Basis B^±_{m=0} (4.12). The constant (zero-energy) solution is excluded and not spanned. | §4.1–4.2 |
| Photon | Mellin gauge representative ϕ^{Δ,±}_{μa} (5.21) = Mellin (5.23) | Maxwell symplectic product (5.19); plane waves (5.20) | "the same argument in Section 4.2 can be directly borrowed": δ-normalizable on d/2 + iℝ, spans (5.18). Pure gauge only at Δ = 1 (5.14). | §5.2–5.3 |
| **Graviton** | h^{Δ,±}_{μ₁μ₂;a₁a₂} (6.4), traceless, Lorenz and radial gauge (6.7), □h = 0 (6.8). Mellin representative ϕ^{Δ,±}_{μ₁μ₂;a₁a₂} = (∓i)^ΔΓ(Δ)P ∂q ∂q/(−q·X∓iε)^Δ (6.18), equal to (6.4) up to residual diffeomorphisms (6.19), (6.14). | Linearized-Einstein symplectic product (6.16); plane waves (6.17) ∝ ωq⁰δ^{(d+1)}(ωq^i − ω′q′^i) | (6.20) ϕ^{Δ,±} = ∫₀^∞dω ω^{Δ−1} g^{±}(X; ω, q). "Following the same argument used in Section 4.2", δ-normalizable (w.r.t. (6.16)) on **Δ ∈ d/2 + iℝ** and spanning the plane waves (6.15). Large diffeomorphisms and zero-energy modes are set aside (footnote 21). | §6.1–6.3 |

Precise reading. Only the **scalar** cases are computed in full (massive: §3; massless: (4.7)–(4.9)). The spin-1 and spin-2 statements are reductions. One picks a gauge/diffeomorphism representative that is literally the Mellin transform of the plane wave, and the ω-integral is then the scalar one. That is a physics-level argument: distributional δ-normalization, no domain statements, modulo gauge, nonzero energy. For d = 2 (4D), the graviton principal series is Δ = 1 + iλ, (h, h̄) = (3/2 + iλ/2, −1/2 + iλ/2) for J = +2 (DPS 1810.05219 §1).

**PSS 1701.00049** gives the massless wavefunction only as the m → 0 limit (2.22), φ^±_{Δ,m=0} = 1/(−X·q ∓ iε)^Δ, "assuming Re(Δ) > 1", from the Bessel form (2.17). Its App. A, (A.1)–(A.5), computes the KG product for **massive, equal-mass scalars only**. It finds δ-normalizability whenever the two weights have the **same real part a** (A.4)–(A.5), with coefficient 1/(Δ₁ + Δ₂* − 2). App. A does **not** single out Re Δ = 1 and does not treat spin. The principal-series selection Re Δ = d/2 is PS 2017's (3.1)–(3.8) for the massive case and (4.7)–(4.9) for the massless one. So the graviton theorem to cite is **PS 2017 §6.3 (eqs. 6.16–6.20), resting on §4.2 (eqs. 4.6–4.9)**.

Our symbolic check `r17b_spin2_sympy.py` (CHECKED) verifies the reduction for the Mellin representative ε⁺_με⁺_ν(−q·X)^{−Δ}, ε⁺ = ∂_z q. It is traceless, transverse and satisfies □ = 0 for every Δ. Under the rotation z → e^{iθ}z the polarization picks up e^{−2iθ} (|J| = 2), and under the boost z → λz the profile scales as λ^{−Δ}. The ω-integral is the scalar Γ(Δ)/a^Δ. **The helicity sits entirely in εε, and Δ only in the scalar Mellin factor.**

**1.2 Special graviton weights (d = 2).**
- **Pure diffeomorphisms (PS 2017 §6.2, Table 4).** h^{Δ=0} (6.9) and h^{Δ=1} (6.10) are pure diffeomorphisms in any d. For d = 2, h^{Δ=1} equals its own shadow (6.11), and the shadow of Δ = 0, which has dimension 2, is also a pure diffeomorphism (6.12). PS note that "Δ = 0 is the shadow dimension of the stress-tensor".
- **Conformally soft (DPS 1810.05219).** "Conformally soft" means h = 0 or h̄ = 0. The celestial stress tensor comes from a (2, 0) primary wavefunction. Supertranslations are generated by a (3/2, 1/2) operator whose OPE shifts (h, h̄) → (h + ½, h̄ + ½), i.e. Δ → Δ + 1. At Δ = 1, the Goldstone mode and the conformally soft (memory) mode form the canonical pair (DPS §4.2).
- **Stress tensor (KMRS 1609.00282).** T_zz is built from the subleading soft graviton.
- **Tower (GHPS 2103.03961).** Hᵏ = lim_{ε→0} ε G⁺_{k+ε}, k = 2, 1, 0, −1, … The Δ = 1 current generates supertranslations, and the Δ = 0 triplet's shadow is the Vir_R stress tensor. In momentum language, a soft factor ω^{n−1} gives a pole of the celestial amplitude at Δ = 1 − n: Δ = 1 (Weinberg leading), 0 (subleading), −1 (sub-subleading).
- **Shadow and spin-shadow.** The 2D shadow is (h, h̄) → (1 − h, 1 − h̄), i.e. (Δ, J) → (2 − Δ, −J). In PS's real-index form it is (6.6): h̃^{Δ} = (−X²)^{d/2−Δ} h^{d−Δ}, with the same SO(d) spin-2 representation. In helicity components, +2 ↔ −2.

**1.3 Daniel's dictionary (draft ch. 7).** The dictionary is Δ = 2s: the critical line goes to Re Δ = 1, s ↔ 1 − s to Δ ↔ 2 − Δ, and the pole s = 1 to Δ = 2. GPPVerify formalizes the same coordinate map as `celestialDelta s = 2s` (CelestialRiemannAffineBridge.lean:23) and as the diagonal lift (h, h̄) = (s, s) (DiagonalConformalLift.lean). That lift is **spin 0**. For the graviton, J = ±2 gives (h, h̄) = (s ± 1, s ∓ 1), with Δ = h + h̄ = 2s unchanged.

## 2. "ζ's s IS the graviton's Δ/2": the map

**2.1 Kinematic lemma (PROVED).** For a massless wavepacket Φ = ∫ d³p/((2π)³2ω) a(ω) g(ŵ) εε e^{ip·X}, PS (4.7)/(6.17) give a KG/symplectic norm of the form C ∫|g|²d²w · ∫₀^∞ |a(ω)|² ω dω. (Here d³p/2ω ∝ ω dω d²w. In PS's dω-density b = ωa, (4.7) reads ∫|b|²ω^{−1}dω.) The overlap with the conformal primary of weight Δ is ∝ g(w) ∫ω^{Δ̄−1} a(ω) dω.
Put y = ω^k. Then ω dω = (1/k) y^{2/k−1} dy and ω^{Δ−1}dω = (1/k) y^{Δ/k−1} dy, so s = Δ/k. The measure is flat (Lebesgue, Plancherel line Re s = ½) **iff k = 2**, and then Δ = 2s and Δ ↔ 2 − Δ becomes s ↔ 1 − s. In general d the same argument gives k = d = D − 2, which reproduces the α(D) = D − 2 table in ch. 7 "Why four dimensions is special". The variable that carries s is **y = ω²** (= |p⃗|² for massless; = |κ|⁴ in spinor-helicity variables, ω = |κ|²).

**2.2 The lattice graviton primary.**
- *Raw form (CHECKED to 14 digits, three Δ).* r_ζ(ω) = Σ_{n≥1} e^{−π n² ω⁴} = Σ_n e^{−π (n y)²}, so the integers scale y = ω². Then ∫₀^∞ ω^{Δ−1} r_ζ dω = ¼ Γ_R(Δ/2) ζ(Δ/2), where Γ_R(s) = π^{−s/2}Γ(s/2).
- *Completed state (PROVED; CHECKED to 1e−22 at nine Δ, including off-line ones, against mpmath closed forms).* a_ζ(ω) = 2 f_ζ(ω²) = 2ω^{−1}Φ(2 log ω), where f_ζ = x^{−1/2}Φ(log x) is Round 17's BPY state and Φ is 2× Pólya's Φ. Then **A_ζ(Δ) = ∫₀^∞ω^{Δ−1}a_ζ(ω)dω = ξ(Δ/2) for all Δ ∈ ℂ**. Proof: ω = e^{v/2} turns the integral into ∫Φ(v)e^{(Δ/2−½)v}dv = ξ(Δ/2) (laws17.py convention).
- *Graviton state.* Φ^±_ζ = ∫ d³p/((2π)³2ω) a_ζ(ω) g(ŵ) ε^±ε^± e^{ip·X}, in the Mellin diffeomorphism representative (6.18). By (6.17) and (6.20), (Φ^±_ζ, ϕ^{Δ,±}(·; w)) = c_± g(w) ξ(Δ̄/2), with c_± ≠ 0 a convention constant (PROVED at PS's level of rigour). Because ξ is real on ℝ, the overlap vanishes iff ξ(Δ̄/2) = 0. The zero set is conjugation-symmetric, so it is {Δ : ξ(Δ/2) = 0}.
- *Normalizability (CHECKED to 1e−26).* Using the dω-density b = ωa (so ∫|b|²ω^{−1}dω = ∫|a|²ω dω), ‖Φ_ζ‖² ∝ ∫|a_ζ|²ω dω = 2∫Φ² = 0.639503623923243, which equals (1/2π)∫|A_ζ(1 + iν)|²dν. This is Mellin–Plancherel exactly on PS's normalizable line.
- *Shadow.* A_ζ(2 − Δ) = A_ζ(Δ) ⇔ ξ(s) = ξ(1 − s). With equal profiles in both helicities, the state is invariant under the spin-shadow (Δ, +2) ↔ (2 − Δ, −2).
- *Unitary equivalence with Round 17 (PROVED).* U: L²(ℝ₊, dy) → L²(ℝ₊, ω dω), (Ub)(ω) = √2 b(ω²), is unitary and carries Round 17's d = 1 primary x^{−s} to the graviton energy profile with Δ = 2s. **Every statement of Round 17 transports verbatim, and no new constraint can come from the radial/energy data.**

**2.3 The s = Δ alternative fails (CHECKED).** An integer-frequency tower (k = 1) gives ∫ω^{Δ−1}/(e^ω − 1)dω = Γ(Δ)ζ(Δ) (Planck) and ∫ω^{Δ−1}Σ_n e^{−πn²ω²}dω = ½Γ_R(Δ)ζ(Δ). Both are CHECKED to 14 digits. Their functional equation is Δ ↔ 1 − Δ, and the critical line lands at Re Δ = ½, where PS (4.7) diverges (c ≠ d/2). So among the energy-Mellin readings, **only Δ = 2s is compatible with the graviton's normalizable line**. This is a real, falsifiable check on Daniel's dictionary, and the dictionary passes.

**2.4 Normalizability versus zeros.**
- (i) PS: h^{Δ,±} is δ-normalizable ⇔ Re Δ = 1 ⇔ Re(Δ/2) = ½. This is a statement about **basis labels**.
- (ii) RH ⇔ every Δ with A_ζ(Δ) = 0 has Re Δ = 1, i.e. every graviton primary that the ζ-graviton state cannot see is a physical (normalizable) one.
- Normalizability of Φ_ζ gives A_ζ ∈ L²(1 + iℝ) and nothing else. An entire function that is L² on a line can vanish anywhere off it (EPS, §3). So the conclusion of the PS theorem is "**the physical Δ live on Re Δ = 1**", not "the zeros are on Re s = ½".
- **GAP (= RH):** (ii) is `GravitonZeroIdentification` (Lean, OPEN), with `rh_of_graviton_zero_identification` PROVED conditionally. Round 17's `XiZeroIsNormalizable1dState` is the same premise after the unitary U of §2.2.

**2.5 Can the graviton-specific structure supply (ii)? Each tested.**

| Structure | What it does under Δ = 2s | Constrains zeros? | Label |
|---|---|---|---|
| Helicity ±2 / spin-2 polarization | Enters only through ε^±ε^± (sympy). Δ-dependence is the scalar Mellin; PS §6.3 reduces to §4.2. | No | PROVED |
| The factor 2 | The measure d³p/2ω is Lebesgue in y = ω² (§2.1); fixes Δ = 2s and Γ_R(Δ/2). | No (kinematic) | PROVED |
| Soft theorems / conformally soft poles | Raw lattice graviton: r_ζ = 1/(2ω²) − ½ + O(ω^{−2}e^{−π/ω⁴}) (CHECKED: residual 6e−22 at ω = 0.5). Poles only at **Δ = 2 (res ½) and Δ = 0 (res −½)**. Finite at Δ = 1, −1, −2, −3 (the Γ(Δ/4) poles at −4, −8, … cancel against trivial zeros). So ζ's poles s = 1, 0 are the stress-tensor and superrotation weights, and Poisson summation truncates the soft expansion. | No: soft data live at integer Δ ≤ 2, and EPS has identical soft data | CHECKED; reading HEURISTIC |
| Weinberg leading soft (Δ = 1 pole) | **Absent** for the lattice graviton: its IR is ω^{−2}, not ω^{−1}. So it is not a Weinberg-soft amplitude. | — | CHECKED |
| Memory / Goldstone, Δ = 1 (self-shadow) | Maps to the centre s = ½: A_ζ(1) = ξ(½) = 0.497120778… ≠ 0. The ζ-graviton couples to the memory channel. | No (one point) | CHECKED |
| Superrotation / stress tensor, Δ = 0, 2 | A_ζ(0) = A_ζ(2) = ξ(0) = ξ(1) = ½ | No | CHECKED |
| BMS supertranslation Ward identity, Δ → Δ + 1 (DPS) | Becomes s → s + ½, sending the critical line to Re s = 1, where ξ ≠ 0 (Hadamard–de la Vallée Poussin). But the shifted object is the overlap of ωΦ_ζ, a different state, with zeros at Δ = 2ρ − 1. | No | PROVED (algebra), HEURISTIC (reading) |
| Poincaré covariance | The zero set of A_Φ is invariant under Lorentz and dilations (a boost multiplies A by λ^{−Δ}) but **not** under translations (e^{ip·a} changes ∫ω^{Δ−1}a(ω)e^{iωq·a}dω). The zero set is not a Poincaré/BMS invariant of the state, so translation or supertranslation Ward identities alone cannot pin it. | No | PROVED |
| c = 0, Virasoro, w₁₊∞ | Act on the celestial weights, not on a fixed state's overlap zeros | No | HEURISTIC |

## 3. Epstein test (x² + 5y², class number 2, no Euler product)

The same construction, a_M(ω) = 2f_M(ω²), is applied to each model (`r17b_graviton.py` → `r17b_graviton.json`).

| Step | ζ | EPS | DH | F_θ |
|---|---|---|---|---|
| S1 radial profile from integer data | ✓ | ✓ (θ(x/√5)θ(√5x), class-group theta) | ✓ | ✓ |
| S2 A_M(Δ) = ξ_M(Δ/2), all Δ (CHECKED vs closed forms; EPS: Z_Q = ζ·L(χ₋₂₀) + L(χ₋₄)L(χ₅)) | ✓ 1e−22 | ✓ 1e−20 | ✓ 1e−20 | ✓ (PROVED; convolution of ζ's) |
| S3 KG-normalizable, KG norm = Plancherel on Re Δ = 1 | 0.6395036 = 0.6395036 | 0.7018925 = 0.7018925 | 2.5216260 = 2.5216260 | ✓ |
| S4 shadow Δ ↔ 2 − Δ | ✓ | ✓ | ✓ | ✓ |
| S5 soft poles exactly at Δ = 2, 0 (raw lattice) | ✓ (½, −½) | ✓ (raw identity CHECKED; r_E = 1/ω² − 1 + exp. small) | none (ξ_DH entire) | ✓ |
| S6 positive profile | ✓ | ✗ (a_EPS(1) = 2Φ_EPS(0) < 0, Round 17 certificate) | ✓ | ✓ |
| **S7 zero ⇒ normalizable primary** | **OPEN (= RH)** | **✗**: ρ = 0.932969697 + 15.668249531i, Δ = 2ρ = 1.865939 + 31.336499i, A_EPS(Δ) = 4e−32 by ω-quadrature (scale at Re Δ = 1, same height: 2.3e−8) | ✗ (off-line zeros, Round 17) | ✗ (PROVED) |

**Verdict.** S1–S5 hold for EPS. S6 fails for EPS, but it is not a physical requirement of a graviton wavepacket, and F_θ passes it anyway. **Only S7 must use the Euler product.** That Euler-product input is GPPVerify's `GravitonZeroIdentification`, still OPEN (not yet located). Where the Euler product sits in the graviton picture (PROVED, formal, Re s > 1) is this. In the raw profile, the n-th harmonic is the dilation y → ny of y = ω². Multiplicativity of n gives Σ_n U_n = Π_p (1 − U_p)^{−1} over commuting prime dilations. EPS's r_Q is not multiplicative, so its sum does not factor. Whether that factorization can force S7 is the open question (GAP-17B). Round 16's two halves (c(n) = Λ(n) ≥ 0, Poisson on ℤ) are the known necessary ingredients. Poisson already appears here as the two-term soft truncation (S5), and ζ shares it with EPS.

## 4. Plain-words answer to Daniel's correction

Using the graviton instead of a generic Δ is the right physics. It explains *why* Δ = 2s: the massless on-shell measure is flat in ω². It tells us *which* ζ object is the graviton's amplitude: ξ(Δ/2), with no freedom. It even rules out the tempting "Planck gas" reading. But it is a change of variables away from Round 17's statement. So the graviton theorem (PS 2017 §6.3) tells us where ζ's critical line has to sit in the graviton's Δ-plane, namely on the physical principal series. It does not tell us that ζ's zeros are on it. What's missing is still one statement: every graviton primary invisible to the ζ-graviton state is physical. That holds for ζ only if RH does, fails for Epstein, and has to come from the Euler product.

## 5. Notes for the draft and the team (no edits made)

- **Ch. 7, Lemma "Intertwining Isometry" (≈ lines 15420–15452).** The celestial Mellin there is taken in κ = √ω, so its Δ is conjugate to the spinor modulus, while the prose says "Δ conjugate to energy ω, s conjugate to |κ|", which would give Δ = s/2. With PS (4.4)/(6.20), where Δ is conjugate to ω, s is carried by y = ω² (= |κ|⁴). The conclusion Δ = 2s stands. A cleaner proof is §2.1: the measure d³p/2ω is Lebesgue in ω². The Haar measures d^×ω and d^×κ used in the lemma have no Plancherel offset, so on their own they cannot fix it.
- **Ch. 7 Ghost Theorem.** It is right that the basis labels lie on Re Δ = 1 (PS (4.7)–(4.9), (6.16)–(6.20)). The last sentence, "Under Δ = 2s, this locus is Re(s) = ½", is about labels. It does not cover ζ's zeros: that is step (ii), §2.4.
- **Ch. 5 ≈ line 8770** pairs O^{+2}_Δ with O^{−2}_{4−Δ}. The 2D shadow partner of (Δ, +2) is (2 − Δ, −2). The pairing 4 − Δ needs its own justification.
- **GPPVerify tension (flag only).** GravitonDeltaIdentification.lean's round-6 transport sOfCelestial(Δ) = (Re Δ/2, Im Δ) is not holomorphic in Δ, so it cannot intertwine Mellin transforms. The energy-Mellin pairing here realizes the whole-variable map Δ = 2s (= CelestialRiemannAffineBridge's `celestialDelta`) exactly, i.e. ν = 2t. The round-6 "Check D" falsification of the whole-variable map should be re-examined against this (OPEN for the team).
- **Lean targets.** (a) `∫ω^{Δ−1} b(ω²) dω = ½ M[b](Δ/2)` and the unitary U of §2.2. (b) The graviton overlap identity A_ζ(Δ) = ξ(Δ/2) via the BPY state (still not formalized, see ROUND17.md §0).

## Files (this folder)
- `r17b_graviton.py` → `r17b_graviton.json`, `r17b_graviton.log`: A identity (all models, nine Δ); B KG vs Plancherel; C power maps and k = 1 Planck/linear-theta; D raw lattice graviton soft poles (ζ, EPS); E special weights; F EPS invisible non-normalizable primary. mpmath at 25 digits, about 3 min.
- `r17b_spin2_sympy.py` → `r17b_spin2_sympy.json`: spin-2 Mellin representative is traceless, transverse and satisfies □ = 0; helicity phase and boost weight; scalar Mellin factor.

## References
S. Pasterski, S.-H. Shao, "Conformal Primary Basis for Flat Space Amplitudes", PRD 96 (2017) 065022, arXiv:1705.01027, §§3.1–3.3, 4.1–4.2, 5.3, 6.1–6.3, eqs. (3.6)–(3.8), (4.4), (4.6)–(4.12), (6.4), (6.6), (6.9)–(6.12), (6.16)–(6.20), Table 4 · S. Pasterski, S.-H. Shao, A. Strominger, "Flat Space Amplitudes and Conformal Symmetry of the Celestial Sphere", PRD 95 (2017) 065026, arXiv:1701.00049, §2.2 eqs. (2.15)–(2.22), App. A eqs. (A.1)–(A.5) · L. Donnay, A. Puhm, A. Strominger, JHEP 01 (2019) 184, arXiv:1810.05219 · D. Kapec, P. Mitra, A.-M. Raclariu, A. Strominger, PRL 119 (2017) 121601, arXiv:1609.00282 · A. Guevara, E. Himwich, M. Pate, A. Strominger, arXiv:2103.03961 · Biane–Pitman–Yor, Bull. AMS 38 (2001) · Round 17 (`3383dd6`) · GPPVerify main `7ebc350`.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).
