# Round 11 literature pass: Pasterski papers vs the AKG target

Heart (Grok / GPPDiscovery5). Sources: Daniel's Drive folder GPP / "The Literature", read with pdftotext.
**Target (from README.md §4):** AKG asks for a Hilbert space H and a map c ↦ Φ_c, built only from the primes and the archimedean term, with ‖Φ_c‖²_H = kform(arch) − kform(prime) − 2Ψ(i/2)²S² on ker ℓ_c for every L. AKG is OPEN and RH-equivalent.
Labels: PROVED / CHECKED / NUMERICAL / CONJECTURE / GAP / OPEN. **Result of the pass:** no paper supplies AKG, and none supplies a new candidate T that evades obstruction O2 (Bochner forcing). Two items sharpen Round 11 (L1, L2 below). Neither is a new theorem, so nothing was posted to PR #11.

## Per paper (two lines each)
- **2108.04801**: Pasterski, *Lectures on Celestial Amplitudes* (SAGEX 2021).
  - What it says: radiative data is the principal series Δ ∈ 1 + iR, continued analytically to Δ ∈ ℂ (§4, eq. 4.19). Eq. 4.22 gives the IR piece of the Mellin integral as ∫₀^Λ ω^{Δ−1+a}dω = Λ^{Δ+a}/(Δ+a), which has poles at Δ = 1, 0, … "on and to the left of the principal series".
  - Relevance: CHECKED support for R11.5. It is the same ∫₀¹ + ∫₁^∞ split as R11 P0, and off-principal modes appear as Mellin poles, not as states. It does not supply AKG.
- **2111.11392**: Pasterski–Pate–Raclariu, *Celestial Holography* (Snowmass white paper).
  - What it says: normalizable finite-energy wavepackets correspond to the spectrum Δ = 1 + iλ (eq. 3.3, citing PS 1705.01027). Soft and collinear structure appears at special Δ off the line.
  - Relevance: CHECKED, but nothing new beyond 1705.01027. It neither supplies nor rules out AKG.
- **2201.01630**: Pasterski–Verlinde, *Chaos in Celestial CFT*.
  - What it says: the CCFT Hilbert space comes from radial quantization. The inner product is *chosen* so that L₀ + L̄₀ is anti-Hermitian on principal-series primaries, which makes dilation U(t) unitary (§4.1). Local operators O(τ) = ∫dλ e^{iλτ}/(2 sinh πλ) O(λ) carry KMS periodicity τ → τ + 2πi (Rindler temperature 1/2π).
  - Relevance: the Hilbert space is explicitly a "guiding assumption", i.e. AKG's analogue is postulated, not built. The resulting thermal/KMS kernel is translation-invariant in τ, so any Q = A*A built from it falls under O2. It gives no new T.
- **2204.02505**: Mizera–Pasterski, *Celestial Geometry*.
  - What it says: App. B lists the unitary irreps of SL(2,ℂ): principal (Δ = 1 + iλ), supplementary (J = 0, Δ ∈ (0,2)∖{1}), "discrete" (J = 0, Δ ∈ ℤ∖{1}, their wording) and trivial (Δ = 1).
  - Relevance: this is the source of L2. Unitarity alone does not force Re Δ = d/2, because the complementary series exists. The trivial representation sits where ζ has its poles (L2c).
- **2310.04932v2**: Pasterski, *A Chapter on Celestial Holography* (already cited in R11).
  - What it says: Mellin residues sit at Δ = 1, 0, −1, …, "on and to the left of the principal series", and translations shift Δ → Δ + 1, so "we are forced to contend with a spectrum analytically continued to complex Δ".
  - Relevance: CHECKED support for R11.5. It does not supply AKG.
- **2411.10527v4**: May–Pasterski–Waddell–Xu, *Cryptographic tests of the python's lunch conjecture*.
  - What it says: it concerns AdS/CFT entanglement wedges and complexity, and contains no Mellin, boost-basis or zeta content.
  - Relevance: none to AKG.

## L1: invariant forms on mode pairs (PROVED, elementary; restates R11 P0 as a uniqueness statement)
Let dilations u ↦ u + a act by e_Δ ↦ e^{(Δ−1/2)a}e_Δ. A sesquilinear form with ⟨e_x, e_y⟩ ≠ 0 is dilation-invariant iff e^{(x + ȳ − 1)a} = 1 for all a, i.e. y = τx = 1 − x̄.
- On span{e_ρ, e_τρ} with τρ ≠ ρ, every invariant Hermitian form is [[0,h],[h̄,0]]. Its eigenvalues are ±|h|, so it is never positive unless h = 0.
- A single off-line mode therefore admits *no* nonzero invariant positive form, and the only invariant positive forms live on the line.
- This is the precise content of "normalizability = Re Δ = d/2" for the abelian dilation group (Daniel's step 2, which is correct as a statement about invariant forms). The open part is unchanged: whether the zero modes carry such a form at all, which is AKG.

## L2: what SL(2) unitarity adds, and the loophole it opens (CHECKED, classical: Bargmann 1947 / Gel'fand–Naimark; quoted in 2204.02505 App. B)
- **(a)** For the spherical principal series π_s of SL(2,ℝ) (d = 1), π_s ≅ π_{1−s} and the Hermitian dual of π_s is π_{1−s̄}. An invariant positive form exists iff π_s is equivalent to its Hermitian dual, i.e. Re s = 1/2 (principal) or s ∈ (0,1) real (complementary series), with the trivial representation at the endpoints s ∈ {0,1}.
- **(b)** So, unlike the abelian case L1, non-abelian unitarity alone allows real Δ off the line: the complementary series. Daniel's note's sentence "normalizability ... excludes every other weight" is right for the KG/Mellin norm (abelian dilations, L1). It is not true of unitarity in general: the complementary-series inner product is the τ-pairing composed with the shadow intertwiner, and it is positive for real Δ ∈ (0,1).
- **(c)** For ζ the loophole is closed classically: ξ has no real zeros in (0,1) because ζ(σ) < 0 there (PROVED, classical). In this language the poles s = 0, 1 are the **trivial representation**, the constant function, i.e. the residue of the Eisenstein series at s = 1. This matches R11.1, where ker ℓ_c is the removal of the trivial representation.
- **Consequence:** if every zero mode generated a unitary SL(2,ℝ) representation, RH would follow (L2a + L2c). The open step is unchanged: realizing the zero modes as vectors of a unitary representation. That is AKG in representation-theoretic form, and it is RH-equivalent. In the one existing SL(2) realization (Lax–Phillips on SL₂(ℤ)\H), the zeros of ζ(2s) appear as poles of the scattering matrix, i.e. resonances, not as unitary constituents (O3).

## Candidates examined, labelled
- **KMS/Rindler kernel from 2201.01630.** The weight 1/(2 sinh πλ) with KMS periodicity 2π does not apply. Our K is even in τ, so K̂ is even, and the KMS relation K̂(−λ) = e^{−2πλ}K̂(λ) forces K̂ = 0 at every λ ≠ 0. With the KMS weight removed, the kernel is translation-invariant, so O2 applies. It is not a new T, and no numerical test is meaningful.
- **Light-sheet / detector Q from 2307.16801** (R11.6): already shown to be the zero measure itself (O2).
- **Complementary-series intertwiner as T** (L2): the intertwiner is positive only for real Δ. For off-line ρ with Im ρ ≠ 0 the Hermitian dual is inequivalent (L2a), so no intertwiner makes the hyperbolic block positive. This is killed as a route to positivity off the line (PROVED via L2a).

## Net for the joint plan
AKG remains the single named target (OPEN, RH-equivalent). The literature confirms the R11 reading, in which off-principal modes appear as Mellin poles or resonances and the Hilbert space is postulated rather than built. It adds the representation-theoretic restatement L2, which shows that the only extra freedom unitarity allows (the complementary series) is already excluded for ζ. Any construction of H must therefore realize the zero modes inside a unitary representation, and must break translation invariance (O2) to do so.

Credit: Daniel Toupin / Golden Physics Project (sources and the Re(Δ) = d/2 argument). This seat: Heart of the Cosmos / Grok (GPPDiscovery5). Audit: Turbo.
