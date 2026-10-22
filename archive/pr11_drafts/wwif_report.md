# *Which Way Is Forward?* v26: read and assessed (Grok, 2026-10-14, read-only)

## 0. Where v26 is, and which versions exist
- **Copies of v26:**
  - Drive: `Which_Way_Is_Forward_v26.tex` (id 1xt48ePgGB8jD1n0ftYhue8w93HdpzoCI) and `Which_Way_Is_Forward_v26.pdf` (id 1LT50dXN7tQe-HvbkD3JzIre89DMBlbNM), both modified 2026-10-08.
  - Repo: `GPPVerify/discovery/which_way_is_forward/Which_Way_Is_Forward_v26.tex` (3,954 lines). It was rebuilt from Drive's text export, and its README says "If it disagrees with the PDF in Drive, the PDF is authoritative".
  - The paper is dated 22 September 2026 (README:3). I read the repo `.tex`, not the PDF.
- **Other versions on Drive** (the first 50 title hits; more pages exist and I did not page through them):
  - v4, v5, v6, v7, v8, v9 (several copies), v10, v11: 2026-09-25;
  - v21, v22, v23 (several copies): both 2026-09-25 and 2026-10-08 copies;
  - v25 and v26: 2026-10-08;
  - also `_CJP`, `_Rigorous`, `_Preprint`, `_Preprint_TikZ`, `_UPDATED`, `_FINAL`.
  - Copies of v21 to v26 sit in `GPP/Newer rh physics stuff`.
- **The paper's own check script** `verify_which_way_is_forward.py` runs green on the box. All nine checks pass ("All exact-core checks passed"); they are exact or floating-point arithmetic, not interval certificates.

**Important finding.** v26 is mainly a paper about *orientation, matter, mass and time direction*, not about ζ.
- Its full title is *"Which Way Is Forward? Orientation, Matter, Mass, and the Two-Lift Structure of Relativistic Physics"*.
- ζ and RH appear in a single subsection (§"One sphere, two orientations", tex:1272-1357).
- The 1D-CFT framing you described (Re Δ = d/2, s as Δ of a 1D CFT) is **not** in v26:
  - no hits for "d/2", "1D CFT", "one-dimensional CFT", Berry–Keating or Connes;
  - v26 still writes "the bookkeeping normalization Δ = 2s" (tex:1308).
- That framing lives in GPPVerify instead (§4) and in Turbo's notes (`MuseResearch/turbo-note-2026-10-12-workbench-state.md:19,68`: "Dimensional dictionary `Δ = d·s` (2D celestial `Δ=2s`, 1D arithmetic `Δ=s`)"; TURBO_HEURISTICS.md:31).
- The later Codex note `work/codex/2026-10-13_cayley_mirror_circle_infinite_length.md:3` takes as its starting point "Daniel Toupin's one-dimensional CFT boundary as an infinite line compactified to a circle".

## 1. Central claims and definitions in v26 (quoted)
- **The question (tex:121-123):** "If every standard used to define an orientation is transformed together with the system, what observable distinguishes the two complete descriptions?"
- **Operational criterion (tex:127-131):** "x and ιx represent the same operational state precisely when O(ιx)=O(x) for every O ∈ 𝒜_op. If one admissible O violates [this], the proposed quotient is physically false." It is called "almost tautological, and that is its advantage" (:134).
- **The core object (abstract, tex:73-79):**
  - q = ±1 is the charge-representation orientation, and t = ±1 is the orientation of the line relative to an observer ("the arrow along which records accumulate").
  - "χ = qt, (q,t) ∼ (−q,−t), Q_obs = Q⋆χ".
  - So (++) and (−−) carry electron charge, while (+−) and (−+) carry positron charge.
  - "t = sgn(−u·k)": the time sign is a relative sign, so no second time coordinate is needed.
- **The finite carrier.** "Its operators I_q, I_t, D generate a group of order 16 with center {1,−1,χ,−χ}"; D is "an involutive real structure whose fixed carrier is {(a,b,b,a)}". It is realized as "the complexified spin module of Cl(2,2), with χ equal to the negative split volume element" (abstract; status, tex:106).
- **The orientation-lift hypothesis (Def., tex:1262-1268).** "A double cover π: X̃ → X, ι² = 1, π∘ι = π, such that every operational observable factors through π." The paper labels the extension honestly: "The claim that charge/time orientation, chirality, duality sector and celestial orientation are all images of the same ι is a conjectural extension and must be demonstrated by commuting maps" (:1270).
- **The time arrow (abstract).** "𝒟_t S_coarse ≥ 0 ... invariant under complete temporal reversal ... thermodynamic irreversibility need not imply an absolute fundamental arrow of time." This answers "which way is forward": a direction is defined only relative to each observer's own record-forming arrow.
- **The critical-line part (Prop. "Anti-linear fixed locus", tex:1311-1316).**
  - "Let ℛ(Δ) = 2−Δ̄. Then ℛ² = 1 and Fix(ℛ) = {Re Δ = 1}."
  - Then: "Under the bookkeeping normalization Δ = 2s ... s ↦ 1−s̄ with fixed locus Re s = 1/2. This is a fixed-locus identity, not a proof of the Riemann hypothesis" (:1308).
  - The comparison table, D ↔ ρ(ν) = −ν̄ and R_t ↔ σ(ν) = −ν (:1340-1352), "is an exact comparison of two involution architectures, not an identification of their physical carriers and not a proof of RH" (:1353).
- **The principal-series Hilbert lift (tex:1358-1415).**
  - "d³p/2p⁰ = 2ω dω d²z ⇒ Δ = 1 + iλ ... the unitary principal-series line in the massless one-particle basis; it is not a dynamical equation for Δ" (citing Pasterski–Shao).
  - Prop.: R̂² = (−1)^{2J}, a real structure for integer spin and a quaternionic structure for half-integer spin.
- **The geometric and unification part (status, tex:108-110):**
  - Grassmannian order-four map, Pin lift, the doubled parent W^{6,2}, and ℝ^{6,2} ⊕ ℝ^{4,0} = ℝ^{10,2}.
  - "Until those are derived, Spin₀(10,2) is a sharply constrained completion candidate, not a completed unified theory."
- **What the paper explicitly does not claim (tex:105):** "No theorem below identifies bare celestial shadow with Wigner time reversal, turns zitterbewegung into a hidden classical trajectory, solves baryogenesis or closes the nonlinear googly problem."
- **Epistemic ledger (tex:3647-3727).**
  - About 50 rows. Most are "theorem; finite model / Clifford classification / prior art". Prior art cited includes Stueckelberg 1941 and Feynman sgn E = t.
  - "Complete orientation conjugation is operationally redundant" is labelled **hypothesis**.
  - "Googly is an orientation reconstruction problem" is labelled **programme**.

## 2. The 1D-CFT picture, as formalized in GPPVerify rather than in v26
- `GppVerify/CelestialHolography/ConformalShadowPrincipalSeries.lean` (last change 2026-10-13, c29290c). Quoted:
  - `def conformalShadow (d : ℝ) (Delta : ℂ) : ℂ := (d : ℂ) - Delta`
  - `theorem conformalShadow_eq_conj_iff (d : ℝ) (Delta : ℂ) : conformalShadow d Delta = (starRingEnd ℂ) Delta ↔ Delta.re = d / 2`
  - `theorem one_dimensional_shadow_eq_arithmetic (s : ℂ) : conformalShadow 1 s = 1 - s`
  - `theorem one_dimensional_shadow_eq_conj_iff (s : ℂ) : conformalShadow 1 s = (starRingEnd ℂ) s ↔ s.re = 1 / 2`
  - `theorem doubled_arithmetic_shadow_eq_two_dimensional (s) : celestialDelta (1 - s) = conformalShadow 2 (celestialDelta s)`
  - The docstring says: "These are representation-theoretic identities only; no zeta-zero statement is used."
- `GppVerify/RiemannHypothesis/ArithmeticConformalKinematics.lean`:
  - the half-density character χ_s(a) = exp(log a · (s − ½)) becomes ψ_s(x) = e^{x(s−½)} on log scale;
  - it is unitary (‖·‖ = 1) on s = ½ + iτ, and shadow sends τ ↦ −τ (that is, the conjugate character);
  - the docstring adds: "It does not by itself construct a full 1d CFT: no local operator algebra, OPE/crossing, stress tensor/Virasoro structure, or full PGL(2) action is asserted here."
- **So the picture, precisely.** In boundary dimension d = 1, a scaling weight Δ = s labels a dilation character x^{s−½} of the multiplicative group. Shadow Δ ↦ 1 − Δ equals conjugation exactly on Re Δ = ½, the principal-series line, where the character is unitary.
  - Mellin transformation turns ξ(s) into a function of that weight. A zero ρ = ½ + iγ then corresponds to a principal-series (unitary, tempered) weight Δ = ½ + iγ. An off-line zero would be a weight off the unitary axis, "complementary-series-like" but complex.
  - The 2D celestial version is the same statement under Δ = 2s, with Re Δ = 1 = d/2 for d = 2.
- **Time direction.** Log-scale reflection x ↦ −x (r ↦ 1/r) "becomes time-like reflection after logarithmic parametrization ... not a proof that radial time is physical proper time" (v26 tex:1288-1291). The arithmetic reflection s ↦ 1 − s̄ is the centered anti-linear map ρ(ν) = −ν̄, which the paper compares with D (tex:1325-1352).

## 3. Standard versus new
**Standard:**
- The principal series of SO(d+1,1) (the conformal group of ℝ^d) is Δ = d/2 + iν, and shadow Δ ↔ d − Δ is the Knapp–Stein intertwiner. For d = 1 this is the SL(2,ℝ) ≅ SO(2,1) principal series s = ½ + iν. Celestial Δ = 1 + iλ is Pasterski–Shao.
- That ρ ↦ 1 − ρ̄ preserves the zeros, and that Re s = ½ is the unitary axis of the dilation characters, is classical (Tate; the Mellin–Plancherel theorem).
- Spectral readings of the zeros as principal-series data are prior art: Hilbert–Pólya, Berry–Keating xp (a dilation generator, i.e. 1D conformal QM), Connes' adelic trace formula, and the Selberg/Eisenstein picture, where ζ zeros appear as scattering resonances (at depth ¼, as Codex noted).
- In the physics part: oriented Wilson lines, Stueckelberg–Feynman reversal, the charged-KG split, the Dirac/Compton phase, CAR algebra, the Cl(p,q) classification, Spin(10)/Pati–Salam branching, and Spin^c obstructions.

**New (Daniel's):**
- The *relational* quotient (q,t) ∼ (−q,−t) as the organizing hypothesis.
- The order-16 carrier and its Cl(2,2) intertwiner.
- The claim that several standard doublings descend from one orientation lift (a hypothesis).
- The Spin₀(10,2) completion candidate.
- The orientation-covariant reading of the second law.
- The architectural comparison of D with s ↦ 1 − s̄.
- The "1D CFT boundary" framing (d = 1, Δ = s), as a dictionary choice.

**Caveat on the 1D-CFT language.** In a unitary 1D CFT, such as de Alfaro–Fubini–Furlan conformal quantum mechanics, physical primaries have real Δ in the discrete series. Complex Δ = ½ + iγ labels principal-series representations or resonances, not unitary-CFT operator dimensions. So "s is the Δ of a 1D CFT" is exact only at the level of representation labels and characters, which is exactly what the Lean asserts. It does not give an actual 1D CFT whose spectrum is the ζ zeros.

## 4. What GPPVerify formalizes
- **Finite core** (`GppVerify/StandardModel/OrientationCliffordCore.lean`, kernel `decide` over ℤ[i]): the Cl(2,2) carrier, the orientation algebra, the four-lift intertwiner, the order-16 group with center {±1,±χ}, Fix(D) = {(a,b,b,a)} with the two half flips agreeing on it, and the Cl(3,1) realization (README).
  - Not yet formalized: the signature scan, the Spin(6,2)/Spin(10,2) BB̄ signs, the determinant cover (README).
- **The paper's Lean inventory is on main.** The README says it has 34 modules, all now on main and building. My regex pulled 30 module paths out of the appendix, and every one of them exists on main. They were ported from Codex's PR #173 branch, with repairs listed in the README. They include:
  - `RiemannHypothesis/OrientationCriticalRealStructureBridge.lean`, with `centered_shadow_eq_conj_iff_critical_real_form (nu) : centeredShadow nu = conj nu ↔ nu.re = 0`;
  - TauDifferential;
  - the Grassmannian, Klein, Pin and Spin-cover modules.
- **Related** StandardModel modules: FourOrientationGaugeProjection, CPTSymmetryDoesNotForceAntimatter, SpinCOrientationQuotient, Spin62BLFiniteCore, Spin102ChiralRealityCore.
- **1D-CFT framing:** ConformalShadowPrincipalSeries and ArithmeticConformalKinematics (§2).
- **Character of the results:** all finite algebra or one-line complex identities. They are correct, but they contain no dynamics, no QFT and no ζ positivity.

## 5. Honest assessment
**What holds:**
- Every listed theorem I checked is true and is stated with its own limits. The paper's three-tier labelling (theorem / dictionary / hypothesis) is well done.
- Re Δ = d/2 as "the locus where shadow equals conjugation" is correct and standard. It is the right correction to "Δ = 2s is a theorem": in d dimensions it is Δ = d·s, which makes d = 1 the natural arithmetic case.

**What is conjecture:**
- The orientation-lift hypothesis: that charge, time, chirality, duality and celestial orientation share one ι.
- That "complete orientation conjugation is operationally redundant" (ledger: hypothesis).
- The Spin₀(10,2) completion. And the zeros-as-1D-CFT-spectrum reading, which no theorem constructs.

**Where the falsifiers bite.** The entire ζ-facing content of v26 and the 1D-CFT modules uses only the involution s ↦ 1 − s̄ and its fixed line.
- **Davenport–Heilbronn.** Its completed form satisfies Λ(s) = Λ(1−s) (I checked numerically to 2.5e-32), so it has exactly the same shadow/conjugation structure. Yet it has a zero at 0.80852 + 85.6993i, off the unitary axis.
- **F_θ = ξ(s−θ)ξ(s+θ).** It also satisfies the functional equation and the reality condition. Its zeros at ½ ± θ + iγ are shadow pairs off the axis.
- **Consequence.** Nothing in the shadow / principal-series / 1D-character dictionary can force zeros onto Re s = ½. The dictionary *describes* what RH would mean (every zero is a unitary weight), but it cannot supply it. The Euler product, and positivity of the Weil form, must enter, and neither appears in v26.
- **The physics side has no such falsifier yet.** The paper's own test is the operational criterion: find an observable that is odd under complete conjugation (tex:3631 onward, "A programme of attempted destruction"; App. "A minimal falsification checklist", tex:3858).

**Bottom line.** v26 is a careful, well-labelled paper on relational orientation in relativistic physics, with an exact finite algebraic core. Its RH content is one honest fixed-locus remark. The 1D-CFT reframing is mathematically correct as kinematics and RH-neutral as a mechanism.
