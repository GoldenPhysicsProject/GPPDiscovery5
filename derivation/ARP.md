# ARP: Arithmetic Reflection Positivity, a first-principles derivation (working document)

Daniel Toupin / Golden Physics Project, with AI assistance from Grok ("Heart of the Cosmos") and the RH sub-team.
Started at S102, 2026-10-27, on Daniel's directive of 2026-10-22 07:48 EDT: "We need a first principles derivation of arithmetic reflection positivity."

**RH is OPEN. Nothing in this file proves or implies RH.**

**Labels.** Every step carries one of these:
- PROVED: literature or an elementary identity, with a citation.
- PROVED-LEAN: kernel-checked, with the GPPVerify name.
- CONJECTURE.
- GAP: the sharp open step.
- NUMERICAL: orientation only, never a proof step.

## 0. Honest framing

By Weil 1952 (and Bombieri 2000), positivity of the Weil form on every window [−L, L] **is** RH. So a valid derivation of arithmetic reflection positivity (ARP) for every L would prove RH.

The deliverable is therefore a **rigorous chain**:
- every step except the named GAPs is PROVED;
- each GAP is stated as sharply as we can;
- for each step we record where every known fake breaks.

A chain whose only GAP is "Weil positivity" in disguise is a restatement. We say so wherever that happens (see §0.3).

### 0.1 Requirements (every candidate step must meet them)
- **R1. Uniform in L.** No bounded window implies RH: by Goddess of PNT's DMV theorem (`team/prime-number-theory/finite_windows_dmv.md`), DMV/P_X systems agree with ζ on every bounded window.
- **R2. Euler product and functional equation together.**
  - Epstein fails: FE and modularity, no Euler product.
  - DMV fails: Euler product, no FE.
  - DH fails: FE, no Euler product.
- **R3. Not through the sign of the prime weights.**
  - F_θ (positive weights) fails; G_θ′ (signed weights 2Λ(n)cos(θ′ log n)) passes (Queen).
  - So positivity cannot come from weight signs.
- **R4. The positive Hilbert space must contain the soft mode / pole state.**
  - Fakes first fail on ζ's nodeless soft mode (S97/S98).
  - Its window ground state has |cos| = 0.995–0.9997 with the pole vector 2cosh(x/2) (Queen S1, NUMERICAL).
- **R5. Each step names where each fake breaks**, and checks that G_θ′ passes. The fakes: F_θ, DH, H₋ₜ, DMV, Epstein, P_X.
- **R6 (Golden Eye; CONJECTURE, well supported).** The closing step should be an exact identity (Q = ‖·‖² or a trace). What is PROVED is narrower: §S no-slack lemma (S(c)4, audited by Golden Eye on 9184613) excludes any MF-relative margin, i.e. no c > 0 with Fl_L ≥ −(1−c)·MF_L uniformly. It does not exclude Q_L ≥ c·N_L for another norm N_L degenerating on dilated nodeless bumps (e.g. weighting only frequencies near the zeros).

### 0.2 Dead ends (cite; do not re-derive)
| dead end | where |
|---|---|
| Log-time OS reflection positivity for ζ is false (Hankel ≠ Toeplitz) | bridge `work/grok/2026-10-21_reflection_positivity.md` (S82-OS); GPPVerify `NoGo/HaarSquareNotReflection.lean`, `NoGo/GaussianHeatOS.lean` |
| Schur / Kreĭn / de Branges / doubling reductions only restate window positivity | GROK_RESEARCH_NOTES S93–S97; Turbo waves 229–237 (S₊, (DOM±)), adjudicated S95–S96 |
| The positive cone is fake-blind (Λ-defects) | S98 |
| Selberg trace formula on the modular surface | scattering, non-arithmetic Hecke G₅ (S-series) |
| Connes' trace formula: the open step is uniformity in the cutoff | Connes 1999; Connes–Consani 2021 |
| η_ζ = 0 is RH-equivalent, not a route | S100 (Queen's quartet identity) |
| Positive prime weights / signed prime split as a forcing | S100–S101 (H100a) |

### 0.3 Inventory (build on these definitions; do not reinvent)
- **Lean OS/RP inventory.** Bridge `work/grok/2026-10-21_reflection_positivity.md` §0a/§0b, GPPVerify e1c1aa9: `reflectionForm_nonneg_iff`, `rh_iff_weil_pairedForm_nonneg`, `osReflection`, `ArithmeticOSGram`, `ArithmeticPrimeLocalOS`, `CriticalRiggedGram`, `SemigroupPosDef`, `PositiveDefiniteGNS`, `HaarSquareNotReflection`, `GaussianHeatOS`, `ShadowCutOSObstruction`, …
- **Lean names for this derivation.** Principal series `team/principal-series/ARP_plan.md` §1. This file serves as the "LEAN_INVENTORY" that ARP_plan.md waits on:
  - Tate local zeta integral, `euler_factor_bridge`;
  - `RationalIdeleGroup`;
  - `completedRiemannZeta_one_sub`;
  - `GppWeilLadder.primeSide_eq_truncation` / `primeSide_eq_zero_of_support_lt_log_two`;
  - `GppNoGhost.weilForm`;
  - new since Oct 22: `GppScalarSchurForm` (#412), `GppArchimedeanHankelGram` (#413), `GppSineSquarePacket` (#416), `GppGroundStateDetector` ((GS) ⇒ RH, #418, conditional).
- **Not citable as unconditional:**
  - `open_*` placeholders (`True := trivial`);
  - ArithmeticTemperedPrincipalSeries.lean (custom axiom `exp_growth_not_tempered`);
  - ArithmeticCoPoissonOrientationTemplate.lean (co-Poisson assumed).
- **Not in GPPVerify:** global Poisson summation on 𝔸/ℚ, Tate's global FE, the explicit formula (placeholder only), Connes' map E.

---

## 1. Sections by member (each owner fills their own; Grok audits and merges)
| § | owner (seat) | assignment |
|---|---|---|
| §P | **Principal series** (adelic / rep-theory skeleton) | Local reflections Θ_v: Tate local FE, normalized intertwiners, temperedness as a filter (not a source of positivity); the polarization lemma for F_θ on L²(C_ℚ). Plan already merged: `team/principal-series/ARP_plan.md` (incorporated by reference; its §2–§5 outline is the skeleton for §G.1–§G.2 below). |
| §Q | **Queen of Mathematics** (arithmetic OS axioms) | A precise axiom list (OS0–OS4 analogue) for an "arithmetic Euclidean field" on C_ℚ: which axiom encodes the Euler product (clustering/factorization over places), which encodes the FE (reflection), which encodes positivity. Show which axioms each fake satisfies. |
| §C | **Goddess Of Celestial CFT** (Bost–Connes / adele class space OS reflection) | OS reflection on the adele class space: the BC system's KMS_β states and the reflection x ↦ x⁻¹ on C_ℚ; what is positive at β = 1 (PROVED, RH-silent) and where the step to the zeros needs the gluing GAP (§G.4). Connect Astra's mirror-OS superselection (Oct 26) and doubled-prime OS test (Oct 27). |
| §S | **God of the Spectrum & Measure Theory** (Bochner / Kreĭn measure) | The Weil distribution as a Kreĭn-space (indefinite) measure: Bochner for the mean field (§G.3) and the signature of the fluctuation; finish the 2.84 replay (in-flight certificate, allowed). |
| §N | **Goddess of Prime Number Theory** (Weil form as a limit of finite Euler-product forms) | Q^{P_X} → Q^ζ as X → ∞; the negative index of Q^{P_X} on [−L, L] for L > (log X)/2; exactly what the FE adds that P_X lacks (R1/R2); the Ward/mean-field split of §G.3 for P_X. |
| §G | **Grok (Heart of the Cosmos)**: synthesis | "Time reflection for the primes": local at p, archimedean at ∞, global gluing by Poisson. Written below. |

### §P Principal series: *to be filled by Principal series* (see `team/principal-series/ARP_plan.md`)
### §Q Queen of Mathematics: *to be filled by Queen of Mathematics*
### §C Goddess Of Celestial CFT: *to be filled by Goddess Of Celestial CFT*
### §S God of the Spectrum & Measure Theory: *to be filled by God of the Spectrum & Measure Theory*
### §N Goddess of Prime Number Theory: *to be filled by Goddess of Prime Number Theory*

---

## §G Grok synthesis: time reflection for the primes

**Notation.**
- Log scale: u = log x on ℝ₊ˣ.
- Test functions: f ∈ C_c^∞([−L, L]), real; g = f ∗ f̃ with f̃(u) = f(−u); FT f̂(t) = ∫ f(u) e^{itu} du.
- **Weil explicit formula (PROVED; Weil 1952; Barner 1981 for the regularity used):**

  Σ_ρ ĝ(γ_ρ) = ĝ(i/2) + ĝ(−i/2) − Σ_{n≥2} Λ(n) n^{−1/2} (g(log n) + g(−log n)) + W_∞(g),  ρ = ½ + iγ_ρ.

- Weil form: Q_L(f) := Σ_ρ ĝ(γ_ρ). Then RH ⇔ Q_L(f) ≥ 0 for all L and f (Weil 1952).
- **Windows see only finitely many primes (PROVED-LEAN `primeSide_eq_truncation`):** window L sees exactly the n < e^{2L}.

### G.1 Local reflection at a finite prime p (Tate's local functional equation)

**G.1a (PROVED; Tate 1950, Thm 2.4.1).** For φ ∈ 𝒮(ℚ_p) and a quasi-character χ,

  Z(φ̂, χ⁻¹|·|, 1−s) = γ_p(χ, s) Z(φ, χ, s).

For χ unramified trivial (with ε_p = 1),

  γ_p(s) = L_p(1−s)/L_p(s),  L_p(s) = (1−p^{−s})⁻¹.

**G.1b (PROVED, elementary).** On Re s = ½, L_p(1−s) is the complex conjugate of L_p(s), so |γ_p(½+it)| = 1.
- Hence the local reflection Θ_p is **unitary** on L²(ℚ_p^*, d^*x) in the half-density normalisation f ↦ |x|^{1/2} f: Θ_p is p-adic Fourier composed with x ↦ x⁻¹, and in Mellin coordinates it is multiplication by γ_p(½+it).
- This uses Plancherel on ℚ_p^* (Tate 1950 §2).
- On spherical vectors it is the normalized intertwiner of the unramified principal series of GL(2) at the centre (Shahidi 1990; Langlands LNM 544).
- Replay: `s102/twist_and_ward_replay.py` (5).

**G.1c (PROVED, elementary): local positivity is not available.**
- The local prime distribution at p is W_p(g) = −Σ_{k≥1} (log p) p^{−k/2} (g(k log p) + g(−k log p)).
- Its spectral density is −2 Re[(log p) z/(1−z)], with z = p^{−1/2−it}, and

  Re[z/(1−z)] = (Re z − |z|²)/|1−z|²,

  which changes sign in t for every p.
- So **no place-by-place positivity statement exists** for the prime part. Weil's form is positive, if at all, only after summing over places.
- This is R3 made local: neither the weight signs nor the local spectral signs carry positivity.

**G.1d (PROVED): what Θ_p does give.**
- Unitarity means the local reflection loses no L² mass.
- The ⊗′ factorization 𝒮(𝔸) = ⊗′𝒮(ℚ_v), with F_𝔸 = ⊗′F_v, is where the **Euler product** enters (R2, first half).
- **R4 at this step.** The soft-mode/pole state is not local. The two boundary values φ(0) and φ̂(0) are global (§G.3).

### G.2 Archimedean reflection at ∞

**G.2a (PROVED; Tate 1950; Mellin).** On ℝ₊ˣ with half-density f ↦ x^{1/2} f, the map x ↦ 1/x is u ↦ −u (time reversal in log scale). In Mellin coordinates it is s ↦ 1−s, i.e. ½+it ↦ ½−it.
- The archimedean γ-factor is

  γ_∞(s) = π^{−(1−s)/2} Γ((1−s)/2) / (π^{−s/2} Γ(s/2)),

  and |γ_∞(½+it)| = 1.
- PROVED-LEAN for the algebra: `ArithmeticOSReflection.osReflection`, `osReflection_eq_self_iff` (fixed set Re s = ½).

**G.2b (PROVED; Yoshida 1992; Connes–Consani 2021, Selecta 27:77): base-case positivity.**
- For supp f ⊂ [−L, L] with L < (log 2)/2, no prime enters (PROVED-LEAN `primeSide_eq_zero_of_support_lt_log_two`).
- In that range Q_L(f) = pole + W_∞ ≥ 0.
- Connes–Consani identify the archimedean form with a Sonin-space trace plus a prolate correction. This is the **only window on which ARP is a theorem.**

**G.2c (PROVED-LEAN conditional; NoGo).** Plain log-time OS positivity of the *full* form is false. Toeplitz positivity does not give Hankel positivity: `HaarSquareNotReflection.lean`; one prime heat amplitude is not OS positive, `GaussianHeatOS.gaussian_hankel_det_neg`.
- So the reflection is the **Fourier-composed** inversion Θ = F ∘ inv (Tate), not bare u ↦ −u.

### G.3 The pole state and the mean field (R4): the positive part, PROVED

**G.3a. Astra's pole/first-prime Ward identity (PROVED; Astra, PR #11 Oct 27 22:05Z; bridge `work/2026-10-27_astra_prime_pnt_os_ward_identity.md`).**
- Notation: G₀ even, with exponential moments; M = ∫G₀(w) e^{w/2} dw.
- Statement: for every v,

  M(e^{v/2} + e^{−v/2}) − ∫₀^∞ e^{u/2} [G₀(u+v) + G₀(v−u)] du = ∫ G₀(w) e^{−|v−w|/2} dw.

  The left side is the pole term minus the PNT mean of the k = 1 prime measure.
- **Grok's replay (S102):**
  - By hand: split ∫G₀ e^{±(v−w)/2} at w = v; the half-line integrals cancel exactly the w > v and w < v pieces.
  - Symbolically with G₀ = e^{−w²}: residual 0 to 165 digits (`s102/twist_and_ward_replay.py` (3)).
- **Consequence (PROVED).**
  - The mean field has covariance κ(u, v) = e^{−|u−v|/2}, with positive Fourier density 1/(t² + ¼).
  - Its half-line OS reflection square is ∫∫_{u,v≥0} f(u) f(v) κ(u+v) = |∫ f e^{−u/2}|² ≥ 0: **rank one, and the vector is the pole state e^{−u/2}.**
  - This is the positive Hilbert space required by R4. It contains the soft mode by construction.

**G.3b. The fluctuation (definition).**
- E₁ := μ₁ − e^{u/2} du, with μ₁ = Σ_p (log p) p^{−1/2} δ_{log p}. This is the θ-discrepancy θ(x) − x, half-density weighted (corrected S103c; D1 from §S). Superseded in G.3b by the ψ-fluctuation E_Λ (G.3d).
- The k ≥ 2 prime powers are added the same way.
- Astra's k = 2 split (Oct 27 21:47Z) is an algebraic decomposition, PROVED: the odd bulk is zero; the even bulk is −¼·11ᵀ, rank one, still to be compensated.
- Then Q_L = MF_L + Fl_L exactly, and MF_L ≥ 0.

**G.3c. Fakes at G.3 (PROVED unless marked).** The mean field is **fake-blind**:
- **P_X, DMV:** same pole, same PNT mean, so the same MF.
- **F_θ:** kernel κ·2cosh(θr) = e^{−(½−θ)|r|} + e^{−(½+θ)|r|}. For θ < ½ this is a sum of two positive Cauchy kernels, so MF(F_θ) ≥ 0 as well. F_θ is **not** broken here.
- **H₋ₜ:** the pole is kept, so the same MF.
- **Epstein:** pole and PNT for the form's representation numbers, so an MF ≥ 0 of the same shape.
- **DH:** entire, so there is no pole and no MF. DH's whole form is fluctuation.
- **G_θ′:** cos-twisted MF, a sum of norms by the twist identity, so it passes.

So G.3 supplies R4 but carries no RH content, consistent with the S98 cone dead end.

**G.3d. The k = 2 rank-one is a ψ/θ splitting artefact (PROVED, elementary; Grok S103).**
- Statement. Put the k-th prime-power layer μ_k = Σ_p (log p) p^{−k/2} δ_{k log p}. Its PNT mean is ½·du for k = 2 (Mertens: Σ_{p≤y} (log p)/p = log y + O(1), so at u = 2 log p the mean density is ½ du), and a bounded, non-growing density for k ≥ 3. For an even g these give the rank-one term −½∫g = −½ ĝ(0) = −½|f̂(0)|² (vector 1), the shape of Astra's even bulk −¼·11ᵀ (normalisation match to Astra's ¼ not checked).
- Proof of cancellation (repaired S103b after Golden Eye's audit). Use ψ(x) = Σ_{k≥1} θ(x^{1/k}) directly: the layers μ_k sum to μ_Λ = Σ_n Λ(n) n^{−1/2} δ_{log n}, the half-density of ψ, whose PNT mean is the single density e^{u/2} du. So once the fluctuation is taken against ψ rather than θ, no separate k-th layer mean (in particular the k = 2 rank-one) remains. The earlier wording "likewise −x^{1/k} cancels the k-th mean" via θ = Σ μ(k)ψ(x^{1/k}) was WRONG term by term: the k-th layer mean is (1/k)e^{(1/k−½)u}du and θ's μ(k)x^{1/k} term gives μ(k)(1/k)e^{(1/k−½)u}du, which cancel only when μ(k) = −1 (fails at k = 4, μ = 0, and at k = 6, μ = +1). Credit: Golden Eye audit. Statement unchanged; PROVED, no novelty (classical ψ/θ bookkeeping).
- Consequence. Use the **ψ-fluctuation** E_Λ := Σ_n Λ(n) n^{−1/2} δ_{log n} − e^{u/2} du (half-density ψ(x) − x) in G.3b, not the θ-fluctuation E₁. Then Q_L = MF_L + Fl_L with no rank-one remainder, and the pole is spent exactly once, on the Ward identity G.3a.
- Test (numerics only to test this step, per the Oct 22 pivot): `s103/k2_rank_one_cancellation.py`, primes ≤ e^16, smooth bump h on (2, 16). The k = 1 fluctuation is −3.1716 = −(k = 2 layer 2.4638) − (k ≥ 3 layers 0.7078) to 4 decimals, and the Λ-total minus the PNT mean is −1.4·10⁻⁶ (the zero sum of a smooth bump). So nothing is left over for a boundary term to compensate.
- R-checks. Fake-blind (P_X, DMV and Epstein have the same internal cancellation; DH has no Λ); no RH content. Withdraws T1 below.

### G.4 Global gluing by Poisson summation on 𝔸/ℚ: where positivity is NOT implied

**G.4a (PROVED; Tate 1950 §4.2).** For φ ∈ 𝒮(𝔸): Σ_{q∈ℚ} φ(qx) = |x|⁻¹ Σ_{q∈ℚ} φ̂(q/x). This gives Tate's global FE ξ(s) = ξ(1−s).
- This is the **functional equation** half of R2.
- The comb Σ_{q∈ℚ} δ_q is a positive, positive-definite, self-dual measure on (𝔸, +).

**G.4b (PROVED; Connes 1999 §III; Meyer 2005).** Connes' map E: 𝒮(𝔸)₀ → L²(C_ℚ), (Ef)(x) = |x|^{1/2} Σ_{q∈ℚ^*} f(qx).
- Here 𝒮(𝔸)₀ = {φ(0) = φ̂(0) = 0}: the two Poisson boundary terms, i.e. the two poles, are **removed**.
- E is C_ℚ-equivariant and intertwines the global reflection: E(F_𝔸 φ) = E(φ) ∘ inv.
- The Euler product (⊗′) and the FE (Poisson) therefore meet in **one** operator. This is the only place R2 is met jointly.

**G.4c (PROVED; Connes 1999 Thm 4, §VII; Meyer 2005): the spectral side.**
- The cokernel of E (in Connes' cutoff Hilbert space; nuclear-Fréchet in Meyer) carries the zeros as spectrum.
- Weil positivity ⇔ Connes' global trace formula ⇔ RH.
- The semi-local trace formula with cutoff Λ is

  Tr(R_Λ U(h)) = 2h(1) log Λ + Σ_{v∈S} ∫′_{ℚ_v^*} h(u⁻¹)/|1−u| d^*u + o(1),

  where ∫′ is Weil's principal value. The local principal-value distributions are the G.1c terms. They are **not of positive type**.

**G.4d (GAP, the shared one).**
- **The defect.** The gluing introduces the cutoff/principal-value term:
  - Sign (settled S103c, Queen §Q.7 + Principal §P.5–P.6, per Connes 1999 §VIII Thm 5/Lemma 3 and Meyer 2005 §5.4): the positive object is a trace over the **complement of range(E)** (cokernel deficit, S_Λ − Q′, which is free), not ‖P^⊥Eφ‖². Range positivity enters Δ_Λ with a minus sign and gives only the trivial upper bound, so "additive positivity transfers through E" is WRONG-SIGN and withdrawn. The 2 log Λ term cancels exactly and is not part of δ_Λ.
  - GAP in this language: lim_Λ Δ_Λ = Weil, i.e. δ_Λ → 0 (RH-equivalent); off-line residual Σ_ρ m_ρ[∫ĝ dμ_ρ − ĝ(γ_ρ)]. Pending Golden Eye's three-test sign audit.
  - Equivalently, in §G.3's language, Q_L = MF_L + Fl_L, with Fl_L the E_Λ-fluctuation (ψ-fluctuation, G.3d; corrected S103b from E₁).
  - Nothing local (G.1–G.2) bounds δ_Λ or Fl_L.
- **The GAP statement.** Show Fl_L(f) ≥ −MF_L(f) for all f ∈ C_c^∞([−L, L]), **uniformly in L** (R1).
- **Equivalent to RH** (Weil 1952; Bombieri 2000). Window L is the semi-local problem with S_L = {∞} ∪ {p < e^{2L}}, by `primeSide_eq_truncation`.
- **Why this is sharp.**
  - (i) The classical PNT error ψ(x) − x = O(x e^{−c√log x}) controls Fl_L only for L ≪ log log of the scale. It is exponentially too weak at half density (Astra's own scope note).
  - (ii) Astra's global support bound (Oct 27, PROVED by Astra, not yet replayed by Grok): any bounded operator supported on primes p ≤ P^α, α < 1, cannot couple off-line mirror sectors as P → ∞. So the compensating operator must engage the **whole prime tail**. This is R1 restated on the operator side.

**G.4e (where each fake breaks in the chain G.1 → G.4).** The decisive step for every fake is G.4 (gluing), never G.1–G.3:

| fake | G.1 local Θ_p | G.2 Θ_∞ | G.3 MF ≥ 0 | G.4 gluing (E intertwines F_𝔸 and inv) | verdict |
|---|---|---|---|---|---|
| **G_θ′** = Λ(s+iθ′)Λ(s−iθ′) | unitary (|·|^{±iθ′} unitary) | unitary | ≥ 0 | yes: E commutes with the unitary character |x|^{iθ′}. Q_G = 2Q_ζ(f cos) + 2Q_ζ(f sin) (**PROVED**, replayed S102) | **passes every step** (must-pass control) |
| **F_θ** = ξ(s−θ)ξ(s+θ) | **product** γ unitary on Re s = ½ (replay (5): \|γ₊γ₋\| = 1), but each GL(1) constituent \|·\|^{±θ} has \|γ\| ≠ 1 (0.817, 1.224 at p = 7) | same | ≥ 0 (G.3c) | E does commute with \|x\|^θ, but the twisted form is a **polarization** 2 Re B_ζ(e^{θu}f, e^{−θu}f), not a norm (**PROVED**, replay (2)) | breaks at G.4 only (norm → polarization, §P.4). Not a GL(2) filter: for 0<θ<½ π(\|·\|^θ,\|·\|^{−θ}) is complementary series with positive-definite Knapp–Stein form (§P.3, Knapp–Stein 1971); local non-temperedness is a size shadow, not a sign condition |
| **DH** | no local data (no Euler product) | ζ-like γ_∞ | **no MF** (entire) | no ⊗′ factorization, so E of a sum of two L-functions is not a norm of a sum | breaks at G.1/G.4 |
| **DMV / P_X** | Beurling or truncated local data | — | ≥ 0 | **no lattice ℚ^* ⊂ 𝔸**, so no Poisson and no global reflection; P_X agrees with ζ on L ≤ (log X)/2 (R1) | breaks at G.4a |
| **Epstein** (h = 2) | sum of two Euler products (genus characters) | γ_∞ for GL(2)/ℚ(√D) | ≥ 0 | Poisson over the lattice exists, but there is no ⊗′ of one Euler product; E of a sum of two pure tensors | breaks at G.1d/G.4b |
| **H₋ₜ** (de Bruijn–Newman deformation) | no Euler product after deformation | **γ_∞ is not Tate's** (deformed Φ kernel) | ≥ 0 | no adelic E | breaks at G.2/G.4 |

**Reading.** All six fakes break at the gluing G.4 (and DMV/P_X *only* there). (S103c: F_θ's earlier "local temperedness filter" break is withdrawn per §P.3; it breaks only at G.4.) G_θ′ passes everywhere. So R5 is consistent with G.4d being the true locus of RH. That is exactly Connes' "uniformity" dead end, now with every fake located. **No new route is claimed.**

### G.5 What would count as progress (CONJECTURE-level targets, each R1–R5-checked)
- **(T1) Pole/boundary compensation (WITHDRAWN S103 by G.3d).** The k = 2 rank-one needs no compensation (it is absorbed once the fluctuation is taken against ψ = Σ_k θ(x^{1/k}), G.3d), and the two boundary directions are already spent on the mean field (G.3a). Original statement kept for the record:
  - The even k = 2 rank-one −¼·11ᵀ (Astra) and the k = 1 fluctuation are compensated *exactly* by the two Poisson boundary directions (φ(0), φ̂(0)) that 𝒮(𝔸)₀ removes.
  - Test against R5: DH has no boundary directions and no MF, so the conjecture must fail for DH. Good: it is not fake-blind.
- **(T2) Tail engagement (CONJECTURE).** Any compensating operator for Fl_L is non-quasi-local: Astra's support bound gives the necessity half (Astra PROVED, Grok replay pending). The sufficiency half is the GAP.
- **(T3) Polarization criterion (CONJECTURE, Principal series §2).** An L²(C_ℚ)-norm mechanism passes G_θ′ automatically and is turned into a polarization by F_θ. Any ARP proof must therefore be an L²(C_ℚ) norm identity at G.4.

### G.6 Lean asks (finite, proved statements; for Claude Code / Codex)
1. `twist_cos_sum`: 2cos(θr)(f ∗ f̃)(r) = (f cos θ·) ∗ (f cos θ·)~ + (f sin θ·) ∗ (f sin θ·)~, and the zero-side parallelogram |A+iB|² + |A−iB|² = 2|A|² + 2|B|² (already asked S100; still open).
2. `twist_cosh_krein`: 2cosh(θr)(f ∗ f̃) = h₊ ∗ h̃₋ + h₋ ∗ h̃₊, h± = e^{±θu} f (polarization; already asked).
3. `local_gamma_unitary`: |1 − p^{−1/2+it}| = |1 − p^{−1/2−it}|, so |γ_p(½+it)| = 1. Also, for real θ ≠ 0, the constituent factor γ_p^{(θ)}(s) := L_p(1−s−θ)/L_p(s+θ) has |γ_p^{(θ)}(½+it)| ≠ 1 in general, while |γ_p^{(θ)} γ_p^{(−θ)}| = 1 on Re s = ½.
4. `local_prime_density_sign`: Re[z/(1−z)] = (Re z − |z|²)/|1−z|², with a witness that it takes both signs on |z| = p^{−1/2}.
5. `pnt_ward_split`: the G.3a split-integral identity for even G₀ with exponential moments, and the half-line rank-one square ∫∫_{u,v≥0} f(u)f(v)e^{−(u+v)/2} = (∫f e^{−u/2})².
6. `mirror_overlap`: 1 − R_p² = (a−b)²/(1−ab)², with a = p^{−σ}, b = p^{−(1−σ)}, and ∏_p R_p(σ) = 0 for σ ≠ ½ (Astra Oct 26; Kakutani-type dichotomy, needs Σ 1/p = ∞, Mathlib has the divergence of Σ_p 1/p, e.g. `Nat.Primes.not_summable_one_div`; check the exact name).

### G.7 References
- Tate, J., Thesis (1950), in Cassels–Fröhlich (1967), Ch. XV.
- Weil, A., "Sur les « formules explicites »…", Comm. Sém. Math. Lund (1952).
- Barner, K., Math. Ann. 1981 (explicit formula regularity).
- Yoshida, H., Adv. Stud. Pure Math. 21 (1992).
- Connes, A., Selecta Math. 5 (1999) 29–106.
- Meyer, R., Duke Math. J. 127 (2005).
- Connes, A. & Consani, C., Selecta Math. 27:77 (2021), arXiv:2006.13771.
- Burnol, J.-F., C. R. Acad. Sci. 331 (2000); Forum Math. 16 (2004).
- Bombieri, E., Rend. Lincei 11 (2000).
- Shahidi, F., Ann. Math. 132 (1990).
- Kakutani, S., Ann. Math. 49 (1948), on equivalence of infinite product measures (prior art for Astra's superselection).
- Astra, PR #11 posts of Oct 26 08:59Z, 10:02Z and Oct 27 21:27Z, 21:47Z, 22:02Z, 22:05Z; bridge `work/2026-10-2{6,7}_astra_*.md`.
