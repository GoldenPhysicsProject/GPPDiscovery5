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
- **R1. Uniform in L.** No bounded window implies RH: by Goddess of PNT's DMV theorem (`team/prime-number-theory/finite_windows_dmv.md`), DMV/P_X systems agree with ζ on every bounded window. Sharp ladder (§N N.2): P_X agrees with ζ exactly for 2L ≤ log p⁺(X), p⁺(X) the least prime > X.
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
- **Bost–Connes / KMS / Astra mirror-OS route (S103c, §C C.2d, C.8a, C.9).** KMS/primon OS ⇔ prime weights ≥ 0 (range-side, wrong sign under the deficit); no covariant KMS identity equals Q_L; every reading of Astra's H1 (zero ↔ coherent mode) is either ⇔ RH/G.4d or an RH-silent unconditional theorem; DMV is a countermodel for any H1 weaker than H1-unitary. Adds nothing beyond G.4d.

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
- **P_X, DMV:** same pole, so the same MF, but "same PNT mean" holds **only by convention** (§N N.0/N.1: keep ζ's pole and W_∞); P_X's true PNT mean is O(log x).
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
  - Sign (settled S103c, Queen §Q.7 + Principal §P.5–P.6, per Connes 1999 §VIII Thm 5/Lemma 3 and Meyer 2005 §5.4): the positive object is a trace over the **complement of range(E)** (cokernel deficit, S_Λ − Q′; exactly free only for function fields. For ℚ, Q′ ≤ S_Λ is not exact with the prolate cutoff, and the exact positive object is Δ^♮_Λ = Tr((S − SQ′S)VV*) (§Q Q.8 correction c2)), not ‖P^⊥Eφ‖². Range positivity enters Δ_Λ with a minus sign and gives only the trivial upper bound, so "additive positivity transfers through E" is WRONG-SIGN and withdrawn. The 2 log Λ term cancels exactly against Tr(S_Λ V(h)) and is not part of δ_Λ; source: Connes 1999 §VIII, eq. (16) and the proof of Thm 5, eqs. (17)–(26) (Queen sign_check.md). Connes–Consani 2021 Thm 1 confirms the deficit sign only: it is the Λ = 1 (Sonin) case, where 2 log Λ = 0, so it cannot test this cancellation (Golden Eye audit).
  - Pole side (§P P8.1, PROVED from Connes App. I (37) = Meyer (31); ℚ as Connes states it): as functionals f(0) ∈ P_Λ (s = 0) and f̂(0) ∈ P̂_Λ (s = 1), swapped by F; Connes' D term uses the same pair B_Λ ⊃ B_{Λ,0} as the deficit trace, and the removed directions are exactly G.3a's OS4 span{e^{±u/2}}, spent once. The s = 0 side is local (survives with S-units, §N N.10); the s = 1 side is global (primes outside S). Projector dictionary S_Λ ⊖ E(B_{Λ,0}) vs CC's Sonin 1 − P₁∨P̂₁ (formerly a minor GAP): **closed negatively at the projector level** (§P P.9.1, REMARK-level beyond (D1)–(D3)). The two act on different spaces (local L²(ℝ)_ev vs L²(C_ℚ) after E), and no choice of Λ identifies them. **The trace-level dictionary is PROVED** (§P P.9.2–P.9.3, Grok audit S103e: SOUND).
  - R6 identity recorded (§P (P9.3), exact for every Λ ≥ 1): L_Λ(f) = W_∞(f) + Tr(ϑ(f)R_Λ) − Tr(ϑ(f)Ŝ_Λ) + χ_Λ(f), with Tr(ϑ(f)Ŝ_Λ) = 2f(1) log Λ exactly.
    - Equivalently L_Λ = L + ε_Λ + χ_Λ, with ε_Λ = −Tr(ϑ(f)Ŝ_ΛP^c_Λ) and χ_Λ = Tr([P^c_Λ, ϑ(f)]P_ΛŜ_Λ), both trace class and uniformly bounded in Λ.
    - Λ = 1 recovers CC 2021. (P9.4), CC §4 at bandwidth 2πΛ², is PROVED given (P9.3) and CC §4 at a = b = Λ; Golden Eye checked (P9.1)–(P9.3) line by line, and CC's (69)–(84) were checked by inspection (§P P.10.0).
  - **Λ-track CLOSED NEGATIVELY (§P P.10, merged S103f, Grok audit: SOUND; RH-silent).** This replaces the former "Λ > 1 analogue of CC Lemma 6.10 / Thm 1 is the [S&M] GAP".
    - Lemma 6.10_Λ, ⟨ξ|N^{(Λ)}_Iξ⟩ ≤ ⟨ξ|N_Iξ⟩ ≤ γ|⟨ξ₀|ξ⟩|², and Thm 6.11_Λ with c_Λ ≤ 4γ/log 2 hold uniformly in Λ. Both are PROVED given CC Lemma 6.10, which is CERTIFIED (computer-assisted).
    - The source is a monotonicity lemma: Tr(ϑ(h)S^{So}_Λϑ(h)^*) ≥ 0 is nonincreasing in Λ and tends to 0. So **Sonin Λ = 1 is optimal**, and the Λ-track adds nothing beyond CC.
    - Exact limits (P10.4) are PROVED. The naive ε^{(Λ)}∘Q bound fails (PROVED given (P9.4)). Fixed-Λ compactness is a GAP, not needed.
    - **The remaining route past CC is the support extension (G.4g), which is where primes enter (R2).**
    - Consequence for Q.8.E (PROVED as stated): a rule L ↦ D_L using only archimedean data cannot satisfy (I1) for both ζ and a datum such as P_X (X < 2) once 2L > log 2. So D_L must contain E through Q′_Λ.
  - Definition of μ_ρ: for a zero ρ = ½ + β + iγ, μ_ρ is the harmonic (Poisson-kernel) measure of ρ with respect to the critical line, a positive probability measure on ℝ, with μ_ρ = δ_γ when β = 0. The identity lim Δ_Λ(f) = Σ_ρ m_ρ ∫ĝ dμ_ρ is PROVED for function fields (Connes 1999 §VIII Lemma 3); for ℚ it is Connes' stated analogue (§VIII, last part), labelled CONJECTURE here until a proof for ℚ is cited.
  - GAP in this language: lim_Λ Δ^♮_Λ = Weil, i.e. δ_Λ → 0 (RH-equivalent). Over ℚ the unsharpened Δ_Λ is not of positive type (§Q c2), so the form "lim Δ_Λ = Weil ⇒ RH" holds only given M1 (Golden Eye §Q.8 audit, fix 3); off-line residual Σ_ρ m_ρ[∫ĝ dμ_ρ − ĝ(γ_ρ)]. Sign audit (Golden Eye, S103c): SOUND on all three tests; shared GAP sign SETTLED, with the ℚ-case of the μ_ρ residual formula still CONJECTURE.
  - **M1 (Queen §Q M1, merged S103e, Grok audit: SOUND).** M1 says Lk_Λ = Δ_Λ − Δ^♮_Λ → 0. Here ℓ_Λ := ‖(1−S_Λ)Q′_Λ‖²_HS.
    - **Thm M1.A (PROVED):** |Lk_Λ(f∗f^*)| ≤ ‖f‖_1²ℓ_Λ + 2‖f‖_1‖f‖_2(2 log Λ)^{1/2}ℓ_Λ^{1/2}.
    - **Cor M1.B (PROVED-conditional):** ℓ_Λ log Λ → 0 ⇒ M1. On Q.8.E's diagonal Λ = e^L, uniformity needs ℓ_Λ = o(1/log²Λ).
    - **Thm M1.C / Cor M1.D (PROVED-conditional)** on Spectrum's prolate inputs (S2), (S3): ℓ_Λ ≤ 2κ_Λ(ε_Λ/Λ + 4ζ(1+σ)²ρ_Λ).
    - **GAP H-κ (Queen):** conditioning κ_Λ of E on B_{Λ,0}. It is a minimum over B_{Λ,0}'s Mellin profiles of the |Mφ|²-average of |ζ(½+is)|².
      - It is arithmetic but uses the critical line only. Fixed-weight mean values do not bound it directly, and it probably needs local zero-spacing input.
      - **H-κ status (§Q ARP_queen_Hkappa.md, merged S103f, Grok audit: SOUND after a low-height fix).**
        - Polynomial bound **Thm H.1, κ_Λ ≤ Λ^B, B absolute: PROVED (sketch)**. It uses the two-constants lemma plus Kovrijkine 2001. It is fake-blind (DH, Epstein), hence RH-silent. The local input is thickness of a superlevel set of |ζ|, not zero spacing.
        - The sharp form Λ^{o(1)} reduces to **(Z-avg), a GAP**. (Z-avg⁺) holds under Lindelöf, hence under RH. (Z-avg⁻) is a CONJECTURE.
        - **Circularity flag:** the sharp form may be RH-conditional, so it cannot be used on any route to RH.
        - **Open: is B < 2?** That is what the untrimmed Connes space needs (ε_Λ/Λ ≍ log Λ/Λ², HEURISTIC). Otherwise M1 needs a trimmed prolate space (CONJECTURE, Spectrum).
        - **Part II update (ARP_queen_Hkappa.md H.8–H.11, merged S103g; Grok audit: SOUND after two wording/arithmetic fixes).**
          - The window exponent improves to O(log(e/a)), down from O(1/a). This is Lemma H.2′, PROVED (citation: Titchmarsh 9.6(B) and 9.2; the page and the "<1" strictness are not re-verified). It gives **Thm H.1′, κ_Λ ≤ Λ^{B′}, PROVED (sketch)**.
          - **Λ^{o(1)} is unreachable by this method.** The O(log t) error and the zero count each cost T^{−c}. The sharp form still needs (Z-avg⁻), which is a GAP.
          - **Trimmed space proposed as primary (OPINION, Queen).** On I^{(K)} = {1 − λ_n ≤ Λ^{−K}} with K ≥ 2B + 8, Thm H.9 gives ℓ^{(K)} ≤ CΛ^{−1}, so **M1 holds on the trimmed space**, diagonal included. Label: PROVED-conditional on the H.1/H.1′ sketch plus the GAP-minor extension of Lemma M1.1. It uses no Spectrum input. The endpoint bound is Lemma H.9.1, which is elementary.
          - **Identifying it with Connes' Q′_Λ is conditional on (T3)**, the trimming-harmless check Tr(SΔQSV(h)) = o(1), which is also M2-compatibility. (T3) splits into (T3a), the plunge-state low-Mellin mass (HEURISTIC, Spectrum), and (T3b), the E-transport (GAP, Queen).
          - **(T3a) at the conjectured rate is expected FALSE** (S&M T3a_plunge_mellin 65bba74, merged S103g; CONJECTURE, HEURISTIC + NUMERICAL). Plunge states oscillate at frequency ≈ Λ at every x, so ω^{arch}_Λ(log Λ) ≍ (log Λ)/Λ², not Λ^{−B−2−δ}. The model |M cos(c·)|² = π/(2c) is PROVED. **So, along Queen's reduction, (T3) needs B < 2** (κ_Λ ≤ CΛ^B), or a (T3b) that avoids the m_Λ loss. It also needs the upper bound ω^{arch} ≤ CS/Λ² (GAP) and (T3b) itself (GAP). B < 2 is not available: H.1/H.1′'s B′ is absolute but not explicit, and the sharp Λ^{o(1)} is barred.
          - **Part III update (Queen 0fb6f9c, merged S103g; Grok audit: SOUND with two label corrections).**
            - H.12 closes the Lemma M1.1 density extension (PROVED). The proof is pointwise symmetric Poisson for truncated prolates off a null set (Dirichlet–Jordan; the Zygmund theorem number is unchecked), the pole terms vanish, and a Mellin approximate identity lands in S_0^ev.
            - Hence **M1 on the trimmed space is PROVED (sketch), given H.1/H.1′ only**.
            - (T3a) is withdrawn and replaced by the trace form (T3a*): 𝒟_rem(s) ≤ Cr_K(1 + |s|)^AΛ^{−2} on |s| ≤ Λ^ε. This is a CONJECTURE.
            - **(T3) reduces to (T3a*) + (PD) + (L), with no B < 2 condition.** Thm H.13 is PROVED-conditional at rate O(K(log Λ)^{A+2}Λ^{−2}). (PD), profile domination of the E-orthonormal complement, is a GAP and must also hold in integrated form for the jump sum. (L), the remainder-leak low-pass, is HEURISTIC.
            - The worst-case route still needs B < 2 plus kept-feedback transport.
            - **(T3) is NOT proved.**
          - **Part IV update (Queen b42c7ef, merged S103h(4) fbabd66; Grok audit: SOUND, three label/citation notes).**
            - **(T3) ⇐ (A_β) + (LD_γ) (Queen, GAP) + (MB) + (F) + (T3a*) + (L) (archimedean, Spectrum; RH-silent), with β + γ < 2.** This is Cor H.15, PROVED-conditional, via Thm H.15 ((PD) pointwise, loss Λ^{β+o(1)}) and H.15.5 ((PD) integrated), using Lemma H.15.3 (PROVED (sketch)). Exact Schur formulas and the spectral-splitting Lemma H.15.2 are PROVED. cond(G_m) ≤ C log Λ is PROVED. The BCHB range check (only |s| ≥ Λ^{4/3+ε}) is PROVED.
            - DMS-type inversion with cond ≍ Λ^B is fatal as a method (PROVED for the bound; "no theorem beats it" is a REMARK). The global Neumann series has no uniform rate (PROVED (sketch); "diverges" is not shown). The e^{−|Δu|/2} decay of G_osc is HEURISTIC. Citations were checked; Gröchenig–Klotz I is quasi-polynomial in cond, which does not change the conclusion.
            - Circularity flag: (A_β) and (LD_γ) must not come from RH-conditional arguments. **(T3) remains OPEN.**
          - H.10's RH-conditional line is flagged and barred.
          - The trimming count is corrected to O(K log²Λ).
    - **M1 overall: CONJECTURE.** RH-silent: no implication to or from RH.
  - Equivalently, in §G.3's language, Q_L = MF_L + Fl_L, with Fl_L the E_Λ-fluctuation (ψ-fluctuation, G.3d; corrected S103b from E₁).
  - Nothing local (G.1–G.2) bounds δ_Λ or Fl_L.
- **The GAP statement.** Show Fl_L(f) ≥ −MF_L(f) for all f ∈ C_c^∞([−L, L]), **uniformly in L** (R1).
- **Equivalent to RH** (Weil 1952; Bombieri 2000). Window L is the semi-local problem with S_L = {∞} ∪ {p < e^{2L}}, by `primeSide_eq_truncation`.
- **Why this is sharp.**
  - (i) The classical PNT error ψ(x) − x = O(x e^{−c√log x}) controls Fl_L only for L ≪ log log of the scale. It is exponentially too weak at half density (Astra's own scope note).
  - (ii) Astra's global support bound (Oct 27, PROVED by Astra; replayed in §C C.6d and transferred to the BC algebra): any bounded operator supported on primes p ≤ P^α, α < 1, cannot couple off-line mirror sectors as P → ∞. So the compensating operator must engage the **whole prime tail**. This is R1 restated on the operator side.

**G.4e (where each fake breaks in the chain G.1 → G.4).** The decisive step for every fake is G.4 (gluing), never G.1–G.3:

| fake | G.1 local Θ_p | G.2 Θ_∞ | G.3 MF ≥ 0 | G.4 gluing (E intertwines F_𝔸 and inv) | verdict |
|---|---|---|---|---|---|
| **G_θ′** = Λ(s+iθ′)Λ(s−iθ′) | unitary (|·|^{±iθ′} unitary) | unitary | ≥ 0 | yes: E commutes with the unitary character |x|^{iθ′}. Q_G = 2Q_ζ(f cos) + 2Q_ζ(f sin) (**PROVED**, replayed S102) | **passes every step** (must-pass control) |
| **F_θ** = ξ(s−θ)ξ(s+θ) | **product** γ unitary on Re s = ½ (replay (5): \|γ₊γ₋\| = 1), but each GL(1) constituent \|·\|^{±θ} has \|γ\| ≠ 1 (0.817, 1.224 at p = 7) | same | ≥ 0 (G.3c) | E does commute with \|x\|^θ, but the twisted form is a **polarization** 2 Re B_ζ(e^{θu}f, e^{−θu}f), not a norm (**PROVED**, replay (2)) | breaks at G.4 only (norm → polarization, §P.4). Not a GL(2) filter: for 0<θ<½ π(\|·\|^θ,\|·\|^{−θ}) is complementary series with positive-definite Knapp–Stein form (§P.3, Knapp–Stein 1971); local non-temperedness is a size shadow, not a sign condition |
| **DH** | no local data (no Euler product) | ζ-like γ_∞ | **no MF** (entire) | no ⊗′ factorization, so E of a sum of two L-functions is not a norm of a sum | breaks at G.1/G.4 |
| **P_X** | truncated local data | — | ≥ 0 (by convention) | no lattice ℚ^* ⊂ 𝔸; agrees with ζ for 2L ≤ log p⁺(X) | **PROVED to fail** past the ladder: ind_− ≥ 1 for L ≥ ½ log X + O(1), odd-pole bump witness (§N, headline) |
| **DMV (glued, full prime tail)** | Beurling local data | — | ≥ 0 | no lattice ℚ^* ⊂ 𝔸, so no Poisson and no global reflection | breaks at G.4a structurally; **PROVED to fail for large L** (§N N.11.3, ineffective; explicit onset GAP N.G1′): ind_− ≥ 1 for some L > ½ log p⁺(X) (onset past the ladder given RH), r → −∞ along a sequence |
| **Epstein** (h = 2) | sum of two Euler products (genus characters) | γ_∞ for GL(2)/ℚ(√D) | ≥ 0 | Poisson over the lattice exists, but there is no ⊗′ of one Euler product; E of a sum of two pure tensors | breaks at G.1d/G.4b |
| **H₋ₜ** (de Bruijn–Newman deformation) | no Euler product after deformation | **γ_∞ is not Tate's** (deformed Φ kernel) | ≥ 0 | no adelic E | breaks at G.2/G.4 |

**Reading.** All six fakes break at the gluing G.4 (and DMV/P_X *only* there). (S103c: F_θ's earlier "local temperedness filter" break is withdrawn per §P.3; it breaks only at G.4.) G_θ′ passes everywhere. So R5 is consistent with G.4d being the true locus of RH. That is exactly Connes' "uniformity" dead end, now with every fake located. **No new route is claimed.**

**G.4f (§N N.11, merged S103d; RH-silent).** Every borrowed-W_∞ system that meets the hypotheses of §N N.11.3 (prime measure with
∫e^{−σu}d|μ| < ∞ for σ > ½, D_μ meromorphic in Re z > a for some a < ½, a pole z₀ of D_μ with Re z₀ > max(a, 0); for glued DMV with θ > ½ this is a real restriction, met by its zeros for large t (Golden Eye N.11 audit)) has ind_− ≥ 1 for large L (PROVED,
ineffective; two-bump identity N.11.2 plus the classical Phragmén–Landau/Ingham Laplace argument). This covers glued DMV (off-line zeros)
and, as the ε = −1 case, plain P_X (uncancelled pole, = N.3.3). It is the easy direction of Weil's criterion, uses no FE, and gives
nothing toward RH. DH and Epstein are outside the hypotheses. F_θ is covered by a sketch only. G_θ′ is not falsely failed (its poles
lie on Re z = 0 under RH). §N N.11.4: the FE adds nothing to the failure direction. The positivity direction needs it, because
"no poles in Re z > 0" does not imply positivity (P_X counterexample, PROVED).
- **Open items from §N after N.11:**
  - GAP N.G1′: explicit onset L₁(DMV), via effective oscillation (Turán/Pintz-type) and explicit DMV zero heights.
  - GAP N.G1″: polynomial bounds for ζ_P′/ζ_P on a line σ = θ₁ ∈ (θ, 1) (Landau lemma; needed only for the zero-side/quartet route).
  - CONJECTURE N.C3: ind_−(Q^{P_DMV}_L) → ∞.
  - (N.G1 itself is CLOSED.)

**G.4g (p = 2 step past CC, Celestial §C.13 + Prime §N N.12, merged S103e; RH-silent on this range).**
- **(IF), the reduced GAP for the p = 2 step.** On the pole-free class, for ½ log 2 < L ≤ ½ log 3, Weil positivity at Λ = 1 is equivalent to **(IF)**: 2c₂B(ξ) ≤ R_L(ξ). Here:
  - c₂ = (log 2)/√2;
  - B(ξ) = ∫[ξ′ξ′(·+log 2) + ¼ξξ(·+log 2)] = Qh(log 2);
  - R_L = 𝒮_L + 2ε′(1⁺)(‖ξ‖² − ⟨ξ, K_Lξ⟩) = W_∞(Qh) is the archimedean form.
  - The reformulation is PROVED (§C C.13.A). Only n = 2 enters, with tail vanishing PROVED-LEAN `GppWeilLadder.primeSide_eq_truncation` / `primeSide_eq_zero_of_support_lt_log_two`.
  - Status: **GAP**, RH-silent on this range (ζ ≡ P_X for X ≥ 2 here).
  - Positivity up to L = 1.42 is CERTIFIED (Arb replay with rigorous remainder, main 2576c89), conditional on Mori's Prop 3.1, Lemma 3.2 and Thm N (unrefereed). So positivity on (½ log 2, ½ log 3] is already CERTIFIED (conditional). The analytic p = 2 step is for mechanism (which step carries positivity), not to establish the result; it is not an analytic L₀. (Label corrected per Golden Eye; was NUMERICAL.)
- **CC's verbatim chain fails at T8 (PROVED).** The p = 2 term is unbounded on L² after Q, of order c₂m². So CC's ε-margin cannot absorb it.
  - This obstructs only the margin route, not (IF). 𝒮_L is also of order m², and (IF) is dominated at high frequency where σ_S ≥ 2c₂.
- **Prime-side reduction (§N N.12, PROVED).**
  - S1: |2c₂B| ≤ c₂‖ξ‖²_{H(E)} on the disjoint edge strips. This is sharp.
  - Exact parity/J splitting: (IF) ⇔ (K2) ∧ (K3), the odd and even harmful-cone inequalities with constant 2c₂.
  - At edge order ℓ = 2L − log 2, the p = 2 term **helps** the even nodeless sector and **hurts** the odd sector.
  - The comparison step S2, R_L ≥ c₂‖ξ‖²_{H(E)}, implies archimedean-only positivity on the pole-free class, and overpays.
  - CC 2106.01715 §2.2 (NUMERICAL, full test space) indicates the bound-first route caps at 2L₀ < log 2.44. That is indicative only, per the Grok audit: the pole-free restriction is not covered by CC's matrices.
- **Named obstruction N.G2 (GAP):** no analytic lower bound for R_L on the odd harmful cone (K2) or the even harmful cone (K3) for 2L > log 2. The required constant is 2c₂ = 0.98026 against the H(S₋)-mass of the harmful J-component.
- **Q.8.E first test now reduces to (IF).** Under (I1), D_L must carry the p = 2 line with diagonal sums −c₂(δ_a + δ_{−a}), e.g. **−(c₂/ℓ)S_Λ(V(2) + V(2)^*)S_Λ**. *(Fixed S103g per Principal P.14.3, Grok-checked: SV(2)S has diagonal sum ℓδ_a, since |S ∩ (S − a)| = 2L − a = ℓ, so the earlier "−c₂S_Λ(V(2) + V(2)^*)S_Λ" was off by the factor ℓ. In the form picture on 𝒯_L the uncompressed −c₂(V(2) + V(2)^*) has weight 1. No result used the block.)*
  - For a given D_L, (I2) on this window ⇒ (IF) (PROVED). Existentially, (IF) ⇔ some D_L ≥ 0 (abstract GNS).
  - An explicit D_L must dominate the harmful cones of N.12.3 with constant 2c₂. Existence of D_L remains CONJECTURE (Q.8.E).
- **Thm N.13.A (§N N.13) is conditional on GAP N.G2-m₀** (the H-margin m₀ at ½ log 2). Mechanism: the p = 2 line carries the even edge traces, and the archimedean margin carries the odd sector. See G.4h for the full split.
- Labels: R1 fixed window | R2 FE only through W_∞ | R3 the sign of B, not weight signs | R4 soft mode helped, odd partner harmed | R5 DH positive here numerically but structurally without an a-line; fake-blind in truth. G.4d is unchanged.

**G.4h (S103g synthesis: N.G2 split into sub-gaps, merged §N N.13 + §C C.13.B + §S SECTOR_COERCIVITY_p2/CELESTIAL_INPUTS_C13; RH-silent on this range).** Grok audit: all three SOUND (constants replayed: Ω(0) = −5.37218, Ω(12.8) = 0.711 ≥ c₂, h* ≤ 4.5, ℓ*-coefficient 6.98·10⁻⁴; Lemma C.13.B bracket at ℓ = 0.119 is +0.0019 with M_ℓ = 4.96; plunge mass D(1) = 0.486326 matches direct quadrature).
- **One symbol (§S, PROVED from CC (102)).** (IF) on C_c^∞(−L, L) is ∫(¼ + t²)[2θ′(t) − 2c₂cos(t log 2)]|ξ̂|² ≥ 0. The Sonin split σ_S + (−ε̂) is a choice of partition. κ̂ > 1 exactly forced on |t| < 6.28984 (where 2θ′ < 0). Stirling: 2θ′ ≥ c₂ for |t| ≥ 10.63, ≥ 2c₂ for |t| ≥ 17.11 (PROVED, sufficient thresholds; exact crossings 10.2615 and 16.7481, NUMERICAL, Golden Eye mpmath replay — ≈0.4 slack in t).
- **N.G2 now splits into four named pieces.**
  1. *Strip-localized cone:* closed for 2L ≤ log 2 + 0.119 (Celestial Lemma C.13.B; **PROVED (computer-assisted)**, S103g: Arb certifies Ω(0), Ω(14), c₂ and M_ℓ ≤ 5.0245, with bracket(0.119) = 0.0019316… > 0. Monotonicity in ℓ is analytic: the bracket is δ − Cℓ³ with M_ℓ nondecreasing. Grok replay identical; team/celestial-cft/arp_c/c13b_lemma_certify.*. Golden Eye audit: SOUND; the method ceiling is ℓ ≈ 0.119; the proof needs the Re ψ monotonicity line, see item 3).
  2. *Edge traces (ramps):* exact rank-two boundary form; its contribution to Q_L (i.e. the ramp part of −2c₂B) is **+2c₂(1/ℓ − ℓ/24)α₋α₊** (Prime N.13.4, PROVED). Even edge pair (α₋ = α₊ = α) is **helped** by +2c₂α²(1/ℓ − ℓ/24); odd pair (α₋ = −α₊) is **charged** exactly 2c₂α²(1/ℓ − ℓ/24) (≤ 2c₂α²/ℓ). *(Sign and ℓ/24 term fixed per Golden Eye G.4h audit; the earlier text quoted 2c₂B's sign.)*
  3. *Coupling 2R(ξ_b, ξ_s):* **GAP N.G2-bulk**, the one irreducible bilinear term (Celestial C.14 asked).
     - **Celestial §C.14 (merged S103g; Grok audit: SOUND).**
       - **(K2) reduces to the symbol inequality of Thm C.14.A:** σ_S(t) + 2ε′(1⁺)(1 − κ^{odd}_L)/(¼ + t²) > c₂ for all t. This is the semilocal {∞, 2} frame D_L^{(w)}, with (I1) exact and the Fejér minimum → −c₂. It is coupling-free and sufficient. It is a pointwise θ = 1 sibling of S&M Thm 2.1, the polynomial route, with no ℓ* restriction. It also gives m₀^o ≥ inf_t[same expression at L₀] (PROVED given CC). Its value is GAP.
       - **(K3) stays GAP.** It goes via (P), the even ĝ(0) direction (Spectrum's), plus **Conjecture C.14.K**, the high-band almost-orthogonality |2P_{>T}(ξ_b, ξ_s)| ≤ γ[P_{>T}(ξ_b) + P_{>T}(ξ_s)] with γ < Ω(T)/(Ω(T) + |Ω(0)|). Given C.14.K and m₀^e > 0, (K3) follows (Cor C.14.D).
       - **Explicit bounds (Thm C.14.B, PROVED).** The local and low-band parts of the coupling are polynomially small: O(ℓ⁴)‖φ‖²_H, ε₂α²/ℓ, and O(T³ℓ³)‖ξ_b‖²_H.
       - **Obstructions (Thm C.14.C, PROVED).**
         - (i) No bound by α alone exists. The off-diagonal kernel sign Qκ < 0 is Arb-certified on [0.001, 1.1], which covers every distance for 2L ≤ log 3.
         - (ii) Cauchy–Schwarz cannot close the high band.
         - (iii) The coupling lives at the inner-edge kink, at frequency ~1/ℓ.
       - **C.14.A indicative check (Celestial §C.14.2, merged S103g; Grok audit: SOUND as NUMERICAL/indicative). Superseded by the actual-K_L value below (Prime §N.14.1).** With 1 − κ^{odd} = 0.58 (Golden Eye sinc-kernel proxy, not proved; it overestimates the gap about 1.85×) and 2ε′(1⁺) = 45.993 (CC, NUMERICAL), the ε-term is 26.68/(¼ + t²). This drops below c₂ = 0.4901 at |t| ≈ 7.36. Under every σ_S lower bound we have (σ_S ≥ 0; σ_S ≥ 2θ′ − B/(¼ + t²) with B = 46 or 66), the symbol inequality fails on roughly 7.37 ≤ |t| ≤ 11.8–13.0. The worst deficit is 0.222 at t ≈ 9.97 (B = 46) and 0.267 at t ≈ 10.92 (B = 66). **NUMERICAL (indicative); σ_S is GAP.**
         - **Actual K_L (Prime §N.14.1, merged S103g; NUMERICAL).** The top odd eigenvalue is 0.68649, so δ_o⁰ = 1 − κ^{odd}_{L₀} = 0.3135 (Nyström N = 200–600, CC kernel (99)/(104), 11 prolate terms). The ε-term is **14.42/(¼ + t²)** (45.993 × 0.3135 = 14.419). It is ≥ c₂ only on |t| ≤ 5.40. The remaining σ_S band is **≈ [5.41, 13.6]**, with σ_S needed ≈ 0.35 at t = 10 and 0.40 at t = 13. The worst deficits under the known σ_S bounds are 0.345 (B = 46) and 0.369 (B = 66). Grok checked all of these numbers.
         - **Frame route to (K2)** needs a **certified σ_S(t) > c₂ − 14.42/(¼ + t²) on 5.41 ≤ |t| ≤ 13.6** (≈ 0.35–0.40) **plus a certified 1 − κ^{odd}**. *This matters less now: S&M's direct odd certificate (item 4 update) closes odd (K2) on 2L ≤ log 2 + 0.028 without it.* (Earlier proxy version: σ_S > c₂ − 26.68/(¼ + t²) on [7.37, 13].) Otherwise use the **band-split route** (S&M Thm 2.1 / C.13.H′), which has no pointwise requirement but a small ℓ-range.
         - **Golden Eye audit of c13b_lemma_certify (item 1): SOUND.** The method ceiling is ℓ ≈ 0.119, with T = 14 optimal, so this route cannot go past 2L ≈ log 2 + 0.119. The written proof needs the line **"Re ψ(x + iy) is increasing in |y|"** (for x > 0), and Celestial has been asked to add it. *(Grok: this is true, since Re ψ(x + iy) − ψ(x) = Σ_{n≥0} y²/((n + x)((n + x)² + y²)), and each term is increasing in |y|.)*
       - **Near-null directions.** They are low-band bulk prolate-type vectors, not ramps or strip vectors (low-band fractions ≤ 8.97ℓ³ and 5.45ℓ³, PROVED). "m₀ at L₀ is not ~10⁻⁸¹" is HEURISTIC.
  4. **Update S103g, S&M CERT_p2_constants (7906287/f6c7257, merged; Grok audit SOUND; Arb replay at 80/128 bits identical).**
     - **Odd (K2) is CERTIFIED on 2L ∈ (log 2, log 2 + 0.028], CC-free.** The bound is Q_L ≥ 0.0126‖ξ‖²_H on odd ξ (≥ 0.0461‖ξ‖²_H for ℓ ≤ 0.025). The inputs are Prime S1 and (Ω5) (both PROVED), Lemma A (sharp pointwise H-bound on |ξ̂|²) and Lemma B (bathtub bound). It uses no CC ε/§5/Thm 6.11, no κ^odd and no Mori. The exact form Q_L = R − 2c₂B is the Weil explicit formula on the window.
     - **This is a small sub-window of the p = 2 step (log 2, log 3], not all of it.** Beyond log 2 + 0.028 the odd entry is open analytically. Positivity there is CERTIFIED only conditional on Mori.
     - **m₀^odd(L₀) ≥ 0.46255 (CERTIFIED, H-margin).** The 24-mode Ritz value ≤ 0.548 (NUMERICAL) is consistent with it.
     - **C.13.H (odd, L₀) is CERTIFIED, CC-free,** at (T, μ_hi, μ_lo) = (12.875, 0.6, 0.25) and (16, 0.7, 0.3).
     - **The odd margin is O(1) in H-norm (0.46), so there is no CCM tiny scale there.** N.13.A's e^{−21/m₀} loss is a method artefact.
     - **The frame route (C.14.A) is off the odd critical path.** δ_o⁰, k_∞, B and 2ε′(1⁺) are not certified (δ_o⁰ and B are GAPs, ε′(1⁺) is NUMERICAL), and none of them is needed for the odd sector. λ₁ = 0.424 (Golden Eye's sinc proxy) does not bound δ_o⁰. C.14.A's pointwise σ_S is not established.
     - **The even sector stays GAP:** (P), m₀^e and C.14.K. The same bathtub method gives ≈ 0 for even at ℓ = 0 (inconclusive).
     - **Independent audit: Golden Eye, SOUND.** Golden Eye confirms the Arb chain is rigorous end to end, monotonicity is analytic, Lemmas A and B hold, and CC-free / no Mori is accurate. For 2L < log 3, Q_L = R − 2c₂B is the classical **Weil explicit formula (Weil 1952 / Barner)**, with the pole terms killed by Q-conjugation, alongside N.13.1. S1 and (Ω5) were re-derived as PROVED.
     - **The 0.028 ceiling comes from the method, not the form (NUMERICAL, Golden Eye Galerkin).** min Q_L/‖ξ‖²_H on odd ξ is ≈ 0.243, 0.073, 0.013 and 0.0005 at ℓ = 0.1, 0.2, 0.3 and 0.4. **Next odd step:** a certified Galerkin lower bound (Temple / Lehmann–Goerisch, plus a Lemma A tail bound) toward ℓ ≈ 0.2–0.3.
     - **m₀ is O(1) on both routes.** The Ritz upper bounds are m₀^odd ≤ 0.549 and m₀^even ≤ 0.073 (Prime §N.14.2, NUMERICAL), so 0.46255 ≤ m₀^odd ≲ 0.549. N.13.A's e^{−21/m₀} smallness is an artefact of its band split.
     - **{2, 3} forms for 2L ∈ (log 3, log 4] (Prime §N.14.3, merged S103g; PROVED, elementary; Grok-checked).**
       - The 3-edge block is +c₃(1/ℓ₃ − ℓ₃/24)σ_x, with c₃ = log 3/√3, so the odd 3-edge charge is exactly 2c₃β²(1/ℓ₃ − ℓ₃/24).
       - B₂(ξ_s) = 0, so the 2-line has no 1/ℓ₃ term.
       - The combined S1 bound is |2c₂B₂ + 2c₃B₃| ≤ (c₂ + c₃)‖ξ‖²_{H(E₂)}, which is sharp by Kronecker.
       - At log 4 the p² = 4 line enters with block +c₄(1/ℓ₄ − ℓ₄/24)σ_x, and ζ = P₃ on the window for 2L ≤ log 5.
       - "The {2, 3} odd charge is paid by the (≤ 2) margin" is a **CONJECTURE**. Caveat: Q^{(≤2)} on the larger window is the P₂ form, which is not known to be ≥ 0.
     - **C.14.K cross-check (Prime, uncommitted; NUMERICAL).** It holds only for ℓT ≲ 0.26 and fails for ℓ ≥ 0.03.
     - **Celestial C.14.4c (merged S103h(3), b02e86f; Grok audit SOUND, NUMERICAL).** **Conjecture C.14.K is refuted numerically** on every useful range: nested-Galerkin lower bounds give even γ ≥ 0.187 / 0.449 / 0.837 at ℓ = 0.03 / 0.06 / 0.1155 (T = 14, γ* = 0.130); it also fails at T = 10 and 20. (Quantifier caveat: this refutes ℓ_K ≳ 0.03; ℓ = 0.01 is unresolved. Promotable to rigorous with one explicit Arb-evaluated profile pair.) Hence **the C.14.K / m₀^e / band-split route to even (K3) is closed beyond ℓ ≈ 7×10⁻⁴**: C.13.H^e (CONJECTURE; NUMERICAL evidence with μ^e_lo < 0.048) gives (K3) only for ℓ ≲ 7.4×10⁻⁴ (Prop C.13.H^e′, PROVED given C.13.H^e); Prop C.13.H^e″, (P) + (HB½) ⇒ C.13.H^e, is PROVED given CC; R alone loses even positivity at ℓ ≈ 0.05 (NUMERICAL). **Adopted route for even (K3) on ℓ ∈ [0.05, 0.1]:** a direct certified even Galerkin lower bound, meaning an interval-certified finite matrix, a Schur tail bound (the complement bound is still a GAP to supply) and a Lipschitz step in ℓ, with the p = 2 term kept signed and no m₀^e, C.14.K or (P) needed. Target m_cert ≈ 0.005 at ℓ = 0.1 (Golden Eye minima 0.046 / 0.015, NUMERICAL). Analytic framework: Celestial; Arb certification: S&M. ℓ = 0.2–0.3 stays GAP. Also merged: the C.13.B Re ψ monotonicity proof, C.13.H recorded as CERTIFIED (S&M), and the C.14.A rerun with 14.42 (window [5.41, 13.6], NUMERICAL).
     - **Prime §N.15 (merged S103h(5); Grok audit SOUND).** The joint {2, 3} edge block has eigenvalues 1.8019 / −0.5333, i.e. 0.9009 / −0.2666 per unit energy (PROVED; matches P.15.1). So the odd high-band payer needs Ω(T) > 0.9009(1 + η), i.e. T ≳ 15.47 (NUMERICAL, mpmath; Stirling agrees). **NUMERICAL (Ritz upper bounds): the near-degenerate direction appears between ℓ₂ ≈ 0.1 and log 1.5 ≈ 0.405, and on (log 3, log 4] the {2, 3} minima are small and falling.** The P₂ form drops from 0.0168 / 0.2452 (even / odd) at 2L = log 2 + 0.1 to 4.5e−4 / 0.0142 at log 2 + 0.3 and 2.5e−4 / 1.34e−3 at log 3 (N = 320, still falling), then goes negative by log 3.5. The full {2, 3} form stays ≥ 0 with minima ≤ 8.5e−4 even and ≤ 2.3e−3 odd (N = 160, falling). Route (A), a P₂ margin, is a dead end. **The CC-free route on (log 3, log 4] reduces to N.15.H** (band coercivity of the P₂ form at log 3 with μ_hi > 0.9009; CONJECTURE) **plus N.15-G1** (a certified non-multiplier margin of size ≲ 1e−3; GAP) **and N.15-G2** (high-band coercivity including the B₂ cross terms; GAP). B2 (N.15.H ⇒ odd entry on ℓ₃ ≤ ℓ₃* ≈ 1e−4) is a PROVED reduction. The even {2, 3} entry is GAP. Positivity over the whole range stays CERTIFIED, conditional on Mori. Also folded in: Prime's C.14.K cross-check (NUMERICAL; agrees with C.14.4c), the P.15.5 fixes, and m₀ quotes now given as Ritz upper bounds (odd ≤ 0.5493, even ≤ 0.0730).
     - **Principal §P.16 (merged S103h(6), 4b4dab0; Grok audit SOUND, with notes).**
       - **The sharp joint S1 constant on 2L ∈ (log 3, log 4) is c₀\* = ½(c₃ + √(c₃² + 4c₂²)) = 0.90093, in both parities (PROVED).**
         - It comes from a fibre decomposition into weighted 4-paths (c₂, c₃, c₂) plus 2-paths (c₂): |2c₂B₂ + 2c₃B₃| ≤ c₀\*‖ξ‖²_{H(F₃)} + c₂‖ξ‖²_{H(G₂)}, with |F₃| = 4ℓ₃.
         - The bound is attained by Perron-weighted bumps on one 4-path.
         - This is a FIX to N.14.3: c₂ + c₃ = 1.1244 is sharp on ℝ but not on 𝒯_L.
       - **High-band payer: Ω(T) ≥ c₀\*(1 + η) for |t| ≥ 15.5 (η = 0, Principal's Arb value Ω(15.5) − c₀\* ∈ [1.86e−3 ± 7e−12]) and for |t| ≥ 16 with η ≤ 0.0373 (S&M's certified Ω(16)). CERTIFIED.**
       - **3-entry (Thm P.16.B, PROVED given (H₃^σ)):** ℓ₃ ≤ πμ/(4c₀\*(1 + 1/η)T).
         - (H₃^σ), the low-band margin of Q^G_L = R − K_G, is a GAP in both parities. It is N.15.H with the joint-block part of G2 discharged, and N.15-G1 for Q^G_L is still needed.
         - **There is no unconditional 3-entry range, odd or even.**
         - The indicative odd size ℓ₃\* ≲ 1.5e−5 is NUMERICAL (Prime Ritz).
         - "S1 + Ω5 alone cannot pay" is a REMARK.
       - **Q.8.E:** D_L = A_LM_ρA_L^\* is the coupled (I2) operator, and (I2) ⇔ s_L ≥ 0 ⇔ Q_L ≥ 0 on the pole-free class.
         - This is PROVED by citation to Krein extension. Audit: apply it to (¼ − ∂²)d, read "invisible" in the weak sense, and cite a distributional version.
         - The explicit n_L and the coupled reach of c₀\* are GAPs. RH OPEN.
  4. *Archimedean margin at ½ log 2:* odd sector needs it. Two conditional entry theorems now exist: Prime's Thm N.13.A (H-margin m₀, rate ℓ* ∝ m₀e^{−21.06/m₀}) and Spectrum & Measure's Thm 2.1 (L²-margin δ_o⁰ of K_L on odd functions, rate ℓ₀ ≈ min{0.006δ_o⁰, (δ_o⁰)²/(8.32k_∞²)}, polynomial). **Grok note:** the second avoids the exponential loss because the low band is controlled at L²-level by 1 − K_L, consistent with Celestial (3b) (K_I is blind only at H¹ level, i.e. at high frequency). So the odd entry is reduced to three certifiable numbers: δ_o⁰, k_∞ = ‖κ‖_{L^∞[0, log 3]}, B = 2ε′(1⁺) + ‖q‖₁ (GAPs; PROLATE_BOUNDS gives the analytic tails, head terms need Arb).
- **Principal §P.14 (merged S103g; Grok audit: SOUND).** P.14 answers the Tate question **negatively (PROVED)**. Piece 2's ramp form is *not* Tate's local reflection Θ₂ = P∘γ₂(V(2)) restricted to the edges. It is the Q-norm edge compression of the m = 1 mode of the 2-adic time delay γ₂^*[U, γ₂] = Σ_m Λ(2^m)2^{−m/2}(V(2)^m + V(2)^{−m}).
  - Four obstructions: (i) Λ(2) needs a derivative; (ii) J winds once while Θ₂ winds zero times; (iii) compressed Θ₂ has diagonal ½ on both parities, so no sign split; (iv) −ℓ/24 is the ¼-mass of Q = s(1 − s).
  - The parity split comes from J∘P (archimedean-global). V(2) + V(2)^* is central in C^*(⟨P, V(2)⟩), and the harmful cones are its V(2) = +1 two-site pieces.
  - **Consequence: piece 4 (the archimedean margin) is the only payer for the odd cost 2c₂α²(1/ℓ − ℓ/24).** It cannot be paid 2-adically.
  - {2, 3}: P is the unique window-preserving reflection commuting with all B_p, and there is no common J_p (PROVED). The N.14 prediction is a CONJECTURE.
- **Principal §P.15 (merged S103h, 2cbbdb1; Grok audit: SOUND).** On 2L ∈ (log 3, log 4] (RH-silent; ζ ≡ P₃ there).
  - **Lemma P.15.A (PROVED):** the orbit-closed node set N_L = {±L, ±(L − log 3), ±(−L + log(3/2)), ±(L − log 2)} is forced. The left 2-strip splits into pieces of lengths (ℓ₃, q, ℓ₃), q = log(9/2) − 2L, and τ_{log 2} maps them onto −I₃, −I₂, −I₁.
  - **Thm P.15.B (PROVED):** exact B₂, B₃ on joint ramps. *Grok replay:* independent exact sympy integration of ⟨u, J₂u⟩_H gives B₂ − hand = 0 (`s103h/p15_joint_ramp_check.py`), and the single-ramp profile recovers N.13.4's c₂(1/ℓ₂ − ℓ₂/24).
  - **p = 3 entry (PROVED on the ramp space):** the order-1/ℓ₃ line block is σ(1/ℓ₃)[[2c₃, 2c₂],[2c₂, 0]], indefinite in **both** P-parities; per unit energy the eigenvalues are σ·{+0.9009, −0.2666} (replayed). So "P-parity decides help/charge per line" (P.14.8(d)) is refuted, the even sector acquires a charged direction at the 3-entry, and the odd worst charge is 0.9009, not c₃ = 0.6343. That 0.9009 is the sharp S1-type constant there is a CONJECTURE.
  - **Reconciliation with Prime N.14.3:** both are correct. "No 1/ℓ₃ term from the 2-line" holds relative to ‖ξ′‖_∞, not relative to H-energy. This sharpens the caveat on Prime's {2, 3} odd-entry CONJECTURE (Q^{(≤2)} on the larger window); it does not refute it.
  - **Q.8.E (PROVED):** the weight-1 block −(c₂/ℓ)S(V(2) + V(2)^*)S has spectrum {0, ±c₂/ℓ}; positive line blocks A_nM_ρA_n^* cost c₂sec(π/(n + 2))·g(0), with infimum c₂ sharp. A **separable (I2) forces R ≥ c₂‖ξ‖²_H on 𝒯_L (global S2)**, whose truth is a GAP (HEURISTIC: false). A working D_L must couple the line block to the archimedean ε-part, band-split.
  - **Piece 4:** still the sole odd payer, now in band-split, ℓ-dependent form beyond log 2 + 0.028 (bookkeeping PROVED; values GAP N.G2-m₀). The C.14.A failing band ≈ [5.41, 13.6] is NUMERICAL (Celestial/Prime). Audit note: P.15's "sympy-checked" claim had no committed script; the Grok replay above now supplies one.
  - **Grok audit S103h(2) (15e5b58):** the {2,3} joint ramp block σ(1/ℓ₃)[[2c₃, 2c₂],[2c₂, 0]] is indefinite in **both** P-parities (PROVED; det = −4c₂²σ²/ℓ₃² < 0; replayed). Hence on 2L ∈ (log 3, log 4] no parity-by-parity payer exists: the payer there must be band-split and ℓ-dependent, covering the joint per-energy charge 0.9009 at frequency ~π/ℓ₃ (size CONJECTURE). Label correction: "separable (I2) ⇔ global S2" is a REMARK (only ⇒ is PROVED, for line blocks with diagonal sums at {0, ±a}); global S2 is a GAP. P.14's G.4h-item-2 complaint and N.13.4 fix request are withdrawn (P.14.1, P.15.3); P.15.5's two wording FIXes to Prime N.14 (endpoint rounding 1.2005/0.6930; "no 1/ℓ₃ term" only relative to ‖ξ′‖_∞) are passed to Prime, not edited here.
- **Even sector has two separate costs (Grok reconciliation, PROVED as bookkeeping).** (i) the CC ĝ(0) eigen-direction (μ₁ ≈ 1.0518, cost 2ε′(μ₁ − 1) ≈ 2.38, O(1), low-frequency, L²-level), which must be paid by 𝒮_L: rank-one Schur condition (P) (§S 4.3, GAP b-1: σ_min on |t| ≲ 6.3); (ii) the edge traces, which the p = 2 line pays (N.13.4). Prime's "even entry is carried by p = 2" and §S's "even payment comes from 𝒮_L" concern (ii) and (i) respectively and do not conflict.
- **Status.** On (½ log 2, ½ log 3] positivity is already CERTIFIED (conditional on Mori; main 992a011). These results answer *which step carries it* and reduce the analytic proof to: N.G2-bulk (GAP), (P) (GAP), and certified values of δ_o⁰, k_∞, B (computational GAPs). None of this touches G.4d (uniformity in L), which remains the RH-equivalent GAP.
- Labels: R1 fixed range | R2 FE via 2θ′, Euler product via the single line at log 2 | R3 the sign of B by parity | R4 even nodeless edge helped, ĝ(0) direction paid by 𝒮_L | R5 fake-blind in truth; DH has no a-line.

**G.4i (S103i synthesis, routine 06:10Z: merged §C C.15 + §N N.16 + Queen H-κ Part V + S&M citation fix; Grok audit: all four SOUND).**
- **Celestial §C.15 (f13efdd).** *Thm C.15.A (PROVED; audited step by step):* Q_L ≥ Q′_T because q − q_T = [Ω(t) − Ω(T)] + 2c₂(1 − cos at) ≥ 0 on |t| ≥ T (Ω = Re ψ(¼ + it/2) − log π is increasing in |t|: ∂_y Re ψ(a + iy) = Σ 2y/((n + a)² + y²)² > 0); the even Dirichlet f-basis is H-orthonormal, so the q_T part has no V_N–W cross term; the cross term is Cauchy–Schwarz in the |a_T|-weighted measure with Lemma C.15.K (|ŵ|² ≤ K_W‖w‖²_H, Σ|c_k|² = ‖w‖²_H) on W; the 2×2 Schur step is exact. This **closes the C.14.4c complement GAP** with a proved bound. *Prop C.15.B (PROVED):* for 2L ≤ log 3 the symbol q = Ω − 2c₂cos(at) does not depend on L, and extension by zero preserves ξ̂, evenness and the H-norm, so λ_e(L) is nonincreasing: **one certificate at ℓ_max covers (0, ℓ_max]** (the 250-point Lipschitz grid is correctly withdrawn; its C_B was unproved). *Thm C.15.D* (even (K3) on ℓ ∈ (0, 0.1], float λ_cert = 0.01385 with analytic majorants) stays CONJECTURE until S&M encloses μ_A, s_V, ε^{abs}, ε^{neg} in Arb. *C.14.K:* CERTIFIED false (Arb, 128 bits) at (ℓ, T) = (a/6, 14), ratio 0.7482 vs γ* = 0.1298, and at (a/13, 14), ratio 0.2932; refuted NUMERICALLY on the useful range ℓ_K ≳ 0.03.
- **Prime §N.16 (5c544bf).** Mechanism for the tiny {2,3} margins (RH not assumed). *PROVED (counting):* the crossover t* = 2πe^{2L} is where the zero density (1/2π)log(t/2π) meets the window's Nyquist density 2L/2π; per parity the dimension surplus below t* is D(t*) = 2Le^{2L} − [e^{2L}(2L − 1) + 7/8] − 1 = e^{2L} − 15/8 (smoothed RvM; Grok re-derived), so for D > 0 a window-limited vector exists that vanishes at every zero below t* and pays only prolate leakage ("zero-dodging"). Even onset is below log 2, odd onset ℓ₂ ≈ 0.17 (matches N.15). *Upper bound* m^σ(L) ≤ C log t*·(1 − λ_k(Lt*)): Grok label **PROVED (sketch), constants untracked, upper bound only**. *CONJECTURE (Landau–Widom plunge):* log m^σ(L) ≈ −2π²e^{2L}/(2L + log 2πL), super-exponentially small in L, which is the analytic content of G.4d: **any step that gives up a fixed fraction of the margin closes at bounded L**. ONE check: the explicit-formula zero sum reproduces the N.15 Ritz values to ≤ 1.1%, and the even near-null vector has |ξ̂(γ₁)|², |ξ̂(γ₂)|² ~ 1e−9 of its peak; 94% of its Ritz value is discretisation, so the true even margin at log 3 is likely ≲ 1e−4 (NUMERICAL). The zero side gives no unconditional lower bound ("all terms ≥ 0" is BARRED).
- **Queen H-κ Part V (0d48b28).** S^{−1} = (G^{−1})_{RR} and the dual variational bound ⟨S^{−1}ρ, ρ⟩ = sup_ψ |⟨ψ, ρ⟩|²/G(ψ, ψ) (PROVED; checked). (LD_γ) retired (no unconditional source; Radziwiłł–Soundararajan controls the upper tail only). Pointwise lower tail ≫ ηT^{5/6−ε} (PROVED sketch; Grok check: |ζ| ≤ sup|ζ′|·|t − γ| near any zero, so simple zeros are not even needed), which corrects Part IV's "super-polynomially rare". New chain (T3) ⇐ (WK), (WK₂) + (NM₂) + (MB)(a–c) + flat density + (T3a*) + (L): PROVED-conditional. **(NM₂), a window-floor negative 4th moment ≤ Λ^{2−δ}/r_K, is the single arithmetic GAP**; unconditional worst case gives exponent ≥ 8/3 (insufficient; checked: Chebyshev with c ≥ 1/6, T ≍ Λ² gives η_W^{−4} ≤ Λ^{16c}); RH ⇒ (NM₂) so it is BARRED as a source; the zero-density route (Euler product, R2) is OPEN. (T3) OPEN.
- **S&M (bad5094, 945d87c).** Citation fix: the Λ > 1 analogue of CC Thm 6.11 is Principal §P.10.3 (theirs); γ in c = 4γ/log 2 is CC Lemma 6.10's 2.94355, not Euler's constant. t3a_test.py committed (NUMERICAL, RH-silent).
- **Connection (Grok, HEURISTIC).** N.16 and Queen's (NM₂) are the same phenomenon at two scales: the Weil form's margin is set by how small ξ̂ can be made on the zeros (N.16, low window) and by how small |ζ| is on unit windows (Queen, height Λ²). Both say the obstruction is quantitative small-value control of ζ, where RH-equivalence enters, not a missing positive structure.
- Labels: R1 fixed range for C.15/N.16, Λ-uniform for Queen | R2 Euler product via the n = 2 line (C.15), zero density (Queen) | R3 p = 2 term kept signed | R4 f_0 carries the ĝ(0) direction inside V_N | R5 DH structurally absent from C.15; N.16 mechanism holds for any zero set, so fake-blind as stated. RH OPEN.

### G.5 What would count as progress (CONJECTURE-level targets, each R1–R5-checked)
- **(T1) Pole/boundary compensation (WITHDRAWN S103 by G.3d).** The k = 2 rank-one needs no compensation (it is absorbed once the fluctuation is taken against ψ = Σ_k θ(x^{1/k}), G.3d), and the two boundary directions are already spent on the mean field (G.3a). Original statement kept for the record:
  - The even k = 2 rank-one −¼·11ᵀ (Astra) and the k = 1 fluctuation are compensated *exactly* by the two Poisson boundary directions (φ(0), φ̂(0)) that 𝒮(𝔸)₀ removes.
  - Test against R5: DH has no boundary directions and no MF, so the conjecture must fail for DH. Good: it is not fake-blind.
- **(T2) Tail engagement (CONJECTURE).** Any compensating operator for Fl_L is non-quasi-local: Astra's support bound gives the necessity half, but only conditional on H1 (Astra's zero ↔ coherent-mode identification, a GAP at least as strong as RH; §C C.7c, replay C.6d). The sufficiency half is the GAP.
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
