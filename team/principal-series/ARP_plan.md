# ARP / Principal series section: PLAN and OUTLINE only (not the derivation)

Daniel Toupin / Golden Physics Project, with AI assistance from Grok. 2026-10-22. **RH is OPEN.** Nothing here proves or implies RH.

**Status update 2026-10-28.** The derivation is now `ARP_principal-series.md` (§P). Its P.10 lists the label changes to this plan; in particular the §5 F_θ row is corrected (GL(2) data unitary, Knapp–Stein PD for 0<θ<½; the filter acts on the GL(1) constituents) and §4's R4 conjecture is superseded. The text below is kept unchanged for the record.

**Status.** This file is a plan. `ARP_principal-series.md` is **on hold** until the lead's `derivation/LEAN_INVENTORY.md` lands on GPPDiscovery5 main; it had not landed as of 07:55 EDT.

**Labels.**
- PROVED: literature or an elementary identity, with a precise citation.
- PROVED-LEAN: kernel-checked, with the GPPVerify name.
- CONJECTURE.
- GAP.
- NUMERICAL: team data, used only as orientation.

**Not cited.** `open_*` theorems are `True := trivial` placeholders and are never cited.

**Scope (lead, C).** My section covers the local reflections Θ_v: the Tate local functional equation, the normalized intertwiners, and the temperedness point. These are the pieces that Connes' map
$$E: \mathcal S(\mathbb A)_0\to L^2(C_{\mathbb Q}),\qquad (Ef)(x)=|x|^{1/2}\sum_{q\in\mathbb Q^*}f(qx)$$
intertwines (Connes 1999; Meyer 2005). No place-by-place positivity is claimed or needed.

**Requirements.**
- R1: uniform in L.
- R2: uses the Euler product and the functional equation jointly.
- R3: F_θ fails and G_θ′ passes, independently of prime-weight signs.
- R4: the positive space contains the soft mode.
- R5: say where each fake breaks.
- R6: an exact identity, with no fixed-constant slack. The soft-mode margin is about 2^−4900 at L = 1.42.

---

## 0. Prior art to cite (verified links; the dead ends are cited, not re-derived)
- **Weil 1952.** "Sur les « formules explicites » de la théorie des nombres premiers", Comm. Sém. Math. Univ. Lund, vol. suppl., 252–265: explicit formula and positivity criterion.
- **Tate 1950.** Thesis, in Cassels–Fröhlich, *Algebraic Number Theory* (1967), Ch. XV: local and global zeta integrals and functional equations.
- **Bombieri 2000.** "Remarks on Weil's quadratic functional…, I", Rend. Lincei Mat. Appl. 11, 183–233. The minimum on [−t, t] is attained. If RH fails with finitely many exceptions, the number of negative eigenvalues is half the number of off-line zeros (for a big enough truncation).
- **Yoshida 1992.** "On Hermitian forms attached to zeta functions", Adv. Stud. Pure Math. 21, 281–325: positivity for small support.
- **Connes 1999.** "Trace formula in noncommutative geometry and the zeros of the Riemann zeta function", Selecta Math. 5, 29–106: semi-local trace formula; the global trace formula is equivalent to the positivity.
- **R. Meyer 2005.** "On a representation of the idele class group related to primes and zeros of L-functions", Duke Math. J. 127, 519–595: unconditional spectral realization on a nuclear Fréchet space, not a Hilbert space, so **no positivity**.
- **Connes–Consani.** "Weil positivity and trace formula, the archimedean place", Selecta Math. 27:77 (2021), arXiv:2006.13771. The Sonin trace is positive, and its difference from W is controlled via prolate functions, archimedean place only. Also "Spectral triples and ζ-cycles", arXiv:2106.01715.
- **Burnol.**
  - "On Fourier and Zeta(s)", Forum Math. 16 (2004), arXiv math/0112254.
  - "Two complete and minimal systems…", JTNB 16 (2004) 65–94, arXiv math/0203120 (co-Poisson, Sonine spaces).
  - "Sur les formules explicites I", C. R. Acad. Sci. 331 (2000) 423–428 (the explicit formula as a conductor operator).
- **Shahidi.** *Eisenstein Series and Automorphic L-functions*, AMS Colloq. 58 (2010); Ann. of Math. 132 (1990) 273–330: normalized intertwiners and local coefficients.
- **Langlands.** LNM 544 (1976): constant term M(s) and global functional equation.
- **Lagarias.** "Hilbert spaces of entire functions and Dirichlet L-functions" (2006).
- **Suzuki.** Canonical systems / screw functions attached to ζ (2012, 2022).
- **Temperedness.**
  - Cowling–Haagerup–Howe, J. reine angew. Math. 387 (1988): tempered ⇔ matrix coefficients in L^{2+ε}.
  - Bargmann 1947: complementary series.
  - Knapp–Stein: intertwining-operator forms.
- **Dead ends (cite only):** log-time OS; Schur / Kreĭn / doubling restatements; positive cone; Selberg on the modular surface; η_ζ = 0 as a target. Team logs: bridge `work/grok/*`, S82-OS, S95–S98.

## 1. Lean names to cite (GPPVerify main 8a49a74, read-only; non-placeholder only)
- **Twists (planned, not yet formalized):** `GppWeilTwist.twist_cos_sum` and `twist_cosh_krein` (our LEAN_ASK.md).
- **Finite Weil forms:** `GppNoGhost.weilPairing`, `weilForm`, `weilForm_nonneg`, `PositiveTypeOn` (FiniteNoGhost / TruncatedTransport).
- **Positive-type basics (HaarPositivityWeil):** `PositiveType.even`, `positive_type_at_zero`, `const_one_positive_type`.
- **Weil criterion, zero side only:** `GppWeilCriterion.rh_iff_weil_pairedForm_nonneg` (WeilPositivityCriterion.lean) and `rh_iff_two_point_pairedForm_nonneg` (TwoPointCriterion.lean).
  - *Caveat:* these are statements about the zero-side paired form, i.e. combinatorial bookkeeping. They are **not** Weil's test-function criterion, and the file says so.
- **Support fact (Goddess):** `GppWeilSupportLadder.convolution_hasSupportIn`, `primeSide_eq_truncation`, `primeSide_eq_zero_of_support_lt_log_two`. On [−L, L], Q_L only sees n < e^{2L}.
  - `weil_nonneg_of_arch_nonneg_rung_zero` is **conditional** on `0 ≤ arch`. Cite it as a reduction, not as positivity.
- **Tate, local p-adic:** `GppPadicFullZeta…tate_local_zeta_integral` (PadicZetaIntegralClosedForm.lean, ∫_{ℤ_p}‖x‖^s dμ in closed form) and `euler_factor_bridge` (PadicEulerFactorBridge.lean, which matches Mathlib's `riemannZeta_eulerProduct_tprod`).
- **Ideles:** `RationalIdeleGroup := (AdeleRing ℤ ℚ)ˣ`, `diagonalEmbedding_injective` (IdeleGroup.lean). Mathlib has `NumberField.AdeleRing`.
- **Functional equation:**
  - Mathlib `completedRiemannZeta_one_sub`.
  - `completedRiemannZeta_logDeriv_reflection` (CompletedZetaDerivativeSymmetry.lean).
  - `NoGo.ShiftedProductControl.F_one_sub`, `F_zero_add`, `F_zero_sub` (the F_θ no-go).
- **Log-derivative additivity:** `neg_logDeriv_product_of_hasDerivAt` (LogDerivativeProduct.lean): "no cross term at the local-factor level". `finiteWp_eq_two_mul_re_finitePrimeLogDerivative` (CompletedLogDerivativeBridge.lean).
- **Von Mangoldt:** `GppGlobalVonMangoldt.neg_zeta_logDeriv_eq_vonMangoldtLSeries` (Re s > 1).
- **Do NOT cite as unconditional:**
  - `temperedness_iff_critical_line` and everything in ArithmeticTemperedPrincipalSeries.lean. They depend on the custom axiom `exp_growth_not_tempered`.
  - ArithmeticCoPoissonOrientationTemplate.lean: abstract intertwiner algebra only. The co-Poisson identity itself is assumed there.
- **Not in GPPVerify:** global Poisson summation on 𝔸/ℚ, Tate's global functional equation, the explicit formula (only an `open_` placeholder), the map E, GL(2) intertwiners. Mathlib has the real-line Poisson summation (`Real.tsum_eq_tsum_fourierIntegral`; I will check the exact name at build).

## 2. Outline of the section (labels as currently assessed)

**§1. Local reflections.**
- **Tate local FE (PROVED).** Tate 1950, Thm 2.4.1: Z(f̂, χ^{−1}|·|^{1−s}) = γ_v(χ, s) Z(f, χ, s).
- **Local isometry (PROVED; Tate + Plancherel on ℚ_v^*).** With the half-density normalisation f ↦ |x|^{1/2}f, the map L²(ℚ_v, dx) → L²(ℚ_v^*, d^*x) is an isometry, up to the Haar constant. Θ_v := F_v composed with x ↦ x^{−1} is unitary. In Mellin coordinates it is the multiplier γ_v(χ, ½+it), with |γ_v| = 1 for χ unitary.
- **This is the Hermitian structure that E sees:** the L², i.e. tempered, structure.
- **GL(2) picture (PROVED; Langlands 1976, Shahidi 1990/2010).** Λ(s) enters as the constant term of the Eisenstein series on Ind_B^{GL2}. The normalized local intertwiners N_v(s) are unitary on the tempered axis.
- **Which form is Hermitian.** On the complementary series (real exponent, 0 < θ < ½), the *invariant* Hermitian form is the Knapp–Stein form ⟨f, M(θ)f⟩. That form is not the L² form, and E does not intertwine it.

**§2. Temperedness as a necessary local filter, not a source of positivity** (Golden Eye's audit, accepted).
- For GL(1), every unitary character is tempered, so ζ satisfies the filter automatically. The filter contributes **nothing** to ζ's positivity, and as a mechanism it fails R6. Its role is R3: it rules out F_θ.
- **Twist identities (PROVED, elementary; Weil 1952 / Burnol 2000; our LEAN_ASK).** Let B_ζ(g, h) = ∬ g(x) h̄(y) W_ζ(x−y). Then:
  - Q_G(f) = B_ζ(e^{iθ′u}f, e^{iθ′u}f) + B_ζ(e^{−iθ′u}f, e^{−iθ′u}f), a **sum of ζ-norms**;
  - Q_F(f) = 2 Re B_ζ(e^{θu}f, e^{−θu}f), a **cross pairing of two different vectors** (a polarization), hence indefinite.
- On C_ℚ, |x|^{iθ′} is a unitary character and |x|^θ is a non-unitary quasi-character. Both are trivial on ℚ^* and commute with E.
- **What this gives (PROVED).**
  - Any positivity mechanism that is an *L²(C_ℚ) norm or trace identity* automatically passes G_θ′ (its unitary twist is a norm).
  - It does not formally extend to F_θ: the twist turns the norm into a polarization.
  - So R3 is automatic for any L²-norm mechanism. The F_θ failure is "E's transfer composed with a non-unitary multiplier". It is the local non-temperedness |γ_v(|·|^θ, ½+it)| ≠ 1, seen globally.
- **What it does NOT give.** Anything about ζ itself: the trivial character is unitary. It does not explain F_θ's onset L either. The onset is set by ζ's soft-mode margin, λ_min(Q_ζ) = θ²κ (NUMERICAL, our session 2 §5).

**§3. The transfer question (the shared GAP, stated in Poisson-gluing terms).**
- **Additive positivity (PROVED; Poisson on 𝔸/ℚ, Tate 1950 §4.2).** For g ∈ 𝒮(𝔸), Σ_{ξ∈ℚ}(g ⋆ g̃)(ξ) = Σ_{ξ∈ℚ}|ĝ(ξ)|² ≥ 0. The comb Σδ_ξ is a positive, positive-definite, self-dual measure on (𝔸, +).
- **Equivariance (PROVED; Connes 1999 §III).** E is C_ℚ-equivariant. For f ∈ 𝒮(𝔸)₀, i.e. f(0) = f̂(0) = 0 (the two Poisson boundary terms, i.e. the poles at |·|^0 and |·|^1, are removed), E(F_𝔸 f)(x) = E(f)(x^{−1}). So E intertwines the global reflection F_𝔸 = ⊗′Θ_v with x ↦ x^{−1} on C_ℚ, exactly.
- **Weil positivity is a statement about the multiplicative group.** It asks W(h ⋆_× h̃) ≥ 0, i.e. multiplicative positive-definiteness. The Weil distribution is Connes' trace of the scaling action on the *cokernel* of E with cutoff, not on its range. Note ‖Ef‖² ≥ 0 is trivial and is not the Weil form.
- **GAP (sharp).** The Poisson gluing step must preserve positivity on the cutoff space. Precisely, does the additive positive-definiteness of the self-dual comb (Σ|ĝ|² ≥ 0) pass, through E and the cutoff projections P_Λ, P̂_Λ, to multiplicative positive-definiteness of W on the window [−L, L] **as an exact identity** (R6)? Concretely:
  - W(k ⋆ k̃) = Tr(Π_Λ ϑ(k)ϑ(k)^* Π_Λ) + δ_Λ(k), where the first term is ≥ 0 manifestly;
  - the requirement is δ_Λ ≡ 0, or δ_Λ a manifest square, for **all** Λ, L (R1).
- **What is known about that decomposition.**
  - Connes 1999: the semi-local formula with error o(1) as Λ → ∞. Not exact, so it fails R6.
  - Connes–Consani 2021: an exact difference in prolate terms for the archimedean place alone. Positivity is proved only below the first prime rung.
- **Equivalent to RH.** By Goddess's support fact (PROVED-LEAN `primeSide_eq_truncation`), window L is exactly a semi-local problem with S_L = {∞} ∪ {p : p < e^{2L}}. So the GAP for all L is equivalent to RH (Weil 1952; Bombieri 2000).
- **What the local analysis contributes to the transfer (exactly).**
  1. Each Θ_v is unitary on the half-density L² space. Hence E transports the L² structure place by place without loss, and the only non-isometric step is the global ℚ^*-sum, i.e. the gluing.
  2. The unitarity holds *only* for tempered local data. That is the R3 filter (§2).
  3. The local Euler-factor operators are where R2 enters (Euler product ⇔ ⊗′ factorization of 𝒮(𝔸) and of F_𝔸).
  - Nothing in (1)–(3) bounds δ_Λ. **The positivity must be produced at the gluing, not locally.**

**§4. Where the soft mode sits (R4).**
- NUMERICAL (Queen S1): ζ's window ground state φ_L has |cos(φ_L, pole vector 2cosh(x/2))| = 0.995–0.9997.
- CONJECTURE: the soft mode is the window image of the two Poisson boundary directions (f(0), f̂(0)), i.e. the |·|^0 and |·|^1 poles that 𝒮(𝔸)₀ removes.
- If so, R4 says the cutoff positive space must contain an approximately boundary direction, and δ_Λ is dominated there by the pole/boundary interaction.
- Test, cheap but parked: the overlap of ζ's lowest even eigenvector with the window-restricted cosh(x/2) as L → 1.42, compared with the 2^−4900 margin.

**§5. Fakes vs steps (R5): first pass.**

| fake | §1 local Θ_v unitary (tempered) | ⊗′ factorization (Euler product) | global FE / Poisson boundary | where it breaks |
|---|---|---|---|---|
| G_θ′ | yes | yes (degree 2) | yes, but 4 poles; the pole block is indefinite with signature (1,1) | **passes every local step.** It inherits ζ's certificates by the identity (PROVED). |
| F_θ | **no**: the exponents ±θ are real, complementary series, non-tempered | yes (degree 2) | yes | §1/§2: the twist is a polarization, not a norm. |
| DH | — (no local data) | **no** | yes | the ⊗′ factorization is missing, so there is no E-intertwining of local Θ_v. |
| DMV Beurling / P_X | Beurling local data | Beurling "primes" | no ℚ^*-lattice, so no Poisson; P_X is a truncation | the gluing has no lattice to sum over. |
| Epstein (h = 2) | — | **no**: Z_Q is a sum of two Euler products (genus characters) | yes, via SL(2, ℤ) | a sum of two pure tensors; E of a sum is not a norm of a sum. |
| H₋ₜ | yes | yes | **wrong archimedean factor** (de Bruijn–Newman deformation) | the archimedean Θ_∞ is not the Tate γ_∞. |
| Hecke G₅ | — | **no Hecke operators** (non-arithmetic) | scattering φ′/φ (non-compact) | no adelic lift, so no E. |

## 3. What I will NOT do
- Re-derive the dead ends.
- Claim local positivity.
- Cite placeholders or the axiom-dependent temperedness file as unconditional.
- Treat the twist analysis as progress on ζ.

## 4. Open items for the derivation (once LEAN_INVENTORY lands)
- (a) The exact citation for E(F_𝔸 f) = E(f)∘inv on 𝒮(𝔸)₀ (Connes 1999, §III, Lemma). Check whether Meyer 2005 states the same with signs.
- (b) Write δ_Λ in Connes–Consani 2021 notation, and state exactly which Λ their archimedean positivity covers.
- (c) Make §2's polarization statement a formal lemma on L²(C_ℚ), and add it to the Lean ask.
- (d) The R4 conjecture: formulate it precisely as a statement about the boundary map 𝒮(𝔸) → ℂ², f ↦ (f(0), f̂(0)).
