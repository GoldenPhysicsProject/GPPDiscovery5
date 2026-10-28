# §P Principal series: local reflections, the polarization lemma, and the E-intertwining (for ARP.md §P)

Daniel Toupin / Golden Physics Project, with AI assistance from Grok. Seat: Principal series. 2026-10-28 (S102/S103 directions).
**RH is OPEN. Nothing in this section proves or implies RH.** The section's result is negative and structural: every local step is an isometry, the only non-isometric step is the global ℚ^*-sum, and nothing local bounds the shared GAP G.4d.

**Labels.** PROVED (literature or elementary, cited) · PROVED-LEAN (GPPVerify name, main 5f54359, checked 2026-10-28) · CONJECTURE · GAP · NUMERICAL (orientation only).
**Conventions (team-binding).** L is the HALF-width: f ∈ C_c^∞([−L, L]) in log scale u = log|x|, g = f ∗ f̃ on [−2L, 2L], window L sees n < e^{2L}. Lean's `HasSupportIn g L′` in `WeilSupportLadder.lean` is stated for the convolution g, so **L′ = 2L**.
**Builds on.** ARP.md §G (G.1–G.4, GAP G.4d) and §0.3 (inventory); `ARP_plan.md` §1–§2 (this file supersedes its §2 outline and corrects its §5 table, see P.12). Exact step checks: `arp/p_checks.py` → `arp/p_checks.log` (sympy; no exploratory numerics).

**Fakes and requirements (each step's R-line).** F_θ = Λ(s−θ)Λ(s+θ) (θ real ≠ 0); G_θ′ = Λ(s+iθ′)Λ(s−iθ′) (must pass); DH (Davenport–Heilbronn); H₋ₜ (de Bruijn–Newman deformation); DMV (Beurling); P_X (finite Euler product); Epstein x²+5y² (h = 2); G₅ (Hecke-group, non-arithmetic). R1 uniform in L; R2 Euler product + FE jointly; R3 not via weight signs; R4 soft mode / pole state; R5 locate each fake; **R6 exact identity, no slack** (every step below is an equality unless stated). R6's scope: what is PROVED is "no MF-relative margin" (§S S(c)4). That the closing step *must* be an exact identity is a CONJECTURE (P.7, Golden Eye).

---

## P.0 Normalizations (the two source texts, verbatim numbering)

Checked against arXiv math/9811068v1 (Connes; = Selecta Math. 5 (1999) 29–106) and arXiv math/0311468v3 (Meyer; = Duke Math. J. 127 (2005) 519–595).

| item | Connes 1999 | Meyer 2005 |
|---|---|---|
| additive character | α nontrivial on 𝔸, α|_k = 1, pairing ⟨x, y⟩ = α(xy): App. I eqs (27)–(28); k^⊥ = k: (29) | ψ ∈ 𝔸̂_K, ψ|_K = 1, "normalized": §5.1; F_ψ f(ξ) = ∫ f(x) ψ(xξ) dx: eq (11) |
| Fourier sign | f̂(β) = ∫ f(a) β(a) da, β(a) = α(βa): App. I eq (33) (no minus sign) | same, eq (11) (no minus sign); F^{−1} = F^* = λ_{−1}F, F^*f(x) = Ff(−x): §3.1 |
| additive Haar | Covol(k) = 1: App. I eq (34) | self-dual; 𝔸_K/K has volume 1: §5.1 (the same measure) |
| multiplicative Haar | ∫_{|g|∈[1,Λ]} d^*g ∼ log Λ: §III eq (14) | ∫_{1≤|x|<|u|} d^×x = ln|u|: Intro (Weil distribution paragraph) |
| half-density | E carries |g|^{1/2}: §III eq (18); norm (7) uses |x| d^*x | no |x|^{1/2}; works in L²(C_S)_{1/2} := L²(C_S, |x| d^×x): §3.2 |
| removed boundary | 𝒮(𝔸)₀ = {f(0) = 0, ∫ f dx = 0}: §III eq (6) = App. I eq (31) (∫ f dx = f̂(0)) | no 𝒮(𝔸)₀; keeps the boundary as the rank-2 map R(f)(x) = f(0) − |x|^{−1}Ff(0): §5.4 |
| summation map | E(f)(g) = |g|^{1/2} Σ_{q∈k^*} f(qg): §III eq (18), App. I Lemma 2 | Σf(x) = Σ_{a∈K^×} f(ax): Intro, §5.3 |
| inversion | g ↦ g^{−1} on C_k (unitary on L²(C_k, d^*g)) | J f(x) = |x|^{−1} f(x^{−1}): eq (12), unitary on L²(C_S)_{1/2} |
| equivariance | E U(a) = |a|^{1/2} V(a) E, U(j)f(x) = f(j^{−1}x): §III eqs (9), (15), (19) | 𝔉 := F∘J is λ-equivariant: eqs (13)–(15) |

Dictionary: **E = M_{|x|^{1/2}} ∘ Σ** (Connes' E is Meyer's Σ in half-density coordinates).

---

## P.1 Local Tate functional equation — PROVED; PROVED-LEAN (local factor only)

**Statement.** For φ ∈ 𝒮(ℚ_v) and a quasi-character χ of ℚ_v^*, Z(φ̂, χ^{−1}|·|, 1−s) = γ_v(χ, s) Z(φ, χ, s); for χ unramified trivial, γ_p(s) = L_p(1−s)/L_p(s), L_p(s) = (1−p^{−s})^{−1} (Tate 1950, Thm 2.4.1; ARP G.1a).
**Lean.** `GppPadicFullZeta.tate_local_zeta_integral` (PadicZetaIntegralClosedForm.lean:131, ∫_{ℤ_p}‖x‖^s dμ in closed form, real s) and `euler_factor_bridge` (PadicEulerFactorBridge.lean:67, real s > 1; matches Mathlib's factor in `riemannZeta_eulerProduct_tprod`): **PROVED-LEAN for the local factor identity at real s > 1**. The FE at complex s is not in GPPVerify.
**R-line.** F_θ: local data exist (two GL(1) constituents |·|^{±θ}) — passes. DH: no Euler product, no local factor — breaks (R2). H₋ₜ: no Euler product after deformation — breaks. DMV: Beurling local factors, no ℚ_p and no local Fourier — breaks formally. P_X: same as ζ for p ≤ X, trivial for p > X — passes locally (fails at P.7). Epstein: each genus-character L-function has local data, but Z is a sum of two Euler products — breaks at the ⊗′ (P.6). G₅: no Hecke operators — breaks. **G_θ′ passes** (unitary twists |·|^{±iθ′}). R6: exact.

## P.2 Local reflections are isometries — PROVED (Lean: fixed-line algebra only)

**Statement.** In half-density coordinates f ↦ |x|^{1/2} f, the local reflection Θ_v := F_v ∘ inv is unitary on L²(ℚ_v^*, d^*x); in Mellin coordinates it is the multiplier γ_v(½+it), and |γ_v(½+it)| = 1 for v = p and v = ∞ (ARP G.1b, G.2a).
**Proof / citations.**
- p: |1 − p^{−(a+it)}|² = 1 − 2p^{−a}cos(t log p) + p^{−2a} is even in t (`p_checks.log` (P1)), so |L_p(½−it)| = |L_p(½+it)|.
- ∞: Connes–Consani 2021, **Lemma 1.3** (arXiv:2006.13771, eqs (17), (24)–(26)): with w(ξ)(λ) = λ^{1/2}ξ(λ), the even Fourier transform is F_{e_ℝ} = I ∘ u_∞, I the unitary inversion on L²(ℝ₊^*, d^*λ) and u_∞(s) = e^{2iϑ(s)} (Riemann–Siegel ϑ), a unimodular multiplier. (Their character is e_ℝ(x) = e^{−2πix}; the sign is immaterial on even functions.)
- S-local: Meyer 2005, **Lemma 3.1**: the equivariant Fourier transform 𝔉 = F∘J and 𝔉^* descend to mutually inverse **unitary** operators on L²(C_S)_{1/2}; and 𝔉 = 𝔉^* there because −1 ∈ K_S^× (§3.2, after Lemma 3.1).
**Lean.** `ArithmeticOSReflection.osReflection_eq_self_iff` (fixed locus of s ↦ 1 − s̄ is Re s = ½): **PROVED-LEAN, algebra only**. |γ_p(½+it)| = 1 is ARP G.6 ask 3 (`local_gamma_unitary`), open.
**R-line.** F_θ: each constituent has |γ_v(½+it ± θ)| ≠ 1 (0.817 / 1.224 at p = 7, t = 1.3, θ = 0.2: S102 replay (5)); the product is unitary (P.3). DH, G₅: no local Θ_p. H₋ₜ: Θ_∞ is not Tate's (deformed Φ-kernel) — breaks at ∞. DMV: no local field. P_X: same as ζ for p ≤ X. Epstein: per genus constituent, unitary. **G_θ′ passes.** R6: exact (equalities |γ| = 1).

## P.3 GL(2) data of F_θ are locally unitary; the Knapp–Stein form is positive definite — PROVED (this is correction (b))

**Statement (b1).** For real θ and s = ½+it, |γ_v(s−θ) γ_v(s+θ)| = 1 at every place, while |γ_v(s+θ)| = |L_v(½+it−θ)/L_v(½+it+θ)| ≠ 1 for θ ≠ 0.
*Proof.* L_v(z̄) = conj L_v(z) (p real; Schwarz reflection for Γ). Then |γ_v(s−θ)| = |L_v(½+it+θ)|/|L_v(½+it−θ)| and |γ_v(s+θ)| = |L_v(½+it−θ)|/|L_v(½+it+θ)|; multiply. ∎ Exact check: S102 replay (5) (|product| = 1.000…, constituents 0.81705…, 1.22391…); `p_checks.log` (P1): |γ_p(s−θ)γ_p(s+θ)|² − 1 = 0 exactly at p = 7, and |γ_∞(s−θ)γ_∞(s+θ)| = 1.0000… with |γ_∞(s+θ)| = 1.3906.
*Representation theory.* By multiplicativity of γ-factors (Shahidi 1990, Ann. of Math. 132, 273–330; Shahidi 2010, AMS Colloq. 58, Ch. 5), the product is γ_v(s, π_v) for the GL(2) principal series π_v = π(|·|^θ, |·|^{−θ}) = Ind_B(|·|^θ ⊗ |·|^{−θ}), the local component of the isobaric sum |·|^θ ⊞ |·|^{−θ} whose standard L-function is F_θ. For any **unitary** π_v, conj γ(s, π_v) = γ(s̄, π̄_v) = γ(s̄, π̃_v) and γ(s, π_v)γ(1−s, π̃_v) = 1 (up to a unimodular ε), so |γ(½+it, π_v)| = 1. (b1) is this identity.

**Statement (b2).** For 0 < θ < ½, π(|·|^θ, |·|^{−θ}) is irreducible and unitarizable (complementary series) at every place v. The invariant Hermitian form is the **Knapp–Stein form** ⟨f, f⟩_θ := ⟨f, A(w, θ) f⟩_{L²(K)}, A(w, θ): I(θ) → I(−θ) the standard intertwining operator, and it is **positive definite**. For θ > ½ (θ ≠ ½ + ℤ) it is indefinite.
*Proof (literature, with the explicit check).* Knapp–Stein, "Intertwining operators for semisimple groups", Ann. of Math. 93 (1971) 489–578; Knapp, *Representation Theory of Semisimple Groups: An Overview Based on Examples* (Princeton 1986), Ch. VII (intertwining operators) and the SL(2, ℝ) complementary series (Ch. II); Bargmann 1947; p-adic: Shahidi 1990 (complementary series ending at the pole of the normalizing L-function). The argument: the normalized operator is unitary at θ = 0; the form stays nondegenerate while I(θ) is irreducible (0 ≤ θ < ½); reducibility first occurs at θ = ½ (trivial ⊕ Steinberg). Explicitly:
- ∞ (K = O(2), weights 2m): Knapp–Stein eigenvalue ratio c_{2m}(θ)/c_0(θ) = ∏_{j=1}^{m} (j−½−θ)/(j−½+θ) (standard, weight vectors normalized so c_0 > 0). All > 0 for 0 < θ < ½ (`p_checks.log` (P2), m ≤ 200 at θ = 0.49); for θ = ¾ the m = 1 ratio is −1/5 (indefinite).
- p (spherical vector): A(w, θ) f° = c_p(θ) f°, c_p(θ) = ζ_p(2θ)/ζ_p(1+2θ) = (1−p^{−1−2θ})/(1−p^{−2θ}) > 0 for θ > 0 (Gindikin–Karpelevich / Langlands LNM 544); `p_checks.log` (P2): c_7(0.2) = 1.7277.

**Consequence (b3): what filters F_θ.** For 0 < θ < ½ (the range where F_θ's mean field is ≥ 0, ARP G.3c), F_θ's local GL(2) data are unitary, the local reflection on them is unitary, and the Knapp–Stein form is positive definite. **No local GL(2) unitarity or Knapp–Stein-form test separates F_θ from ζ ⊞ ζ.** The filter acts on the **GL(1) isobaric constituents** |·|^{±θ}, which are non-unitary quasi-characters of C_ℚ; it is global, and it acts through the polarization lemma P.4. Reason: the Weil distribution is additive over isobaric constituents, W_{F_θ} = |·|^θ·W_ζ + |·|^{−θ}·W_ζ, because −(Λ₁Λ₂)′/(Λ₁Λ₂) = −Λ₁′/Λ₁ − Λ₂′/Λ₂ (**PROVED-LEAN** `neg_logDeriv_product_of_hasDerivAt`, LogDerivativeProduct.lean:54; functional equation of the product **PROVED-LEAN** `GppNoGo.ShiftedProductControl.F_one_sub`). The explicit formula sees the constituents, not the GL(2) representation.

**Correction of earlier wording.** `ARP_plan.md` (§1 last bullet, §3 item 2, §5 F_θ row) said in effect "GL(2) temperedness is the filter / local Θ_v not unitary for F_θ". Corrected:
- *unitary* ≠ *tempered*. F_θ's local GL(2) data are unitary (complementary series) but non-tempered.
- Local temperedness (Ramanujan, |α_p| = |β_p| = 1, i.e. |a(p^k)| ≤ k+1) does separate F_θ from G_θ′. But it is a bound on the **sizes** of the local coefficients, i.e. on weight magnitudes. It is fake-blind for DH, DMV and P_X (tempered-size weights; P_X has exactly ζ's), and it contributes nothing to ζ's positivity. So it is a filter, never a source (Golden Eye's audit stands).
- For an isobaric sum of idele-class quasi-characters, "tempered at one place" ⇔ "every constituent unitary on C_ℚ", because |·|^θ has the same exponent at every place. So local non-temperedness is the **shadow** of the global non-unitary constituent, and the global statement is the one that does the work (P.4).
- **R3 relevance.** R3 forbids sign-of-weight arguments: F_θ has positive weights 2cosh(θ log n)Λ(n) and fails; G_θ′ has signed weights 2cos(θ′ log n)Λ(n) and passes. P.3 + P.4 locate the distinction exactly: neither the sign of the weights nor local GL(2) unitarity, but whether the twisting quasi-character is unitary on C_ℚ (norm vs polarization). Any ARP mechanism that is an L²(C_ℚ) norm identity therefore satisfies R3 automatically (T3, ARP G.5).
- *Observation (PROVED that both thresholds are θ = ½; a structural identification is CONJECTURE):* F_θ's mean-field kernel e^{−(½−θ)|r|} + e^{−(½+θ)|r|} (G.3c) stops being positive definite exactly where the complementary series ends.

**R-line.** F_θ: **passes** P.3 for 0 < θ < ½ (this is the correction); for θ > ½ the Knapp–Stein form is indefinite, but a filter must work for all θ ≠ 0, so P.3 is not the filter. DH: a combination of two Dirichlet L-functions with complex coefficients, not the L-function of one representation — breaks. H₋ₜ, P_X, DMV, G₅: no automorphic representation — break. Epstein: a sum of L-functions of two dihedral GL(2) representations, not one — breaks. **G_θ′ passes** (tempered unitary principal series |·|^{iθ′} ⊞ |·|^{−iθ′}). R6: exact.

## P.4 Polarization lemma on L²(C_ℚ) — PROVED (elementary); finite parts PROVED-LEAN; Lean ask open (this is answer (a))

**Setting.** C_ℚ = 𝔸^*/ℚ^* ≅ Ẑ^* × ℝ₊^*, d^*x Haar (Connes §III (14)); H := L²(C_ℚ, d^*x) (the half-density space, Connes §III (12) with δ = 0). For L > 0, 𝒟_L := Bruhat–Schwartz functions on C_ℚ supported in {|log|x|| ≤ L}, and H_L its closure. On Ẑ^*-invariants, u = log|x| identifies 𝒟_L with C_c^∞([−L, L]). Involution ξ^*(x) := conj ξ(x^{−1}); convolution (ξ ⋆ η^*)(x) = ∫ ξ(y) conj η(x^{−1}y) d^*y (in log scale ∫ ξ(u) conj η(u−r) du).
**Hypotheses.**
- (H1) W is a Hermitian distribution on C_ℚ (W(h^*) = conj W(h)), continuous on 𝒟_{2L}. Define B_W(ξ, η) := W(ξ ⋆ η^*) on 𝒟_L. (For ζ: W = W_ζ, the Weil distribution in half-density, and B_W(f, f) = Q_L(f) of ARP §G.)
- (H2) ω: C_ℚ → ℂ^× is a continuous quasi-character, ω = χ|·|^σ with χ unitary and σ = Re ω ∈ ℝ (Tate 1950, §1.2 / §4.3). ω·W denotes h ↦ W(ωh), and ω^∨ := ω̄^{−1} = χ|·|^{−σ}.

**Lemma P.4.**
1. *(Exact twist identity.)* For all ξ, η ∈ 𝒟_L: (ω·W)(ξ ⋆ η^*) = B_W(M_ω ξ, M_{ω^∨} η).
2. *(Support.)* M_ω maps 𝒟_L bijectively onto 𝒟_L: multiplication by e^{±θu} preserves support in [−L, L] (no window inflation). On H_L, ‖M_ω‖ = ‖M_ω^{−1}‖ = e^{|σ|L}, and M_ω|_{H_L} is unitary ⇔ σ = 0 ⇔ ω^∨ = ω.
3. *(Unitary case: norm.)* If σ = 0, then (ω·W)(ξ ⋆ ξ^*) = B_W(M_ω ξ, M_ω ξ). Hence B_W ≥ 0 on 𝒟_L ⇔ B_{ω·W} ≥ 0 on 𝒟_L, for the **same** L. Symmetrized, with χ = |·|^{iθ′}, W real and even, and real ξ: (χW + χ̄W)(ξ ⋆ ξ^*) = 2B_W(ξ cos θ′u, ·) + 2B_W(ξ sin θ′u, ·), a sum of norms. This is G_θ′'s identity Q_G = 2Q_ζ(f cos) + 2Q_ζ(f sin).
4. *(Non-unitary case: polarization.)* If χ = 1 and σ = θ ≠ 0, then with h_± := |·|^{±θ}ξ:
   (ω·W + ω^{−1}·W)(ξ ⋆ ξ^*) = B_W(h₊, h₋) + B_W(h₋, h₊) = 2 Re B_W(h₊, h₋) = ½[B_W(h₊+h₋, ·) − B_W(h₊−h₋, ·)],
   and for real ξ and real even W this equals 2B_W(ξ cosh θu, ·) − 2B_W(ξ sinh θu, ·). This is F_θ's identity Q_F = 2Re B_ζ(e^{θu}f, e^{−θu}f), a cross-pairing of two different vectors. When B_W ≥ 0 the only bound is Cauchy–Schwarz, |Q_F| ≤ 2‖h₊‖_B‖h₋‖_B, which has no sign.
5. *(Sharpness: 2-atom witness.)* B_W ≥ 0 on 𝒟_L does **not** imply that the symmetrized non-unitary twist is ≥ 0 on 𝒟_L. Witness:
   - W = P(r) = 1 + cos r, which is positive type (it is the Fourier transform of δ₀ + ½(δ₁ + δ₋₁)).
   - Atoms at u = 0, ½ (inside [−L, L] for L ≥ ½), c = (1, −1), θ = 10.
   - Twisted Gram M_ij = 2cosh(θ(x_i−x_j)) P(x_i−x_j) = [[4, a], [a, 4]] with a = 2cosh 5 · (1 + cos ½) = 278.67.
   - Eigenvalues 4 ± a = 282.67 and −274.67; c^T M c = 8 − 2a = −549.34 < 0.
   - Unitary twist 2cos(θ r)P: eigenvalues 5.065 and 2.935 (≥ 0); untwisted: 3.878 and 0.122 (`p_checks.log` (P4)).
   - Smooth version: replace the atoms by bumps ξ_ε of masses ±1. By continuity of P the twisted form tends to c^T M c < 0, so some ξ ∈ 𝒟_L violates positivity.

**Proof of 1–4.**
- (1) ω is a homomorphism, so ω(x) = ω(y)·ω(x^{−1}y)^{−1}, and conj ω^∨(z) = ω(z)^{−1}. Hence ω(x) ξ(y) conj η(x^{−1}y) = (ωξ)(y) · conj((ω^∨η)(x^{−1}y)). Integrate in y, then apply W.
- (2) |ω(x)| = |x|^σ = e^{σu}, with sup and inf over |u| ≤ L equal to e^{±|σ|L}.
- (3) ω^∨ = ω when |ω| = 1. Bijectivity of M_ω on 𝒟_L gives ⇔. The cos/sin form: 2cos(θ′r) = e^{iθ′r} + e^{−iθ′r}, then expand B_W(ξe^{±iθ′u}, ·) into real and imaginary parts (S102 replay (1)).
- (4) Apply (1) with ω and ω^{−1} = ω^∨ (σ real, χ = 1). The polarization algebra is checked in `p_checks.log` (P3), and the integrand 2cosh(θr) = e^{θu}e^{−θ(u−r)} + e^{−θu}e^{θ(u−r)} in S102 replay (2). ∎

**Lean status.**
- PROVED-LEAN (finite parts):
  - `GppNoGhost.weilForm_nonneg` (FiniteNoGhost.lean:36): a positive-type kernel gives a nonnegative finite form.
  - `GppCosinePositiveType.cos_mul_positiveType` (CosinePositiveType.lean:65): cos(ωt) is positive type. With `positiveType_add` (FinitePrimePoissonPositiveType.lean:30) and `GppHaarPositivityWeil.const_one_positive_type` (HaarPositivityWeil.lean:163), this makes the witness kernel P = 1 + cos positive type.
  - `positiveType_weighted_gram` (SchurWeilClass.lean:46) is exactly the **norm** case: real weights d_i d_j are absorbed into the test vector.
  - The F_θ twist is the weighting d_i d_j^{−1} (d_i = e^{θx_i}), which cannot be absorbed. That is (4) in finite form.
- Open: `twist_cos_sum`, `twist_cosh_krein` and the 2-atom witness (our LEAN_ASK.md; ARP G.6 asks 1–2). Suggested addition: `twist_unitary_norm` (P·cos(θ′·) positive type from `hP n x (c·e^{∓iθ′x})`).

**What P.4 does not give.** Nothing about ζ itself, whose constituent (the trivial character) is unitary. Nor does it explain F_θ's onset L: that is set by ζ's soft-mode margin, λ_min = θ²κ (NUMERICAL, session 2 §5). For ζ's own F_θ, non-positivity on some window is PROVED independently by Weil's criterion, since F_θ has off-line zeros ρ ± θ.
**R-line.** F_θ is case (4): **breaks** here (polarization). **G_θ′ is case (3): passes** on the same window. DH, H₋ₜ, DMV, Epstein, P_X and G₅ are not twists of ζ, so the lemma is silent on them (it is an R3 tool). R4: (3) carries ζ's soft mode to its twist unchanged (M_ω is unitary). R6: exact identities; the witness is exact in closed form.

## P.5 The E-intertwining: exact citation and sign conventions — PROVED (this is answer (c); closes ARP_plan open item (a))

**Connes 1999.** Appendix I, **Lemma 2** (arXiv math/9811068v1, p. 56): for f ∈ 𝒮(𝔸)₀ the series E(f)(g) = |g|^{1/2} Σ_{q∈k^*} f(qg) converges absolutely, satisfies ∀n ∃c: |E(f)(g)| ≤ c e^{−n|log|g||}, and **E(f̂)(g) = E(f)(g^{−1})**. Proof via Poisson (App. I eq (36)): |x| Σ_{q∈k} f(xq) = Σ_{q∈k} f̂(x^{−1}q), i.e. (37) |x| Σ_{k^*} f(xq) = Σ_{k^*} f̂(x^{−1}q) + δ with δ = −|x|f(0) + ∫ f(y) dy, which on 𝒮(𝔸)₀ gives
  **(38)  E(f)(x) = E(f̂)(1/x),  f ∈ 𝒮(𝔸)₀.**
The same identity is quoted as §VIII eq (22). Here 𝒮(𝔸)₀ is §III eq (6) = App. I eq (31). The Fourier transform is App. I eq (33) with the character conventions (27)–(29) and dx normalized by (34). Equivariance is §III eq (19), E U(a) = |a|^{1/2} V(a) E; the isometry property is §III eq (18) with norms (7), (12) (δ-weighted).
**Meyer 2005.** §5.4 "The Poisson Summation Formula":
- Σ_{a∈K} f(ax) = |x|^{−1} Σ_{a∈K} Ff(ax^{−1}) for f ∈ 𝒮(𝔸_K), with ψ normalized (§5.1) so that (x·K)^⊥ = x^{−1}·K.
- R(f)(x) := f(0) − |x|^{−1}·Ff(0), and **eq (31): J Σ F = Σ + R**, with J f(x) = |x|^{−1}f(x^{−1}) (eq (12)).
- Lemma 5.5: the image of i₊ + i₋ is X = {(f₀, f₁) : f₀ − f₁ = c₀ − c₁|x|^{−1}}, whose two extra dimensions are the poles, spec π₊ = {1, |x|} (Intro).
- Mellin form, eq (33): Σf^(χ) = ΣFf^(|x|χ^{−1}) for χ ∉ {1, |x|} (Tate's FE).
- Symmetry: §5.6 (F on H₊ and J on H₋ agree on H₊ ∩ H₋).
**Reconciliation (exact).** Put f ∈ 𝒮(𝔸)₀, i.e. f(0) = 0 = ∫f dx = Ff(0), so R(f) = 0. Then (31) reads Σf(x) = |x|^{−1}ΣFf(x^{−1}). Multiplying by |x|^{1/2} gives E(f)(x) = |x|^{−1/2} ΣFf(x^{−1}) = E(Ff)(x^{−1}), which is Connes (38). So:
- *Signs.* Both texts use the character without a minus sign (Connes (33), Meyer (11)). The other sign convention gives F^* f = Ff(−·), and E(F^*f) = E(Ff) because −1 ∈ ℚ^* and Σ_{q∈ℚ^*} is invariant under q ↦ −q (Meyer: "since −1 ∈ K_S^×, F = F^*", §3.2). Equivalently F² = λ_{−1} and E∘λ_{−1} = E, so **E(F_𝔸φ) = E(φ)∘inv and E(φ) = E(F_𝔸φ)∘inv are the same statement.** ARP G.4b's form is correct as written.
- *|x|^{1/2}.* Connes puts it in E (half-density, unitary inversion g ↦ g^{−1} on L²(C, d^*g)). Meyer keeps it in the measure, L²(C_S)_{1/2} = L²(|x|d^×x), and in J (the |x|^{−1} in (12)). These are the same structure.
- *Haar.* Both use the self-dual dx with vol(𝔸/ℚ) = 1, so no constant appears in (38)/(31). The multiplicative normalizations (Connes (14), Meyer's Intro) agree up to the log-scale convention, and they only rescale d^*x.
- *𝒮(𝔸)₀.* Connes removes the two boundary values; Meyer keeps them as the rank-2 operator R, whose spectrum is {1, |x|} = the two poles. 𝒮(𝔸)₀ = ker R.
**Lean.** Not in GPPVerify: global Poisson, Tate's global FE and E (ARP §0.3). No PROVED-LEAN tag is possible here.
**R-line.** F_θ: M_ω acts on the C_ℚ side and is C_ℚ-covariant (V(a)M_ωV(a)^{−1} = ω(a)^{−1}M_ω), so F_θ passes the intertwining; it breaks at P.4. DH: there is no single ⊗′ Schwartz datum, so E of a combination of two isotypic pieces is not a norm — breaks. DMV, P_X: no lattice ℚ ⊂ 𝔸 and no Poisson (31) — break. Epstein: Poisson over ℤ² exists, but there is no ℚ^*-sum on 𝔸_ℚ — breaks. H₋ₜ, G₅: no adelic E — break. **G_θ′ passes.** R2: (38)/(31) is the one place where ⊗′ (Euler product, Meyer Lemma 5.3) and Poisson (FE) meet. R6: exact.

## P.6 The only non-isometric step is the global ℚ^*-sum — PROVED

**Statement.** In the chain 𝒮(𝔸) → (local Θ_v) → Σ / E → L²(C_ℚ):
- every local operator is unitary on the half-density space: Θ_v (P.2), and J and 𝔉 on L²(C_S)_{1/2} (Meyer Lemma 3.1, eq (12));
- the boundary map R has rank 2 (the poles);
- the global summation is the only non-unitary operator. On O^×_{∁S}-invariants, Σ ≅ **L_S := ∏_{v∈∁S}(1 − λ_{p_v}^{−1})^{−1}: T(C_S) → 𝒮(C_S)_>** (Meyer **Lemma 5.3**).
- ‖λ_{p_v}^{−1}‖_{L²(C_S)_α} = q_v^{−α}, so the product converges in operator norm on L²(C_S)_α **only for α > 1**. Meyer: "equivalent to the convergence of the Euler product defining the ζ-function in the region Re s > 1".
- On the critical-line space α = ½, L_S is the Mellin multiplier t ↦ ∏_{v∈∁S} L_v(½+it). This does not converge absolutely. By Kronecker–Weyl the partial products have operator norms ∏(1 − p^{−1/2})^{−1} → ∞, so there is no operator-norm limit. (Pointwise convergence on the line is itself a deep question: K. Conrad, Canad. J. Math. 57 (2005).)
- In Connes' Hilbert version, E is made an isometry by **defining** the norm (7) through the sum. The price is the global mismatch dx = lim_{ε→0} ε|x|^{1+ε}d^*x (§III eq (8)) and the δ-weight; Meyer §5.1 states that L²(𝔸_K, dx) and L²(𝔸_K^×, |x|d^×x) become different globally.
Either way, the non-unitary content sits entirely in the ℚ^*-sum, i.e. in the Euler product on the critical line.
**R-line.**
- P_X: L_S is a finite product, bounded at α = ½. It has no Poisson R-term and agrees with ζ on L ≤ (log X)/2 (R1, DMV theorem) — this is exactly why it is not a counterexample.
- F_θ: L_S = ∏(1 − p^θλ)^{−1}(1 − p^{−θ}λ)^{−1}, bounded only for α > 1 + θ.
- DH, H₋ₜ, G₅: no Euler operator. Epstein: a sum of two. DMV: Beurling operator without a lattice.
- **G_θ′ passes** (α > 1, as for ζ).
- R6: (31) is exact; the boundedness statement is a theorem, with no slack introduced.

## P.7 Link to the shared GAP G.4d — GAP (unchanged; equivalent to RH)

**What §P contributes to G.4d (exactly).**
1. The local reflections are isometries (P.2, P.3), even for F_θ's GL(2) data when 0 < θ < ½. So no local step loses or creates positivity.
2. The polarization lemma (P.4) explains F_θ completely at the level of mechanism: a non-unitary GL(1) constituent turns the norm into a polarization on the same window. It also shows that any L²(C_ℚ)-norm identity passes G_θ′ automatically (T3).
3. The only non-isometric step is the global ℚ^*-sum, i.e. the Euler operator L_S on the critical-line space (P.6), glued to Poisson by (38)/(31) (P.5).
4. **Nothing local bounds the GAP.** Fl_L ≥ −MF_L uniformly in L (G.4d) is a statement about L_S at α = ½ composed with the cutoff. It is equivalent to RH (Weil 1952; Bombieri 2000), and window L is the semi-local problem S_L = {∞} ∪ {p < e^{2L}} (PROVED-LEAN `GppWeilLadder.primeSide_eq_truncation`, WeilSupportLadder.lean:90, with Lean L′ = 2L).
**Base case.** For 2L < log 2, no prime enters (PROVED-LEAN `primeSide_eq_zero_of_support_lt_log_two`, line 106). There Q_L = pole + W_∞. Connes–Consani (Selecta 27:77) Thm 1 gives W_∞(g⋆g^*) ≥ Tr(ϑ(g)Sϑ(g)^*) ≥ 0 on the pole-free class ĝ(±i/2) = ĝ(0) = 0 (Cor. 2.3 gives the L − D decomposition; see P.8c). Full small-support positivity with the pole is as cited in ARP G.2. In Lean it is only conditional: `weil_nonneg_of_arch_nonneg_rung_zero` (named hypothesis 0 ≤ arch).

**Scope caveat (Golden Eye audit, 2026-10-28).**
- **PROVED:** (1) nothing local bounds G.4d (items 1–4 above); (2) §S S(c)4, *no MF-relative margin*: there is no c > 0 with Q_L ≥ c·MF_L for all L.
- **CONJECTURE (well supported, not proved):** that the closing step for G.4d *must be an exact identity*. S(c)4 does not exclude Q_L ≥ c·N_L for a norm N_L that degenerates on the dilated nodeless bumps (for example one weighting only frequencies near the zeros).
- So R6 is a working requirement at CONJECTURE level, not a theorem.
**R-line.** The GAP is R1 by construction. All fakes break before or at P.6, and G_θ′ passes every step (ARP G.4e, refined in P.12). **No new route is claimed.**

## P.8 Closing open item (b): the pole in the cutoff pair, the deficit form, and Connes–Consani 2021 Thm 1 (S103d)

Sources (verbatim numbering): Connes 1999 = arXiv math/9811068v1 §VII eqs (15)–(16), (31), (33), Thm 3; §VIII eqs (16)–(28), Cor. 2, Thm 5, Lemma 3; App. I eqs (36)–(38), Lemma 2; App. II eqs (5)–(6). Meyer 2005 = arXiv math/0311468v3 eq (12), §5.4 eq (31), Lemma 5.5. CC 2021 = Connes–Consani, Selecta 27:77, arXiv:2006.13771, Thm 1 (eq (4)), eq (5) (= Thm 6.11), Thm 3 (eqs (13)–(14) = Thm 4.7), eqs (8)–(11), Cor. 2.3 (eq (51)), eqs (48)–(49), Lemma 4.1 (eqs (60), (62)–(63)), Def. 4.4 (eq (72)). The 2 log Λ cancellation and Thm 5's positivity argument are written out in Queen's `team/queen/arp/sign_check.md` and cited on main (ARP G.4d, 7328e39). They are **not** repeated here.

### P.8a Queen's Q3: where the pole sits — PROVED (literature + elementary algebra)

**Cutoff pair.** P_Λ is multiplication by ρ_Λ = 1_{|x|≤Λ} (Connes §VII eq (15), "infrared cutoff"). P̂_Λ := F P_Λ F^{−1} is the "ultraviolet cutoff", and R_Λ := P̂_Λ P_Λ (eq (16)). Globally B_Λ := Im P_Λ ∩ Im P̂_Λ (§VIII, before Cor. 2), and B_{Λ,0} ⊂ B_Λ is the span of f ∈ 𝒮(𝔸)₀ with f, f̂ vanishing for |x| > Λ (§VIII, after (21)). Since F² = λ_{−1} and λ_{−1} commutes with P_Λ, **F P_Λ F^{−1} = P̂_Λ and F P̂_Λ F^{−1} = P_Λ: F swaps the pair.**

**The two boundary terms with their exact coefficients.** Write Ẽf(g) := |g|^{1/2} Σ_{q∈ℚ^*} f(qg) for any f ∈ 𝒮(𝔸), so Ẽ = E on 𝒮(𝔸)₀. Multiplying Connes App. I eq (37), |x| Σ_{k^*} f(xq) = Σ_{k^*} f̂(x^{−1}q) + δ with δ = −|x| f(0) + ∫ f dx, by |x|^{−1/2} gives the exact identity

  **(P8.1)  Ẽf(g) = Ẽf̂(g^{−1}) + f̂(0)|g|^{−1/2} − f(0)|g|^{1/2},  f ∈ 𝒮(𝔸), g ∈ C_ℚ.**

On 𝒮(𝔸)₀ this is (38) / §VIII (22). It is also Meyer's eq (31) in half density: JΣF = Σ + R with R(f)(x) = f(0) − |x|^{−1}Ff(0) (§5.4) and J from eq (12). Multiplying Σf = JΣFf − Rf by |x|^{1/2} reproduces (P8.1) term by term, including the sign.

**Mellin (Tate).** With d^*g normalized as in Connes §III eq (14), the boundary part of (P8.1) on |g| < 1 has Mellin transform ∫_{|g|<1} (f̂(0)|g|^{−1/2} − f(0)|g|^{1/2}) |g|^{s−1/2} d^*g = f̂(0)/(s−1) − f(0)/s. So:
- **Res_{s=1} = f̂(0) and Res_{s=0} = −f(0)** (Tate 1950 Main Thm; Golden Eye's lead is correct).
- F (f ↦ f̂, f̂ ↦ f(−·)) exchanges the two coefficients, and with them the vectors |g|^{−1/2} ↔ |g|^{1/2}. **F swaps the poles** (PROVED).

**Answer to Q3, in two senses, both exact.**
1. *As functionals on 𝒮(𝔸) — Golden Eye's lead holds.*
   - f ↦ f(0) is the P_Λ-side (position) functional: (P_Λ f)(0) = f(0) for every Λ > 0, since x = 0 lies in {|x| ≤ Λ}. It is the residue at s = 0.
   - f ↦ f̂(0) is the P̂_Λ-side (momentum) functional: (P̂_Λ f)^(0) = f̂(0). It is the residue at s = 1.
   - F swaps them, together with P_Λ ↔ P̂_Λ.
   - Neither functional depends on Λ.
2. *As vectors in L²(C_ℚ) after E (the sense relevant to S_Λ and the deficit) — both poles sit in the P̂_Λ half.* Take f ∈ B_Λ.
   - For |g| > Λ: Ẽf(g) = 0 for **every** f ∈ 𝒮(𝔸) vanishing on |x| > Λ, because |qg| = |g| (product formula). The P_Λ half of the inclusion E(B_{Λ,0}) ⊂ S_Λ needs no boundary condition.
   - For |g| < Λ^{−1}: Ẽf̂(g^{−1}) = 0 by the P̂_Λ condition, so (P8.1) gives **Ẽf(g) = f̂(0)|g|^{−1/2} − f(0)|g|^{1/2} exactly.** Both pole vectors appear together at the end that Connes controls with P̂_Λ through (22). This is exactly why Connes must restrict to 𝒮(𝔸)₀ (§VIII, after (21): "the equality (22) shows that E(f)(g) vanishes for |g| < Λ^{−1}").
   - In the conjugate coordinate Ẽ∘F the same statement holds with P_Λ and P̂_Λ exchanged.
   - So the honest answer is **P̂_Λ (in Connes' coordinate E), carrying both poles, not one each.** The "one each" picture is true only at the level of functionals (sense 1).
   - Meyer has no cutoff projections. Both poles sit in the single rank-2 map R (eq (31)): Lemma 5.5's two extra dimensions, spectrum {1, |x|}, i.e. |x|^{±1/2} in half density.

**Semi-local check (§N N.10, PROVED for pure tensors; checked against Connes §VIII).** With only the S_X-units O_S^* in place of ℚ^*, the Mellin transform of the sum is Z_S(f, s) = ∏_{v∈S} Z_v(f_v, s) (N.10(i)).
- It is holomorphic at s = 1: **the f̂(0) direction is absent.**
- −f(0) survives as a pole of order ≤ |S| at s = 0. It comes from each local Tate integral near x_v = 0.
- The s = 1 residue f̂(0) reappears only with the full ℤ[1/S]-sum, as ζ^{(S)}(s)·Z_S(f, s) (N.10(ii)). So it is carried entirely by the primes outside S.

This matches Connes:
- (a) His semi-local cutoff formulas, §VII Thm 4 and §VIII Cor. 2 (Q_Λ on L²(X_S)), have right-hand side 2h(1) log′Λ + Σ_{v∈S}∫′ with **no D term**.
- (b) D appears only globally, in the proof of Thm 5, as the B_Λ / B_{Λ,0} nuance in L²(X) with the global E of (18) and the global window S_Λ of (21).
- (c) The trivial-character case, the only one that carries the poles, is where Connes needs "more care since 1_R ∉ 𝒮₀(𝔸)" (p. 41). Here 1_R = ⊗_{v∉S} 1_{R_v} is exactly the tail outside S.
- It also matches P.6. The s = 1 pole of ζ is the pole of the cofinite Euler operator L_S (Meyer Lemma 5.3); finite products have none. This is why P_X is holomorphic at 1 and has no FE (N.10(i)).

**How this fits Golden Eye's one-pole-per-side picture.** Sense 1 survives, refined into an *asymmetry*:
- The **P_Λ side (position, f(0), s = 0) is local.** It is present at every finite S, degraded to order ≤ |S|.
- The **P̂_Λ side (momentum, f̂(0), s = 1) is global.** It exists only after Poisson (the ∫ f dx term of App. I (37)) and the infinite prime tail.
- **F swaps the two only after the full ℚ^*-sum**, through Tate's global FE. Semi-locally there is no FE (N.10(i)), so there is no swap.
- Sense 2 (both vectors at the |g| < Λ^{−1} end in the coordinate E) is the global statement. E needs ℚ^*.
- Golden Eye's picture "P_Λ carries s = 0, P̂_Λ carries s = 1, F swaps them" is therefore **correct globally as a statement about functionals**. It is not a statement about where the E-images sit, and it is not semi-local. R2 is visible here: the s = 1 pole needs the Euler tail and the FE together.

**Is it the same projection in OS4 / G.3a and in the deficit trace? Yes: the same pair B_Λ ⊃ B_{Λ,0}.** PROVED for function fields (Connes Thm 5 proof, p. 43). For ℚ it is only *as stated by Connes* (prolate B_Λ, §VIII end) and not re-derived.
- **The pole term.** The pole distribution is D(f) = ∫_{C_k} f(w)(|w|^{1/2} + |w|^{−1/2}) d^*w (App. II eq (6)), i.e. D(f) = f̂(i/2) + f̂(−i/2) with CC's Fourier convention, CC eq (11). This is exactly the OS4 / G.3a pole term ĝ(i/2) + ĝ(−i/2). Its two characters are the two vectors |g|^{±1/2} of (P8.1).
- **The deficit is built on B_{Λ,0}, not B_Λ.** Δ_Λ(f) = Tr((S_Λ − Q′_{Λ,0}) V(f)) (§VIII eq (24)), with Q′_{Λ,0} = E Q_{Λ,0} E^{−1}. Hypothesis a) of Thm 5 is stated for Q_Λ on B_Λ. Transport by (19) with f(g) = |g|^{−1/2}h(g^{−1}), insert Q_Λ, and use (26):
  Δ_Λ(f) = [Tr(S_Λ V(f)) − Tr(Q_Λ U(h))] + Tr((Q_Λ − Q_{Λ,0}) U(h)).
  - The first bracket tends to log|d^{−1}| f(1) − Σ_v D_v(f). This uses a), Tr(S_Λ V(f)) = 2f(1) log′Λ, and (26). The 2 log Λ cancels here.
  - Comparing with (17), Δ = log|d^{−1}|δ₁ + D − Σ_v D_v, **lim Tr((Q_Λ − Q_{Λ,0})U(h)) = D(f)**. This is Connes' "the term D … comes from the nuance between the subspaces B_Λ and B_{Λ,0}".
- **Rank and sector.** B_{Λ,0} is cut out of B_Λ by the two functionals f(0), f̂(0). These vanish off the trivial-character sector, so Q_Λ − Q_{Λ,0} has rank ≤ 2, supported in χ = 1 (finite-dimensional B_{Λ,χ}, Lemma 1; Connes p. 41: "the proof for χ = 1 … requires more care since 1_R ∉ 𝒮₀(𝔸)").
- **Consequence for OS4 / G.3a.** The two directions removed by 𝒮(𝔸)₀ are **exactly** the two pole directions that G.3a spends on the mean field (OS4). They enter Δ_Λ once, on the complement (deficit) side, consistent with "spent exactly once" (G.3a, §Q OS4).
- **No sign claim for D alone.** On a convolution square, D(g⋆g^*) = 2 Re(ĝ(i/2) ĝ(−i/2)‾) has signature (1,1), the cosh/sinh pair of §Q. Positivity needs the PNT mean (G.3a: pole − PNT mean = MF_L ≥ 0).
- **Old ARP_plan §4 CONJECTURE** ("soft mode = window image of (f(0), f̂(0))"):
  - The window image of the two boundary directions is now **identified, PROVED**: span{|g|^{1/2}, |g|^{−1/2}} restricted to the window, i.e. span{e^{u/2}, e^{−u/2}} = span{cosh(u/2), sinh(u/2)} in log scale. The even sector is Queen's pole vector 2cosh(u/2).
  - The alignment of ζ's window ground state with it remains **NUMERICAL** (Queen S1: 0.995–0.9997).
  - As a *route* it stays SUPERSEDED. G.3a already spends these directions, and T1 is withdrawn.

**R-line (P.8a).**
- F_θ: the two characters become |g|^{±1/2±θ}. There are four boundary vectors, and F still swaps them pairwise (G.3c) — passes here; breaks at P.4.
- DH: entire, so no boundary terms and no D. Q3 is vacuous; there is no E for a sum of two L-functions — breaks (P.5).
- DMV, P_X: no lattice, so (P8.1) does not exist — break.
- Epstein: (P8.1) holds per genus character, not for the sum — breaks.
- H₋ₜ, G₅: no E — break.
- **G_θ′ passes** (unitary characters |g|^{±1/2±iθ′}, a sum of two ζ pictures).
- R4: this is R4's location (the pole pair). R6: (P8.1) is exact.

### P.8b The ARP_plan §3 formula rewritten as a deficit — definition + PROVED sign; the limit is the GAP

The old formula, W(k⋆k̃) = Tr(Π_Λ ϑ(k)ϑ(k)^* Π_Λ) + δ_Λ(k) with a range-type positive term, is **withdrawn** (WRONG-SIGN, Queen sign_check.md items 1–2; ARP G.4d). Replacement, consistent with Queen §Q.7 and ARP G.4d:

  **(P8.2)  W(k⋆k^*) = Tr( V(k)^* (S_Λ − Q′_{Λ,0}) V(k) ) − δ_Λ(k),  δ_Λ(k) := Δ_Λ(k⋆k^*) − W(k⋆k^*).**

  *Corrected in P.9.0 (Queen §Q.8 (c2)).* For ℚ, read this as (P9.0), Tr((S_Λ − S_ΛQ′S_Λ)V) − δ^♮_Λ. The two forms agree exactly when Q′_{Λ,0} ≤ S_Λ, which holds for function fields (Connes (23)).

- **Terms.**
  - S_Λ is the window projection (§VIII eq (21)). Q′_{Λ,0} = E Q_{Λ,0} E^{−1} is the projection onto E(B_{Λ,0}) (eq (23)).
  - S_Λ − Q′_{Λ,0} is the projection onto **S_Λ ⊖ E(B_{Λ,0}), the complement of range(E) inside the window.**
  - V is the regular representation on L²(C_ℚ), eq (20). k is the test function on C_ℚ (log scale u = log|g|). The transport to h is Connes' f(g) = |g|^{−1/2}h(g^{−1}).
- **Labels.**
  - (i) Tr(V(k)^*(S_Λ − Q′_{Λ,0})V(k)) ≥ 0 for every Λ: PROVED (eqs (23)–(25); a projection inequality, "free").
  - (ii) The 2 log Λ term appears in both Tr(S_Λ V) and Tr(Q′V) and cancels; it is **not in δ_Λ** (Connes §VIII (16), proof of Thm 5 (17)–(26); Queen sign_check.md).
  - (iii) δ_Λ → 0 for each k ⇔ (16) ⇔ RH: PROVED for function fields (Thm 5), where RH is Weil's theorem, so δ_Λ → 0 there. **For ℚ it is the GAP** (≡ G.4d; RH-equivalent).
  - (iv) The limit residual lim δ_Λ = −Σ_{ρ off-line} m_ρ[∫ĝ dμ_ρ − ĝ(γ_ρ)] (μ_ρ as defined in ARP G.4d) is PROVED for function fields only (Lemma 3) and is a **CONJECTURE for ℚ.** In function fields the off-line set is empty.
  - No identification of (MF_L, Fl_L) with (Δ_Λ, δ_Λ) is claimed (Queen item 4).
- **Sign convention.** W = (trace over the complement, ≥ 0) − δ_Λ. Weil positivity follows from δ_Λ ≤ Tr(V^*(S_Λ − Q′)V), in particular from δ_Λ → 0.

### P.8c Reduction to Connes–Consani 2021 Thm 1 (Golden Eye's test 1) — sign PROVED to agree; 2 log Λ: test vacuous at Λ = 1

**CC Thm 1 in our notation.**
- g ∈ C_c^∞(ℝ_+^*), supp g ⊂ [2^{−1/2}, 2^{1/2}], i.e. f(u) := g(e^u) ∈ C_c^∞([−½ log 2, ½ log 2]). This is **team half-width L = ½ log 2**; the convolution g⋆g^* is supported in |u| ≤ log 2 = 2L.
- So no prime enters: p = 2 sits at the closed boundary, where g⋆g^* vanishes to all orders (CC eq (149) and our P.7 base case; Lean L′ = 2L = log 2).
- With ĝ(i/2) = 0 = ĝ(0) (Fourier convention ĝ(z) = ∫ g(ρ)ρ^{−iz} d^*ρ, CC eq (11)), CC eq (4) reads
  W_∞(g⋆g^*) ≥ Tr(ϑ(g) S ϑ(g)^*),
  where W_∞ := −W_ℝ (CC p. 3), ϑ(λ)ξ(v) = λ^{−1/2}ξ(λ^{−1}v) on L²(ℝ)_ev (eq (61)), and S = the projection onto **Sonin's space S(1,1) = {ξ : ξ(q) = 0 for |q| ≤ 1, (F_{e_ℝ}ξ)(p) = 0 for |p| ≤ 1}** (Def. 4.4, eq (72)). That is, **S = 1 − P₁ ∨ P̂₁ = (1 − P₁) ∧ (1 − P̂₁), the orthogonal complement of range(P₁ ∨ P̂₁)** (CC p. 6).
- The vanishing ĝ(±i/2) = 0 is exactly D(g) = 0 (App. II (6)). **CC's class is the archimedean analogue of 𝒮(𝔸)₀ / B_{Λ,0}: the two boundary (pole) directions of P.8a are removed.** ĝ(0) = 0 is an extra condition, "irrelevant for RH" (CC after Cor. 2).

**Exact correspondences (PROVED, each cited):**
1. *Cutoff pair.* CC Lemma 4.1 works with Connes' pair at Λ = 1: P = 1 − P₁ and P̂ = 1 − P̂₁ (eq (62)), with R₁ = P̂₁P₁ (Connes §VII (16)).
2. *Remainder = cutoff trace.* CC's trace-remainder δ(ρ) (eq (8)) satisfies δ(ρ) = Tr(ϑ(ρ^{−1}) P̂₁ P₁) for ρ ≥ 1 (eq (63)), with δ(ρ^{−1}) = δ(ρ). So D_CC(f) := ∫ f(ρ^{−1})δ(ρ) d^*ρ (eq (9)) equals **Tr(ϑ(f) R₁)** for supp f ⊂ (0, 1], and is the Λ = 1 cutoff trace of Connes §VII.
   - Shape check: Connes §VII eq (31), Tr(R_Λ T) = ∫_{|u|≤Λ²} ĝ(u)(2 log′Λ − log|u|) du, at Λ = 1 is ∫_{|u|≤1} ĝ(u)(−log|u|) du. That is the shape of CC eq (48), δ(ρ) = 4ρ^{1/2}∫₀¹ cos(2πρt) cos(2πt)(−log t) dt.
3. *Decomposition and sign.* CC Cor. 2.3 / eq (9) gives W_∞ = L − D_CC, with L(f) = Tr(ϑ(f) P P̂ P) = Tr(ϑ(f)(1−P₁)(1−P̂₁)(1−P₁)) ≥ 0 on f = g⋆g^* (Prop. 2.2(iii)). CC Thm 3 / eq (13) gives the refinement Tr(ϑ(f)S) = W_∞(f) + ∫ f ε d^*ρ. Hence
   **W_∞(g⋆g^*) = Tr(ϑ(g) S ϑ(g)^*) − δ^{CC}(g),  δ^{CC}(g) := ∫ (g⋆g^*) ε d^*ρ,**
   which is (P8.2) term for term: a trace over a *complement* (≥ 0, free) minus a remainder. **CC Thm 1 is exactly the statement δ^{CC}(g) ≤ 0 on its class.** The sign is the same as Connes Thm 5's deficit: PROVED, in two independent proved instances (function fields, Connes; archimedean small support, CC).
4. *2 log Λ.* At Λ = 1, 2 log′1 = ∫_{|λ|∈[1,1]} d^*λ = 0 on ℝ^*. The term is **identically absent**, so CC Thm 1 **cannot test** the 2 log Λ cancellation (Golden Eye's test 1; already on main, ARP G.4d).
   - It also cannot be reached by rescaling. ϑ(a)P_Λϑ(a)^{−1} = P_{aΛ} and ϑ(a)P̂_Λϑ(a)^{−1} = P̂_{Λ/a} (since F ϑ(a) = ϑ(a^{−1}) F), so (P_Λ, P̂_Λ) ≅ (P₁, P̂_{Λ²}). The genuine parameter is Λ², matching |u| ≤ Λ² in (31), and Λ = 1 is a single point of it.
   - The cancellation is verified only in Connes' global argument (cited in P.8b(ii)).

**What is NOT claimed (scope).**
- The two positive projectors are **different objects.** Connes' is S_Λ − Q′_{Λ,0}: the window minus the E-range of the doubly-cut space; global; any Λ. CC's is 1 − P₁ ∨ P̂₁: Sonin, the complement of the *join*; local at ∞; Λ = 1; no window, no E.
  - Locally at ∞, P₁ ∧ P̂₁ = 0 (no nonzero even L² function is both time- and band-limited), so the literal archimedean analogue of E(B_{Λ,0}) is {0}.
  - The global identity (22) has no local counterpart: F = I∘u_∞ (CC Lemma 1.3) carries the phase u_∞.
  - So "reduces exactly" holds at the level of the cutoff pair, the remainder, the sign and the (absent) 2 log Λ, **not** as an equality of projectors. A projector-level dictionary (Sonin ↔ S_Λ ⊖ E(B_{Λ,0})) is **open (GAP, minor)**.
- CC Thm 1 is positivity for W_∞ on the pole-free class only. Full small-support Weil positivity with the pole is as cited in ARP G.2, not re-checked here.

**R-line (P.8b–c).**
- F_θ: (P8.2) becomes Tr(V(k₊)^*(S_Λ − Q′)V(k₋)), a polarization (P.4) — breaks. Its archimedean factor is shifted, so CC Thm 1 does not apply verbatim.
- **G_θ′ passes**: two copies of (P8.2) with k cos, k sin (G.4e identity).
- DH, H₋ₜ, G₅: no E, so (P8.2) is undefined — break. DMV, P_X: no lattice — break. Epstein: sector-wise only — breaks.
- CC Thm 1 itself is fake-blind for every fake with ζ's γ_∞ (DH, P_X, DMV), as a PROVED local statement must be.
- R1: CC is a single window (L = ½ log 2); P.8b(iii) is the all-Λ statement (GAP).
- R4: CC's class removes the pole pair (P.8a), so Thm 1 is silent on R4.
- R6: the identities (P8.1)–(P8.2) and CC (13) are exact; Thm 1 is an inequality on one window.

**Verdict on "under review".**
- The **sign can drop "under review"**: the positive object is a trace over a complement, and the remainder enters with a minus sign. This is PROVED in Connes Thm 5 (function fields) and in CC Thm 3/Thm 1 (archimedean, Λ = 1). Meyer makes no positivity claim, so there is no conflict.
- Main already records "sign SETTLED per Golden Eye" (7328e39). Recommend that §S S(c)5 and ARP G.4d drop "under review / pending" **for the sign only.**
- Still open: G.4d itself (δ_Λ → 0 for ℚ, ⇔ RH); the ℚ residual formula (CONJECTURE); the projector dictionary (GAP, minor).

## P.9 The projector dictionary, why Λ = 1 is special, and the exact Λ ≥ 1 identity (S103e; analytic only)

Daniel Toupin / Golden Physics Project, with AI assistance from Grok. **RH is OPEN.** No numerics are used in this section.

**Track labels (Lead's split).**
- **[S&M]** = Spectrum & Measure: the Sonin/prolate side, and the cutoff Λ > 1 at 2L < log 2. This is RH-silent.
- **[Cel]** = Celestial §C.13: support past p = 2 at L = ½ log 2. This is the main line, with primes (R2).
- **[Q]** = Queen: M1, the prolate leak (1 − S_Λ)Q′_Λ. It is not duplicated here.

**Notation.** All of this lives in L²(ℝ)_ev ≅ L²(ℝ_+^*, d^*λ), via CC (17) w(ξ)(λ) = λ^{1/2}ξ(λ), unless stated otherwise.
- ϑ(λ)ξ(v) = λ^{−1/2}ξ(λ^{−1}v) (CC (61)) and ϑ(f) = ∫ f(ρ)ϑ(ρ) d^*ρ.
- F = F_{e_ℝ} (CC (24)). F² = 1 on L²(ℝ)_ev.
- P_a = 1_{|x|≤a}, P^c_a := 1 − P_a, and P̂_b = F P_b F^{−1} (Connes VII (15)), with P̂^c_b := 1 − P̂_b.
- CC's pair is P = P^c_1, P̂ = P̂^c_1 (CC (62)).
- Local windows: S_Λ := 1_{Λ^{−1} ≤ |x| ≤ Λ} (the local trace of Connes VIII (21)) and Ŝ_Λ := F^{−1}S_ΛF.
- R_Λ := P̂_ΛP_Λ (Connes VII (16)).
- L(f) := Tr(ϑ(f) P P̂ P) (CC (51), Prop 2.2(iii)).
- L_Λ(f) := Tr(ϑ(f) P^c_Λ P̂^c_Λ P^c_Λ).

### P.9.0 Correction to (P8.2) (Queen §Q.8 (c2); ARP G.4d 1999c2a) — PROVED [Q]
(P8.2) wrote the positive part as Tr(V^*(S_Λ − Q′_{Λ,0})V). That form silently uses Connes VIII (23), Q′_{Λ,0} ≤ S_Λ.
- (23) is exact for function fields.
- For ℚ it is **not** exact with the prolate B_Λ, because at ∞ no nonzero f has f and f̂ both vanishing beyond Λ (Connes §VIII, before (31); Queen (c2)).

The corrected form is:

  **(P9.0)  W(k⋆k^*) = Δ^♮_Λ(k⋆k^*) − δ^♮_Λ(k),  where Δ^♮_Λ(h) := Tr((S_Λ − S_ΛQ′_ΛS_Λ)V(h)) ≥ 0 on h = k⋆k^*  (Queen Lemma Q.8.4),**

with δ^♮_Λ := Δ^♮_Λ − W.
- **(P8.2) = (P9.0) exactly when S_ΛQ′S_Λ = Q′**, i.e. when E(B_{Λ,0}) ⊂ S_Λ. That is Connes (23), so it holds for function fields.
- For ℚ, Δ_Λ − Δ^♮_Λ = −Tr((Q′ − SQ′S)V(h)). This is controlled by the leak (1 − S_Λ)Q′_Λ only, which is Queen's M1 [Q].
- By Queen's Prop Q.8.3 (PROVED), the P_Λ and P̂_Λ placements give the same Δ^♮ and δ^♮ on real f. So §P.8a's pole analysis does not change δ^♮.
- The sign verdict of P.8c (a trace over a complement minus a remainder) is unchanged.
- From here on, every "(P8.2)" in §P is to be read as (P9.0).

### P.9.1 Question (1): the projector dictionary — PROVED; answer: they do NOT agree, at Λ = 1 or at any Λ [S&M]
Compare CC's Sonin projection S^{So}_1 (onto S(1,1), CC Def 4.4 (72), with S^{So}_1 = 1 − P₁∨P̂₁ = P∧P̂, CC (81)/(90)) with Connes' S_Λ ⊖ E(B_{Λ,0}) (VIII (21)–(24)). They differ in six exact ways.

**(D1) Opposite Halmos corners.** By the two-projection theorem (CC Lemma 4.2, (64)–(65)), L²(ℝ)_ev = (generic part ⊕_n E_n) ⊕ (four corners). The generic part is the same for (P₁, P̂₁) and (P, P̂), with the same angle operator (CC Prop 4.5 proof, the "Π_n" paragraph: U = 1 − 2P₁, V = 1 − 2P̂₁).
- Sonin is the corner **ker P₁ ∩ ker P̂₁** (outside both cutoffs).
- Connes' B_Λ is the corner **ran P_Λ ∩ ran P̂_Λ** (inside both). At ∞ this corner is **{0}** (Connes §VIII before (31)).
- So locally, Connes' cutoff space is replaced by the prolate span (Connes VIII (31)–(33); Queen Q.8.3 (P)). That span is generic-part, not corner, data.

**(D2) Setting.** S^{So} acts on the local L²(ℝ)_ev, with no ℚ^*-sum. Connes' S_Λ ⊖ E(B_{Λ,0}) acts on L²(C_ℚ) *after* E (VIII (18)).
- In the ζ-sector, E(φ⊗1_Ẑ)(λ) = λ^{1/2} Σ_{n≠0} φ(nλ). This is w plus the integer sum, which is the only non-isometric step (P.6, Meyer Lemma 5.3).
- No choice of Λ identifies an operator on one space with one on the other. Any dictionary can only be at the level of trace formulas (P.9.2). *(Grok audit: this sentence is a REMARK. It records that (D1)–(D3) leave no identification. No theorem excluding every possible operator-level identification is proved here.)*

**(D3) At Λ = 1, Connes' object is zero.** S_1 = 1_{|g|=1} is the indicator of a null set in C_ℚ = ℝ_+^* × Ẑ^*, so S_1 = 0 and S_1 ⊖ E(B_{1,0}) = {0}. Then Δ^♮_1 ≡ 0 and Tr(S_1V) = 2f(1) log 1 = 0.
- By contrast, S(1,1) ≠ {0}: it is the eigenvalue-1 eigenspace of P P̂ P (CC (81)) and is infinite-dimensional (Sonin).
- So **"agree at Λ = 1" is FALSE (PROVED).** CC's Λ = 1 is not a degenerate instance of Connes' window. It sits in the opposite corner (D1).

**(D4) Even/odd: consistent.** Since −1 ∈ ℚ^*, E only sees the even part at ∞ in the trivial character sector. CC work in L²(ℝ)_ev, using even prolates only (CC footnote 7: the even λ(n)² sum to 2.237…, and to 4 with the odd ones included). So CC's even sector matches the ζ-sector. The odd sector belongs to odd characters and is irrelevant here.

**(D5) Half-density: consistent.** CC's w (17) has the same |λ|^{1/2} as E (18). CC's ϑ (61) is the unitary normalization of Connes' U, as in Connes VIII (19), EU(a) = |a|^{1/2}V(a)E.

**(D6) The pole pair.**
- Connes removes it by B_{Λ,0} ⊂ 𝒮(𝔸)₀ (f(0) = f̂(0) = 0; P.8a).
- CC remove it through ĝ(±i/2) = 0, (11) and (54), implemented by Q = −(ρ∂_ρ)² + ¼ (55). Its local shadow is CC (93): the range of 𝒟 = D_u² + D_u satisfies f(0) = f̂(0) = 0, the archimedean trace of 𝒮(𝔸)₀.
- CC Thm 1 removes one **extra** direction, ĝ(0) = 0 (s = ½, the constant vector ξ₀ of Lemma 6.10). This is not a pole condition and has no counterpart in Connes' B_{Λ,0}.

**What does hold (PROVED).**
- **(i) Cutoff traces agree at Λ = 1.** CC (63), δ(ρ) = Tr(ϑ(ρ^{−1})P̂₁P₁), is literally Connes' R_1 trace (VII (16)), restricted to the even sector. CC (48), with weight (−log t) on [0, 1], is Connes VII (31) at Λ = 1, where the weight 2log′1 − log|u| = −log|u| on |u| ≤ 1.
- **(ii) Windows correspond.** u^*S_Λu = Ŝ_Λ, because F = I∘u (CC Lemma 1.3, (25)) and IS_ΛI = S_Λ. So the band window Ŝ_Λ is the local, ultraviolet image of Connes' window (21), and Tr(ϑ(f)Ŝ_Λ) = Tr(ϑ(f^♭)S_Λ) = 2f(1) log Λ (f^♭(ρ) = f(ρ^{−1})). This is the local analog of Connes VIII p. 43, Tr(S_ΛV(f)) = 2f(1) log′Λ.
  - Proof: ϑ(f) has kernel f(λ/μ) for d^*μ, so its diagonal is f(1), and ∫_{Λ^{−1}}^{Λ} d^*λ = 2 log Λ. Trace class: ϑ(f)χ with smooth compact χ ≡ 1 on the window has a rapidly decreasing kernel.
- **(iii) The Sonin family is indexed by ab only.** ϑ(λ)S(a, b) = S(λa, b/λ). So S(Λ, Λ) ≅ S(1, Λ²), and every ϑ(f)-trace depends only on c = ab.

**Conclusion (1).** The minor GAP is closed **negatively**. There is no projector-level identification, because of (D1)–(D3). The correct dictionary is at the trace-formula level, (i)–(ii), and P.9.2 makes it exact at every Λ ≥ 1.

### P.9.2 Question (2): why Λ = 1 is special, and the exact Λ ≥ 1 identity — PROVED [S&M]
**(a) Dilations — the prior is confirmed, with one correction.**
- Fϑ(λ) = ϑ(λ^{−1})F (direct computation; this is CC's ϑ(ρ) = F^{−1}ϑ(ρ^{−1})F in the proof of Prop 4.5(iv)).
- Hence ϑ(λ)P_aϑ(λ)^* = P_{λa} and ϑ(λ)P̂_bϑ(λ)^* = P̂_{b/λ}. The pair (P_a, P̂_b) moves along the orbit ab = const, and **c = ab = Λ² is the complete dilation invariant** (PROVED).
- **Correction to the prior:** F-exchange (F P_a F^{−1} = P̂_b and back) holds iff a = b. Each orbit contains exactly one such pair, so (P_Λ, P̂_Λ) *is* F-self-dual for every Λ.
- What is special at Λ = 1 is **I-self-duality**, I P I = 1 − P, for the inversion λ ↦ λ^{−1}. I is the non-scaling part of F = I∘u, and u commutes with all ϑ(λ). In general I P^c_a I = P_{1/a}. So the pair (P^c, P̂^c) is I-adapted iff c = 1.
- The obstruction is the ratio c = Λ², which no dilation changes.

**(b) The exact steps of CC that use Λ = 1.**
- **I-self-duality, c = 1:**
  - Prop 2.2(iii), (43): P̂ = u^*(1 − P)u, from (25) and I P I = 1 − P. This turns P P̂ P into −P(½u^*du)P and gives L = ∫ f(δ − τ) = W_∞ + D (Cor 2.3 (51)).
  - Lemma 4.1, (60): the same identity, (u^g)^*P u^g = 1 − P̂, which gives (63).
  - Through (85), the proof of Thm 4.7 inherits this.
- **F-exchange only (a = b), so valid at every Λ:**
  - Prop 2.2(ii), the symmetry δ(ρ) = δ(ρ^{−1}).
  - Prop 4.5 (69)–(80): prolates at bandwidth 2πΛ², with P̂_Λξ_n = λ_n η_n.
  - (77), (82)–(84): the ρ ≥ 1 support argument, with ϑ(ρ^{−1})ξ_n supported in |v| ≤ Λ/ρ.
  - (81)/(90): the eigenvalue-1 eigenspace of P^c_ΛP̂^c_ΛP^c_Λ is S(Λ, Λ).
- **The support condition supp g ⊂ [2^{−1/2}, 2^{1/2}]** (Golden Eye: the second, separate restriction):
  - Prop 5.5 (I ⊂ [−log 2, log 2], length ≤ log 2), K_I (Thm 3.6, §6), Lemma 6.10 (140) and Thm 6.11 (141). These use the support and the ε-kernel, not I-self-duality.
  - Cor 2 (6): Weil = W_∞, i.e. there are no primes because supp g⋆g^* ⊂ [½, 2].
- **So the self-duality step uses Λ only and never the support. The support step uses the support only.** These are independent.

**(c) The exact identity (R6).** Fix Λ ≥ 1. Two facts:
- Dilation by ϑ(Λ) carries the I-adapted pair to (P^c_Λ, P̂^c_{1/Λ}).
- The actual band projection differs from it by exactly the band window: **P̂^c_Λ = P̂^c_{1/Λ} − Ŝ_Λ** (since 1_{|ξ|>Λ} = 1_{|ξ|>1/Λ} − 1_{Λ^{−1}≤|ξ|≤Λ}).

Since Tr(ϑ(f)ϑ(Λ)Xϑ(Λ)^*) = Tr(ϑ(f)X) (ℝ_+^* is abelian):

  **(P9.1)  L_Λ(f) = L(f) − Tr(ϑ(f) P^c_Λ Ŝ_Λ P^c_Λ).**

Split Ŝ_ΛP^c_Λ = P^c_ΛŜ_ΛP^c_Λ + P_ΛŜ_ΛP^c_Λ, and evaluate Tr(ϑ(f)Ŝ_ΛP^c_Λ) with Connes' symbol calculus, VII (21)–(30). The product of a Fourier multiplier and a multiplication operator has symbol region {|x| > Λ, Λ^{−1} ≤ |ξ| ≤ Λ}, and its log-measure at u = xξ is log min(|u|, Λ²)·1_{|u|>1}. This gives

  **(P9.2)  L_Λ(f) = L(f) + ε_Λ(f) + χ_Λ(f),**

where
- ε_Λ(f) := −Tr(ϑ(f)Ŝ_ΛP^c_Λ) = −∫_{|u|>1} ĝ(u) log min(|u|, Λ²) du, with g from f by Connes VII (28) in the even-sector normalization of CC (42)/(48);
- χ_Λ(f) := Tr(ϑ(f)P_ΛŜ_ΛP^c_Λ) = Tr([P^c_Λ, ϑ(f)] P_ΛŜ_Λ).
- *Grok audit (uniform bound made explicit).*
  - Here ĝ(u) = G(u) := ∫ f(ρ) ρ^{1/2} e^{2πiu(1−ρ)} d^*ρ. For f ∈ C_c^∞(ℝ_+^*), G is the Fourier transform of a C_c^∞ function of ρ, so it is Schwartz.
  - Hence |ε_Λ(f)| ≤ ∫_{|u|>1} |G(u)| log|u| du < ∞, uniformly in Λ, and ε_Λ(f) → ε_∞(f) as Λ → ∞ by dominated convergence.
  - |χ_Λ(f)| ≤ ‖[P^c_Λ, ϑ(f)]‖_1 is uniform as well, since the commutator norm is dilation-invariant (c = ab). So the uniform bound is verified.
  - A parity factor of 2 (even-sector versus full line) may be absorbed in the normalization of CC (42)/(48). It does not affect boundedness.

Equivalently, using Connes VII (31) for R_Λ and R_1, and Tr(ϑ(f)(R_Λ − R_1)) = Tr(ϑ(f)Ŝ_ΛP_Λ) = 2f(1) log Λ + ε_Λ(f):

  **(P9.3)  L_Λ(f) = W_∞(f) + Tr(ϑ(f)R_Λ) − Tr(ϑ(f)Ŝ_Λ) + χ_Λ(f),  with Tr(ϑ(f)Ŝ_Λ) = 2f(1) log Λ exactly.**

At Λ = 1, Ŝ_1 = 0 and χ_1 = 0, and (P9.3) is CC (85)/(43). Running the proof of CC Thm 4.7 verbatim at bandwidth 2πΛ², with (85) replaced by (P9.3), gives the **Λ-form of CC Thm 4.7**:

  **(P9.4)  Tr(ϑ(f) S^{So}_Λ) = W_∞(f) + ∫ f(ρ^{−1}) ε^{(Λ)}(ρ) d^*ρ − 2f(1) log Λ + χ_Λ(f),**

where
- ε^{(Λ)}(ρ) = Σ_n λ_n(Λ) (1 − λ_n(Λ)²)^{−1/2} ⟨ξ_n^Λ | ϑ(ρ^{−1})ζ_n^Λ⟩ for ρ ≥ 1, extended by ρ ↦ ρ^{−1} (CC (84) at bandwidth 2πΛ²);
- S^{So}_Λ is the projection onto S(Λ, Λ).

At Λ = 1 this is CC (83) ((13) in Thm 3).

**Justification of the trace-class steps.**
- [P^c_Λ, ϑ(f)] is ϑ(Λ)-conjugate to [P, ϑ(f)], which is trace class (CC proof of Prop 2.2(iii)).
- ϑ(f)R_Λ is trace class (Connes VII Thm 3).
- ϑ(f)Ŝ_Λ is trace class by (ii) of P.9.1.
- For f = g⋆g^*, 0 ≤ P^cP̂^c_{Λ²}P^c ≤ P P̂ P, so 0 ≤ L_Λ(f) ≤ L(f) is finite.

Labels:
- (P9.1)–(P9.3): **PROVED**. These are elementary operator identities plus Connes VII (21)–(31) and CC (43)/(51), valid for every f ∈ C_c^∞(ℝ_+^*) and every Λ ≥ 1. **No support condition is used.**
- (P9.4): **PROVED**, given (P9.3), by CC's §4 argument at a = b = Λ *(Grok audit: proved by running CC §4 at bandwidth 2πΛ². This was checked by inspection, not replayed line by line; Golden Eye audit pending)*. Every step used is in the "F-exchange only" list of (b).

**(d) Answers to (2).**
1. **Where 2 log Λ lives.** It is exactly the window trace Tr(ϑ(f)Ŝ_Λ) = f(1)·log(Λ²), the log of the dilation invariant c. Ŝ_Λ is precisely the defect between the actual pair and the dilated I-self-dual pair.
   - Formally, Ŝ_Λ = ∫_{−log Λ}^{log Λ} (d/dt)P̂_{e^{−t}} dt, the dilation-generator commutator −i[𝔇, P̂] integrated over a log-interval of length log Λ². So the prior's "trace of the commutator with the dilation generator" is right per unit log-scale. This integral form is a REMARK (formal, since [𝔇, P̂] is a boundary distribution); the window form (P9.3) is the exact statement.
2. **The positive object does not pay 2 log Λ.** In (P9.3) the 2 log Λ inside Tr(ϑR_Λ) cancels exactly against the window. This is the local, archimedean mirror of Connes VIII (16) and the Thm 5 cancellation p. 43 (Queen sign_check.md), now an exact finite-Λ identity.
   - The true cost of moving from Λ = 1 to Λ is **ε_Λ + χ_Λ**, which is **bounded uniformly in Λ**: |ε_Λ(f)| ≤ ∫_{|u|>1} |ĝ| log|u| du, and |χ_Λ(f)| ≤ ‖[P, ϑ(f)]‖₁, a dilation-invariant bound (PROVED).
   - So "restores self-duality at the cost of exactly 2 log Λ" is **FALSE as stated**. Dilation cannot restore it, since c is invariant. The band window restores it at a cost of exactly 2f(1) log Λ in the *cutoff* trace, while the *positive* trace changes only by the bounded ε_Λ + χ_Λ.
   - Monotonicity: **0 ≤ L_Λ ≤ L on f = g⋆g^*** (PROVED). The inequality is the operator ordering P̂^c_{Λ²} ≤ P̂^c_1; (P9.2) is its exact form.
3. **Link to (1).** Ŝ_Λ is the u-image of Connes' window (P.9.1 (ii)). χ_Λ is the archimedean (pre-E) cross term between the inside and outside of P_Λ across the band window.
   - It is the local analog of a leak term. Whether it matches Queen's M1 leak (1 − S_Λ)Q′_Λ after E is **not claimed** [Q].
4. **Λ versus support.** Λ is an archimedean phase-space cutoff. Changing Λ conjugates by ϑ(Λ), which commutes with ϑ(f), so it **never moves the support of f** and **brings in no primes** (PROVED: the trace identity above). The support of g is a separate parameter.

**(e) Consequence for R2 (stated plainly).**
- **Every identity in P.9.2 is archimedean, local and prime-free.** It is RH-silent and cannot engage R2.
- Pushing CC from Λ = 1 to Λ > 1 at 2L < log 2 [S&M] changes only the Sonin/prolate bookkeeping (P9.4). A Λ > 1 analog of CC Thm 1 would need a new compact-operator bound on the ε^{(Λ)} ∘ Q functional on the support class (the analog of K_I, Lemma 6.10). That is a **GAP [S&M]**, and RH-silent either way. *(Superseded by P.10.3: closed NEGATIVELY. Lemma 6.10_Λ / Thm 6.11_Λ hold uniformly given CERTIFIED CC Lemma 6.10, Λ = 1 is optimal in the Sonin family, and no K_I^{(Λ)} bound is needed. Grok audit S103f.)*
- **Only pushing the support past p = 2, at L = ½ log 2 [Cel], brings in primes, and with them R2 and the ℚ^*-sum (P.6).** The Λ-identity does not transfer to that step. It supplies only the archimedean Λ-bookkeeping for whatever [Cel] builds.

**(f) Queen's Conj Q.8.E (an explicit D_L ≥ 0 with Q_L = Tr(D_L VV^*)).**
- **Candidate from this section:** the Sonin projection, D^{loc}_Λ := S^{So}_Λ. By (P9.4), W_∞ = Tr(ϑ S^{So}_Λ) − [∫ f ε^{(Λ)} − 2f(1) log Λ + χ_Λ]. The bracket is not zero, so this is not an exact identity of the Q.8.E type.
- **Obstruction (R5; PROVED as stated):** suppose the rule L ↦ D_L uses only the archimedean data (P_Λ, P̂_Λ, S_Λ, Ŝ_Λ, ϑ), with no arithmetic input.
  - Then Tr(D_L VV^*) is the same for ζ and for any fake with ζ's gamma factor.
  - For 2L > log 2 the Q_L of ζ and of such a fake differ by the prime side (PROVED-LEAN `GppWeilLadder.primeSide_eq_truncation`). So (I1) can hold for at most one of them, and the rule cannot be what singles out ζ.
  - *Grok audit (wording fix).*
    - The fixed inequality is 2L **>** log 2: at 2L = log 2, compact support inside (−L, L) kills the n = 2 atom.
    - The witness fake should be read as a **datum**, not an L-function: e.g. P_X with X < 2, whose Weil form on the window is W_∞ alone. DH does not qualify, because its gamma factor has conductor 5.
    - No Dirichlet series with ζ's gamma factor, functional equation and finitely many poles is a fake: Hamburger's theorem makes it cζ.
    - The conclusion stands, with "archimedean data" read as excluding the atoms ϑ(n) at prime powers.
  - So every admissible D_L must contain E, through Q′_Λ. This is the operator form of P.7's "nothing local bounds G.4d" and of Queen's (N3).
- What P.9 does give Q.8.E is a constraint. In any D_L built from S_Λ, P_Λ, P̂_Λ, the 2 log Λ must be spent exactly against the window (P9.3), and an archimedean cross term of the χ_Λ type must be accounted for. CONJECTURE (structural; consistent with (N1)–(N3)).

**R-line (P.9).**
- R1: every Λ ≥ 1 and every f; nothing is uniform in L.
- R2: **not engaged**, since everything is archimedean and prime-free (e). R2 lives on the [Cel] support step only.
- R3: the positivity is the operator statement 0 ≤ P^cP̂^cP^c, not a weight sign.
- R4: n/a. The pole pair is removed as in (D6), and the extra ĝ(0) direction is CC's, not OS4.
- R5: **G_θ′ passes every step**, as do all fakes sharing W_∞, so the section is fake-blind by construction and cannot be a closing step.
- R6: (P9.0)–(P9.4) are exact identities; the only inequality is the monotonicity 0 ≤ L_Λ ≤ L.

**Status table (P.9).**

| item | label | track |
|---|---|---|
| (P9.0) correction S − SQ′S; agrees with (P8.2) iff E(B_{Λ,0}) ⊂ S_Λ (function fields) | PROVED (Queen Q.8.4, (c2)) | [Q] |
| (D1)–(D3) the projectors do not agree; at Λ = 1 Connes' object is 0 | PROVED | [S&M] |
| (D4)–(D5) even sector / half-density consistent; (D6) extra ĝ(0) direction in CC | PROVED | [S&M] |
| P.9.1 (i)–(iii) trace-level dictionary | PROVED | [S&M] |
| c = ab is the complete dilation invariant; F-exchange holds iff a = b; I-self-duality iff c = 1 | PROVED | [S&M] |
| CC steps using c = 1: Prop 2.2(iii)/(43), Lemma 4.1/(60); support only in §5–6 and Cor 2 | PROVED (by inspection) | [S&M]/[Cel] |
| (P9.1)–(P9.3) exact Λ ≥ 1 identity, 2 log Λ = window trace, cancels exactly | PROVED | [S&M] |
| (P9.4) CC Thm 4.7 at every Λ ≥ 1 | PROVED | [S&M] |
| ε_Λ, χ_Λ bounded uniformly in Λ; 0 ≤ L_Λ ≤ L | PROVED | [S&M] |
| Λ > 1 analog of CC Thm 1 (a K_I^{(Λ)} bound) | GAP (RH-silent) | [S&M] |
| χ_Λ ↔ M1 leak after E | not claimed | [Q] |
| an arithmetic-free D_L cannot single out ζ (so D_L must contain E) | PROVED (R5, as stated) | [Q] |
| Λ-extension engages R2 | NO (PROVED: Λ never moves the support) | [Cel] |

## P.10 The Λ > 1 analogue of CC 2021 Lemma 6.10 / Thm 6.11 (S103f; analytic only; [S&M] Λ-track)
Daniel Toupin / Golden Physics Project, with AI assistance from Grok. **RH is OPEN.**
- **Scope.** The support is kept fixed at CC's supp g ⊂ [2^{−1/2}, 2^{1/2}] (I = [−½ log 2, ½ log 2] in additive form, supp g⋆g^* ⊂ [½, 2]). So no primes enter (CC Cor 2), and every statement here is **RH-silent**.
- **Track.** Spectrum & Measure owns the Λ-track; this Lemma 6.10 analogue is ours by agreement with them.
- **Prolate inputs.** These are cited from S&M's `team/spectrum-measure/PROLATE_BOUNDS.md` (95ea14a, credit Spectrum & Measure), in their normalization: prolate bandwidth c_pr = 2πΛ², and dilation invariant ab = Λ².
- **Named inputs still pending from S&M** (the trace-norm and HS bounds for [P^c_Λ, ϑ(f)] and P_ΛŜ_Λ) are labelled **[S&M-pending]**. Nothing below depends on them except where it says so.

### P.10.0 Addendum and erratum to §P.9.2 (Golden Eye audit of main 25f8f83: (P9.1)–(P9.3) PASS line by line; credit Golden Eye)
1. **(P9.4) relabelled:** **PROVED given (P9.3) and CC §4 at a = b = Λ. The internals (69)–(84) were checked by inspection, not replayed.**
2. **Dilation convention, pinned once against CC (61):** (ϑ(λ)ξ)(v) = λ^{−1/2} ξ(λ^{−1}v), i.e. **ϑ(λ) dilates by λ** (support ×λ). This is the convention of all of §P.9, with ϑ(f) = ∫ f(λ)ϑ(λ) d^*λ.
   - **Kernel.** In this convention ϑ(f) has kernel k(y, x) = f(y/x)(x/y)^{1/2}/x on ℝ_+ (even extension). Put y = ρx; then the Connes VII (21)-type symbol is σ(x, ξ) = ∫ f(ρ)ρ^{1/2} e^{2πi xξ(ρ−1)} d^*ρ = G(−xξ), with G(u) = ∫ f(ρ)ρ^{1/2}e^{2πiu(1−ρ)} d^*ρ as in the §P.9 audit note.
   - The sign of u is immaterial in ε_Λ, because the region {|u| > 1} and the weight log min(|u|, Λ²) are even. The f(ρ^{−1}) in (P9.4) is CC's (83) with this same ϑ.
   - **Opposite convention.** ϑ(λ) ↦ ϑ(λ^{−1}) replaces f by f^∨(ρ) = f(ρ^{−1}). For f = f^* (the positivity case, f^*(ρ) = conj f(ρ^{−1})), f^∨ = f̄. Then G_{f̄}(u) = conj G_f(−u), and since the weight is real and even, ε_Λ(f̄) = conj ε_Λ(f). So for complex f the two conventions swap ε and conj ε (Golden Eye's point), and they agree for real-valued f = f^*. **PROVED.**
3. **Simpler trace route for Tr(ϑ(f)Ŝ_Λ). PROVED.**
   - Ŝ_Λ = F^{−1}S_ΛF = u^* I S_Λ I u = u^*S_Λu, using F = I∘u (CC Lemma 1.3, (25)) and I S_Λ I = S_Λ (the window is inversion-symmetric).
   - Since [u, ϑ(λ)] = 0, Tr(ϑ(f)Ŝ_Λ) = Tr(u ϑ(f) u^* S_Λ) = Tr(ϑ(f)S_Λ) = f(1)·∫_{Λ^{−1}}^{Λ} d^*x = **2f(1) log Λ**. (Trace class: write f = f₁⋆f₂ by Dixmier–Malliavin; S'ϑ(f₁) and ϑ(f₂)S_Λ are HS for a window S' ⊃ supp f₁·[Λ^{−1}, Λ].) This is the same value as the symbol route in (P9.3).
   - Recalled: for f = g⋆g^*, 0 ≤ L_Λ(f) ≤ L(f).
4. **Sanity constraint (analytic). PROVED, given (P9.4).**
   - Set ϑ(Λ^{−1}). Then S(Λ, Λ) ≅ S(1, Λ²) ⊂ S(1, 1): ξ vanishes on |x| ≤ Λ with ξ̂ vanishing on |p| ≤ Λ ⇒ ϑ(Λ^{−1})ξ vanishes on |x| ≤ 1 with transform vanishing on |p| ≤ Λ².
   - So for f = g⋆g^*, the LHS of (P9.4), Tr(ϑ(g)S^{So}_Λϑ(g)^*), is ≥ 0 and nonincreasing in Λ. This forces
     **∫ f ε^{(Λ)} + χ_Λ(f) ≥ 2f(1) log Λ − W_∞(f)** and **∫ f ε^{(Λ)} + χ_Λ(f) − 2f(1) log Λ ≤ ∫ f ε^{(1)}**.
   - With |χ_Λ(f)| ≤ ‖[P^c_1, ϑ(f)]‖_1 uniform in Λ (§P.9 audit note), this gives **∫ f ε^{(Λ)} = 2f(1) log Λ + O_f(1), two-sided.** It holds for every f ∈ C_c^∞(ℝ_+^*), by polarization over finite sums of g⋆h^* (Dixmier–Malliavin).
   - So ε^{(Λ)} carries mass log Λ² at ρ = 1 (+ bounded). That is the "growth like log Λ²".
   - **Consistency with P.10.3:** the 6.10-type bound of P.10.3 is for the total remainder, where this growth cancels exactly against −2f(1) log Λ. In P.10.4 the same growth reappears as the explicit leading term 2 log Λ·⟨ξ|Q_+ξ⟩, with the same sign. **No inconsistency, so no sign error.**
   - **Plunge-region mechanism. CONJECTURE (heuristic, not used).** The log Λ² mass is carried by the plunge-region modes. Their number is ≤ (2/π²) log(200Λ² + 25)·log(5/(ε(1 − ε))) + 7 (KRD Thm 3, as cited in PROLATE_BOUNDS §5), i.e. ∝ log c with the Landau–Widom constant (2/π²) log c·log(1/ε). This matches log Λ² in order.
   - **R1 is not claimed for the ε^{(Λ)} series.** It converges for each Λ (CC App. F at Λ = 1; the Λ > 1 analogue is S&M-pending), but not uniformly. P.10.3 does not use the series.

### P.10.1 What CC actually prove and use (citations verified in arXiv:2006.13771)
- **Thm 6.11 (141).** For g ∈ C_c^∞([2^{−1/2}, 2^{1/2}]) with ĝ(−i/2) = 0: W_∞(g⋆g^*) ≥ Tr(ϑ(g)Sϑ(g)^*) − c|ĝ(0)|², with c = 4γ/log 2.
- **Its proof uses exactly four inputs:**
  - (i) Thm 4.7 in the form (142): Tr(ϑ(f)S) = W_∞(f) + E(f), E(f) = ∫ f(ρ^{−1})ε(ρ) d^*ρ (91).
  - (ii) The factorization (143): k := Y⋆g has support in [2^{−1/2}, 2^{1/2}] (by ĝ(−i/2) = 0), Q(k⋆k^*) = g⋆g^* with Q = −(ρ∂_ρ)² + ¼ (55), and k̂(0) = −2ĝ(0).
  - (iii) Prop 5.5 (103): E∘Q(k⋆k^*) = ⟨ξ|N_Iξ⟩ for ξ = k∘exp ∈ H = L²(I).
  - (iv) Lemma 6.10 (140): ⟨ξ|N_Iξ⟩ ≤ γ|⟨ξ₀|ξ⟩|², γ ≈ 2.94355, ⟨ξ₀|ξ⟩ = (log 2)^{−1/2}ξ̂(0).
- **Lemma 6.10's proof is computer-assisted.**
  - Lemma 6.9 is the abstract rank-one restoration (138)–(139).
  - It is applied to the finite-rank approximant T of K_I (§6.4, Lemma 6.3/6.8: λ_max = 1.05158 > 1, λ₂ < 0.773), with ‖K_I − T‖ ≤ ε₁ ≈ 0.00122 (134).
  - Then (140) follows from N_I = −½ε′(1⁺)(Id − K_I) (Prop 5.5(ii), (104); ε′(1⁺) from Lemma 5.4).
  - **Label of CC Lemma 6.10 / Thm 6.11 for us: CERTIFIED** (published, with computer-assisted error control (134); not replayed here).
  - K_I itself is not positive or negative: λ₁(K_I) > 1 (Remark 6.12), so the rank-one term is necessary.
- **Inputs (ii) and (iii) are Λ-independent** (support only). Input (i) has the Λ-form (P9.4). **The only Λ-dependent input is (iv).**

### P.10.2 Definitions at Λ
- The Λ-remainder, defined series-free: **E^{tot}_Λ(f) := Tr(ϑ(f)S^{So}_Λ) − W_∞(f).**
  - By (P9.4), E^{tot}_Λ(f) = ∫ f(ρ^{−1})ε^{(Λ)}(ρ) d^*ρ − 2f(1) log Λ + χ_Λ(f).
  - At Λ = 1, E^{tot}_1 = E (CC (142)).
- The Λ-form of the CC operator: **⟨ξ|N^{(Λ)}_Iξ⟩ := E^{tot}_Λ(Q(k⋆k^*)), with ξ = k∘exp ∈ C_c^∞(int I).** At Λ = 1 this is (103).
- E^{tot}_Λ, not ∫ f ε^{(Λ)} alone, is the right Λ-analogue of CC's E. It is what enters W_∞ = Tr(ϑ(f)S^{So}_Λ) − E^{tot}_Λ, which is the role (142) plays in the proof of Thm 6.11. P.10.4 treats the naive ε^{(Λ)}∘Q.

### P.10.3 Main result: the Λ-analogue holds uniformly, with the same constant. PROVED, given CC Lemma 6.10 (CERTIFIED). [S&M]
**Monotonicity lemma. PROVED (elementary).**
- Let h ∈ C_c^∞(ℝ_+^*). Then Λ ↦ Tr(ϑ(h)S^{So}_Λϑ(h)^*) is ≥ 0 and nonincreasing on [1, ∞), and → 0 as Λ → ∞.
- **Proof.**
  - ϑ(h) commutes with ϑ(Λ), and S^{So}_Λ = ϑ(Λ)S^{So}_{(1,Λ²)}ϑ(Λ)^*. So the trace equals Tr(ϑ(h)S^{So}_{(1,Λ²)}ϑ(h)^*), with S(1, Λ²) decreasing in Λ (P.10.0 item 4).
  - The operators are positive and trace class, dominated by ϑ(h)S^{So}_1ϑ(h)^*, which is trace class by CC Thm 4.7.
  - ∩_Λ S(1, Λ²) = {0}, since ξ̂ ≡ 0. So S^{So}_{(1,Λ²)} → 0 strongly, and normality of the trace gives → 0.
- **Positivity class.** Q(k⋆k^*) = h⋆h^* with h = (½ + ρ∂_ρ)k (this is the factorization Q = (½ − ρ∂_ρ)(½ + ρ∂_ρ) used in CC (143) and Remark 6.12). Hence:

  **(P10.1)  E^{tot}_Λ(Q(k⋆k^*)) ≤ E^{tot}_1(Q(k⋆k^*)) = E∘Q(k⋆k^*), for every Λ ≥ 1** (W_∞ is Λ-free).

**Lemma 6.10_Λ.** For I = [−½ log 2, ½ log 2], every Λ ≥ 1 and ξ ∈ C_c^∞(int I):

  **(P10.2)  ⟨ξ|N^{(Λ)}_Iξ⟩ ≤ ⟨ξ|N_Iξ⟩ ≤ γ|⟨ξ₀|ξ⟩|²,  γ ≈ 2.94355 (CC (140)).**

**Thm 6.11_Λ.** For g as in CC Thm 6.11 and every Λ ≥ 1:

  **(P10.3)  W_∞(g⋆g^*) ≥ Tr(ϑ(g)S^{So}_Λϑ(g)^*) − c_Λ|ĝ(0)|²,  c_Λ ≤ c = 4γ/log 2.**

**Explicit Λ-dependence.**
- The optimal γ_Λ and c_Λ (sup of E^{tot}_Λ over the class) are **nonincreasing in Λ**, bounded by CC's constants, and have a limit c_∞ = sup over the class of −W_∞(g⋆g^*)/|ĝ(0)|².
  *(Grok audit: as proved, lim c_Λ exists (monotone) and lim c_Λ ≥ c_∞ := sup(−W_∞/|ĝ(0)|²), since E^{tot}_Λ ≥ −W_∞. Equality needs the sup and the Λ-limit to be interchanged, which is not shown: **GAP-minor**, not used.)*
- The Λ → ∞ limit is exact:

  **(P10.4)  lim_{Λ→∞} [∫ f ε^{(Λ)} − 2f(1) log Λ + χ_Λ(f)] = −W_∞(f)** (from E^{tot}_Λ(f) → −W_∞(f) by the monotonicity lemma; for every f by polarization).

- Also lim L_Λ(f) = 0 and lim Σ_n λ_n(Λ)²⟨ζ^Λ_n|ϑ(f)ζ^Λ_n⟩ = 0. These are positive and nonincreasing, since P̂^c_{Λ²} ↓ 0 after dilation. Together with (P9.2) this gives **L(f) = −ε_∞(f) − lim χ_Λ(f)**, an exact identity. **PROVED.**

**What this does and does not give.**
- (P10.2) and (P10.3) are **uniform in Λ** (R1 in Λ).
- But (P10.3) is **weaker than CC Thm 6.11**, since Tr(ϑ(g)S^{So}_Λϑ(g)^*) ≤ Tr(ϑ(g)S^{So}_1ϑ(g)^*). **Among the Sonin family, Λ = 1 is optimal.**
- The Λ-track adds no positivity beyond Λ = 1, and its Λ → ∞ end, W_∞(g⋆g^*) ≥ −c_∞|ĝ(0)|², is implied by CC Thm 6.11.
- **No K_I^{(Λ)}-type input is needed**, beyond CC's Λ = 1 computation. The last [S&M] GAP on the Λ-track is therefore **closed (PROVED given CERTIFIED CC Lemma 6.10).** It is closed in this sense, with no new information on RH.
- *Grok audit (S103f): SOUND; wording.*
  - The Λ-track is **closed NEGATIVELY**. Thm 6.11_Λ is an immediate corollary of CC Thm 6.11 plus the monotonicity lemma, so it carries no positivity beyond CC.
  - The only route left past CC is the support extension, which is where primes enter (R2; Celestial §C.13, Prime §N.12).
  - Checks:
    - N^{(Λ)} ≤ N_I is a **quadratic-form (operator) inequality**, not a domination of kernels. It comes from the projection inequality S^{So}_{(1,Λ²)} ≤ S^{So}_1 (S(1, Λ²) ⊂ S(1, 1)), conjugated by ϑ(h), which preserves ≤ (X P₁ X^* ≤ X P₂ X^* for P₁ ≤ P₂). W_∞ is Λ-free and cancels.
    - Unitary dilation invariance of the trace uses only [ϑ(h), ϑ(Λ)] = 0.
    - CC numbers checked against arXiv:2006.13771: γ ≈ 2.94355 (140); λ_max = 1.05158; ε₁ ≈ 0.00122 (134); c > 0.227784 (Lemma 6.8, ⇒ λ₂ < 0.773); N_I = −½ε′(1⁺)(Id − K_I); c = 4γ/log 2 (141).
    - CC Lemma 6.10 is labelled CERTIFIED (computer-assisted), not PROVED. ✓

### P.10.4 Where the naive analogue fails: ε^{(Λ)}∘Q alone. PROVED, given (P9.4). [S&M]
- **Notation.** Write ⟨ξ|Q_+ξ⟩ = ‖ξ′‖² + ¼‖ξ‖² (Dirichlet form of Q_+ = −∂_x² + ¼ on I, CC proof of Thm 3.6). Then Q(k⋆k^*)(1) = ⟨ξ|Q_+ξ⟩, since (ξ⋆ξ^*)″(0) = −‖ξ′‖² in the additive variable.
- **Exact identity.** From (P9.4) and P.10.2:

  **(P10.5)  ∫ Q(k⋆k^*)·ε^{(Λ)} = ⟨ξ|N^{(Λ)}_Iξ⟩ + 2 log Λ·⟨ξ|Q_+ξ⟩ − χ_Λ(Q(k⋆k^*)).**

- **Two-sided growth.** The lower bound uses P.10.0 item 4. So with f_ξ := Q(k⋆k^*) and sup_Λ |χ_Λ(f_ξ)| ≤ ‖[P^c_1, ϑ(f_ξ)]‖_1 < ∞ for each fixed ξ:
  2 log Λ⟨ξ|Q_+ξ⟩ − W_∞(f_ξ) − ‖[P^c_1, ϑ(f_ξ)]‖_1 ≤ ∫ f_ξ ε^{(Λ)} ≤ γ|⟨ξ₀|ξ⟩|² + 2 log Λ⟨ξ|Q_+ξ⟩ + ‖[P^c_1, ϑ(f_ξ)]‖_1.
- **Failure point (explicit).** Take any ξ ≠ 0 in C_c^∞(int I) with ⟨ξ₀|ξ⟩ = 0, i.e. ξ̂(0) = 0. A CC-type bound ∫ f_ξ ε^{(Λ)} ≤ γ′|⟨ξ₀|ξ⟩|² would need the left side to be ≤ 0, but it is > 0 as soon as

  **log Λ > (W_∞(f_ξ) + ‖[P^c_1, ϑ(f_ξ)]‖_1) / (2⟨ξ|Q_+ξ⟩).**

  So **no Λ-uniform rank-one (Lemma 6.10-type) bound holds for ε^{(Λ)}∘Q alone.** The obstruction is exactly the term 2 log Λ·Q_+, which is unbounded on H. This is the window trace of (P9.3) pushed through Q.
- **For fixed Λ > 1, GAP:**
  - whether ε^{(Λ)}∘Q alone admits some rank-one bound with Λ-dependent γ′_Λ;
  - whether N^{(Λ)}_I keeps CC's form −a_Λ(Id − K^{(Λ)}) with K^{(Λ)} compact.
- **Both GAPs need two inputs:**
  - (a) C² convergence on (1, 2] of the ε^{(Λ)} series, i.e. the App. F analogue: edge values ξ^Λ_n(1) (CC (100)-type) and the KRD tail (PROLATE_BOUNDS §5). [S&M-pending]
  - (b) A bound for χ_Λ∘Q relative to Q_+, from the trace-norm bound for [P^c_Λ, ϑ(f)] and the HS bound for P_ΛŜ_Λ. [S&M-pending]
- **Not needed for P.10.3.** If (a) holds, ∫ Q(k⋆k^*)ε^{(Λ)} is a bounded form (boundary term −2ε^{(Λ)}′(1⁺)‖ξ‖² plus HS), and then (P10.5) forces χ_Λ∘Q to carry the −2 log Λ·Q_+ (CONJECTURE, conditional on (a)).

### P.10.5 Status, labels, R-line
| item | label | track |
|---|---|---|
| P.10.0 (1) relabel (P9.4) | PROVED given (P9.3) + CC §4 at a = b = Λ; (69)–(84) by inspection | [S&M] |
| P.10.0 (2) convention ϑ(λ) dilates by λ (CC (61)); ε ↔ conj ε for complex f | PROVED | [S&M] |
| P.10.0 (3) Tr(ϑŜ_Λ) = Tr(ϑS_Λ) = 2f(1) log Λ via F = I∘u | PROVED | [S&M] |
| P.10.0 (4) ∫ f ε^{(Λ)} = 2f(1) log Λ + O_f(1); consistent, no sign error | PROVED (given (P9.4)) | [S&M] |
| plunge-mode carrier of the log Λ² mass | CONJECTURE (heuristic) | [S&M] |
| CC Lemma 6.10 / Thm 6.11 at Λ = 1 | CERTIFIED (CC, computer-assisted; not replayed) | — |
| monotonicity lemma; Λ → ∞ limits (P10.4), L = −ε_∞ − lim χ_Λ | PROVED | [S&M] |
| Lemma 6.10_Λ (P10.2), Thm 6.11_Λ (P10.3), γ_Λ, c_Λ nonincreasing, ≤ CC's | PROVED given CERTIFIED CC Lemma 6.10 | [S&M] |
| no Λ-uniform rank-one bound for ε^{(Λ)}∘Q alone (P10.5), explicit threshold | PROVED (given (P9.4)) | [S&M] |
| fixed-Λ form of ε^{(Λ)}∘Q; N^{(Λ)} = −a(Id − K) | GAP (inputs (a), (b) S&M-pending) | [S&M] |
| support past p = 2 | not touched (Celestial §C.13) | [Cel] |
| M1 leak | not touched (Queen) | [Q] |
| RH | OPEN | — |

**R-line:**
- **R1:** P.10.3 holds for every Λ ≥ 1, uniformly. This is not claimed for the ε^{(Λ)} series.
- **R2:** not engaged (support fixed at [2^{−1/2}, 2^{1/2}], no primes).
- **R3:** operator positivity / monotonicity of Sonin projections.
- **R4:** n/a.
- **R5:** fake-blind. Every object with ζ's archimedean place, and G_θ′ in particular, passes every step, since only archimedean data enter.
- **R6:** (P10.4) and (P10.5) are exact; the bounds are (P10.1)–(P10.3).

## P.11 Not used
- HiddenZeros (GPPVerify #454, 5f54359): 4–6-point amplitude identities; no celestial amplitude enters §P. **Not cited.**
- `temperedness_iff_critical_line` / ArithmeticTemperedPrincipalSeries.lean (custom axiom `exp_growth_not_tempered`) and ArithmeticCoPoissonOrientationTemplate.lean (co-Poisson assumed): **not cited**. P.3 replaces them with literature statements.
- `open_*` placeholders: never cited.

## P.12 Corrected fakes table (replaces ARP_plan §5; proposed refinement of ARP G.4e's F_θ and H₋ₜ rows)

| fake | P.1 local data | P.2 local Θ_v unitary | P.3 GL(2) unitary / Knapp–Stein PD | GL(1) constituents unitary | P.4 twist | P.5 E / Poisson | P.6 Euler operator | breaks at |
|---|---|---|---|---|---|---|---|---|
| **G_θ′** | yes | yes | yes (tempered) | **yes** | norm | yes | α > 1 | **passes every step** |
| **F_θ** | yes | product yes; constituents no | **yes for 0<θ<½** (complementary series, non-tempered) | **no** (\|·\|^{±θ}) | **polarization** | yes | α > 1+θ | **P.4** (global constituent), and G.4 |
| DH | no | — | no single rep | — | n/a | no ⊗′ | none | P.1 / P.5 |
| H₋ₜ | no (deformed) | Θ_∞ not Tate's | no rep | — | n/a | no E | none | P.1 / P.2 |
| DMV | Beurling | no local field | no rep | — | n/a | no lattice | Beurling | P.5 |
| P_X | = ζ for p ≤ X | = ζ | no rep (no FE) | — | n/a | no Poisson | finite (bounded) | P.5 / P.6 (and R1) |
| Epstein | two Euler products | per constituent | sum of two dihedral reps | per constituent | n/a | no single ⊗′ | sum of two | P.5 |
| G₅ | no Hecke operators | — | no adelic rep | — | n/a | no E | none | P.1 |

## P.13 Label changes relative to ARP_plan.md (2026-10-22)
1. **§5 F_θ row:** "local Θ_v unitary (tempered): no (non-tempered)" → **"GL(2) local Θ_v unitary: YES (product γ, P.3 b1); Knapp–Stein form PD for 0<θ<½ (b2); tempered: no; GL(1) constituents unitary: NO; breaks at P.4 (polarization)."** The column header "unitary (tempered)" conflated two different notions; it is now split.
2. **§1 "Which form is Hermitian"** (Knapp–Stein form "not the L² form, E does not intertwine it"): superseded by P.3. The Knapp–Stein form is positive definite on the complementary series. The relevant point is that the explicit formula and E are linear in the GL(1) constituents.
3. **§3 item 2** ("unitarity holds only for tempered local data; that is the R3 filter"): corrected to "the R3 filter is unitarity of the GL(1) isobaric constituents on C_ℚ; local GL(2) unitarity holds for complementary series."
4. **§2 temperedness:** "filter, not source" is kept. The filter is now located as global (P.3 b3); local non-temperedness is its shadow.
5. **§5 H₋ₜ row:** "local Θ_v yes, ⊗′ yes" → **no Euler product after deformation; Θ_∞ not Tate's** (agrees with ARP G.4e).
6. **§4 R4 CONJECTURE** (soft mode = window image of the two Poisson boundary directions): **SUPERSEDED**, not tested. ARP G.3a supplies R4 through the mean field, and S103/G.3d withdrew T1 (OS4 is fully used by the mean field).
7. **Open items:** (a) E citation → **closed** (P.5). (c) polarization lemma → **closed** (P.4). (b) δ_Λ in Connes–Consani notation → **closed** (P.8). The ARP_plan §3 formula W = Tr(Π ϑϑ^* Π) + δ_Λ is withdrawn (range-type, WRONG-SIGN) and replaced by the deficit form (P8.2). CC 2021 Thm 1 matches its sign; the 2 log Λ test is vacuous at Λ = 1 (P.8c). (d) is superseded with item 6.
8. **New PROVED-LEAN citations:** `weilForm_nonneg`, `cos_mul_positiveType`, `positiveType_add`, `const_one_positive_type`, `positiveType_weighted_gram`, `neg_logDeriv_product_of_hasDerivAt`, `F_one_sub`, `osReflection_eq_self_iff`, `primeSide_eq_truncation`, `primeSide_eq_zero_of_support_lt_log_two`, `tate_local_zeta_integral`, `euler_factor_bridge` (all checked on GPPVerify main 5f54359). Lean `HasSupportIn` uses L′ = 2L.
9. **R6 / exact identity (S103d, Golden Eye):** ARP_plan §3's "as an exact identity (R6)" and §P's R6 wording → PROVED part = nothing local bounds G.4d + §S S(c)4 no MF-relative margin. "The closing step must be an exact identity" → **CONJECTURE**, since margins in norms degenerating on the dilated nodeless bumps are not excluded (P.7 scope caveat).
10. **§4 boundary-direction identification (S103d):** the window image of (f(0), f̂(0)) is span{e^{±u/2}}, the G.3a/OS4 pole pair → **PROVED** (P.8a). Soft-mode alignment stays **NUMERICAL**, and the route stays SUPERSEDED (item 6).
11. **Base case wording (P.7):** "positivity there is CC Thm 1" → CC Thm 1 is W_∞ ≥ Sonin trace on the pole-free class only.
12. **S103e (P.9):**
   - (P8.2) → (P9.0), the S − SQ′S form (Queen (c2)); the two agree iff Q′ ≤ S_Λ.
   - The minor GAP "projector dictionary" is closed NEGATIVELY: the projectors do not agree (D1–D3), and the dictionary holds at the trace level.
   - CC Thm 4.7 holds at every Λ ≥ 1 (P9.4), and the 2 log Λ is the window trace, which cancels exactly (P9.3).
   - The Λ > 1 CC Thm 1 is a GAP [S&M], RH-silent (closed in P.10.3; see item 13).

13. **S103f (P.10):**
   - (P9.4) relabelled (Golden Eye).
   - The Λ > 1 CC Thm 1 / Lemma 6.10 GAP [S&M] is **closed**: PROVED given CERTIFIED CC Lemma 6.10, uniformly in Λ, with the same γ and c. It is weaker than Λ = 1, since Λ = 1 is optimal in the Sonin family.
   - The naive ε^{(Λ)}∘Q-only analogue fails (P10.5); its fixed-Λ form is a GAP (S&M-pending inputs).

14. **S103g (P.14):**
   - The p = 2 edge form is **not** the 2-adic Tate reflection restricted. It is the edge compression of the m = 1 time delay γ₂^*[U, γ₂] in the Q-norm (PROVED).
   - The even-help / odd-charge split comes from J = V(2)P (archimedean-global).
   - Sign fix for the N.13.4 D_L edge block (+c₂(1/ℓ − ℓ/24)σ_x). It was applied on main by Prime's erratum N.14.0(a) (656de12, merged 77730c4), so P.14 no longer asks for it. The exact constant is (1/ℓ − ℓ/24). *(Reworded S103h: the G.4h item-2 display was already correct in 5354fde, and the complaint about it was written against 3470f1f and is withdrawn.)*

15. **S103h (P.15):**
   - The {2,3} prediction of P.14.8 is refined: per-line constants are confirmed, and "P-parity decides per line" is refuted at the p = 3 entry (PROVED). The orbit-closed node set is forced.
   - The Q.8.E weight-1 block has spectrum {0, ±c₂/ℓ}. Positive line blocks cost exactly c₂·g(0) (sharp). A separable (I2) implies a global S2 (⇒ only; REMARK); S2 stays a GAP. *(Grok audit S103h(2): "equivalent … (PROVED)" overstates P.15.2. What is proved is the ⇒ direction (necessary condition) for separable choices whose line block has diagonal sums only at {0, ±a} on (−2L, 2L). The converse needs an operator D^∞ ≥ 0 realizing the given diagonal sums and is not shown. So read "(I2)-separable ⇔ global S2" as a REMARK; global S2 itself is a GAP.)*
   - Status: C.14 is merged; C.14.K is a CONJECTURE *(S103j: CERTIFIED false at (log 2/6, 14))*; the C.14.A indicative failure is NUMERICAL (Celestial's). Piece 4 is the sole odd payer, in band-split form.
   - N.14.3 cross-check: the identities agree. "No 1/ℓ₃ term from the 2-line" holds in ‖ξ′‖_∞ but not in H-energy. (L4) is consistent with P.14.3. The N.13.4 D_L sign was fixed by N.14.0(a), so nothing is left open from P.14.

16. **S103i (P.16):**
   - The sharp joint S1 constant on (log 3, log 4) is c₀\* = ½(c₃ + √(c₃² + 4c₂²)) = 0.90093 (fibre decomposition into 4-paths; PROVED). The N.14.3 c₂ + c₃ sharpness claim gets a FIX.
   - The high band pays the joint charge for T ≥ 16 (CERTIFIED Ω). The 3-entry follows for ℓ₃ ≤ πμ/(4c₀\*(1 + 1/η)T), given the low-band GAP (H₃^σ). No unconditional range in either sector.
   - Q.8.E (I2): D_L = A_LM_ρA_L^\* with an invisible correction n_L. (I2) ⇔ s_L ≥ 0 ⇔ Q_L ≥ 0 (Krein, cited). The explicit n_L is the GAP.

17. **S103j (P.16 fixes, P.17):**
   - P.16 wording: T = 16 (η ≤ 0.0373) is the usable certified point; Krein step restated for (¼ − ∂²)d with weak invisibility, a verified distributional citation and a measure identity; "(H₃^σ) = N.15.H" downgraded to REMARK; C.14.K status lines updated to CERTIFIED false at (log 2/6, 14).
   - P.17: (H₃^σ) at T = 16 is NUMERICALLY false; Q^G_L is not monotone in L, the full Weil form is; Temple–Lehmann lemma for Q_L (PROVED) with a C.15-type ν₂ lemma (PROVED) and an exact H\*-residual formula; odd 3-entry range via one certificate at the right endpoint; Celestial certification list; coupled line cost reduced to a strip-supported moment problem on ℤ² (duality PROVED, value → c₀\* NUMERICAL, exact attainment GAP).

## P.14 The p = 2 line as local Tate data: is Prime's edge form the 2-adic Tate reflection? (S103g; analytic only; [Cel]/[Q], support track L, Λ = 1)
Daniel Toupin / Golden Physics Project, with AI assistance from Grok. **RH is OPEN. Nothing in P.14 proves or implies RH.** On the range treated here (2L ∈ (log 2, log 4)), ζ ≡ P_X for X ≥ 2, so everything is **RH-silent**. Positivity on the range is already CERTIFIED (Mori-conditional, main 992a011); P.14 is about **mechanism**.
- **Sources** (read this session): main 3470f1f (ARP G.4g, G.4h; §N N.12–N.13 merged), team/celestial-cft 9f0473b (§C.13, C.13.B, Lemma C.13.B, Conj C.13.H), team/queen ARP_queen_M1.md, ARP_queen_Q8.md (Q.8.E (I1)–(I2)). *Status update (S103h, main 9cd3187):* Celestial §C.14 has been merged since c7c94ee (merge a4153a8) and is no longer pending. Its coupling bound Thm C.14.B is PROVED, Conjecture C.14.K is a CONJECTURE, and Celestial's indicative C.14.A check (ff3a9dc, merged 220f19c) is NUMERICAL/indicative. Prime §N.14 is merged (77730c4), and S&M's CERT_p2_constants is merged (9d72bbb). Main is now b104678. See P.15.3.
- **Notation** (ξ-coordinates, as N.12.0/N.13): a = log 2, c₂ = Λ(2)·2^{−1/2} = a/√2, 2c₂ = 0.98026, ℓ = 2L − a, S₋ = [−L, −L + ℓ], S₊ = S₋ + a, J x = −a − x on S₋, P x = −x (global inversion; the I of F = I∘u, CC (25)), τ_a f(x) = f(x + a) =: V(2) (dilation by 2 in u = log|x|), ‖u‖²_H = ∫(u′² + ¼u²).
  - Q_L(Qh) = R(ξ) − 2c₂B(ξ), B(ξ) = ⟨ξ, τ_aξ⟩_H = (1/2π)∫(¼ + t²)cos(ta)|ξ̂|² (N.12.1, C.13.A).
- **Credits.** Prime (N.12.3 parity/J split, N.13.4 ramp form), Celestial (C.13 T7–T9, C.13.B (3a)–(3d)), Golden Eye (J-commutation; the G.4h sign and −ℓ/24 points), Queen (Q.8.E).

### P.14.1 The core answer
**No: Prime's rank-two edge form is not the local Tate reflection at p = 2 restricted to edge traces.**
- **What it is (PROVED, P.14.2–P.14.3):** it is the edge-trace compression of the **m = 1 Fourier mode of the phase derivative (Wigner–Smith time delay) of the local Tate unitary γ₂ at p = 2**, measured in the Q-norm H = H¹_¼.
  - The 2-adic input is (a, c₂) only.
  - The reflection that produces the even-help / odd-cost split is **J = V(2)∘P**, an archimedean-global involution: the global inversion P composed with one 2-step. J is not the 2-adic Tate reflection Θ₂.
  - **What breaks the identification** (P.14.4):
    - (i) a derivative, shown by the factor Λ(2) = log 2;
    - (ii) winding number (J ~ V(2), Θ₂ ~ P);
    - (iii) rank/normalization (Θ₂ has a diagonal ½ and asymmetric hops);
    - (iv) the −ℓ/24 term is the ¼ of Q, from the global pole pair, not 2-adic.
- **Sign (Golden Eye (1)), PROVED from first principles (P.14.5):** even edge pairs are **helped** and odd pairs are **charged** exactly 2c₂α²(1/ℓ − ℓ/24).
  - As a contribution to Q_L, or as a block of D_L, the sign is **+**: Q_L gains +2c₂(1/ℓ − ℓ/24)α₋α₊. G.4h item 2 on main has read exactly this since 5354fde, and it is correct. *(Reworded S103h. The first version of this bullet called the G.4h display wrong; it was written against 3470f1f, was out of date, and is withdrawn.)*
  - N.13.4's D_L bullet "−c₂(1/ℓ − ℓ/24)·[[0,1],[1,0]]" had the wrong sign. It must be +c₂(1/ℓ − ℓ/24)σ_x, with eigenvalue + on the even pair, as Celestial C.14.1 also records. *(S103h: fixed on main by Prime's erratum N.14.0(a), 656de12. Nothing is left to fix.)*
  - The text "even helped, odd charged" (N.13.4, G.4g, G.4h) is **right**.

### P.14.2 Local Tate data at p = 2 in the ζ-sector — PROVED (Tate thesis §2; elementary) [Q]
- **Schwartz–Bruhat side.** On S(ℚ₂), with ψ₂ standard and vol(ℤ₂) = 1, F₂1_{ℤ₂} = 1_{ℤ₂} and F₂1_{2^kℤ₂} = 2^{−k}1_{2^{−k}ℤ₂}.
  - Z(1_{2^kℤ₂}, s) = 2^{−ks}L₂(s), with L₂(s) = (1 − 2^{−s})^{−1}. The local FE is Z(F₂φ, 1 − s) = γ₂(s)Z(φ, s), with **γ₂(s) = L₂(1 − s)/L₂(s)**, ε₂ = 1 (unramified).
- **Unitarity.** On Re s = ½, with c := 2^{−1/2} and V := e^{−iat} (the Mellin symbol of V(2)): **γ₂(½ + it) = (1 − cV)(1 − cV̄)^{−1}, |γ₂| = 1.** This is G.6 Lean ask 3, `local_gamma_unitary`.
- **Operator form on the ζ-sector.** In the ζ-sector L²(ℝ_u) the ℤ₂^*-invariant 2-adic data become functions of V(2), via E's q = 2^{±1} summands, the only non-isometric step (P.6).
  - So the local Tate unitary is the bounded translation-invariant operator γ₂(V(2)) = (1 − cV(2))(1 − cV(2)^*)^{−1}.
  - The **local Tate reflection** is **Θ₂ := P∘γ₂(V(2))**. It is a unitary involution, because Pγ₂(V)P = γ₂(V^*) = γ₂(V)^*.
  - Its archimedean counterpart is Θ_∞ = P∘γ_∞, which is CC's F = I∘u (Lemma 1.3, (25)) in Mellin form.
  - Note that the global P is the product of the local inversions. On C_ℚ only the product exists, so the local I₂ cannot be separated from I_∞ in the ζ-sector.
- **Action on 1_{ℤ₂} and on the q = 2^{±} shifts.** In ℓ²(valuations) ≅ L²(2^ℤ ⊂ ℝ_u), the expansion is
  γ₂(V) = (1 − c²) − cV + c(1 − c²)V^* + Σ_{n≥2} c^n(1 − c²)V^{*n}.
  - So on the unit-shell vector δ₀ (the image of 1_{ℤ₂^*}), γ₂ has stay-amplitude ½, hop +1 amplitude −2^{−1/2}, and hop −1 amplitude 2^{−3/2}.
  - On the window (2L < 2a), only the n = 0, ±1 terms survive the compression S·S.

### P.14.3 The p = 2 Weil operator is the time delay of γ₂ — PROVED (elementary; normalization checked against (IF-F)) [Q]/[Cel]
- **Phase.** Let φ₂(t) := arg γ₂(½ + it) = −2 arg L₂(½ + it) = 2Σ_{m≥1}(c^m/m) sin(mat). Then
  **φ₂′(t) = Σ_{m≥1} 2Λ(2^m)2^{−m/2} cos(m t log 2),**
  and with U = multiplication by u (Mellin: U ↔ −i∂_t),
  **(P14.1)  γ₂(V(2))^*[U, γ₂(V(2))] = φ₂′(D) = Σ_{m≥1} Λ(2^m)2^{−m/2}(V(2)^m + V(2)^{−m}).**
- **Archimedean analogue.** γ_∞(½ + it) = Γ_ℝ(½ − it)/Γ_ℝ(½ + it) = e^{−2iθ(t)}, so −(arg γ_∞)′ = 2θ′ = Ω.
- **The Weil multiplier as time delays.** The Weil multiplier on this range is therefore −Σ_v (arg γ_v)′ truncated by the support:
  **Q_L(Qh) = (1/2π)∫(¼ + t²)[Ω(t) − φ₂′(t)]_{trunc}|ξ̂|²dt = R(ξ) − ⟨ξ, φ₂′(D)ξ⟩_H|_{window}.**
  - For 2L < log 4 only m = 1 survives, and the m = 1 part of φ₂′ is **2c₂cos(ta)**, exactly the (IF-F) multiplier. In operator form it is **c₂(V(2) + V(2)^*)**. So −c₂S(V(2) + V(2)^*)S is the window compression of the m = 1 mode of (P14.1).
  - This is the 2-adic analogue of CC's local trace formula T1 (W_ℝ = phase derivative of u_∞), and Connes 1999 §VII's local term ∫′ h(u^{−1})/|1 − u| d^*u at v = 2.
- **Normalization (REMARK; PROVED by a diagonal-sum computation).**
  - In the **form picture** (ξ ↦ ⟨ξ, ·ξ⟩_H on 𝒯_L), the line is the uncompressed c₂(V(2) + V(2)^*), with weight 1.
  - In Q.8.E's **trace picture** (I1), D_L enters only through its diagonal sums d(b) = ∫D(x, x + b)dx. For S = 1_{[−L,L]} (= S_Λ at Λ = e^L), SV(2)S has kernel 1_S(x)δ(y − x − a)1_S(y), so d(b) = |S ∩ (S − a)|δ_a(b) = **ℓ·δ_a**.
  - So the block "−c₂S_Λ(V(2) + V(2)^*)S_Λ" (C.13.4, N.12.5, G.4g) carries −c₂ℓ(δ_a + δ_{−a}), **off by the factor ℓ**. A weight-1 block is −(c₂/ℓ)S(V + V^*)S, or −c₂(V(2)M_ρ + M_ρV(2)^*) with ρ ≥ 0, ∫ρ = 1 on S ∩ (S − a).
  - *Grok audit (S103g): checked. S ∩ (S − a) = [−L, L − a] has length 2L − a = ℓ, so d(b) = ℓδ_a for SV(2)S. In the form picture, ⟨ξ, SV(2)Sξ⟩ = h(a) has weight 1, so the factor ℓ appears only in the trace picture. ARP.md G.4g and Q.8.E have been corrected accordingly.*
  - Since (I1) fixes only diagonal sums, the block is not unique, and (I2) depends on the choice. Nothing in (IF), N.12 or N.13 uses the block, so no result changes.

### P.14.4 What breaks the "local Tate reflection" identification — PROVED (each item elementary) [Q]
1. **Derivative, not the reflection: the Λ(2) factor.**
   - Every matrix coefficient of γ₂(V(2)) and Θ₂ lies in ℚ(√2): powers of 2^{−1/2} and (1 − c²) = ½.
   - The edge form's weight is c₂ = (log 2)·2^{−1/2}, which contains Λ(2) = log 2. Only ∂_t brings it down, i.e. the commutator with U in (P14.1).
   - So the edge form is (a compression of) **γ₂^*[U, γ₂]**, the infinitesimal generator of the 2-adic scattering phase, and not a restriction of Θ₂ or γ₂.
2. **Winding number (parity of the reflection).**
   - J = V(2)P = Pτ_{−a} is the reflection about −a/2. Its "Tate symbol" is V(2)^{±1} = e^{∓iat}, of winding number ±1 on the valuation circle t ∈ ℝ/(2π/a)ℤ.
   - Θ₂ = Pγ₂(V(2)), and γ₂(e^{−iat}) = (1 − cV)/(1 − cV̄) has **winding number 0**, since c = 2^{−1/2} < 1. So Θ₂ is homotopic, through unitary involutions P∘γ, to P, and **not** to J.
   - Formally γ₂ → −V as c → 1, i.e. as if the local factor sat on the line Re s = 0: then Θ₂ → −PV = −(reflection about +a/2). That is J's mirror with a sign, and it is reached only at the edge of unitarity.
3. **Rank and normalization on edge traces.**
   - On the two-site edge space (u₋ on S₋, u₊ = τ_aξ|_{S₋}), the compression of γ₂ is [[½, −c], [c/2, ½]] (stay ½; hops −2^{−1/2} and 2^{−3/2}). It is not symmetric, has nonzero diagonal, and does not act as ±1 on the even/odd edge pairs.
   - The edge form is the symmetric, traceless c₂·σ_x (from c₂(V + V^*)). Its ± eigenvectors are exactly the J-parity pairs.
4. **The −ℓ/24 is the ¼ of Q (global pole pair), not 2-adic.** On the unit ramp r on S₋ (r(−L) = 0, r(−L + ℓ) = 1), τ_ar₊ = Jr. Hence
   ⟨r, Jr⟩_H = ‖P₊r‖²_H − ‖P₋r‖²_H with **P₊r = ½·1_{S₋}** (because r + Jr ≡ 1) and P₋r = r − ½. So
   **‖P₊r‖²_H = ℓ/16,  ‖P₋r‖²_H = 1/ℓ + ℓ/48,  ⟨r, Jr⟩_H = ℓ/24 − 1/ℓ  (exact).**
   - The 1/ℓ is the H¹ seminorm of the J-odd linear part P₋r.
   - The −ℓ/24 in (1/ℓ − ℓ/24) = ‖P₋r‖²_H − ‖P₊r‖²_H is exactly the ¼-mass term: ¼(‖P₋r‖²_{L²} − ‖P₊r‖²_{L²}) = ¼(ℓ/12 − ℓ/4) = −ℓ/24.
   - The ¼ is the constant in ¼ + t² = s(1 − s) at s = ½ + it, i.e. the pole-removal Q = −∂² + ¼ from the global poles at s = 0, 1 (OS4, spent once at C.13 T4; P.8a).
   - So any "Tate identification" reproduces −ℓ/24 only through the **global** factor s(1 − s), never from γ₂. It is not a non-local correction to a local identity. It is the Q-geometry, where H¹_¼ is the pole-free norm.

### P.14.5 The sign, from first principles (Golden Eye (1)); which display is wrong — PROVED [Cel]/[Q]
- **Even ξ.** ξ(x + a) = ξ(−a − x), so τ_aξ|_{S₋} = J(ξ|_{S₋}) and B(ξ_e) = ⟨u, Ju⟩_H.
- **Odd ξ.** B(ξ_o) = −⟨u, Ju⟩_H. Equivalently τ_a = J∘P on the strip, the two unitary involutions of the problem: J on H(S₋) (an isometry of the Q-norm) and P (global inversion).
- **On the edge ramp** u = αr: ⟨r, Jr⟩_H = ℓ/24 − 1/ℓ < 0 (P.14.4 item 4), because a monotone edge trace is dominantly J-odd. So:
  - **even:** B(ξ_b) = α²(ℓ/24 − 1/ℓ) < 0, and the Q_L-contribution is −2c₂B = **+2c₂α²(1/ℓ − ℓ/24)** (help);
  - **odd:** B(ξ_b) = α²(1/ℓ − ℓ/24) > 0, and the Q_L-contribution is **−2c₂α²(1/ℓ − ℓ/24)** (charge).
- **Rule.** The sign of the p = 2 contribution on an edge trace = −(P-eigenvalue)·(J-eigenvalue of the dominant component) = −σ·(−1) = σ. The local unitarity behind the split is that of **J (archimedean, on H¹_¼(S₋)) together with P**. The 2-adic Θ₂ gives no sign split (item 3: diagonal +½ on both parities).
- **Which side has the error.**
  - N.13.4's identity 2c₂B(ξ) = 2c₂B(ξ_s) − 2c₂α₋α₊(1/ℓ − ℓ/24) + (c₂/2)X is **correct**. So is its interface form [R(ξ_b) + 2σc₂α²(1/ℓ − ℓ/24)].
  - **Wrong sign:** N.13.4's "D_L … is the matrix −c₂(1/ℓ − ℓ/24)·[[0,1],[1,0]] in (α₋, α₊)" should be **+c₂(1/ℓ − ℓ/24)σ_x**, since D_L carries −c₂(V + V^*) and ⟨ramp₋, τ_a ramp₊⟩_H < 0.
  - *Grok audit (S103g): main 5354fde had already fixed the G.4h item-2 display. It now reads "contribution to Q_L … +2c₂(1/ℓ − ℓ/24)α₋α₊". This section was written against 3470f1f. The N.13.4 D_L-bullet finding is confirmed: team/prime-number-theory/ARP_prime-number-theory.md, §N.13.4 "Interface form", second-bracket sub-bullet.*
  - *Withdrawn (S103h):* the earlier claim that G.4h item 2's display had the (IF)-side sign. On main 5354fde and later it reads "contribution to Q_L … +2c₂(1/ℓ − ℓ/24)α₋α₊", which is correct. The complaint referred to 3470f1f.
  - "Charged exactly 2c₂α²(1/ℓ − ℓ/24)" (Golden Eye (2)) is now on main in G.4h item 2 and the P.14 record. N.13.4's two "exactly 2c₂α²/ℓ" occurrences were corrected by erratum N.14.0(a). The parity verdict is unchanged.
  - **Consistency check on the identification:** a genuine Θ₂-restriction would give the same sign on both parities, contradicting the PROVED help/charge split. This independently confirms P.14.4 item 3.

### P.14.6 Support: isotypic decomposition of −c₂S(V(2) + V(2)^*)S — PROVED (elementary) [Cel]
- **Dihedral structure.** The group generated by P and V(2) is the infinite dihedral group D_∞. J = V(2)P is a reflection in it.
  - V(2) + V(2)^* is central in C^*(D_∞): every reflection conjugates V to V^*, which is C.13.B (3a) in group form.
  - It acts on the irreducible π_θ (θ = ta mod 2π) by the scalar 2cos θ. So the cosine multiplier 2c₂cos(ta) is literally the central character.
- **Isotypic piece hit.** On the two-site edge space, (P, J, V(2) = JP) act by commuting signs. The harmful pieces are exactly the **V(2) = +1** pieces:
  - (P, J) = (+, +), the even harmful cone of N.12.3 (K3);
  - (P, J) = (−, −), the odd harmful cone (K2).
  - The helped pieces are V(2) = −1, i.e. (+, −) and (−, +). These are the θ = 0 and θ = π endpoints of the D_∞ dual.
  - On the Fourier side, Celestial's blind frequencies t_j = (2j + 1)π/a, for symmetric pairs, are exactly the θ = π frequencies. The archimedean interference factor 1 + cos ta vanishes on the frequencies whose global character is V(2) = −1, while the harmful edge piece has two-site label V(2) = +1. That is C.13.B (3c) restated.
- **Can it beat the harmful cones?**
  - Even (K3): on edge traces the p = 2 line itself pays, through the sign above. What remains is the ĝ(0) cost paid by 𝒮_L (G.4h, (P), GAP b-1) and the coupling (N.G2-bulk).
  - Odd (K2): the p = 2 line is a pure cost, 2c₂α²(1/ℓ − ℓ/24). No 2-adic structure can pay it (P.14.4: Θ₂ has no sign split), so it must come from R's margin (N.G2-m₀; S&M Thm 2.1 / Prime Thm N.13.A).
  - For ξ odd, ĝ(0) = 0, so (K2) at ℓ = 0 is exactly CC Thm 1 (≥ 0, no margin).
  - PROVED as bookkeeping; the closures are GAPs as in G.4h.
- **Golden Eye's commutation.** B and R_L commute with J on strip-supported vectors (C.13.B (3a)). J swaps the two pole functionals ∫k e^{±x/2} with weights 2^{∓1/2}: ∫(Jk)e^{x/2} = 2^{−1/2}∫k e^{−x/2}.
  - The 2^{−1/2} is the same half-density factor as in c₂ = Λ(2)·2^{−1/2} (PROVED, one line). So J intertwines the two pole conditions, and the pole pair is still spent once (T4).
- **Sonin/prolate.** V(2) = ϑ(2) moves the Sonin family along the orbit ab = 1 (P.9.2(a)): ϑ(2)S^{So}_{(1,1)}ϑ(2)^* = S^{So}_{(2,1/2)} ≠ S^{So}_{(1,1)}.
  - So V(2) has no isotypic component inside the Sonin/prolate decomposition. It is an intertwiner between members of the family, not an endomorphism. This is the operator form of C.13 T9 (σ_S is smooth, so it carries no cos(ta)). REMARK, PROVED as stated.
- **χ_Λ and Queen's leak.** Nothing exact relates χ_Λ (P.9) to the p = 2 line or to M1 here. Not claimed.

### P.14.7 Placement in G.4h's four-part split of N.G2
| G.4h piece | what P.14 adds | label |
|---|---|---|
| 1 strip-localized cone (Lemma C.13.B, ℓ ≤ 0.119) | the harm there is the V(2) = +1 piece; absorbed by Ω's high band. No 2-adic input beyond (a, c₂) | PROVED (bookkeeping) |
| 2 edge traces (ramps) | the form is the edge compression of the m = 1 time delay of γ₂ (P14.1) in the Q-norm. The help/charge split comes from J∘P (archimedean-global), not Θ₂. Sign: Q_L gains +2c₂(1/ℓ − ℓ/24)α₋α₊ (as in G.4h item 2 since 5354fde); the N.13.4 D_L edge block is +c₂(1/ℓ − ℓ/24)σ_x (fixed by N.14.0(a), 656de12). Exact constant (1/ℓ − ℓ/24), with −ℓ/24 = ¼-mass term of Q | PROVED |
| 3 coupling 2R(ξ_b, ξ_s) | Celestial C.14 merged (c7c94ee): Thm C.14.B PROVED, C.14.C obstructions PROVED, closure reduced to Conj C.14.K (CONJECTURE; S103j: CERTIFIED false at (log 2/6, 14)) + m₀^e. P.14 adds nothing 2-adic here | GAP N.G2-bulk |
| 4 archimedean margin at ½ log 2 | the odd cost cannot be paid 2-adically (no sign split in Θ₂), so it must be R's margin. The ℓ-independent frame version (Thm C.14.A) indicatively fails on 7.4 ≲ \|t\| ≲ 13 (Celestial ff3a9dc, NUMERICAL/indicative); the band-split routes remain (P.15.3) | PROVED (as stated); value GAP N.G2-m₀ |

### P.14.8 Remark for Prime's §N.14 ({2, 3} lines, 2L ∈ (log 3, log 4]) — REMARK / CONJECTURE [Cel]
*(Tested in P.15.1, S103h: (a)–(c) confirmed; "P-parity decides help/charge per line" refuted at the p = 3 entry, PROVED.)*
- **Every reflection.** Every reflection x ↦ c − x conjugates V(p) to V(p)^*, so each V(p) + V(p)^* commutes with every reflection on L²(ℝ). The window [−L, L] is preserved only by c = 0. So **P is the unique window-preserving reflection commuting with every B_p** (PROVED).
- **Strip reflections.** The strip reflections J_p: x ↦ −log p − x are p-dependent. J₂ and J₃ generate a D_∞ whose translation is log(3/2), so there is **no common J** (PROVED). The local Tate data γ₂(V(2)) and γ₃(V(3)) commute with each other and with nothing reflection-like except through P.
- **CONJECTURE (prediction for N.14).**
  - The {2, 3} parity split is by P, with per-line J_p on the nested strip pairs (lengths 2L − log 2 ⊃ 2L − log 3). The two lines' edge forms interact through the overlap of their strips.
  - The rank-two ramp structure persists per line, with per-line constants c_p(1/ℓ_p − ℓ_p/24).
  - Since J is not 2-adic, no 2-adic/3-adic local reflection will appear as the commuting symmetry.

### P.14.9 Labels and R-line
| item | label | track |
|---|---|---|
| γ₂ unitary, Θ₂ = Pγ₂(V(2)) involution, action on 1_{ℤ₂}/shifts | PROVED | [Q] |
| (P14.1) p = 2 Weil operator = time delay of γ₂; m = 1 mode = c₂(V + V^*) | PROVED | [Q]/[Cel] |
| core: edge form is not Θ₂ restricted (items 1–4) | PROVED | [Q] |
| sign: even helped, odd charged 2c₂α²(1/ℓ − ℓ/24); D_L edge block sign fix | PROVED | [Cel] |
| −ℓ/24 = ¼-mass term of Q (global pole pair) | PROVED | [Cel] |
| D_∞ isotypic decomposition; harmful = V(2) = +1 | PROVED | [Cel] |
| (I1) pairing normalization of −c₂(V + V^*) | REMARK | [Q] |
| Sonin/prolate: V(2) intertwines S(1,1) → S(2,½) | REMARK (PROVED as stated) | [Cel] |
| {2, 3}: P unique common reflection, no common J (PROVED); N.14 prediction | REMARK / CONJECTURE | [Cel] |
| N.G2-bulk, N.G2-m₀, (P) | GAP (unchanged; C.14 merged; C.14.K CERTIFIED false at (log 2/6, 14), S103j; see P.15.3) | [Cel] |
| RH | OPEN | — |

**R-line:**
- **R1:** one rule on 2L ∈ (log 2, log 4). (P14.1) holds on every window, with more primes added.
- **R2:** the Euler product enters through the local γ₂ at p = 2 (its time delay); the FE through γ_∞ (Ω = −(arg γ_∞)′).
- **R3:** the sign of B by J∘P parity, not weight signs.
- **R4:** the pole pair is spent once at T4. Its ¼ reappears as −ℓ/24, and J intertwines the two pole functionals.
- **R5, where each fake breaks:**
  - **F_θ:** its constituents γ₂^{(±θ)} are not unitary (G.6 ask 3), so (P14.1) is not a self-adjoint time delay. The p = 2 weight becomes c₂·cosh(θa), a polarization: it breaks at the time-delay step (norm → polarization, P.4).
  - **DH:** no Euler product, no γ₂, no a-line. Breaks at P.14.2.
  - **H₋ₜ:** no Euler product after deformation (no γ₂), and γ_∞ is not Tate's (Ω deformed). Breaks at P.14.2/P.14.3.
  - **DMV (glued) / P_X:** coincide with ζ on this window (N.12.1), so they pass every step here. They are fake-blind on this range, and break only past the ladder (G.4e/G.4f).
  - **G_θ′:** passes. Its local factors are unitary twists, the time delay is shifted t → t ± θ′, and Q_G is a sum of two ζ-forms (G.4e).
- **R6:** (P14.1) and the ramp identities are exact. The closures stay GAPs.

## P.15 The {2, 3} window, the Q.8.E weight-1 block and its (I2) choice, and the status of piece 4 after C.14 (S103h; analytic only; [Cel]/[Q]/[S&M])
Daniel Toupin / Golden Physics Project, with AI assistance from Grok. **RH is OPEN. Nothing in P.15 proves or implies RH.** On 2L ≤ log 4 we have ζ ≡ P_X for X ≥ 3, so everything here is **RH-silent**. Positivity on the range is already CERTIFIED (Mori-conditional, main 992a011); P.15 is about **mechanism**.
- **Sources** (read this session; written on main 9cd3187 and revised on main b104678):
  - ARP G.4g/G.4h (with the lead's Q.8.E factor-ℓ fix);
  - team/queen/ARP_queen_Q8.md, Conjecture Q.8.E (I1)–(I3) and Lemma Q.8.5, whose notation is used in P.15.2;
  - Celestial §C.14 (merged a4153a8, audited c7c94ee);
  - Celestial's indicative C.14.A evaluation at team/celestial-cft ff3a9dc (**not yet on main**; cited as hers, NUMERICAL/indicative);
  - Prime §N.13 (611e39a).
- **Prime §N.14 is merged** (source for P.15.1 and P.15.5: team/prime-number-theory cf4ee37, at Prime's direct request; errata 656de12; merge 77730c4, audit b104678). It was not on main when the first version of this section (641b63c) was written. P.15.1 now cross-checks N.14.3, answering its interface request N.14.4 "Principal: cross-check the per-line ramp constants and the (L4) chain". S&M CERT_p2_constants (merged 9d72bbb) is used in P.15.3.
- **Notation** as in P.14: ‖u‖²_H = ∫(u′² + ¼u²), B_p(ξ) := ⟨ξ, τ_{log p}ξ⟩_H, c_p := Λ(p)p^{−1/2}, so c₂ = log 2/√2 and **c₃ = log 3/√3**. a₂ = log 2, a₃ = log 3, ℓ_p := 2L − log p.
- **Method.** The algebraic identities in P.15.1 were also checked symbolically and exactly (sympy, piecewise integration with symbolic L; `B2 − hand = 0`, `B3 − hand = 0`). This is a check of the algebra, not exploratory numerics. Decimal values below are one-line evaluations of closed forms.

### P.15.1 The two lines on 2L ∈ (log 3, log 4]: per-line constants and the strip-overlap interaction (with Prime's §N.14.3, merged 77730c4) — PROVED [Cel]
On this range, Q_L(Qh) = R(ξ) − 2c₂B₂(ξ) − 2c₃B₃(ξ). At 2L = log 4 the p = 4 term sits at the boundary of supp h and vanishes.

**Lemma P.15.A (orbit-closed node set).** Put x₁ := L − log 3, w := −L + log(3/2) and x₂ := L − log 2. Then
- −L < x₁ < w < x₂ ≤ 0, and N_L := {±L, ±x₁, ±w, ±x₂} is the smallest set containing ±L that is closed under x ↦ x ± log 2 and x ↦ x ± log 3 inside [−L, L].
- The left p = 2 strip S₋^{(2)} = [−L, x₂] splits into the pieces I₁ = [−L, x₁], I₂ = [x₁, w] and I₃ = [w, x₂], of lengths **ℓ₃, q, ℓ₃**, where **q := log(9/2) − 2L** ∈ [log(9/8), log(3/2)). So ℓ₂ = 2ℓ₃ + q, and the p = 3 strip is S₋^{(3)} = I₁.
- τ_{log 3} maps I₁ onto −I₁. τ_{log 2} maps I₁, I₂, I₃ onto −I₃, −I₂, −I₁.

*Proof.*
- Ordering: x₁ < w ⟺ 2L < log(9/2), and w < x₂ ⟺ 2L > log 3.
- Closure: −L + log 2 = −x₂, −L + log 3 = −x₁, x₁ + log 2 = −w, w + log 2 = −x₁, x₁ + log 3 = x₂ + log 2 = L. Every other shift leaves the window, e.g. w + log 3 > L ⟺ 2L < log(9/2). ∎

So the strip-overlap interaction is concrete. **The p = 2 line couples the p = 3 edge strip I₁ to the central piece −I₃** (a strip of length ℓ₃ at distance |x₂| ≤ ½ log(4/3) from the centre), and I₃ to the p = 3 strip −I₁. It does not couple I₁ to −I₁; that is the p = 3 line's job. The p = 3 strips are nested in the p = 2 strips (ℓ₃ < ℓ₂).

**Joint ramp decomposition (PROVED).** Let ξ_b agree with ξ on the bulk [x₂, −x₂] and be piecewise linear between consecutive nodes of N_L on [−L, x₂] ∪ [−x₂, L]. Put ξ_s := ξ − ξ_b, which vanishes at every node.
- For p = 2, 3, **the derivative cross terms of B_p between ξ_b and ξ_s vanish identically.** τ_{log p} maps each piece of S₋^{(p)} onto a piece (Lemma P.15.A), ξ_b′ is constant on pieces, and ∫ over a piece of ξ_s′ (shifted or not) is a difference of ξ_s at two nodes, which is 0.
- Hence B_p(ξ) = B_p(ξ_b) + B_p(ξ_s) + ¼X_p, where X_p is the ¼-mass cross term, exactly as in N.13.4.
- With N.13.4's two nodes {−L, x₂} only, the p = 2 derivative cross term would **not** vanish. ξ_s would need to vanish at x₁ + log 2 = −w, which is not a p = 2 node. So the orbit-closed node set is forced, and the ramp space has **three edge values per side**: v₁ = ξ(x₁), v₂ = ξ(w), v₃ = ξ(x₂).

**Theorem P.15.B (exact joint ramp forms; PROVED, sympy-checked).** For ξ of P-parity σ (ξ(−x) = σξ(x)):
> **B₃(ξ_b) = σ·v₁²(ℓ₃/24 − 1/ℓ₃),**
> **B₂(ξ_b) = σ·[ 2v₁(v₂ − v₃)/ℓ₃ − (v₁ − v₂)²/q + ℓ₃(v₁v₃ + 2v₁v₂)/12 + q(v₁² + v₂² + 4v₁v₂)/24 ].**
- *Proof.* τ_{log p}ξ = σ·J_p ξ on S₋^{(p)}, with J_p x = −log p − x.
  - J₃ reverses I₁, which reproduces N.13.4 at ℓ = ℓ₃.
  - J₂ reverses [−L, x₂] and maps the pieces (ℓ₃, q, ℓ₃) to themselves in reverse order, so J₂ preserves N_L ∩ S₋^{(2)}. Then B₂(ξ_b) = σ⟨u, J₂u⟩_H, with u piecewise linear with values (0, v₁, v₂, v₃) and J₂u with values (v₃, v₂, v₁, 0).
  - The slopes give 2v₁(v₂ − v₃)/ℓ₃ − (v₁ − v₂)²/q. The ¼-mass part comes from ∫_0^h (linear)(linear) = h(2ac + ad + bc + 2bd)/6. ∎
- **Checks.**
  - On the linear profile v₁ = (ℓ₃/ℓ₂)v₃, v₂ = ((ℓ₃ + q)/ℓ₂)v₃ (a single p = 2 ramp), B₂(ξ_b) = σv₃²(ℓ₂/24 − 1/ℓ₂): **N.13.4's per-line constant c₂(1/ℓ₂ − ℓ₂/24) is recovered on its own ramp.**
  - As 2L ↓ log 3 (ℓ₃ → 0, q → ℓ₂), smooth ξ has v₁ → 0 and v₂ − v₃ → 0, which gives N.13.4 at p = 2.

**The p = 3 entry: both parities become indefinite (PROVED).** Let ℓ₃ ↓ 0 (q → log(3/2), bounded). Write u₁ := v₁ (the jump across the p = 3 strip I₁) and u₃ := v₃ − v₂ (the jump across the central piece I₃). The line part of Q_L, −2c₂B₂ − 2c₃B₃, is at order 1/ℓ₃
> **σ·(1/ℓ₃)·(2c₃u₁² + 4c₂u₁u₃) + O(1)·|v|²,  matrix σ(1/ℓ₃)[[2c₃, 2c₂], [2c₂, 0]], determinant −4c₂²/ℓ₃² < 0.**
- The H-energy of ξ_b on I₁ ∪ I₃ and their mirrors is 2(u₁² + u₃²)/ℓ₃ + O(ℓ₃). So per unit energy, the two lines' contribution has the eigenvalues
  > **σ·½(c₃ ± √(c₃² + 4c₂²)) = σ·(+0.9009, −0.2666).**
- **Even sector (σ = +):**
  - The p = 3 line alone only helps (+c₃ = 0.6343 per unit energy).
  - With the p = 2 coupling to the central strip, there is a **charged even direction**: per unit energy −0.2666, with u₃/u₁ = −(c₃ + √(c₃² + 4c₂²))/(2c₂) = −1.838. The outer p = 3 jump and the central jump have opposite signs.
- **Odd sector (σ = −):**
  - The worst charge is **0.9009 per unit energy**, versus c₃ = 0.6343 for the p = 3 line alone (factor 1.420).
  - There is also a helped odd direction (+0.2666).
- The indefiniteness holds for every c₂ ≠ 0. The (I₃, I₃) entry is zero: the p = 2 line maps I₃ to −I₁, off-diagonal, and the p = 3 line does not see I₃.

**Verdict on the P.14.8 prediction.**
- (a) **P is the unique window-preserving reflection commuting with every B_p, and there is no common J_p.** Confirmed and used (PROVED).
- (b) **Per-line J_p exist and act on the nested strips.** Confirmed: B_p(ξ_b) = σ⟨u, J_pu⟩_H on S₋^{(p)}, and J₂ preserves the refined nodes (PROVED).
- (c) **Per-line ramp constants c_p(1/ℓ_p − ℓ_p/24).** Exact for p = 3, and exact for p = 2 on its own single ramp (PROVED). The joint decomposition that both lines need is rank three per parity, not rank two per line.
- (d) **"P-parity decides help/charge for each line" is REFUTED at the p = 3 entry (PROVED).** In each P-parity the order-1/ℓ₃ block is indefinite.
  - Mechanism: this is the quantitative form of "no common J". The J₃-odd p = 3 ramp on I₁ has no definite J₂-parity, because J₂ sends I₁ to I₃, where J₃ does not act. The J₂-even and J₂-odd combinations of the two short pieces (u₁ ± u₃) enter B₂ with opposite signs (the hyperbolic term −2u₁u₃/ℓ₃).
- **Consequences for N.14 (REMARK).**
  - At the p = 3 entry, the even sector is no longer paid by the new line alone, so it also needs an archimedean payer.
  - The odd p = 3 entry cost is larger than c₃ per unit energy.
  - Both costs sit at frequency ~π/ℓ₃, where Ω is large. So an analogue of Lemma C.13.B (strip-localized absorption by the high band of Ω) is the natural tool. **The constant to beat per unit energy is ½(c₃ + √(c₃² + 4c₂²)) = 0.9009, not c₃** (CONJECTURE that this is the sharp S1-type constant at the p = 3 entry; the upper bound is PROVED on the ramp space only).
  - No 2-adic or 3-adic local reflection appears as the commuting symmetry. This is consistent with P.14.4 and confirms that part of the prediction.

**Cross-check with Prime §N.14.3 (merged 77730c4) — PROVED [Cel].**
- **Agreement.** The following agree with N.14.3, as identities:
  - the strip geometry (G), including 2L < log 4.5;
  - the 3-edge block +c₃(1/ℓ₃ − ℓ₃/24)σ_x and the odd charge 2c₃β²(1/ℓ₃ − ℓ₃/24) (= B₃ in Thm P.15.B, with β₋ = v₁);
  - B₂(ξ_s) = 0 for Prime's ξ_s ∈ H¹₀(3-strips) (τ_{log 2} maps I₁ onto −I₃, not onto a 3-strip);
  - the symmetry statement (S).
- **Prime's 2-line derivative term is the hyperbolic term of Thm P.15.B.**
  - In our nodes, a − L = −x₂ and L + a − b = −w, so Prime's (2β₋/ℓ₃)[ξ(L + a − b) − ξ(a − L)] = σ(2v₁/ℓ₃)(v₂ − v₃).
  - "No 1/ℓ₃ term from the 2-line" is right relative to ‖ξ′‖_{L∞}: the term is 2β₋ times a mean slope.
  - It is **not** lower order relative to the H-energy. A profile with a kink at scale ℓ₃ on the central piece I₃ (|v₃ − v₂| ~ |v₁|) has mean slope ~|v₁|/ℓ₃ there. The term is then of the same order 1/ℓ₃ as the 3-edge block, and the combined order-1/ℓ₃ block is the indefinite σ(1/ℓ₃)[[2c₃, 2c₂], [2c₂, 0]] above.
  - Prime's interface form is exact and puts this term inside Q^{(≤2)}(ξ_b), since his ξ_b contains the central pieces. So both statements are correct: they measure the 2-line term in different norms.
- **Consequence for Prime's {2, 3} odd-entry CONJECTURE.**
  - The odd 3-charge must be paid together with the 2-line's 1/ℓ₃-coupling to the central pieces. The worst per-energy odd charge of the two lines on the joint ramp space is **0.9009, not c₃ = 0.6343**.
  - The even sector at the 3-entry has a charged direction (−0.2666 per unit energy).
  - Both live at frequency ~π/ℓ₃, so they are high-band objects: a C.13.B-type absorption, not the low-band (≤ 2) margin, is the natural payer.
  - This sharpens the Grok caveat that Q^{(≤2)} on the larger window is the P₂ form, not known to be ≥ 0. It does not refute the conjecture. REMARK.
- **(L4) chain against the P.14 normalization (PROVED).**
  - The p² = 4 line is the m = 2 mode of (P14.1): Λ(4)·4^{−1/2}(V(2)² + V(2)^{−2}) = c₄(V(2)² + V(2)^{*2}), with c₄ = (log 2)/2.
  - Its weight-1 trace block is −(c₄/ℓ₄)S_Λ(V(2)² + V(2)^{*2})S_Λ, with ℓ₄ = |S ∩ (S − 2a)| = 2L − log 4. In the form picture its edge block is +c₄(1/ℓ₄ − ℓ₄/24)σ_x. Both agree with N.14.3 (L4).
  - Prime's KMS note: Σ_{m≠0} r^{|m|}e^{imat}·log 2 = φ₂′(t) exactly (P14.1), with r = 2^{−1/2}. So the "full 2-adic multiplier −log 2(P_r(at) − 1)" is −φ₂′, the 2-adic time delay, and its Toeplitz sections are the KMS matrices. This is consistent with P.14.3, and no Tate reflection is involved (P.14.4).

### P.15.2 Q.8.E: the weight-1 (I1) block and the (I2) choice, in Queen's notation — PROVED [Q]
Queen's notation (Q.8.E): Λ = e^L, S_Λ = 1_{[−L, L]}, V(u) is the translation group on L²(ℝ_u), and V(f) = ∫f(u)V(u)du. (I1) reads Q_L(f) = Tr(D_L V(f)V(f)^*) for real f ∈ 𝒯_L, with D_L S_Λ V(f) Hilbert–Schmidt. (I2) reads D_L ≥ 0. Let g := f ∗ f^*. For real f, g is even and supported in [−2L, 2L]. For bounded D with kernel D(x, y),
> Tr(D V(f)V(f)^*) = ∫ d_D(b)g(b)db,  d_D(b) := ∫D(x, x + b)dx (the diagonal sums; g even, so the orientation of b drops out).
The p = 2 term of Q_L is −c₂(g(a) + g(−a)), with a = log 2.

**(I1): the weight-1 block (PROVED).** For ℓ = 2L − log 2 > 0:
> **D_L^{(2),1} := −(c₂/ℓ)·S_Λ(V(2) + V(2)^*)S_Λ,  d = −c₂(δ_a + δ_{−a}),  Tr(D_L^{(2),1}V(f)V(f)^*) = −c₂(g(a) + g(−a)).**
- *Diagonal sums.* S_ΛV(2)S_Λ has d = |S ∩ (S − a)|δ_a = ℓδ_a (P.14.3). The trace identity holds for every real f, since the trace pairing is translation-invariant. Support enters only through g.
- *Spectrum.* On L²(S₋) ⊕ L²(S₊), with S₋ = [−L, L − a] and S₊ = [a − L, L] (disjoint, since 2L < 2a), S_Λ(V(2) + V(2)^*)S_Λ = [[0, U], [U^*, 0]] with U unitary, and it is 0 on the bulk and off the window. So
  > **spec D_L^{(2),1} = {−c₂/ℓ, 0, +c₂/ℓ},  ‖D_L^{(2),1}‖ = c₂/ℓ.**
  - The eigenvalue −c₂/ℓ sits on the pairs (φ, U^*φ). These are the V(2) = +1 two-site pieces, the harmful label of P.14.6.
- *HS condition.* ‖S_ΛV(f)‖²_HS = 2L‖f‖²₂, so D_L^{(2),1}S_ΛV(f) is HS.
- *R1 caveat (REMARK).* ‖D_L^{(2),1}‖ = c₂/ℓ blows up as the line enters (ℓ ↓ 0), although its trace does not. The equivalent block −c₂(V(2)M_ρ + M_ρV(2)^*) with ∫ρ = 1 has the same trace for **every** probability density ρ, since localization in S ∩ (S − a) is optional for the trace pairing. With a fixed L-independent ρ it is uniformly bounded, by 2c₂‖ρ‖_∞. **Recommendation: use a fixed ρ for "one rule for all L".**

**(I2): positive line blocks and their exact cost (PROVED).** Fix n ≥ 1 and b_j := sin(π(j + 1)/(n + 2)), j = 0..n. Scale β_j := κ_n b_j so that Σ_j β_jβ_{j+1} = c₂. Let ρ ≥ 0 be bounded with compact support and ∫ρ = 1. Define
> **A_n := Σ_{j=0}^{n} (−1)^jβ_j V(2)^j,  D_n^{(2)} := A_n M_ρ A_n^* ≥ 0.**
- *Diagonal sums.* V(2)^jM_ρV(2)^{−k} has d = (∫ρ)δ_{(j−k)a}, so
  > d_{D_n^{(2)}} = c₂ sec(π/(n + 2))·δ_0 − c₂(δ_a + δ_{−a}) + Σ_{2≤|m|≤n}(·)δ_{ma}.
  - Σ_jb_j² / Σ_jb_jb_{j+1} = sec(π/(n + 2)): b is the Perron vector of the path graph on n + 1 vertices, with eigenvalue 2cos(π/(n + 2)).
  - The atoms at |m| ≥ 2 sit at |m|a ≥ 2a > 2L, so they are invisible on 𝒯_L.
- *Trace.* For real f ∈ 𝒯_L,
  > Tr(D_n^{(2)}V(f)V(f)^*) = ‖M_ρ^{1/2}A_n^*V(f)‖²_HS = c₂sec(π/(n + 2))·g(0) − c₂(g(a) + g(−a)) ≥ 0.
  - n = 1 is the two-site block, at cost 2c₂·g(0).
- *Sharpness: the infimum cost is exactly c₂·g(0), and it is not attained by these blocks.*
  - Any D ≥ 0 whose diagonal sums agree on (−2L, 2L) with c₀δ_0 − c₂(δ_a + δ_{−a}) has c₀ ≥ c₂. Test f = 1_{[−L, −L+ε]} + 1_{[−L+a, −L+a+ε]} with ε ≤ ℓ: then g(0) = 2ε and g(±a) = ε.
  - c₀ = c₂ is the limit n → ∞.
  - This is the trace-picture form of Celestial's Fejér completions (C.14.1: min p_w = −c₂sec(π/(n + 2)) ↓ −c₂). Same constant, same mechanism.
- **The (I2) choice.** Take **D_L = D_L^{∞,(n)} + D_n^{(2)}**, where D_L^{∞,(n)} carries the archimedean and pole diagonal sums **minus c₂sec(π/(n + 2))δ_0**. Then (I1) holds exactly, and (I2) follows from D_L^{∞,(n)} ≥ 0.
- **What this choice costs (PROVED necessary condition).**
  - On the pole-free subclass (f = (½ + ∂)ξ-type, g = Qh), g(0) = ‖ξ‖²_H. So D_L^{∞,(n)} ≥ 0 forces
    > **R(ξ) ≥ c₂sec(π/(n + 2))·‖ξ‖²_H for all ξ ∈ 𝒯_L.**
  - This is an S2-type bound in the **global** H-norm, strictly stronger than S2, whose norm is ‖·‖_{H(E)}. If it held, the p = 2 step would be trivial: S1 gives |2c₂B| ≤ c₂‖ξ‖²_{H(E)}.
  - Its truth is a **GAP**. Even: C.14.3, the rank-one κ^{ev} > 1 direction. Odd: its pointwise-symbol sufficient version is Thm C.14.A, which indicatively fails (P.15.3). Expected to fail: HEURISTIC (C.14.4b near-null directions).
- **Why no trace-picture block can do better (PROVED).** Tr(M_w V(f)V(f)^*) = g(0)∫w for every w. The trace pairing sees only diagonal sums, which are translation-invariant. So a block can be made positive only by borrowing g(0)-mass (global ‖ξ‖²_H), never strip-localized mass ‖ξ‖²_{H(E)}, the S1 localization that makes the form picture work.
  - **So for Queen:** a separable (I2), meaning a positive line block plus a positive archimedean remainder, implies the multiplier-frame bound (global S2). Only ⇒ is proved, and only for line blocks whose diagonal sums sit at {0, ±a} (REMARK; S2 stays GAP; adopted S103h(2)). A working D_L must couple the line block to the archimedean block through non-multiplier pieces: the ε-part 2ε′(1⁺)(1 − K_L) of C.14.1, band-split as in S&M Thm 2.1. This is a REMARK, PROVED as a statement about separable choices.
- *Label summary:* weight-1 block, spectrum, HS: PROVED. Positive blocks D_n^{(2)} and the sharp cost c₂: PROVED (elementary). Separable (I2) ⇒ global S2: PROVED. Its truth: GAP (HEURISTIC: false). *(Grok audit S103h(2): the ⇒ is PROVED as stated in the preceding bullets; any "⇔" reading, e.g. in P.15.0 or the bridge note, is a REMARK only, and S2 (global or strip) remains GAP.)*

### P.15.3 Status update for §P.14 / G.4h after C.14, and what the C.14.A failure means for piece 4 — labels as marked [S&M]/[Cel]
- **C.14 is merged** (a4153a8, audit c7c94ee) and is no longer pending.
  - Thm C.14.A is PROVED as a sufficient criterion ("iff" corrected to sufficient in the audit).
  - Thm C.14.B (coupling bound) and Thm C.14.C (obstructions) are PROVED.
  - **Conjecture C.14.K (high-band almost-orthogonality) is a CONJECTURE.** *(S103j status: CERTIFIED false at (ℓ, T) = (log 2/6, 14); ARP G.4i.)*
  - Cor C.14.D is PROVED given its inputs.
  - So piece 3 (N.G2-bulk) is a GAP reduced to C.14.K + m₀^e. P.14.7 and P.14.9 are updated accordingly.
- **Celestial's indicative check (team/celestial-cft ff3a9dc; NUMERICAL/indicative, hers).**
  - Inputs: 1 − κ^{odd}_{L₀} = 0.58 (Golden Eye sinc-kernel proxy, not proved), so 2ε′(1⁺)·0.58 = 26.68.
  - With σ_S ≥ max(0, 2θ′ − B/(¼ + t²)), the C.14.A symbol inequality fails on |t| ∈ [7.37, 11.78] (B = 46, worst deficit 0.222) and on [7.37, 12.96] (B = 66, deficit 0.267).
  - That is, it **fails around 7.4 ≲ |t| ≲ 13**: the ε-term 26.68/(¼ + t²) drops below c₂ beyond 7.37, and 2θ′ < c₂ below 10.26.
  - *Superseded (main b104678; NUMERICAL).* With the actual K_L, δ_o⁰ = 0.3135 (Prime N.14.1), so the ε-term is 14.42/(¼ + t²). The failing band widens to **≈ [5.41, 13.6]**, with worst deficits 0.345 (B = 46) and 0.369 (B = 66). The 0.58 proxy overestimated the gap about 1.85×.
  - *New certified input (S&M CERT_p2_constants, merged 9d72bbb; CERTIFIED, Arb):*
    - odd (K2) holds on 2L ∈ (log 2, log 2 + 0.028], CC-free, with Q_L ≥ 0.0126‖ξ‖²_H;
    - m₀^odd(L₀) ≥ 0.46255 (H-margin, O(1)).
    - The frame route is off the odd critical path on that sub-window.
- **What it means for piece 4 as the only payer of the odd cost 2c₂α²(1/ℓ − ℓ/24).**
  1. *What fails is the ℓ-independent, globally charged route, not piece 4 itself (PROVED, bookkeeping).*
     - C.14.A charges the line by its Fejér minimum −c₂ at **every** frequency. By P.15.2 this global charge is all a translation-invariant (multiplier or trace-picture) block can do.
     - The indicative failure says the margin available from known σ_S bounds plus the L²-gap is short of c₂ on 7.4 ≲ |t| ≲ 13.
     - It does not say the odd sector is unpaid. C.14.A is sufficient only, and the true σ_S may exceed c₂ − 26.68/(¼ + t²) there (≈ 0.22 near t ≈ 10; ≈ 0.35 with the actual K_L, 14.42). That is Celestial's precise S&M ask. *(Grok audit S103h(2): c₂ − 26.68/(¼ + t²) at t = 10 is ≈ 0.22, not 0.27 (0.267 is Celestial's B = 66 deficit). With the actual K_L (14.42) the needed σ_S there is ≈ 0.35.)*
  2. *The cost itself lives above the deficit band (PROVED, from C.14.4b's bound).*
     - The odd charge sits on the ramps. By C.14.4b, the fraction of a ramp's energy in |t| < T is ≤ (T/4 + T³/3)ℓ³/(3π). At T = 13 this is ≤ 78.1ℓ³: ≤ 0.132 at ℓ = 0.119 and ≤ 0.0098 at ℓ = 0.05.
     - The ramp kinks sit at |t| ~ π/ℓ ≥ 26, where 2θ′ > 2c₂.
     - So for small ℓ the odd cost can be charged in the high band, and only the bulk low band needs the L²-gap δ_o⁰.
  3. *Consequence (PROVED as bookkeeping; values GAP).*
     - Piece 4 remains the **only** payer of the odd cost. P.14 shows there is no 2-adic payer; P.15.1 shows none appears from p = 3 either, since the {2,3} odd entry is worse; and Cor C.14.D shows the odd sector needs no coupling bound, so piece 3 is not on the odd critical path.
     - On 2L ∈ (log 2, log 2 + 0.028] piece 4 **is paid**, CERTIFIED (S&M CERT_p2_constants: S1 + (Ω5) + Lemmas A/B, a band-type argument, CC-free). That is a band-split payment, consistent with this item.
     - Beyond log 2 + 0.028 the odd entry is open analytically (CERTIFIED only conditional on Mori). Piece 4 must there be paid in **band-split, ℓ-dependent form**: S&M Thm 2.1/(3.3) with δ_o⁰ (polynomial range), C.13.H′, or Thm N.13.A with m₀^o (exponential range). It cannot be paid in the ℓ-independent frame form, unless a certified σ_S clears c₂ − 2ε′(1⁺)(1 − κ^{odd}_{L₀})/(¼ + t²) on the band (≈ [5.41, 13.6] with δ_o⁰ = 0.3135, NUMERICAL).
     - For the band-split route, the threshold must be T ≳ 13 (where σ_S ≥ c₂(1 + η) is plausible). This shrinks the range through ν_T(ℓ) ∝ T³ℓ³. The range is a GAP until δ_o⁰ and σ_S are supplied.
  4. *For Queen (Q.8.E):* by P.15.2, a separable (I2) is the global-charging route, so the same indicative failure applies to its pointwise version. Q.8.E's D_L near the p = 2 entry should be built band-split, coupling the line block to the ε-part. REMARK.
- **Display fixes: none open.** Prime's erratum N.14.0(a) (656de12) fixed N.13.4's D_L bullet. The G.4h item-2 complaint is withdrawn (P.14.1, P.14.5).

### P.15.4 Labels and R-line
| item | label | track |
|---|---|---|
| Lemma P.15.A: orbit-closed nodes N_L; pieces (ℓ₃, q, ℓ₃); τ_{log 2}: I₁ → −I₃ | PROVED | [Cel] |
| joint ramp split: derivative cross terms vanish only with N_L | PROVED | [Cel] |
| Thm P.15.B: exact B₂, B₃ on joint ramps; per-line constants recovered | PROVED (sympy-checked algebra) | [Cel] |
| p = 3 entry: both parities indefinite; per-energy eigenvalues σ·½(c₃ ± √(c₃² + 4c₂²)) | PROVED | [Cel] |
| P.14.8 prediction: (a)–(c) confirmed; (d) "parity decides per line" refuted | PROVED | [Cel] |
| 0.9009 as the sharp S1-type constant at the p = 3 entry | CONJECTURE → **PROVED** in P.16.1 (c₀\*, both parities) | [Cel] |
| Q.8.E weight-1 block −(c₂/ℓ)S_Λ(V(2) + V(2)^*)S_Λ: diagonal sums, spectrum {0, ±c₂/ℓ}, HS | PROVED | [Q] |
| fixed-ρ variant uniformly bounded (R1) | REMARK (PROVED) | [Q] |
| D_n^{(2)} = A_nM_ρA_n^* ≥ 0, cost c₂sec(π/(n + 2))g(0); infimum c₂ sharp (not attained by A_nM_ρA_n^*) | PROVED | [Q] |
| separable (I2) ⇒ R ≥ c₂‖·‖²_H on 𝒯_L (global S2); ⇒ only, line blocks with diagonal sums at {0, ±a} | REMARK; S2 GAP (HEURISTIC false) | [Q]/[S&M] |
| trace picture cannot localize to H(E) | PROVED | [Q] |
| C.14 merged; C.14.K | ~~CONJECTURE~~ **CERTIFIED false** at (ℓ, T) = (log 2/6, 14) (S103j status update, ARP G.4i) | [Cel] |
| C.14.A symbol inequality fails on ≈ 7.4 ≤ \|t\| ≤ 13 (proxy 0.58); ≈ [5.41, 13.6] with actual K_L | NUMERICAL (indicative; Celestial ff3a9dc / Prime N.14.1) | [S&M] |
| cross-check of N.14.3: identities agree; "no 1/ℓ₃ term" holds in ‖ξ′‖_∞, not in H-energy | PROVED | [Cel] |
| (L4): c₄ line = m = 2 time-delay mode, trace block −(c₄/ℓ₄)S(V(2)² + V(2)^{*2})S | PROVED | [Q] |
| odd (K2) on (log 2, log 2 + 0.028]; m₀^odd ≥ 0.46255 | CERTIFIED (S&M, cited) | [S&M] |
| piece 4 sole odd payer, band-split form needed | PROVED as bookkeeping; values GAP N.G2-m₀ | [S&M] |
| RH | OPEN | — |

**R-line:**
- **R1:** one rule on (log 3, log 4]. The combinatorics of N_L are constant on the range. The (I2) blocks are L-uniform with fixed ρ, and D_L^{(2),1} is not uniform as ℓ ↓ 0.
- **R2:** the Euler product enters through the two lines (c₂, c₃), the m = 1 time-delay modes of γ₂ and γ₃ (P.14.3); the FE enters through Ω.
- **R3:** the signs come from P-parity and the indefinite 2×2 block at the p = 3 entry, not from weight signs.
- **R4:** the pole pair is spent once. Its ¼ appears as the mass terms (ℓ₃/24, q/24, …), as in P.14.4, and g(0) = ‖ξ‖²_H is the pole-free norm in P.15.2.
- **R5, where each fake breaks:**
  - **F_θ:** P.15.1–P.15.2 are linear in the line weights, so they go through with c_p → c_p cosh(θ log p). It breaks where P.14 R5 says: the line is not the self-adjoint time delay of a unitary γ_p (polarization, P.4).
  - **DH:** no Euler product and no γ_p. −D′/D has sign-indefinite weights, not supported only on prime powers. Breaks at P.14.2.
  - **H₋ₜ:** no Euler product after deformation, and Ω is deformed. Breaks at P.14.2/P.14.3.
  - **DMV / P_X (X ≥ 3):** coincide with ζ on this window. Fake-blind here.
  - **G_θ′:** passes. Q_G is a sum of two ζ-forms (G.4e), and every identity above is linear in the form.
- **R6:** Lemma P.15.A, Thm P.15.B, the trace identities and the sharp cost c₂ are exact. The closures (N.G2-m₀, N.G2-bulk, even (P)) stay GAPs.

### P.15.5 PASS/FIX verdicts on Prime §N.14 (cf4ee37), as requested by Prime

Each item was checked analytically against Thm P.15.B and P.14. Numbers are closed-form evaluations, not exploration.

| # | N.14 claim (cf4ee37) | Verdict | Label |
|---|---|---|---|
| 1 | 3-edge constant c₃(1/ℓ₃ − ℓ₃/24), with ℓ₃ = 2L − log 3 | **PASS**. It equals B₃(ξ_b) = σv₁²(ℓ₃/24 − 1/ℓ₃) in Thm P.15.B; the contribution to Q_L is +2σc₃v₁²(1/ℓ₃ − ℓ₃/24) | PROVED |
| 2 | 2-line c₂(1/ℓ₂ − ℓ₂/24), with ℓ₂ = 2L − log 2 ∈ (log 1.5, log 2] | **PASS, with a cosmetic FIX**. The formula is the single-ramp reduction of Thm P.15.B. The endpoint values are c₂(1/ℓ − ℓ/24) = **1.20053** at ℓ₂ = log 1.5 = 0.40547 and **0.69295** at ℓ₂ = log 2. N.14.3 prints 1.2003 and 0.6929; they should read 1.2005 and 0.6930 (4 d.p.). "Comparison only" is right: at these ℓ₂ the 2-ramp is not an edge object | PROVED |
| 3a | B₂(ξ_s) = 0 for ξ_s ∈ H¹₀(3-strips), 2L < log 4 | **PASS**. τ_{log 2} carries each 3-strip into the bulk (I₁ ↦ −I₃), never onto a 3-strip | PROVED |
| 3b | the 2-line couples the 3-ramps only to the bulk | **PASS** as an identity (the interface form is exact) | PROVED |
| 3c | "at order 1/ℓ₃ the {2, 3} edge form is the 3-line block alone" | **FIX (wording)**. The derivative term (2β₋/ℓ₃)[ξ(L + a − b) − ξ(a − L)] = σ(2v₁/ℓ₃)(v₂ − v₃) is O(\|β\|‖ξ′‖_∞). Relative to the H-energy it is the same order 1/ℓ₃ whenever ξ has an ℓ₃-scale kink on the central piece. The joint block is then σ(1/ℓ₃)[[2c₃, 2c₂], [2c₂, 0]], with eigenvalues per unit energy σ(+0.9009, −0.2666). Suggested wording: "at order 1/ℓ₃ in ‖ξ′‖_∞ (bounded bulk slope)". Note too that ‖ξ′‖_∞ is not finite on all of H¹ | PROVED |
| 4a | only P survives (citing P.14) | **PASS** (P.14 Prop on uniqueness of P and no common J_p) | PROVED |
| 4b | the {2, 3} odd-entry prediction stays CONJECTURE | **PASS**. Add: the odd payer must cover the joint per-energy charge 0.9009 (not c₃ = 0.6343). It sits at frequency ~π/ℓ₃, i.e. high band | CONJECTURE (REMARK on the size) |
| 5a | −log 2·Σ_{m≠0} r^{\|m\|}e^{imθ} = −log 2·(P_r(θ) − 1), r = 2^{−1/2} | **PASS**. Σ_{m∈ℤ} r^{\|m\|}e^{imθ} = (1 − r²)/(1 − 2r cos θ + r²), the unnormalized Poisson kernel. The coefficients log 2·r^{\|m\|} = Λ(2^m)2^{−\|m\|/2} are the m-modes of (P14.1) | PROVED |
| 5b | KMS inverse tridiagonal, symbol ∝ \|1 − re^{iθ}\|² | **PASS**. KMS_n^{−1} = (1 − r²)^{−1}·tridiag(−r; 1, 1 + r², …, 1 + r², 1), with interior symbol 1 − 2r cos θ + r² = \|1 − re^{iθ}\|² (classical) | PROVED |
| 5c | 2×2 section − I = rσ_x; −log 2·r = −c₂; det KMS₂ = 1 − r² = 1 − 1/p | **PASS** (p = 2: 1 − r² = ½). REMARK: with r_p = p^{−1/2} this is the local factor ζ_p(1)^{−1}. That is an identity only, and no Tate reflection follows (P.14.4). The "edge coefficient" match is at the multiplier level; the ramp factor (1/ℓ − ℓ/24) and the sign of the Q_L contribution come from the ramp form (P.14.3), not from the multiplier | PROVED |
| 5d | at log 4 the three-site chain is the 3×3 KMS section, with c₄ = log 2·r² | **PASS**. The corner entry r² gives log 2·r² = (log 2)/2 = c₄ = Λ(4)/√4, the m = 2 mode. The edge ℓ₄ = 2L − log 4 agrees with P.15.1 (L4) | PROVED |

**Summary.** 9 PASS. 2 FIX, both wording only: the rounding in item 2, and the norm qualifier in 3c. Neither changes a theorem of N.14. The interface form, S1 (as an inequality), (3c), (S), (L4), (ζ) and the identities in item 5 are confirmed. *(S103i: S1's sharpness claim "c₂ + c₃ cannot be improved" is a further FIX. On 𝒯_L the sharp constant is c₀\* = 0.90093, by P.16.1.)* R1–R6 as in P.15.4: every check here uses only the {2, 3} line weights and the strip geometry. F_θ, DH, H₋ₜ and DMV fail the identities that use log 2/log 3 ∉ ℚ and Λ(4) = log 2, and G_θ′ passes. RH OPEN.

## P.16 The band-split payer at the p = 3 entry, sector ranges, and the coupled (I2) operator for Q.8.E (S103i; analytic first; [Cel]/[S&M]/[Q])

Daniel Toupin / Golden Physics Project, with AI assistance from Grok. **RH is OPEN. Nothing in P.16 proves or implies RH.**

**Notation.** As in P.15:
- a = log 2, b = log 3, c₂ = a/√2, c₃ = b/√3.
- 2L = b + ℓ₃ with 0 < ℓ₃ < log(4/3), so ℓ₂ = 2L − a = (b − a) + ℓ₃.
- Q_L = R − 2c₂B₂ − 2c₃B₃ on 𝒯_L, with B_p(ξ) = ⟨ξ, τ_{log p}ξ⟩_H and band parts ξ^{<T}, ξ^{≥T}.
- **c₀\* := ½(c₃ + √(c₃² + 4c₂²)) = 0.900927616…** (the P.15.1 joint value).

**Sources** (main 9cdc498):
- S&M SECTOR_COERCIVITY_p2 Thm 2.1, and CERT_p2_constants Thms 1–3 with their certified Ω(12.875) and Ω(16).
- Prime §N.12–N.14: S1, (Ω1)–(Ω5), N.13.A.
- Queen ARP_queen_Q8 Q.8.E/F.

**Verification.** `team/principal-series/scripts/p15_p16_check.py` (log `p15_p16_check.log`) covers:
- the P.15 joint block (Thm P.15.B, exact, sympy);
- the P4 path identity and the fibre combinatorics of P.16.1;
- Arb enclosures (python-flint 0.9.0, 128 bits) of Ω at the thresholds used below.

**Lead corrections adopted (S103h(2)).**
- "Separable (I2) ⇒ global S2" is proved only in the ⇒ direction, for line blocks whose diagonal sums sit at {0, ±a}. It is a REMARK, and S2 stays a GAP. P.15.2, P.13 item 15 and the P.15.4 row are reworded accordingly.
- In P.15.3, c₂ − 26.68/(¼ + t²) at t = 10 is ≈ 0.22, not 0.27. With the actual K_L (14.42) it is ≈ 0.35.

### P.16.1 The sharp joint S1 constant on (log 3, log 4) is c₀\*, not c₂ + c₃ — PROVED (elementary; algebra sympy-checked) [Cel]
**Prop P.16.1 (fibre decomposition).** Shift the window to [0, w), with w = 2L.
- **Graph.** Join x ~ x + a when x ∈ [0, w − a), and x ~ x + b when x ∈ [0, w − b) = [0, ℓ₃). For b < w < 2a, every point has at most one a-neighbour and at most one b-neighbour (w < 2a < 2b). The components are therefore paths. Explicitly (sympy check, part 3):
  - (i) For x ∈ [0, ℓ₃), the component is the 4-path **(x + a) —a— x —b— (x + b) —a— (x + b − a)**, with weights (c₂, c₃, c₂). The constraints used are ℓ₃ < b − a, a + ℓ₃ < b and b − a + ℓ₃ = w − a.
  - (ii) The remaining a-pairs, x ∈ G := [ℓ₃, b − a), with x ~ x + a, are 2-paths of weight c₂. (The points [b − a, ℓ₂) are 4-path ends.)
  - (iii) All other points are isolated.
- **Node sets in window coordinates.** The 4-path nodes form F₃ = (3-strips) ∪ (the two central pieces [a − L, a − L + ℓ₃] and its mirror), with |F₃| = 4ℓ₃. These are exactly the orbit-closed pieces of P.15.1. The 2-path nodes form G₂ = E₂ ∖ F₃.
- **Consequence.** Let K := c₂(τ_a + τ_{−a}) + c₃(τ_b + τ_{−b}), so that 2c₂B₂ + 2c₃B₃ = ⟨ξ, Kξ⟩_H = ⟨ξ′, Kξ′⟩_{L²} + ¼⟨ξ, Kξ⟩_{L²}. Its compression to L²[0, w) is the direct integral of the weighted adjacency matrices of the components. The 4-path spectrum is {±0.90093, ±0.26664} (part 2), and the 2-path spectrum is {±c₂}. Hence, for every ξ ∈ 𝒯_L:
  > **|2c₂B₂(ξ) + 2c₃B₃(ξ)| ≤ c₀\*‖ξ‖²_{H(F₃)} + c₂‖ξ‖²_{H(G₂)} ≤ c₀\*‖ξ‖²_{H(E₂)}.**  (P16.1)
- **Sharpness, in both parities.** Put one bump φ ∈ H¹₀ of width < ℓ₃ at each node of one 4-path, with Perron weights (1, 1.838, 1.838, 1). The Rayleigh quotient is then exactly c₀\* (part 4). P maps the 4-path of x to the 4-path of ℓ₃ − x, a disjoint component. So ξ ± Pξ attains c₀\* in each parity, for both signs (charge or help).

**Remarks.**
- **FIX to N.14.3 and to my P.15.5.** N.14.3 says that c₂ + c₃ = 1.1244 "cannot be improved by phase alignment". The Kronecker argument concerns the symbol on ℝ. On 𝒯_L, the sharp constant against ‖ξ‖²_{H(E₂)} (and against ‖ξ‖²_H) is c₀\* = 0.90093. The inequality S1 ≤ (c₂ + c₃) is still true (PASS); only the sharpness claim needs fixing. P.15.5's "S1 … confirmed" is amended accordingly.
- **P.15.1 is a special case.** The ramp-space eigenvalues σ(0.9009, −0.2666) are the 4-path spectrum. The ramp ansatz ties the 4-path of x to its mirror, which is why P.15.1 saw only 0.2666 as the even charge there. The global per-energy charge is **c₀\* in both parities** (PROVED). The P.15.4 row "0.9009 sharp S1-type constant" is upgraded from CONJECTURE to PROVED.
- **Range.** The constant is constant on all of (log 3, log 4). At 2L = log 4 the a-gap [w − a, a) closes, paths grow, and (L4) begins.
- *(Grok audit S103h(6): P.16.1 checked by hand: the fibre graph (x+a)–x–(x+b)–(x+b−a) and the endpoint inequalities hold, and the symmetric weighted 4-path (c₂, c₃, c₂) has spectrum ±½(c₃ ± √(c₃² + 4c₂²)). Script part 3 replayed, so PROVED. The sup is **attained**, not just approached, by disjoint translates of one H¹₀ bump at the four nodes with Perron weights (1, 1.838, 1.838, 1), provided that combination lies in 𝒯_L. If 𝒯_L carries extra linear constraints, the value is approached by adding small corrections and the sharpness is unchanged. P sends the path of x to the path of ℓ₃ − x, so ξ ± Pξ gives both parities, away from x = ℓ₃/2. The FIX to N.14.3 is right: c₂ + c₃ is sharp for the symbol on ℝ (Kronecker), but not on the finite-window fibres.)*

### P.16.2 The band-split payer — PROVED (conditional where stated); thresholds CERTIFIED [S&M]
Write K = K_F + K_G (exact): K_F is the 4-path part and K_G the 2-path part. Define the **3-entry reduced form**
> Q^G_L(ξ) := R(ξ) − ⟨ξ, K_Gξ⟩_H, so that Q_L = Q^G_L − ⟨ξ, K_Fξ⟩_H exactly, and |⟨ξ, K_Fξ⟩_H| ≤ c₀\*‖ξ‖²_{H(F₃)}.

As ℓ₃ ↓ 0, F₃ shrinks to measure 0 and Q^G_L tends to the full {2}-form at 2L = log 3.

**(Ω5′) (PROVED, as in N.13.2).** For measurable F and the band u = ξ^{<T}: ‖u‖²_{H(F)} ≤ |F|(T/π)‖u‖²_H. The pointwise Cauchy–Schwarz bound for |u′|² + ¼|u|² is (T/π)‖u‖²_H. With |F₃| = 4ℓ₃ this gives (4ℓ₃T/π).

**Thm P.16.A (high-band payer for the joint charge; PROVED; constants CERTIFIED).**
- Suppose Ω(T) ≥ c₀\*(1 + η). Then the high-band part of the joint 4-path charge is covered by R alone:
  > R(ξ^{≥T}) − c₀\*(1 + η)‖ξ^{≥T}‖²_H ≥ [Ω(T) − c₀\*(1 + η)]‖ξ^{≥T}‖²_H ≥ 0.
- Arb (part 5):
  - Ω(15.5) − c₀\* ∈ [0.0018618 ± 7e−12] > 0 (Principal's own Arb value, p15_p16_check part 5; not an S&M certificate; at T = 15.5 it allows only η ≤ 2.07e−3);
  - Ω(16) = 0.934548783922(±4e−13), the same value as S&M's certified Ω(16) = 0.93455, so **η ≤ 0.0373183 at T = 16**;
  - Ω(19.35) = 1.124704 (η ≤ 0.2484);
  - Ω(38.1) = 1.802309 (η ≤ 1.0005).
- By (Ω2), PROVED, each inequality holds for all |t| ≥ T. **CERTIFIED.**
- **Usable certified point (S103j wording fix): T = 16 with η ≤ 0.0373** (S&M's certified Ω(16)). The T = 15.5 statement rests on Principal's own Arb value only and gives η ≤ 2.07e−3; it is not used below.
- *(S103j: at T = 16 the high-band part of Thm P.16.A is CERTIFIED, but the companion low-band hypothesis (H₃^σ)(16, 0.0373, μ) is NUMERICALLY false for every μ ≥ −1, in both parities (P.17.0). So the T = 16 and T = 19.35 rows of the P.16.B table are vacuous as they stand.)*
- By contrast, the old constant c₂ + c₃ needs T ≥ 19.35 even at η = 0, since Ω(19.35) − (c₂ + c₃) = 2.9e−4 and Ω(16) < c₂ + c₃. The sharp c₀\* lowers the high-band threshold to S&M's certified T = 16.
- *(Grok audit S103h(6): the script was re-run and its output is byte-identical to the log. The enclosures are genuine Arb: python-flint acb.digamma at 128 bits, with T given as exact decimals. Ω(15.5) − c₀\* ∈ [0.0018617842 ± 7e−12] and Ω(16) = [0.934548783922 ± 3.4e−13], which agrees with S&M's certified 0.93455. Precision: T = 15.5 lies **below** S&M's Ω(16), so the T ≥ 15.5 statement rests on Principal's own Arb value, not on S&M's. At T = 15.5 it allows only η ≤ 2.07e−3, and Thm P.16.B needs η > 0. The usable certified point is T = 16 with η ≤ 0.0373.)*
- Ramp-scale content, at frequency ~π/ℓ₃, lies above T = 16 when ℓ₃ ≲ π/16 = 0.196 (REMARK).

**Thm P.16.B (3-entry, polynomial range; PROVED given (H₃^σ)).**
- **Hypothesis (H₃^σ)(T, η, μ)** [GAP]: for all ξ ∈ 𝒯_L of parity σ,
  > Q^G_L(ξ) ≥ c₀\*(1 + η)‖ξ^{≥T}‖²_H + μ‖ξ^{<T}‖²_H.
- **Conclusion.** Then Q_L(ξ) ≥ [μ − c₀\*(1 + 1/η)(4ℓ₃T/π)]‖ξ^{<T}‖²_H ≥ 0, provided
  > **ℓ₃ ≤ ℓ₃\*(σ) := πμ/(4c₀\*(1 + 1/η)T).**
- *Proof.* Apply (P16.1) to K_F, then ‖x + y‖² ≤ (1 + η)‖y‖² + (1 + 1/η)‖x‖² on H(F₃), then ‖ξ^{≥T}‖_{H(F₃)} ≤ ‖ξ^{≥T}‖_H and (Ω5′). ∎ This is S&M Thm 3's chain, with R replaced by Q^G_L, 2ℓ by 4ℓ₃ and c₂ by c₀\*.
- **Evaluations** (the Ω values are CERTIFIED; the rest is eval.):

  | T | η (max with Ω(T) ≥ c₀\*(1 + η)) | ℓ₃\* / μ |
  |---|---|---|
  | 16 | 0.0373 | 1/510 |
  | 19.35 | 0.248 | 1/112 |
  | 38.1 | 1.0 | 1/87 |

  The high-band part of (H₃^σ) is consistent with R exactly when Ω(T) ≥ c₀\*(1 + η) (Thm P.16.A). The content of (H₃^σ) is that the low band of Q^G_L has a margin μ after the 2-path line is subtracted.

**Thm P.16.C (3-entry from a plain margin; PROVED given (M₃^σ)).**
- **Hypothesis (M₃^σ)** [GAP]: Q^G_L(ξ) ≥ m‖ξ‖²_H on 𝒯_L^σ, with 0 < m ≤ 1.
- **Setup.** Let 1 − θ := m/(2(|Ω₀| + c₂) + m), where |Ω₀| = 5.37218 (Ω1). Choose T with Ω(T) ≥ c₂ + c₀\*(1 + η)/(1 − θ). The crude (Ω1) bound Ω ≥ Ω₀ + log 2T allows T = ½exp(|Ω₀| + c₂ + c₀\*(1 + η)(2(|Ω₀| + c₂) + m)/m).
- **Conclusion.** Q_L ≥ 0 on parity σ for
  > **ℓ₃ ≤ πθm/(8c₀\*(1 + 1/η)T).**
- *Proof.* Write Q^G = θQ^G + (1 − θ)Q^G.
  - (1 − θ)Q^G ≥ (1 − θ)[(Ω(T) − c₂)‖ξ^{≥}‖² − (|Ω₀| + c₂)‖ξ^{<}‖²], using |⟨ξ, K_Gξ⟩_H| ≤ c₂‖ξ‖²_H (2-path fibres) and the fact that R is diagonal in t.
  - θQ^G ≥ θm‖ξ‖².
  - Bound K_F as in Thm P.16.B. With this θ the low bracket is θm/2 − c₀\*(1 + 1/η)(4ℓ₃T/π), and the high bracket is ≥ 0 by the choice of T. ∎
- The range is exponentially small in 1/m, as in N.13.A: the e^{−c/m} loss comes from the crude low-band Ω₀. Thm P.16.B is the useful form.

**Why S1 + (Ω5) alone cannot pay the 3-entry (REMARK; the logic is PROVED).**
- The 2-path part lives on G₂, with |G₂| = 2(ℓ₂ − 2ℓ₃) ≥ 2(b − a) − 2ℓ₃ > 0.23, which is O(1). (Ω5′) gives no smallness there.
- So the low band of the 2-line must be paid by its signed structure, which is the p = 2 mechanism carried to log 3. That is exactly what (H₃^σ) and (M₃^σ) assume. This agrees with N.12.4 (S2 cannot reach log 3).
- The band-split payer therefore pays **only the new 3-line increment**, which is the joint 4-path charge c₀\*. It covers that charge fully on the high band (Thm P.16.A, CERTIFIED). The low band costs O(ℓ₃) (Thm P.16.B).

### P.16.3 Which ℓ₃-range each sector is paid on — labels as marked [S&M]
| sector | unconditional range | conditional range | what is missing |
|---|---|---|---|
| odd | **none** (the S&M odd certificate stops at 2L = log 2 + 0.028, far below log 3; indicative ℓ₃\* ≲ 1.5e−5, NUMERICAL, P.16.5) | ℓ₃ ≤ πμ_o/(4c₀\*(1 + 1/η)T), given (H₃^o)(T, η, μ_o): Thm P.16.B, PROVED given input. High band CERTIFIED for T ≥ 16 (η ≤ 0.0373). Also ℓ₃ ≤ πθm_o/(8c₀\*(1 + 1/η)T) given (M₃^o): Thm P.16.C | μ_o or m_o for Q^G_L near 2L = log 3: **GAP**. Q at 2L = log 3 is ≥ 0, CERTIFIED conditional on Mori (main G.4); no H-margin is known |
| even | **none** | the same formulas with (H₃^e)/(M₃^e) | μ_e/m_e: **GAP**, plus the even ĝ(0)/pole direction (C.14.3, rank-one κ^{ev} > 1) and S&M GAP m₀^e inherited from the p = 2 step |

- **Parity does not change the line side at the 3-entry.** The sup of the joint charge is c₀\* in both sectors (P.16.1, sharp). Every sector difference sits in the low-band hypothesis on Q^G_L. **REMARK (PROVED).**
- **Nothing here is NUMERICAL.** The Ω enclosures are Arb (CERTIFIED), the constant c₀\* is exact, and the ℓ₃\*/μ ratios are evaluations of closed forms. No margin for Q^G_L has been computed. A future replay of μ_σ would be NUMERICAL unless it carries an Arb remainder. It would need Lemma A/B for a form that is **not** a multiplier (Q^G_L contains K_G), so S&M's bathtub does not apply verbatim. **This is the S&M ask.**

### P.16.4 The coupled (I2) operator for Q.8.E, in Queen's notation — PROVED where stated [Q]
Notation as in P.15.2: Λ = e^L, S_Λ = 1_{[−L, L]}, V the translation group, V(f) = ∫f(u)V(u)du, g = f ∗ f^*, and d_D(b) = ∫D(x, x + b)dx. We work on the pole-free class (f = Qξ, |f̂|² = (¼ + t²)|ξ̂|², g(0) = ‖ξ‖²_H), where (N3) forces Res ≡ 0.

- **(a) The (I1) family (PROVED).** Let n_L be any real even distribution with supp n_L ∩ (−2L, 2L) = ∅ ("invisible"). Define the symbol
  > **s_L(t) := σ_S(t) + 2ε′(1⁺)(1 − κ̂_L(t))/(¼ + t²) − 2c₂cos(at) − 2c₃cos(bt) + n̂_L(t),**
  - Here σ_S(t) + 2ε′(1⁺)(1 − κ̂_L)/(¼ + t²) = Ω(t) is the S&M (5.2) split of the archimedean symbol into its 𝒮-part and its ε-part. K_L is the compression of the convolution by κ_L = q/(2ε′(1⁺)).
  - Let A_L be the Fourier multiplier with |â_L|² = s_L (taken where s_L ≥ 0), ρ ≥ 0 a fixed L-independent bounded density with compact support and ∫ρ = 1, and
  > **D_L := A_L M_ρ A_L^\*.**
  - Then d_{D_L} = (∫ρ)·š_L, and **Tr(D_L V(f)V(f)^\*) = (1/2π)∫s_L|f̂|² = Q_L(f)** for real f ∈ 𝒯_L, since n_L is invisible on supp g ⊂ [−2L, 2L]. So (I1) holds.
  - *(S103j wording fix.)* "Invisible" is used in the **weak sense**: (¼ − ∂²)n_L is supported off (−2L, 2L). Then n_L pairs to zero with every pole-free g (g = ξ∗ξ\*, ξ ∈ H¹₀(−L, L)), which is all (I1) needs; supp n_L itself may meet (−2L, 2L), since (¼ − ∂²)⁻¹ is non-local. In the weak sense n_L = αe^{x/2} + βe^{−x/2} on (−2L, 2L) (PROVED: the kernel of ¼ − ∂² on the interval).
  - The earlier sentence "if s_L is bounded, D_L is bounded" is **vacuous** and is withdrawn: s_L is never bounded in this ansatz (REMARK P.17.5(e), PROVED). So D_L = A_LM_ρA_L^\* is an unbounded operator, and Queen's boundedness requirement in (I1) needs a non-multiplier D (GAP).
- **(b) The (I2) criterion (PROVED).**
  - The construction is available, and then D_L ≥ 0, exactly when s_L ≥ 0 a.e. (|â_L|² = s_L).
  - Conversely, every D ≥ 0 has positive-definite diagonal sums: d_D is a sum of autocorrelations Σλ_k e_k ⋆ ẽ_k. So any (I1) + (I2) solution has d̂_D = s_L for some invisible n_L.
  - Existence of an invisible n_L with s_L ≥ 0 is **equivalent to Q_L ≥ 0 on the pole-free 𝒯_L**, with s_L read as a measure (below). The extension is not explicit, so R6 explicitness remains the GAP of Q.8.E.
  - **(S103j restatement of the ⇐ direction; PROVED by citation.)**
    - Apply the extension theorem to d′ := (¼ − ∂²)d, not to d. d′ is a positive-definite distribution on (−2L, 2L): ⟨d′, ξ∗ξ\*⟩ = Q_L(ξ) ≥ 0 for all ξ ∈ C_c^∞(−L, L), and this class is dense in the pole-free 𝒯_L.
    - d′ has atoms (from the line terms) and a log-singular archimedean part, so the continuous (Krein 1940) and measurable (Sasvári 1994, Thm 4.2.3, Artjomenko) versions do not apply. The distributional version used is: *a positive-definite (tempered) distribution on a finite open interval has a positive-definite extension to ℝ* — Jorgensen, Pedersen & Tian, "Positive definite (p.d.) functions vs p.d. distributions", arXiv:1712.07264, **Corollary 2.2** (statement read and verified this session; its proof rests on their monograph *Extensions of Positive Definite Functions*, LNM 2160, Springer 2016, theorem number there unverified). Schwartz's Bochner theorem then writes the extension as the Fourier transform of a **positive tempered measure** σ_L. Alternatives: Krein–Langer, IEOT 78 (2014) 1–69, Thm 4.4 is the *continuous* case (verified via the MaRDI review); its §5 screw-function continuation fits W = −g″ with g continuous (theorem number unverified). Gel'fand–Vilenkin vol. 4, Ch. II (p. 197 attributes the even-function representation to Krein; exact theorem number unverified).
    - **Measure form of the identity.** The conclusion is σ_L ≥ 0 with σ̌_L = d′ on (−2L, 2L). Then n_L is any tempered solution of (¼ − ∂²)n_L = σ̌_L − d′ (with d′ taken on all of ℝ). The right side is supported off (−2L, 2L), so n_L is weakly invisible, and it is unique up to the kernel αe^{x/2} + βe^{−x/2} on the interval. The old "|â_L|² = s_L" is the **density** statement 2π dσ_L = (¼ + t²)s_L dt, which holds only when σ_L is absolutely continuous. In general (¼ + t²)s_L dt must be read as the measure 2π dσ_L, and Q_L(ξ) = ∫|ξ̂(t)|² dσ_L(t) for ξ ∈ C_c^∞(−L, L) (same normalisation as the trace identity in (a)) is the identity of record.
    - **REMARK (PROVED, elementary).** If D ≥ 0 has Σλ_k‖e_k‖²₂ < ∞, then d̂_D = Σλ_k|ê_k|² ∈ L¹, which is absolutely continuous. So a singular part of σ_L cannot be realized by such a D, and the A_LM_ρA_L^\* realization needs an absolutely continuous extension. Whether one exists (extensions are not unique) is a GAP.
  - *(Grok audit S103h(6), Krein step:*
    - *⇒ is trivial: D ≥ 0 gives Tr(DV(f)V(f)\*) ≥ 0.*
    - *For ⇐, apply Krein to d′ := (¼ − ∂²)d, i.e. the symbol (¼ + t²)s, which is positive-definite on (−2L, 2L) against all ξ∗ξ\* with ξ ∈ H¹₀(−L, L). The resulting n is then invisible in the weaker sense that (¼ − ∂²)n is supported off (−2L, 2L), meaning it pairs to zero with every pole-free g. Its support need not avoid (−2L, 2L), since (¼ − ∂²)⁻¹ is non-local. That suffices for (I1) on the pole-free class, but the definition in (a) should say so.*
    - *Since d has atoms at ±log 2 and ±log 3 and a log-singular archimedean part, the distributional form of the extension theorem is needed. Krein's 1940 theorem is for continuous functions. Sasvári 1994, Ch. 4 §4.2 "The one dimensional case", has it for measurable functions (Thm 4.2.3, Artjomenko). Cite a distributional version explicitly, e.g. Gel'fand–Vilenkin vol. 4 or Krein–Langer.*
    - *The extension is in general a measure, so |â_L|² = s_L needs s_L to be a function. Otherwise read D_L ≥ 0 via a spectral measure.*
    - *Label kept: PROVED by citation, with the citation to be sharpened.)* *(S103j: sharpened as above.)*
- **(c) The line part and its minimal cost (PROVED).**
  - Put the line atoms with an invisible line correction n^F: λ(t) := −2c₂cos(at) − 2c₃cos(bt) + n̂^F(t). Any such choice has
  > **sup_t(−λ(t)) ≥ c₀\* = 0.90093** on (log 3, log 4).
  - *Proof.* Test f = Perron bumps on one 4-path (P16.1): (1/2π)∫λ|f̂|² = −c₀\*g(0). ∎
  - The separable Fejér choice D_n^{(2)} + D_n^{(3)} (P.15.2, both with fixed ρ) costs (c₂ + c₃)sec(π/(n + 2)) → 1.1244.
  - Whether some coupled invisible n^F reaches c₀\* is open. Atoms at ma + nb are invisible iff |ma + nb| ≥ 2L, and the visible set {(m, n) : |ma + nb| < 2L} is an infinite strip in ℤ², so this is a 2-D Fejér–Riesz problem with a strip constraint: **GAP / CONJECTURE (yes)**.
- **(d) The coupling Queen needs (PROVED as a necessary-structure statement; REMARK otherwise).**
  - Split n_L = n^F + n^{lo} and ask each part to fix its own term: line part ≥ −c₀ and archimedean + ε part ≥ c₀. That is a separable choice, and (since n^{lo} is invisible) it forces R(ξ) ≥ c₀‖ξ‖²_H with c₀ ≥ c₀\*, i.e. a global S2 at the 3-entry. Per the lead's correction this is the ⇒ direction only, and S2 stays a GAP. Since R has Ω(0) < 0 directions, it is HEURISTIC-false.
  - A working D_L must therefore take **n^{lo} depending on the line coefficients and frequency-adapted**. The ε-symbol 2ε′(1⁺)(1 − κ̂_L)/(¼ + t²) is not ≥ 0 either, since sup κ̂_L ≥ sup spec K_L > 1 (CC's even eigenvalue).
  - The band-split form of the requirement is the trace-side image of Thm P.16.B: s_L(t) ≥ 0 must hold
    - on |t| ≥ T from σ_S − B/(¼ + t²) ≥ c₀\*(1 + η) (σ_S route; B is S&M GAP B), or directly from Ω(T) ≥ c₀\*(1 + η) (CERTIFIED, T ≥ 16);
    - on |t| < T only after adding an invisible n^{lo} whose existence is equivalent, by (b), to the low-band hypothesis (H₃^σ).
  - So, for Q.8.E: **D_L = A_LM_ρA_L^\* with s_L as in (a) is the coupled (I2) operator.** (I1) is exact for every invisible n_L, and (I2) ⇔ s_L ≥ 0. On every window where Q_L ≥ 0 on the pole-free class, some n_L exists (PROVED by citation), non-explicitly. The explicit rule is the GAP.

### P.16.5 Relation to Prime §N.15 (cd94c27, merged 7953486; read after the first draft of P.16) — REMARK, with labels as marked [S&M]/[Cel]
- **N.15.1** reproduces P.15.1's (0.9009, −0.2666) in the same normalization (PASS). P.16.1 shows this is the **global** sharp joint constant on 𝒯_L in both parities, not only the edge-class value. This answers N.15's flag: the per-line constants and the kink normalization are unchanged, so the 0.9009 and the 15.47 threshold stand. The Arb value is Ω(15.5) − c₀\* = 1.86e−3 > 0 (Principal's own Arb, 128 bits; the usable certified point is T = 16, η ≤ 0.0373, S103j). N.15's 15.471 is the mpmath root.
- **N.15.3 B1** (the S1-only template fails at log 3; NUMERICAL, Ritz) agrees with the REMARK in P.16.2. The sharper (P16.1) does not rescue it, because at ℓ₃ → 0 the F₃ part vanishes and the form is R − c₂‖ξ‖²_{H(E₂)}, which N.15 finds ≤ −0.194 (odd) at log 3.
- **N.15.H vs (H₃^σ).**
  - Both are PROVED sufficient conditions (each implies its 3-entry statement).
  - N.15.H keeps Q^{(≤2)} = R − 2c₂B₂ and needs μ_hi > c₃(1 + η). The 0.9009 enters there only heuristically ("folded in"): the 2-line's kink coupling stays inside Q^{(≤2)}, which is N.15-G2's "B₂ high-band cross terms".
  - (H₃^σ) moves exactly those 4-path a-edges into K_F, which is bounded fibrewise by c₀\* on F₃ (|F₃| = 4ℓ₃). This makes the 0.9009 requirement **exact**, and the cross terms of the joint block are handled by (P16.1) with no t-diagonality needed. The remaining 2-path part K_G is O(1)-supported and stays inside the hypothesis.
  - ~~So (H₃^σ) is N.15.H with N.15-G2's joint-block part discharged (PROVED).~~ **S103j: downgraded to a REMARK.** (H₃^σ) and N.15.H are different forms (Q^G_L vs Q^{(≤2)} = Q^G_L − ⟨K_F^{(a)}⟩) and neither implies the other. **What is PROVED** is narrower: Thm P.16.B needs no G2-type control of the 4-path a-edges, because those are bounded exactly by (P16.1). The non-multiplier certificate N.15-G1 is still needed, now for Q^G_L; the two forms agree as ℓ₃ ↓ 0.
  - *(Grok audit S103h(6): what is PROVED is that Thm P.16.B needs no G2-type control of the 4-path a-edges, because those are bounded exactly by (P16.1). (H₃^σ) and N.15.H are different hypotheses: they concern Q^G_L and Q^{(≤2)} = Q^G_L − (a-edges of K_F), and neither implies the other in general. So "=" is a REMARK identification. N.15-G1 for Q^G_L is still required.)*
- **Size (NUMERICAL, Prime's Ritz upper bounds).** At 2L = log 3 the margins are m ≲ 1.34e−3 (odd) and ≲ 2.5e−4 (even), not converged. With μ ≤ m, Thm P.16.B gives at best ℓ₃\* ≲ 1.34e−3/87 ≈ **1.5e−5** (odd; T = 38.1, η = 1), against N.15 B4's ≈ 1e−4. The difference comes from the 4ℓ₃ (four node pieces) and the c₀\* in place of c₃. Both are indicative only.
- **N.15.2:** the P₂ form becomes negative by log 3.5 (NUMERICAL). So (H₃^σ) and (M₃^σ) can hold only near ℓ₃ = 0, and P.16.B/C are perturbative entries, like S&M Thm 3. A non-perturbative certificate is needed beyond that.
- **Status, for the record (main 5acc5e5, ARP G.4i).** C.14.K is **CERTIFIED false** (Arb, 128 bits) at (ℓ, T) = (log 2/6, 14) and at (log 2/13, 14), and refuted NUMERICALLY on the rest of the useful range ℓ_K ≳ 0.03 (Celestial C.15; credit Prime cd94c27). This does not affect P.16, which does not use C.14.K.

### P.16.6 Labels and R-line
| item | label | track |
|---|---|---|
| Prop P.16.1: 4-path/2-path fibres; ‖compressed K‖ = c₀\* on (log 3, log 4); (P16.1); sharp in both parities | PROVED (sympy-checked) | [Cel] |
| N.14.3 "c₂ + c₃ cannot be improved": FIX (sharp is c₀\* on 𝒯_L); P.15.5 S1 row amended; P.15.4 "0.9009 sharp" CONJECTURE → PROVED | PROVED | [Cel] |
| Thm P.16.A: high band pays the joint charge for T ≥ 16 (η ≤ 0.0373): **the usable certified point** (S&M Ω(16)). T ≥ 15.5 (η ≤ 2.07e−3) rests on Principal's own Arb only | CERTIFIED (Arb Ω; (Ω2) PROVED) | [S&M] |
| Thm P.16.B: ℓ₃\* = πμ/(4c₀\*(1 + 1/η)T) | PROVED given (H₃^σ) | [S&M] |
| Thm P.16.C: exponential range from a plain margin m | PROVED given (M₃^σ) | [S&M] |
| (H₃^σ), (M₃^σ): low-band margins of Q^G_L near log 3, both sectors | GAP (= N.15-G1 for Q^G; sizes ≲ 1.3e−3 odd, 2.5e−4 even, NUMERICAL Prime) | [S&M] |
| (H₃^σ) "=" N.15.H: different forms, neither implies the other | REMARK (S103j downgrade) | [S&M]/[Cel] |
| Thm P.16.B needs no G2 control of the 4-path a-edges ((P16.1) bounds them exactly) | PROVED | [S&M]/[Cel] |
| (H₃^σ)(16, 0.0373, μ) at 2L = log 3, both parities | NUMERICALLY false for every μ ≥ −1 (P.17.0) | [S&M] |
| unconditional 3-entry range, odd and even | none | — |
| P.16.4(a)–(b): (I1) family D_L = A_LM_ρA_L^\*; D_L ≥ 0 ⇔ s_L ≥ 0 | PROVED | [Q] |
| (I2) solvable ⇔ Q_L ≥ 0 on pole-free 𝒯_L, as a measure identity (σ_L ≥ 0, σ̌_L = (¼ − ∂²)d on (−2L, 2L); weak invisibility) | PROVED by citation (Jorgensen–Pedersen–Tian arXiv:1712.07264 Cor 2.2, verified; LNM 2160 number unverified) + Schwartz–Bochner. A.c. realization: GAP | [Q] |
| line cost ≥ c₀\*; separable Fejér cost c₂ + c₃; coupled reach of c₀\* | PROVED / PROVED / GAP (CONJECTURE yes) | [Q] |
| separable split ⇒ global S2 (c₀ ≥ c₀\*) | REMARK (⇒ only); S2 GAP | [Q] |
| RH | OPEN | — |

**R-line:**
- **R1:** one rule on (log 3, log 4). The fibre combinatorics are constant there, ρ is fixed, and the thresholds T are L-independent.
- **R2:** the Euler product enters through the two lines and the fibres of {log 2, log 3}; the FE enters through Ω (and σ_S/ε).
- **R3:** the signs come from the bipartite 4-path spectrum (±c₀\*) and P-mirroring, not from weight signs.
- **R4:** the pole pair is spent once (pole-free class, g(0) = ‖ξ‖²_H; Res ≡ 0 by (N3)).
- **R5, where each fake breaks:**
  - **F_θ:** c_p → c_p cosh(θ log p) changes c₀\* but nothing structural. It breaks where P.14 R5 says (polarization, P.4).
  - **DH:** no Euler product, so the weights on log 2 and log 3 are not Λ(n)/√n and there is no fibre-line identification with γ_p. Breaks at P.14.2.
  - **H₋ₜ:** Ω is deformed, so the Arb thresholds of P.16.A do not transfer. Breaks at P.14.2/P.14.3.
  - **DMV / P_X (X ≥ 3):** agree with ζ on this window. Fake-blind here.
  - **G_θ′:** passes. Q_G is a sum of two ζ-forms, and every statement above is linear in the form or a bound on it.
- **R6:** P.16.1 and the (I1) family are exact. The closures (H₃^σ)/(M₃^σ) and the explicit n_L are GAPs. RH OPEN.

## P.17 The odd entry at log 3 and the first odd 3-entry range: Temple–Lehmann lemma, monotonicity, Celestial certification list, and the coupled line cost (S103j; analytic first; [Cel]/[S&M]/[Q])

Daniel Toupin / Golden Physics Project, with AI assistance from Grok. **RH is OPEN. Nothing in P.17 proves or implies RH.**

**Sources** (main 5acc5e5): Celestial §C.15 (Thm C.15.A, Lemma C.15.K, Prop C.15.B; merged, audited SOUND), Golden Eye's extension of C.15.B to the full Weil form at every L and its Temple–Lehmann PROPOSAL (room, S103j), Prime §N.16 (merged), S&M CERT_p2_constants, P.15–P.16 above.

**Input discipline (steering S103j).**
- **No zero-side input.** Every quantity below is on the prime/archimedean side: Ω, the line atoms Λ(n)/√n for n ≤ 4, and window geometry. Per N.16, the zero side gives only *upper* bounds on margins. Lower bounds come from the prime side (G1) only.
- Any picture of the form Q_L = Σ_γ(¼ + γ²)|ξ̂(γ)|² assumes RH. It is **heuristic, assuming RH**, and is never an input here.
- Numbers from Golden Eye or from my float probes are **NUMERICAL**. They are reference values for the certifiers, not results.

**Verification.**
- `scripts/p15_p16_check.py` parts 6–7 (new): sympy identities for the {2,3} comparison symbol and the 2×2 Schur step, and Arb values of q_T = Ω(T) − 2(c₂ + c₃). Parts 1–5 reproduce the merged log byte for byte.
- `scripts/p17_h3_probe.py`, `p17_probe.py`, `p17_temple_probe.py`, `p17_strip_sdp.py` (+ logs): float64 probes, NUMERICAL.

### P.17.0 What changed: (H₃^σ) at T = 16 fails, and only the full form is monotone
**(i) (H₃^σ)(16, 0.0373, μ) is NUMERICALLY false at 2L = log 3, in both parities, for every μ ≥ −1** (`p17_h3_probe.log`; sine basis N = 20/40/60, quadrature to t = 4000).
- At ℓ₃ = 0, Q^G_L is the {2}-form Q₂. (H₃^σ) at T = 16 asks Q₂ − Ω(16)Π_{≥16} ≥ μΠ_{<16}.
- λ_min(Q₂ − Ω(16)Π_{≥16}) ≤ −0.1734 (odd) and ≤ −0.2156 (even), stable in N. These are Ritz *upper* bounds, so the failure is robust up to quadrature.
- **Mechanism (REMARK).** The 2-line charges 2c₂cos(at) on the high band near t ≈ 4π/a = 18.13, just above T = 16. There R's slack Ω(t) − Ω(16) is still small, so handing the whole high band to Ω(16) over-spends it. The first T where the high-band part is compatible is ≈ 30 (NUMERICAL, from the S103j h3_scan). The best indicative P.16.B ranges are then ℓ₃\* ≈ 5e−6 (odd, T = 55, η = 1) and ≈ 1e−6 (even).
- **Consequence.** The lead's route (ii), "P.16 at T = 16, η ≤ 0.0373, plus G1", has no admissible input. The requested μ_σ for Q^G_L at T = 16 does not exist. P.16.B/C stay PROVED as conditional statements, but they are superseded by P.17.3.

**(ii) Monotonicity.**
- **Full Weil form (PROVED; C.15.B as extended by Golden Eye).**
  - Write W(ξ) := R(ξ) − Σ_n (Λ(n)/√n)·2⟨ξ, τ_{log n}ξ⟩_H. This is one L-independent form on C_c^∞(ℝ), and for ξ ∈ 𝒯_L only the terms with log n < 2L survive.
  - Extension by zero maps 𝒯_L^σ into 𝒯_{L′}^σ isometrically for L < L′, preserving parity (symmetric windows) and ξ̂. So λ₁^σ(L) := inf W/‖·‖²_H over 𝒯_L^σ is **nonincreasing in L**.
  - **Endpoint atoms vanish.** For ξ ∈ H¹₀(−L, L), g = ξ∗ξ\* is continuous and g(±2L) = 0. So the n = 3 atom at 2L = log 3 and the n = 4 atom (Λ(4)/√4 = a/2) at 2L = log 4 contribute 0. The closed endpoints are fine, and no δ-shift or extra c₄ term is needed (PROVED).
  - **Hence one certificate λ₁^σ(L\*) ≥ μ > 0 covers every 2L ≤ 2L\*, in that parity.**
- **The partial form Q^G_L is not monotone (PROVED, elementary). The reason:**
  - Q^G_L = R − ⟨ξ, K_G(L)ξ⟩_H, where K_G(L) keeps only the a-edges whose endpoints have *no* b-neighbour inside the window (the 2-path edges, P.16.1).
  - Take ξ ∈ 𝒯_L ⊂ 𝒯_{L′}. Then Q^G_{L′}(ξ) − Q^G_L(ξ) = 2c₂·(the a-edge form of ξ restricted to the edges that are 2-path edges for L but 4-path edges for L′). Those are the a-edges near the window ends whose b-partner lies in (−L′, L′) ∖ (−L, L).
  - That restricted edge form takes both signs (two bumps on one such edge, with equal or opposite signs). So neither Q^G_{L′} ≥ Q^G_L nor ≤ holds on 𝒯_L.
  - Q^G is not the restriction of an L-independent form: its line part depends on L through G(L) = [ℓ₃, b − a).
- **Answer to Golden Eye (1).** Yes, §P.17's target is phrased for the **full Weil form**, so C.15.B applies.
  - A single certificate at 2L = log 4 would cover all of (log 3, log 4]. But the margin there is ≲ 1.0e−7 (odd) and ≲ 1.1e−5 (even) (NUMERICAL, `p17_probe.log` Part A), consistent with N.16's super-exponential decay.
  - **Certification across the whole of (log 3, log 4] is not viable (HEURISTIC; the lead's record).** The obstruction is the size of the margin, not coverage.
  - REMARK on "block Lehmann on the cluster would cover log 4 only": by the monotonicity just stated, a positive λ₁ certificate at 2L = log 4 would cover every 2L ≤ log 4 in that parity. The difficulty is getting a positive certificate there at all.
- **Golden Eye's consistency test** (NUMERICAL, `p17_probe.log` Part A, N = 60):
  - Odd: Ritz upper bounds 9.48e−4 (log 3) → 1.19e−5 (log 3.5) → 1.0e−7 (log 4). Nonincreasing: **PASS**.
  - Even: 1.246e−4 (log 3) but 1.315e−4 (log 3.5). This **FAILS** at N = 60, which shows the even Ritz value at log 3 is not converged; N.16 already says the true even margin there is likely ≲ 1e−4. Celestial's N = 64 floats are to be trusted over N = 24, and even N = 60 is not enough for even.

### P.17.1 The Temple–Lehmann lemma for Q_L in the H¹_¼ / dual-norm setting — PROVED (by citation + elementary)
**Setting.**
- Fix 2L ∈ [log 2, log 4] and a parity σ. Let 𝓗 := 𝒯_L^σ, the real parity-σ subspace of H¹₀(−L, L), with ⟨u, v⟩_H = ∫(u′v′ + ¼uv) = (1/π)∫_0^∞(¼ + t²)Re(û v̂‾)dt.
- **Multiplier identity (PROVED).** For ξ ∈ 𝓗:
  > Q_L(ξ) = (1/π)∫_0^∞(¼ + t²)s_L(t)|ξ̂(t)|²dt, with s_L = Ω − 2c₂cos(at) − 2c₃cos(bt)·1[2L > b],
  - The n = 4 term vanishes by the endpoint remark in P.17.0. This is the prime/archimedean side of the classical Weil form, with no zero-side input.
- **Operator (PROVED).**
  - Ω(t) ≥ Ω(0) (Ω is increasing in |t|) and the line part is bounded by 2(c₂ + c₃).
  - So Q_L is a closed form on 𝓗, bounded below by Ω(0) − 2(c₂ + c₃). Its form domain is {ξ ∈ 𝓗 : ∫(¼ + t²)log(2 + t)|ξ̂|² < ∞}. It is the restriction to a closed subspace of a closed form plus a bounded one.
  - By Kato's first representation theorem it defines a self-adjoint A = A_L^σ on 𝓗, bounded below, with Q_L(u, η) = ⟨Au, η⟩_H for u ∈ D(A) and η in the form domain.
- **Residual in the dual norm.** For v in the form domain define the functional r_v(η) := Q_L(v, η) − ρ⟨v, η⟩_H.
  - v ∈ D(A) iff r_v is bounded on 𝓗. In that case ‖r_v‖_{𝓗\*} = ‖(A − ρ)v‖_H (Riesz).
  - A sufficient condition: ∫(¼ + t²)s_L²|v̂|² < ∞, by Cauchy–Schwarz. Finite combinations of Dirichlet modes satisfy it, since v̂ = O(t⁻²).

**Lemma P.17.A (Temple 1928; Kato 1949; Lehmann 1949/50).**
- **Hypotheses.**
  - (T1) v ∈ D(A), ‖v‖_H = 1, ρ := Q_L(v).
  - (T2) ν is a lower bound for the second spectral point: σ(A) ∩ (−∞, ν) consists of at most one point counted with multiplicity. Equivalently, the second min-max value λ₂(A) ≥ ν.
  - (T3) ρ < ν.
- **Conclusion.**
  > **λ₁(A) ≥ ρ − ‖r_v‖²_{𝓗\*}/(ν − ρ).**
- *Proof.* On σ(A) ⊂ {λ₁} ∪ [ν, ∞) one has (λ − λ₁)(λ − ν) ≥ 0, so ⟨(A − λ₁)(A − ν)v, v⟩ ≥ 0 by the spectral theorem. Expanding gives ‖Av‖² − (λ₁ + ν)ρ + λ₁ν ≥ 0. Then use ‖Av‖² − ρ² = ‖(A − ρ)v‖² = ‖r_v‖²_{𝓗\*}. ∎
- **Citations** (bibliographic data verified this session):
  - G. Temple, Proc. Roy. Soc. A 119 (1928) 276–293, doi 10.1098/rspa.1928.0098;
  - T. Kato, J. Phys. Soc. Japan 4 (1949) 334–339, doi 10.1143/JPSJ.4.334 (two-sided "Kato–Temple" form);
  - N. J. Lehmann, ZAMM 29 (1949) 341–356 and 30 (1950) 1–16 (the block / interval-enclosure generalisation).
  - The exact statement above, for lower-semibounded self-adjoint A with ψ ∈ D(A) and ⟨ψ, Aψ⟩ below a lower bound for the second eigenvalue, is quoted e.g. in arXiv:1912.08408, eqs. (4.1)–(4.2) (read this session).
  - The theorem numbers inside the 1928/1949 papers are not used.
- **Monotonicity in the inputs (PROVED).** Let R := ‖r_v‖_{𝓗\*}. The bound g(ρ, R, ν) = ρ − R²/(ν − ρ) is increasing in ν, decreasing in R, and increasing in ρ when R < ν − ρ. So interval enclosures are evaluated at (ρ̲, R̄, ν̲), with R computed from the same ρ-interval.
- **Loss is quadratic in the residual (REMARK).** For the bound to stay positive one needs R² < ρ(ν − ρ).

**Lemma P.17.B (coarse lower bound ν₂ for the second eigenvalue; PROVED).**
- **Data**, as in Thm C.15.A, extended to {2,3} when 2L > log 3 (P.17 sympy part 6):
  - q_T = Ω(T) − 2c₂ − 2c₃·1[2L > b];
  - a_T = s_L − q_T on [0, T];
  - m_abs = Ω(T) − Ω(t) + 2c₂(1 − cos at) + 2c₃(1 − cos bt)·1[2L > b];
  - m_neg = Ω(T) − Ω(t);
  - A_N = the matrix of q_T‖·‖² + 𝒜_T on V_N = span{f_0, …, f_{N−1}}, the H-orthonormal parity-σ Dirichlet modes;
  - s_V ≥ λ_max(m_abs-Gram on V_N);
  - ε^{abs}, ε^{neg} ≥ (1/π)∫_0^T(¼ + t²)m·K_{W_N} (Lemma C.15.K);
  - β² := s_Vε^{abs};
  - c_tail := q_T − ε^{neg}.
- **Conclusion.** For every δ > 0:
  > **λ₂(A) ≥ ν₂ := min{λ₂(A_N) − δ, c_tail − β²/δ}.**
- *Proof.*
  - Q_L ≥ Q′_T as forms (C.15.A step, {2,3} majorant identity of part 6), so λ_k(A) ≥ λ_k(Q′_T) by min-max.
  - For ξ = v + w with v ∈ V_N and w ⊥ V_N: Q′_T(ξ) ≥ ⟨A_N c, c⟩ − 2β‖v‖‖w‖ + c_tail‖w‖² ≥ ⟨(A_N − δ)c, c⟩ + (c_tail − β²/δ)‖w‖². The cross term is as in C.15.A; then use 2β‖v‖‖w‖ ≤ δ‖v‖² + (β²/δ)‖w‖².
  - The right side is the form of (A_N − δ) ⊕ (c_tail − β²/δ)I. The second-smallest point of its spectrum is ≥ the stated min. ∎
- This lemma *is* the non-multiplier Lemma A/B the lead asked for:
  - Lemma C.15.K's tail kernel replaces the bathtub Lemma A, and it sees the window constraint, which is the only non-multiplier feature of Q_L.
  - The 2×2 split replaces Lemma B.
  - Only O(β²/δ) and O(δ) losses enter, i.e. coarse (~1e−4 to 1e−2), as Golden Eye proposed.

**Exact residual formula (PROVED).**
- Let v = Σ_{j<n} c_j f_j. Parseval in the H-orthonormal basis {f_k} gives:
  > ‖r_v‖²_{𝓗\*} = Σ_{k<M} (Σ_{j<n} Q_{kj}c_j − ρc_k)² + τ_M², with Q_{kj} := Q_L(f_k, f_j), M ≥ n,
  >
  > τ_M := sup{|Q_L(v, w)| : w ⊥ V_M, ‖w‖_H = 1}, since ⟨v, w⟩_H = 0 there.
- **Tail bound (PROVED).** For any T₁ > 0:
  > **τ_M ≤ (1/π)∫_0^{T₁}(¼ + t²)|s_L||v̂|·√K_{W_M} dt + [(1/π)∫_{T₁}^∞(¼ + t²)s_L²|v̂|²dt]^{1/2}.**
  - The first term uses |ŵ(t)|² ≤ K_{W_M}(t)‖w‖²_H (C.15.K).
  - The second term uses Cauchy–Schwarz against (1/π)∫(¼ + t²)|ŵ|² ≤ 1.
- The second term is O(log T₁/√T₁) for v̂ = O(t⁻²). Imposing the single linear condition v′(±L) = 0 on the coefficients makes v̂ = O(t⁻⁴): v″ = −Σc_jω_j²f_j already vanishes at ±L. That term then becomes O(log T₁·T₁^{−5/2}). **Recommended.**

**Reference values (NUMERICAL; `p17_temple_probe.log`; Löwdin-orthonormalised quadrature, M = 160, the τ_M tail not included).**

| sector, 2L | n | ρ | λ₂ (float) | Σ_{k<160} part of ‖r‖ | Temple with ν = λ₂ | Temple with ν = 0.44 |
|---|---|---|---|---|---|---|
| odd, log 3 | 40 | 9.541e−4 | 0.6635 | 5.96e−3 | 9.01e−4 | 8.73e−4 |
| odd, log 3 + 0.03 | 40 | 4.742e−4 | 0.5335 | 4.55e−3 | 4.35e−4 | 4.27e−4 |
| even, log 3 | 40 | 1.247e−4 | **0.0763** | 8.26e−4 | 1.16e−4 | (ν = 0.44 invalid: λ₂ < 0.44) |

- **Odd at log 3** agrees with Golden Eye's numbers (λ₂ ≈ 0.66, λ₁ ≈ 9.5e−4; NUMERICAL). The positivity budget is R < √(ρ(ν − ρ)) ≈ 2.5e−2. For μ ≥ ρ/2 the budget is ≈ 1.8e−2 (ν = 0.66) or 1.4e−2 (ν = 0.44).
- **Even at log 3 (stretch).** Here λ₂ ≈ 0.076 (NUMERICAL), much smaller than odd. The positivity budget is ≈ 3.1e−3, and ≈ 2.2e−3 for μ ≥ ρ/2. *Flag:* this is looser than the room's "≲ 7e−4". The difference is in the assumed ν₂ and target, and should be reconciled with Golden Eye before Celestial sizes the even run.

### P.17.2 What the certificate gives — PROVED given the certified inputs
**Cor P.17.C.** Suppose Celestial certifies (T1)–(T3) for A_{L\*}^σ with enclosures giving μ_σ(L\*) := ρ̲ − R̄²/(ν̲₂ − ρ̄) > 0. Then:
> **Q_L(ξ) ≥ μ_σ(L\*)‖ξ‖²_H for all ξ ∈ 𝒯_L^σ and every 2L ≤ 2L\*.**

*Proof.* Lemma P.17.A, then the monotonicity of P.17.0(ii). ∎ No CC input is used. This is CC-free and classical on 2L ≤ log 4.

### P.17.3 The odd 3-entry range in terms of μ_odd
- **(Direct; the deliverable; PROVED given the certificate.)** Let 2L\* = log 3 + ℓ₃\* and μ_odd := μ_odd(L\*) > 0. Then:
  > **the odd sector is paid on all of (log 3, log 3 + ℓ₃\*], with margin ≥ μ_odd,**
  - and also on (0, log 3] (monotonicity).
  - The range formula is just "ℓ₃ ≤ ℓ₃\*, the largest certified right endpoint". μ_odd enters only through its positivity at L\*. No (1 + 1/η), c₀\* or T factor is lost.
  - **Recommended targets** (NUMERICAL margins, `p17_probe.log` Part A and Temple table):
    - ℓ₃\* = 0.03 primary (λ₁ ≈ 4.7e−4; Temple float 4.35e−4);
    - 0.01 fallback (≈ 7.6e−4);
    - 0.05 stretch (≈ 2.8e−4).
  - Because the full form is monotone, the first certificate needed is **odd at log 3** (target (i)). The second is the same pipeline at log 3 + ℓ₃\* (target (ii)).
- **(Perturbative, via P.16; kept for the record; PROVED given its hypothesis.)** Given (H₃^o)(T, η, μ_H) at **each** 2L in the range (Q^G is not monotone, so a certificate at one L does not propagate):
  > ℓ₃ ≤ πμ_H/(4c₀\*(1 + 1/η)T).
  - μ_H is a band-split margin of Q^G_L. It is **not** λ₁(Q_{log 3}): a plain λ₁ ≥ μ_odd does not imply (H₃^o), and at T = 16 (H₃^o) is false (P.17.0(i)).
  - Using P.16.C with a plain margin of Q^G_L gives ℓ₃ ≤ πθm/(8c₀\*(1 + 1/η)T) with T exponential in 1/m. For m ~ 1e−3 that is vacuous.
  - So the lead's (ii) "P.16 (T = 16, η ≤ 0.0373) + G1" should be replaced by "Cor P.17.C at L\* = ½(log 3 + ℓ₃\*)" (my recommendation; the lead decides).
- **Even.** At log 3 this is a stretch goal (P.17.1 table). An even 3-entry range would come from the same Cor P.17.C. Its float margin stays ~1.3e−4 for ℓ₃ ≤ 0.1, but it is not converged (consistency FAIL above): NUMERICAL, not a target yet.
- **Comparison-form-only route (C.15.A extended to {2,3}; NUMERICAL, `p17_probe.log` Part B).**
  - At T = 100 and N = 300, the float λ_cert is +4.36e−4 (ℓ₃ = 0.01), +2.41e−4 (0.03) and +9.1e−5 (0.05).
  - It is negative at N ≤ 128, so it needs N ≈ 200–300. For T ≈ 60 it is negative, since q_T(60) = 0.0076 (Arb part 7).
  - It works without Temple, but at roughly 5× the matrix size. Temple only needs ν₂ coarse.

### P.17.4 What Celestial must certify (Arb, 128 bits; prime/archimedean side only)
For each target:
- (i) odd, 2L = log 3;
- (ii) odd, 2L\* = log 3 + ℓ₃\* with ℓ₃\* ∈ {0.03; fallback 0.01; stretch 0.05};
- (iii) stretch: even, 2L = log 3.

Here s_L = Ω − 2c₂cos(at) − 2c₃cos(bt)·1[2L > log 3]. The n = 3 and n = 4 endpoint atoms contribute 0.

1. **Closed-form trial vector v.**
   - v = Σ_{j<n} c_j f_j with n = 40 (odd) or n = 40–60 (even).
   - f_j are the H-normalised parity-σ Dirichlet modes on (−L, L): odd sin(jπx/L), even cos((j + ½)πx/L), with ‖·‖²_H = L(ω_j² + ¼).
   - c_j are **exact dyadic rationals** (rounded float Ritz coefficients, ≥ 40 bits), published with the certificate.
   - Recommended: impose Σ_j c_j f_j′(L) = 0, so that v̂ = O(t⁻⁴).
2. **ρ.** Enclose ρ = Σ_{jk} c_jc_kQ_{jk}/Σc_j², with Q_{jk} = R_{jk} − 2c₂B²_{jk} − 2c₃B³_{jk}, in the C.15 closed forms plus Arb quadrature with a certified remainder for R_{jk}. Absolute precision 1e−10.
3. **Residual, finite part.** Enclose r_k = Σ_{j<n}Q_{kj}c_j − ρc_k for k < M (M = 160–400) at 1e−12 per entry, and output an upper bound for Σ_{k<M}r_k².
4. **Residual tail τ_M.**
   - Upper-bound the two-term bound of P.17.1: ∫_0^{T₁}(¼ + t²)|s_L||v̂|√K_{W_M} plus [∫_{T₁}^∞(¼ + t²)s_L²|v̂|²]^{1/2}, each /π.
   - v̂ in closed form. K_{W_M} = K_σ − Σ_{k<M}|f̂_k|², with K_σ in the C.15.K closed form (odd analogue to be derived; S&M Lemma A's dens_σ gives K_σ = π·dens_σ/(¼ + t²); please confirm the normalisation). The [T₁, ∞) term uses |s_L| ≤ |Ω| + 2(c₂ + c₃) with Ω ≤ log(t/2π) + O(1/t²) (S&M (Ω2)) and v̂'s explicit decay.
   - Target: **R̄ := (Σ_{k<M}r_k² + τ̄_M²)^{1/2} ≤ 1.2e−2** (odd, ensures μ ≥ ½ρ at ν₂ = 0.44), or ≤ 2.2e−3 (even).
5. **ν₂ (Lemma P.17.B).**
   - T = 60–100, N = 64 (Celestial's N = 64 floats are trusted over N = 24); choose δ.
   - Certify the **inertia**: A_N − (ν₂ + δ)I has at most one negative eigenvalue, by interval LDLᵀ with every pivot sign-determined. This is not an enclosed float λ₂.
   - Certify s_V from above by Cholesky of sI − M_abs (Frobenius is affordable as a fallback). Certify ε^{abs} and ε^{neg} from above (C.15.K, relative 1e−3), and c_tail − β²/δ ≥ ν₂.
   - Targets: ν₂ ≥ 0.44 (odd, log 3), ν₂ ≥ 0.40 (odd, L\*), ν₂ ≥ 0.05 (even). At 2L = log 3, q_T = Ω(T) − 2c₂ (no c₃), so c_tail is large and ν₂ is limited by λ₂(A_N) − δ.
6. **Output.** μ = ρ̲ − R̄²/(ν̲₂ − ρ̄), with ρ̄ < ν̲₂ checked.
7. **Consistency checks.**
   - μ_odd(L\*) ≤ ρ(log 3) (monotonicity).
   - ρ(L\*) ≤ ρ(log 3) up to Ritz error.
   - The full {2,3} Ritz values at log 3.5 and log 4 must be ≤ the log 3 value (Golden Eye; odd PASS, even FAIL at N = 60, P.17.0).
8. **Not allowed:** any zero-side input, open_\* placeholders, or custom-axiom Lean theorems as unconditional inputs.

### P.17.5 Question (c): can a coupled invisible n^F reach the line cost c₀\*? — duality PROVED; value → c₀\* NUMERICAL; exact attainment GAP
Setting as in P.16.4(c):
- λ(t) = −F(t) + n̂^F(t), with F := 2c₂cos(at) + 2c₃cos(bt), and cost c(n^F) := sup_t(−λ).
- The class 𝒩: finite **discrete** real even measures supported in |x| ≥ 2L.
- c_min := inf_{𝒩} c(n^F). Already PROVED: c_min ≥ c₀\* (P.16.4(c)), and the separable Fejér split costs c₂ + c₃ = 1.1244.

**(a) Prop P.17.D (dual problem; PROVED).**
> **c_min = max{2c₂G(1, 0) + 2c₃G(0, 1) : G positive definite on ℤ², G(0, 0) = 1, G(m, n) = 0 whenever |ma + nb| ≥ 2L}.**
- *Proof.*
  - All symbols are almost periodic, so sup_t = max over the Bohr compactification bℝ. Then c_min = inf_{h∈𝒩̂} p(F − h), with p(u) = max u on C(bℝ) sublinear.
  - Hahn–Banach gives c_min = max{ℓ(F) : ℓ ≤ p linear, ℓ|_{𝒩̂} = 0}. Here ℓ ≤ p means ℓ is a probability measure on bℝ, and G(x) := ℓ(e^{ixt}) is positive definite on ℝ_discrete with G(0) = 1 and G = 0 on |x| ≥ 2L.
  - Restrict to the subgroup Γ = aℤ + bℤ ≅ ℤ² (a/b irrational). Conversely, a PD G on Γ extended by zero is PD on ℝ_d (coset splitting) and satisfies the support constraint. ∎
- The visible set S = {(m, n) : |ma + nb| < 2L} is an infinite strip in ℤ². So (c) is exactly a **2-D trigonometric moment problem with strip-supported moments**, the problem P.16.4(c) called a 2-D Fejér–Riesz problem.
- The window feasible point G(m, n) = (χ∗χ̃)(ma + nb)/‖χ‖²₂, with χ the Perron bumps on one 4-path, gives c₀\*.

**(b) Primal certificates are explicit (PROVED, structure).**
- Truncate to the box |m|, |n| ≤ K and keep only the finite Toeplitz section [G(p − q)]_{p,q∈box} ⪰ 0. This is a relaxation with value v_K ≥ c_min, nonincreasing in K.
- Its SDP dual is: c − F(θ, φ) − h(θ, φ) = Σ_i|P_i(θ, φ)|², with P_i supported in the box and h having Fourier support off S.
- That **is** a coupled invisible n^F: atoms at ma + nb with |ma + nb| ≥ 2L, |m|, |n| ≤ 2K, of cost v_K. The Gram matrix Z ⪰ 0 of the SOS is an explicit positive line block for Queen's D (R6-explicit at finite K).

**(c) Values (NUMERICAL; `p17_strip_sdp.log`, CLARABEL).**

| K | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|
| v_K − c₀\* (2L = log 3 + 0.05) | 1.08e−1 | 5.69e−2 | 3.53e−2 | 2.42e−2 | 1.76e−2 | 1.34e−2 | 1.05e−2 |

- The same values hold to 1e−4 at 2L = log 3 + 0.2 and log 4 − 0.01. The fitted decay is ≈ K^{−1.8}.
- So **coupled invisible n^F beat the separable cost decisively**: v_8 = 0.9115 vs c₂ + c₃ = 1.1244 (NUMERICAL; an exact rational rounding of the K = 8 dual would make it CERTIFIED).
- The data are consistent with c_min = c₀\*.

**(d) Answer: GAP, with the exact obstruction.**
- c_min = c₀\* holds iff every PD G on ℤ² with moments supported in the strip S satisfies 2c₂G(1, 0) + 2c₃G(0, 1) ≤ c₀\*.
- The window points G = (χ∗χ̃)∘(a, b) achieve c₀\*. The gap would be closed by showing these are extremal, i.e. a 1-D Krein/Boas–Kac-type factorisation for strip-supported PD functions on ℤ² (pulled back along (m, n) ↦ ma + nb).
- For a general box truncation the relaxation gap is positive at every K tested. So **c₀\* is at best an infimum over finite n^F (NUMERICAL)**, and attaining it exactly needs infinitely many atoms. This is the 2-D Fejér–Riesz degree obstruction: in two variables, nonnegative trigonometric polynomials need not be sums of squares at the same degree (Rudin, Illinois J. Math. 7 (1963) 532–539, the extension problem on a square; bibliographic data unverified).
- **Label: duality PROVED; lower bound c₀\* PROVED; c_min = c₀\* CONJECTURE (yes), NUMERICAL support; GAP = extremality of window moments in the strip moment problem.**

**(e) REMARK (PROVED): s_L is unbounded in the multiplier ansatz.**
- Weak invisibility forces n = αe^{x/2} + βe^{−x/2} on (−2L, 2L), so n is smooth there, and local averages of n̂ against φ ≥ 0 with φ̌ ∈ C_c^∞(−2L, 2L) vanish.
- The line atoms average to O(1), while the Ω part averages to → ∞ like log|t|. So s_L is unbounded, and so is D_L = A_LM_ρA_L^\*.
- Hence P.16.4(a)'s "if s_L is bounded" was vacuous (fixed there). Queen's Q.8.E with bounded D needs a non-multiplier D: GAP.
- This does **not** affect (b): the SOS line block is bounded, and only the archimedean part is unbounded.

### P.17.6 Labels and R-line
| item | label | track |
|---|---|---|
| (H₃^σ)(16, 0.0373, μ) at 2L = log 3, both parities, every μ ≥ −1 | NUMERICALLY false (λ_min ≤ −0.173 odd, −0.216 even) | [S&M] |
| full Weil form monotone in L; endpoint atoms n = 3, 4 vanish; one certificate at L\* covers 2L ≤ 2L\* | PROVED (C.15.B + Golden Eye extension) | [Cel] |
| Q^G_L not monotone (2-path/4-path edge reassignment, signed) | PROVED (elementary) | [Cel] |
| certification across all of (log 3, log 4] | not viable (HEURISTIC; margin ≲ 1e−7 odd at log 4, NUMERICAL) | [Cel] |
| Lemma P.17.A (Temple–Kato–Lehmann for Q_L on 𝓗, dual-norm residual) | PROVED (by citation; proof given) | [Cel] |
| Lemma P.17.B (ν₂ via C.15-type split) | PROVED | [Cel] |
| exact H\*-residual formula and tail bound τ_M | PROVED | [Cel] |
| Cor P.17.C (odd margin on 2L ≤ 2L\* from one Temple certificate) | PROVED given the certified inputs | [Cel] |
| μ_odd(log 3) ≈ 9.0e−4, μ_odd(log 3 + 0.03) ≈ 4.35e−4, μ_even(log 3) ≈ 1.2e−4; λ₂ ≈ 0.66 / 0.53 / 0.076 | NUMERICAL (float; τ_M excluded; partly Golden Eye's numbers) | [Cel] |
| first unconditional odd 3-entry range | none yet: pending Celestial's Arb run (target ℓ₃\* = 0.03) | [Cel] |
| {2,3} C.15.A-only route: λ_cert > 0 at T = 100, N = 300 for ℓ₃ ≤ 0.05 (odd) | NUMERICAL | [S&M]/[Cel] |
| Prop P.17.D (coupled line cost = strip moment problem on ℤ²) | PROVED | [Q] |
| coupled n^F with cost 0.9115 at K = 8 (< c₂ + c₃) | NUMERICAL | [Q] |
| c_min = c₀\* | CONJECTURE (yes); GAP = extremality of window moments | [Q] |
| s_L unbounded in the multiplier ansatz; bounded D needs non-multiplier | REMARK (PROVED) / GAP | [Q] |
| Σ_γ(¼ + γ²)\|ξ̂(γ)\|² picture | heuristic, assuming RH; not used | — |
| RH | OPEN | — |

**R-line:**
- **R1:** one rule per certificate. Monotonicity turns the certificate at L\* into the whole range 2L ≤ 2L\*, and T, N and δ are fixed per target.
- **R2:** the Euler product enters through the atoms Λ(n)/√n, n ≤ 4, and the fibres of {log 2, log 3}; the FE enters through Ω.
- **R3:** the signs are kept: s_L is signed, Temple uses the signed form, and the line cost comes from the strip moment problem, not from weight signs.
- **R4:** the pole pair is spent once (pole-free 𝒯_L, g(0) = ‖ξ‖²_H).
- **R5, where each fake breaks:**
  - **F_θ:** c_p → c_p cosh(θ log p) changes s_L. The Temple certificate is form-specific, so a fake needs its own certificate, and the margins ~1e−3 are not stable under θ ≠ 0 of size ≳ 1e−3 (HEURISTIC). It breaks structurally at P.4.
  - **DH:** no Euler product, so the atom weights are not Λ(n)/√n. Breaks at P.14.2.
  - **H₋ₜ:** Ω is deformed, so the Arb inputs do not transfer.
  - **DMV / P_X (X ≥ 5):** these agree with ζ on 2L ≤ log 4. They are fake-blind here, as every finite-window certificate is.
  - **G_θ′:** passes. All statements are linear in the form or are bounds on it.
- **R6:** Lemmas P.17.A/B, Cor P.17.C and Prop P.17.D are explicit. The certificates are pending Celestial's Arb run, and c_min = c₀\* is a GAP. RH OPEN.

## References (added to ARP G.7)
- Temple, G., "The theory of Rayleigh's principle as applied to continuous systems", Proc. Roy. Soc. A 119 (1928) 276–293.
- Kato, T., "On the upper and lower bounds of eigenvalues", J. Phys. Soc. Japan 4 (1949) 334–339.
- Lehmann, N. J., "Beiträge zur numerischen Lösung linearer Eigenwertprobleme" I, II, ZAMM 29 (1949) 341–356; 30 (1950) 1–16.
- Jorgensen, P., Pedersen, S. & Tian, F., "Positive definite (p.d.) functions vs p.d. distributions", arXiv:1712.07264 (Cor 2.2); *Extensions of Positive Definite Functions*, LNM 2160, Springer (2016).
- Rudin, W., "The extension problem for positive-definite functions", Illinois J. Math. 7 (1963) 532–539 (bibliographic data unverified).
- Knapp, A. W. & Stein, E. M., "Intertwining operators for semisimple groups", Ann. of Math. 93 (1971) 489–578.
- Knapp, A. W., *Representation Theory of Semisimple Groups: An Overview Based on Examples*, Princeton (1986), Ch. II, VII.
- Shahidi, F., Ann. of Math. 132 (1990) 273–330; *Eisenstein Series and Automorphic L-functions*, AMS Colloq. Publ. 58 (2010).
- Bargmann, V., Ann. of Math. 48 (1947) 568–640.
- Langlands, R. P., LNM 544 (1976).
- Connes, A., arXiv math/9811068 (Selecta 5 (1999)): §III eqs (6)–(19); App. I Lemma 2, eqs (27)–(38); §VIII eq (22).
- Meyer, R., arXiv math/0311468v3 (Duke 127 (2005)): eqs (1)–(4), (11)–(15), (31), (33); Lemmas 3.1, 5.3, 5.5; §§5.1, 5.4, 5.6.
- Connes, A. & Consani, C., arXiv:2006.13771 (Selecta 27:77 (2021)): Thm 1 (4), (5)/Thm 6.11, Thm 3 (13)–(14)/Thm 4.7, eqs (8)–(11), (17), (24)–(26), (48)–(49), (51)/Cor. 2.3, Lemma 1.3, Lemma 4.1 (60), (62)–(63), Def. 4.4 (72).
- Connes 1999 additionally: §VII (15)–(16), (31), (33), Thm 3; §VIII (16)–(28), Cor. 2, Thm 5, Lemma 3; App. II (5)–(6).
- Queen, `team/queen/arp/sign_check.md` (2 log Λ cancellation, Thm 5 write-up). §N N.10 (`team/prime-number-theory/ARP_prime-number-theory.md`). §S S(c)4–S(c)5.
- Tate, J., Thesis (1950), Cassels–Fröhlich Ch. XV.
