# §P Principal series: local reflections, the polarization lemma, and the E-intertwining (for ARP.md §P)

Daniel Toupin / Golden Physics Project, with AI assistance from Grok. Seat: Principal series. 2026-10-28 (S102/S103 directions).
**RH is OPEN. Nothing in this section proves or implies RH.** The section's result is negative and structural: every local step is an isometry, the only non-isometric step is the global ℚ^*-sum, and nothing local bounds the shared GAP G.4d.

**Labels.** PROVED (literature or elementary, cited) · PROVED-LEAN (GPPVerify name, main 5f54359, checked 2026-10-28) · CONJECTURE · GAP · NUMERICAL (orientation only).
**Conventions (team-binding).** L is the HALF-width: f ∈ C_c^∞([−L, L]) in log scale u = log|x|, g = f ∗ f̃ on [−2L, 2L], window L sees n < e^{2L}. Lean's `HasSupportIn g L′` in `WeilSupportLadder.lean` is stated for the convolution g, so **L′ = 2L**.
**Builds on.** ARP.md §G (G.1–G.4, GAP G.4d) and §0.3 (inventory); `ARP_plan.md` §1–§2 (this file supersedes its §2 outline and corrects its §5 table, see P.11). Exact step checks: `arp/p_checks.py` → `arp/p_checks.log` (sympy; no exploratory numerics).

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
**R-line.** The GAP is R1 by construction. All fakes break before or at P.6, and G_θ′ passes every step (ARP G.4e, refined in P.11). **No new route is claimed.**

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
- Pushing CC from Λ = 1 to Λ > 1 at 2L < log 2 [S&M] changes only the Sonin/prolate bookkeeping (P9.4). A Λ > 1 analog of CC Thm 1 would need a new compact-operator bound on the ε^{(Λ)} ∘ Q functional on the support class (the analog of K_I, Lemma 6.10). That is a **GAP [S&M]**, and RH-silent either way.
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

## P.10 Not used
- HiddenZeros (GPPVerify #454, 5f54359): 4–6-point amplitude identities; no celestial amplitude enters §P. **Not cited.**
- `temperedness_iff_critical_line` / ArithmeticTemperedPrincipalSeries.lean (custom axiom `exp_growth_not_tempered`) and ArithmeticCoPoissonOrientationTemplate.lean (co-Poisson assumed): **not cited**. P.3 replaces them with literature statements.
- `open_*` placeholders: never cited.

## P.11 Corrected fakes table (replaces ARP_plan §5; proposed refinement of ARP G.4e's F_θ and H₋ₜ rows)

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

## P.12 Label changes relative to ARP_plan.md (2026-10-22)
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
   - The Λ > 1 CC Thm 1 is a GAP [S&M], RH-silent.

## References (added to ARP G.7)
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
