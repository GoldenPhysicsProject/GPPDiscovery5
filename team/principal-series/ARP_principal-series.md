# §P Principal series: local reflections, the polarization lemma, and the E-intertwining (for ARP.md §P)

Daniel Toupin / Golden Physics Project, with AI assistance from Grok. Seat: Principal series. 2026-10-28 (S102/S103 directions).
**RH is OPEN. Nothing in this section proves or implies RH.** The section's result is negative and structural: every local step is an isometry, the only non-isometric step is the global ℚ^*-sum, and nothing local bounds the shared GAP G.4d.

**Labels.** PROVED (literature or elementary, cited) · PROVED-LEAN (GPPVerify name, main 5f54359, checked 2026-10-28) · CONJECTURE · GAP · NUMERICAL (orientation only).
**Conventions (team-binding).** L is the HALF-width: f ∈ C_c^∞([−L, L]) in log scale u = log|x|, g = f ∗ f̃ on [−2L, 2L], window L sees n < e^{2L}. Lean's `HasSupportIn g L′` in `WeilSupportLadder.lean` is stated for the convolution g, so **L′ = 2L**.
**Builds on.** ARP.md §G (G.1–G.4, GAP G.4d) and §0.3 (inventory); `ARP_plan.md` §1–§2 (this file supersedes its §2 outline and corrects its §5 table, see P.9). Exact step checks: `arp/p_checks.py` → `arp/p_checks.log` (sympy; no exploratory numerics).

**Fakes and requirements (each step's R-line).** F_θ = Λ(s−θ)Λ(s+θ) (θ real ≠ 0); G_θ′ = Λ(s+iθ′)Λ(s−iθ′) (must pass); DH (Davenport–Heilbronn); H₋ₜ (de Bruijn–Newman deformation); DMV (Beurling); P_X (finite Euler product); Epstein x²+5y² (h = 2); G₅ (Hecke-group, non-arithmetic). R1 uniform in L; R2 Euler product + FE jointly; R3 not via weight signs; R4 soft mode / pole state; R5 locate each fake; **R6 exact identity, no slack** (every step below is an equality unless stated).

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
4. **Nothing local bounds the GAP.** Fl_L ≥ −MF_L uniformly in L (G.4d) is a statement about L_S at α = ½ composed with the cutoff. It is equivalent to RH (Weil 1952; Bombieri 2000), and window L is the semi-local problem S_L = {∞} ∪ {p < e^{2L}} (PROVED-LEAN `GppWeilSupportLadder.primeSide_eq_truncation`, WeilSupportLadder.lean:90, with Lean L′ = 2L).
**Base case.** For 2L < log 2, no prime enters (PROVED-LEAN `primeSide_eq_zero_of_support_lt_log_two`, line 106). Positivity there is Connes–Consani's archimedean theorem (Selecta 27:77, Thm 1 / Cor. 2.3); in Lean it is only conditional, `weil_nonneg_of_arch_nonneg_rung_zero` (named hypothesis 0 ≤ arch).
**R-line.** The GAP is R1 by construction. All fakes break before or at P.6, and G_θ′ passes every step (ARP G.4e, refined in P.9). **No new route is claimed.**

## P.8 Not used
- HiddenZeros (GPPVerify #454, 5f54359): 4–6-point amplitude identities; no celestial amplitude enters §P. **Not cited.**
- `temperedness_iff_critical_line` / ArithmeticTemperedPrincipalSeries.lean (custom axiom `exp_growth_not_tempered`) and ArithmeticCoPoissonOrientationTemplate.lean (co-Poisson assumed): **not cited**. P.3 replaces them with literature statements.
- `open_*` placeholders: never cited.

## P.9 Corrected fakes table (replaces ARP_plan §5; proposed refinement of ARP G.4e's F_θ and H₋ₜ rows)

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

## P.10 Label changes relative to ARP_plan.md (2026-10-22)
1. **§5 F_θ row:** "local Θ_v unitary (tempered): no (non-tempered)" → **"GL(2) local Θ_v unitary: YES (product γ, P.3 b1); Knapp–Stein form PD for 0<θ<½ (b2); tempered: no; GL(1) constituents unitary: NO; breaks at P.4 (polarization)."** The column header "unitary (tempered)" conflated two different notions; it is now split.
2. **§1 "Which form is Hermitian"** (Knapp–Stein form "not the L² form, E does not intertwine it"): superseded by P.3. The Knapp–Stein form is positive definite on the complementary series. The relevant point is that the explicit formula and E are linear in the GL(1) constituents.
3. **§3 item 2** ("unitarity holds only for tempered local data; that is the R3 filter"): corrected to "the R3 filter is unitarity of the GL(1) isobaric constituents on C_ℚ; local GL(2) unitarity holds for complementary series."
4. **§2 temperedness:** "filter, not source" is kept. The filter is now located as global (P.3 b3); local non-temperedness is its shadow.
5. **§5 H₋ₜ row:** "local Θ_v yes, ⊗′ yes" → **no Euler product after deformation; Θ_∞ not Tate's** (agrees with ARP G.4e).
6. **§4 R4 CONJECTURE** (soft mode = window image of the two Poisson boundary directions): **SUPERSEDED**, not tested. ARP G.3a supplies R4 through the mean field, and S103/G.3d withdrew T1 (OS4 is fully used by the mean field).
7. **Open items:** (a) E citation → **closed** (P.5). (c) polarization lemma → **closed** (P.4). (b) δ_Λ in Connes–Consani notation remains open. (d) is superseded with item 6.
8. **New PROVED-LEAN citations:** `weilForm_nonneg`, `cos_mul_positiveType`, `positiveType_add`, `const_one_positive_type`, `positiveType_weighted_gram`, `neg_logDeriv_product_of_hasDerivAt`, `F_one_sub`, `osReflection_eq_self_iff`, `primeSide_eq_truncation`, `primeSide_eq_zero_of_support_lt_log_two`, `tate_local_zeta_integral`, `euler_factor_bridge` (all checked on GPPVerify main 5f54359). Lean `HasSupportIn` uses L′ = 2L.

## References (added to ARP G.7)
- Knapp, A. W. & Stein, E. M., "Intertwining operators for semisimple groups", Ann. of Math. 93 (1971) 489–578.
- Knapp, A. W., *Representation Theory of Semisimple Groups: An Overview Based on Examples*, Princeton (1986), Ch. II, VII.
- Shahidi, F., Ann. of Math. 132 (1990) 273–330; *Eisenstein Series and Automorphic L-functions*, AMS Colloq. Publ. 58 (2010).
- Bargmann, V., Ann. of Math. 48 (1947) 568–640.
- Langlands, R. P., LNM 544 (1976).
- Connes, A., arXiv math/9811068 (Selecta 5 (1999)): §III eqs (6)–(19); App. I Lemma 2, eqs (27)–(38); §VIII eq (22).
- Meyer, R., arXiv math/0311468v3 (Duke 127 (2005)): eqs (1)–(4), (11)–(15), (31), (33); Lemmas 3.1, 5.3, 5.5; §§5.1, 5.4, 5.6.
- Connes, A. & Consani, C., arXiv:2006.13771 (Selecta 27:77 (2021)): eq (17), Lemma 1.3, eqs (24)–(26).
- Tate, J., Thesis (1950), Cassels–Fröhlich Ch. XV.
