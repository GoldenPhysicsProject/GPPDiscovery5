# §Q Arithmetic Osterwalder–Schrader axioms for ζ: OS0–OS4, the fakes, and the one GAP (Queen of Mathematics, for derivation/ARP.md)

Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
Queen of Mathematics (RH sub-team), session 2 (2026-10-28). Written for inclusion in `derivation/ARP.md`, against ARP at GPPDiscovery5 main `36aa6be` (§G incl. the S103b G.3d repair) and LEAD_DIRECTIONS S102/S103 (GPP-bridge main `d47c653`). Lean names are checked against GPPVerify main `5f54359` (read-only clone) and ARP §0.3 (the inventory; `derivation/LEAN_INVENTORY.md` does not exist on any branch, and S102 says the inventory now lives in §0.3).

**RH is OPEN. Nothing in §Q proves or implies RH.** §Q packages ARP §G as five axioms. Four of them (OS0, OS1, OS3, OS4) and half of the fifth (OS2a) are PROVED for ζ. The other half (OS2b) is Weil positivity, so it is RH-equivalent: it is THE GAP, and it is the same GAP as G.4d. Every step carries a label and one R-line (R1 uniform in L; R2 Euler product and FE together; R3 not via the sign of the weights; R4 the positive space contains the soft mode/pole state; R5 where each fake breaks, and G_θ′ passes; R6 an exact identity, not a lossy inequality: ζ's margin is about 2^−4900, so a lossy step cannot survive).

**Labels.** PROVED (with citation), PROVED-LEAN (with GPPVerify name), CONJECTURE, GAP. NUMERICAL is orientation only. "PROVED in sketch" is used once (Turán Thm A), as the lead asked.

---

## Q.0 Setting (§G notation)

- C = C_ℚ = 𝔸^*/ℚ^* is the idele class group, with module |·| : C → ℝ₊^*. The ζ-sector (Ẑ^*-invariant, even at ∞) is L²(C)^{Ẑ^*} ≅ L²(ℝ, du) with u = log|g|.
- **Window.** Let 𝒯_L = C_c^∞((−L, L)), real. For f ∈ 𝒯_L put g = f ∗ f̃ and f̂(t) = ∫ f e^{itu}. Connes' cutoff S_Λ corresponds to Λ = e^{L}, with the window g supported in [−2L, 2L].
- **Explicit formula and Weil form (PROVED; Weil 1952; Barner 1981).** Q_L(f) = Σ_ρ ĝ(γ_ρ) = ĝ(i/2) + ĝ(−i/2) − Σ_n Λ(n) n^{−1/2}(g(log n) + g(−log n)) + W_∞(g).
- **Reflection.** Θ := F_𝔸 ∘ inv. Transported by Connes' E it is Θ_C ξ := ξ ∘ inv, i.e. u ↦ −u: E(F_𝔸 φ) = E(φ) ∘ inv (Connes 1999 VIII (22); ARP G.4b).
- **The reflection pairing.** For real f, f^*(g) := f(g^{−1}) = (Θ_C f)(g). So Q_L(f) = Δ(f ∗ Θ_C f), where Δ is Weil's distribution on C. **This is the OS pairing used in OS2.** It is a group-algebra (Toeplitz) pairing: Θ_C acts as the involution of the convolution algebra. It is *not* the half-line Hankel pairing; see OS2.
- **Mean field and fluctuation (ARP G.3a, G.3d).**
  - MF_L(f) := ∫∫ f(u) f(v) e^{−|u−v|/2} du dv. By the Ward identity at v = 0 for even g, this is exactly pole − PNT-mean = ĝ(i/2) + ĝ(−i/2) − ∫_0^∞ (g(u) + g(−u)) e^{u/2} du.
  - E_Λ := Σ_n Λ(n) n^{−1/2} δ_{log n} − e^{u/2} du, the ψ-fluctuation.
  - Fl_L(f) := W_∞(g) − ∫_0^∞ (g(u) + g(−u)) dE_Λ(u).
  - Then **Q_L = MF_L + Fl_L exactly** (R6).
  - Convention: W_∞ is put in Fl_L. ARP §G.3 leaves this implicit. Moving W_∞ to the other side changes nothing below except which side carries it.

**A datum for a candidate D** (a fake or G_θ′) is the corresponding list: local factors, lattice, reflection, boundary terms, Weil distribution Δ_D, window form Q_{D,L}. "D satisfies OSk" means the OSk statement holds with D's data. "Fails" means the data the statement needs does not exist, or violates it.

---

## Q.1 The axioms

### OS0 Regularity (PROVED for ζ)
For every L:
- Q_L is a real symmetric form on 𝒯_L.
- The zero sum Σ_ρ ĝ(γ_ρ) converges absolutely.
- The explicit formula holds for g = f ∗ f̃ (Lipschitz suffices).
- Δ is a tempered distribution, even under inv (the FE).
- **Window locality.** Q_L depends only on Δ restricted to |u| ≤ 2L. The prime side sees exactly the n < e^{2L}.

Proof and citations:
- Weil 1952; Barner 1981.
- PROVED-LEAN `primeSide_eq_truncation` and `primeSide_eq_zero_of_support_lt_log_two` (RiemannHypothesis/WeilSupportLadder).

R-line: R1, OS0 is per-window and so is necessary, not sufficient (P_X satisfies it on L ≤ (log X)/2) | R2, not engaged | R3, not engaged | R4, pole term present | R5, all six fakes and G_θ′ satisfy OS0 (H₋ₜ on the zero side) | R6, exact.

### OS1 C_ℚ-covariance with the ℚ^* lattice (PROVED for ζ)
- There is a representation U of C on 𝓗 = L²(𝔸/ℚ^*) (Connes' L²(X), δ = 0). There is an isometry E : L²(X)₀ → L²(C) with (Ef)(g) = |g|^{1/2} Σ_{q∈ℚ^*} f(qg) and **E U(a) = |a|^{1/2} V(a) E** (V is the regular representation).
- The window forms are translation-covariant (Toeplitz in u).
- The off-diagonal singular support of Δ (the prime atoms) comes from the q-sum.
- In OS language: the "Euclidean group" is C and the "lattice of the field" is ℚ^* ⊂ 𝔸^*.

Citations: Connes 1999 §III (18)–(19) and §VIII (18)–(19); Meyer 2005 (the virtual representation whose character is Weil's distribution). Lean: only the type `RationalIdeleGroup` (RiemannHypothesis/IdeleGroup) exists. The equivariance itself is not in Lean.

R-line: R1, uniform (a group action, not a window) | R2, half: the lattice ℚ^* is the Poisson/FE half's prerequisite | R3, not via weights | R4, E is defined on 𝒮(𝔸)₀, i.e. with the pole removed, and the pole is returned in OS4 | R5, **DMV and P_X fail here** (Q.2) | R6, exact.

### OS2 Reflection positivity with respect to Θ = F_𝔸 ∘ inv
OS2 has two parts.

**OS2a (the reflection is a unitary involution, intertwined by E; PROVED for ζ).**
- Θ_C = inv is a unitary involution of L²(C), and E(F_𝔸 φ) = E(φ) ∘ inv.
- Locally, Θ_v acts in Mellin coordinates by γ_v(½+it) with |γ_v| = 1 on each GL(1) constituent. This is Tate's local FE: ARP G.1b for p, G.2a for ∞.
- Citations: Tate 1950; Connes 1999 VIII (22).
- PROVED-LEAN for the algebra of the reflection: `ArithmeticOSReflection.osReflection`, `osReflection_eq_self_iff` (the fixed set is Re s = ½).
- **Constituent filter.** OS2a is imposed on every GL(1) constituent, not only on the product (lead S102(b); ARP G.4e).

**OS2b (positivity; GAP, RH-equivalent).** For every L > 0 and every f ∈ 𝒯_L,

  Q_L(f) = Δ(f ∗ Θ_C f) ≥ 0.

- This is the Toeplitz (group-involution) form.
- The half-line Hankel form "∫∫_{u,v≥0} f(u) f(v) W(u+v) ≥ 0" is **false** for ζ's full W. Citations: PROVED-LEAN `HaarSquareNotReflection.hankel_neg`, `GaussianHeatOS.gaussian_hankel_det_neg`; ARP G.2c and §0.2.
- So OS2b must be read in the Toeplitz sense. Only the mean field satisfies the Hankel form too (OS2′ below).

R-line for OS2a: R1, uniform | R2, the global part is Poisson (FE half) | R3, γ_v-unitarity is not a weight sign | R4, not engaged | R5, **F_θ fails on GL(1) constituents only** (Q.2) | R6, exact.

R-line for OS2b: R1, required uniformly in L (bounded windows never suffice, by DMV/P_X) | R2, it must use OS1 and OS3 jointly | R3, F_θ has positive weights and fails, G_θ′ has signed weights and passes | R4, contains the pole state through OS4 | R5, fails for every fake with off-line zeros | R6, an exact identity is needed. This is the GAP.

**Sector bookkeeping note for OS2b on the first prime window (Queen, 2026-10-28; source: main 8ecee4f, derivation/ARP.md G.4h; labels copied from there, not upgraded).** This concerns only the window 2L ∈ (log 2, log 3], i.e. with one prime line (p = 2) present.
- *Who pays for what.* The p = 2 line carries the **even edge**. In G.4h item 2 (Prime N.13.4, PROVED, exact) the rank-two boundary form helps the even edge pair and charges the odd one exactly 2c₂α²/ℓ. The **archimedean margin at ½ log 2 carries the odd sector**; this is G.4h item 4, with conditional entry theorems only (Prime Thm N.13.A, conditional on GAP N.G2-m₀; Spectrum Thm 2.1, conditional on δ_o⁰, k_∞, B, which are computational GAPs). **𝒮_L pays for the ĝ(0) direction**: the CC eigen-direction μ₁ ≈ 1.0518 has cost ≈ 2.38 and is paid through the rank-one Schur condition (P) (§S 4.3, **GAP** b-1). The split itself is labelled in 8ecee4f as "Grok reconciliation, PROVED as bookkeeping".
- *What is still open on this window.* The coupling 2R(ξ_b, ξ_s) is **GAP N.G2-bulk** (Celestial C.14 was *asked*; it is not proved in 8ecee4f). (P) is a **GAP**. The constants δ_o⁰, k_∞, B are computational **GAP**s. Positivity on (½ log 2, ½ log 3] is separately **CERTIFIED, conditional on Mori** (main 992a011). The bookkeeping only says *which term carries it* analytically.
- *How this sits against OS2b.* This is an R1 statement about a bounded window. OS2b is required uniformly in L, and its GAP is **G.4d** (uniformity in L, RH-equivalent). G.4d is **untouched** by 8ecee4f, and by the DMV/P_X remark in the R-line above no bounded-window result can close it. OS2b keeps its label (**GAP, RH-equivalent**). The even/odd split is G.4h's parity split on this window. I do not identify it with any decomposition of Θ in OS2a (OPINION: a natural thing to check, but it is not claimed here).
R-line: R1 one fixed window only | R2 FE through W_∞ / 2θ′, Euler product through the single p = 2 line | R3 sign of B by parity | R4 the ĝ(0) direction is paid by 𝒮_L | R5 per 8ecee4f: fake-blind on this window, and DH has no a-line | R6 bookkeeping is exact; the analytic closure still needs N.G2-bulk + (P) + three certified constants.

### OS3 Euler clustering (restricted tensor product; PROVED for ζ)
- 𝒮(𝔸) = ⊗′_v 𝒮(ℚ_v), F_𝔸 = ⊗′_v F_v. The vacuum φ⁰ = ⊗_v φ_v⁰ with φ_p⁰ = 1_{ℤ_p} and φ_∞⁰ = e^{−πx²}.
- ∫_C E(φ⁰)(g) |g|^{s−½} d^*g = Π_v Z_v(φ_v⁰, s) = ξ(s) (Tate).
- **Clustering in window form.** Δ = pole + Δ_∞ + Σ_p Δ_p.
  - Each Δ_p depends only on the local factor L_p.
  - Each Δ_p is supported on {±k log p} (ARP G.1c).
  - Contributions of distinct places add, as the logarithm of a restricted tensor product. Equivalently, the Dirichlet coefficients are multiplicative.

Citations: Tate 1950 §§3–4; PROVED-LEAN `euler_factor_bridge` (RiemannHypothesis/PadicEulerFactorBridge: Tate's local ℤ_p integral equals Mathlib's Euler factor).

R-line: R1, uniform (all places) | R2, the Euler half | R3, local Δ_p has sign-changing spectral density (G.1c), so no local positivity is claimed | R4, not engaged | R5, **DH, Epstein and H₋ₜ fail here** (Q.2) | R6, exact.

### OS4 Pole/boundary condition (PROVED for ζ; fully consumed by the mean field)
- The two boundary functionals φ ↦ φ(0) and φ ↦ φ̂(0) span the supplement ℂ ⊕ ℂ(1) of 𝒮(𝔸)₀ ⊂ 𝒮(𝔸) (Connes 1999 III (31)). They produce exactly the pole term ĝ(i/2) + ĝ(−i/2).
- Through Astra's Ward identity (ARP G.3a) they are **spent exactly once**: pole − PNT-mean = MF_L.
- **S103: no boundary direction is left over.**
  - Target T1 is withdrawn.
  - G.3d is PROVED. In the S103b repaired form: ψ = Σ_k θ(x^{1/k}), so the layers μ_k sum to the half-density of ψ, whose PNT mean is the single e^{u/2} du. The k = 2 rank-one even bulk is therefore a ψ/θ bookkeeping artefact.
  - No extra mean-field layer exists.
  - Note: the lead's wording "cancels against θ(x)'s −√x term" is the pre-repair proof. Golden Eye showed it is wrong term by term (k = 4, 6). The statement stands via ψ. ARP G.5 (T1) still carries the old phrase "it cancels against θ's −√x" and should be updated.

R-line: R1, uniform | R2, the poles come from Poisson (FE half) | R3, not via weights | R4, **this is R4**: the pole state e^{−u/2} spans the rank-one OS square | R5, DH has no pole, so OS4 holds only vacuously and R4 is unavailable; the others keep an MF ≥ 0 (G.3c) | R6, exact.

---

## Q.2 Fakes table, with proofs

The fakes are F_θ (0 < θ < ½), DH, DMV, P_X, Epstein x²+5y², H₋ₜ (t > 0), and the control G_θ′. ✓ = satisfies, ✗ = fails.

| D | OS0 | OS1 | OS2a | OS2b | OS3 | OS4 | first failure | expectation |
|---|---|---|---|---|---|---|---|---|
| **G_θ′** = Λ(s+iθ′)Λ(s−iθ′) | ✓ | ✓ | ✓ | ⇔ OS2b(ζ) | ✓ | ✓ | none | **verified** |
| **F_θ** = ξ(s−θ)ξ(s+θ) | ✓ | ✓ (non-unitary quasi-characters) | ✗ GL(1) constituents; ✓ product | ✗ | ✓ | ✓ (MF ≥ 0) | OS2a constituent filter, then OS2b | **verified, with a convention caveat** |
| **DH** | ✓ | ✓ (χ, χ̄ sectors) | ✓ | ✗ | **✗** | vacuous (no pole) | OS3 | **verified**; also no R4 |
| **DMV** | ✓ | **✗** | ✗ (no FE; consequence) | ✗ | (✓ Beurling) | ✓ (MF ≥ 0) | OS1 | **verified**, with a refinement |
| **P_X** | ✓ | **✗** past X; ✓ on L ≤ (log X)/2 | ✗ (no FE) | ✗ | ✓ up to X | ✓ | OS1 | **verified**; invisible on bounded windows (R1) |
| **Epstein** x²+5y² | ✓ | ✓ over K = ℚ(√−5) (two C_K sectors) | ✓ | ✗ | **✗** | ✓ | OS3 | determined: **OS3** |
| **H₋ₜ** | ✓ | ✓ (formally) | ✓ | ✗ | **✗** (natural realization) | convention | OS3 | determined: **OS3** (intrinsic form CONJECTURE) |

### Proofs

**G_θ′ (PROVED).**
- The constituents |·|^{±iθ′} are unitary characters of C. So OS1 and OS2a hold, and E commutes with |x|^{iθ′}.
- OS3: G_θ′ has a degree-2 Euler product Π_p L_p(s+iθ′) L_p(s−iθ′).
- OS4: four poles, all with unitary boundary characters. The twisted MF is a sum of norms (ARP G.3c).
- OS2b: Q_G(f) = 2Q_ζ(f cos θ′u) + 2Q_ζ(f sin θ′u) (twist identity; ARP G.4e, replayed S102). So **OS2b(G_θ′) ⇔ OS2b(ζ)**, and every PROVED part transfers.
- The must-pass control passes. R3: G_θ′'s weights 2Λ(n)cos(θ′ log n) are signed.

**F_θ (PROVED).**
- OS1 holds: E commutes with multiplication by |x|^{±θ}. But the C-action on each constituent is through the non-unitary quasi-characters |·|^{±θ}.
- **OS2a fails on the GL(1) constituents.**
  - The constituent factor is γ_p^{(θ)}(s) = L_p(1−s−θ)/L_p(s+θ) = (1 − p^{−(s+θ)})/(1 − p^{−(1−s−θ)}).
  - At s = ½ + it with p^{−it} = 1 its modulus is (1 − p^{−½−θ})/(1 − p^{−½+θ}) < 1 for θ > 0. So |γ_p^{(θ)}| ≢ 1.
  - The product of the ±θ factors is unitary on Re s = ½ (ARP G.4e, replay (5); NUMERICAL 0.817 and 1.224 at p = 7 for orientation).
- **OS2b fails.**
  - F_θ has zeros at ½ ± θ + iγ, off the line.
  - By Lemma 5.2 (Q.5) with Θ_F = ½ + θ, λ_min(Q_{F,L}) ≤ −e^{(2θ−ε)L} along a sequence of L.
  - Structurally, Q_F = 2Q_ζ(f cosh θu) − 2Q_ζ(f sinh θu) = 2Re B_ζ(e^{θu}f, e^{−θu}f) is a Kreĭn polarization, not a norm (ARP G.4e, replay (2)).
- OS3 holds (degree-2 Euler product). OS4 holds in the MF sense: κ·2cosh(θr) = e^{−(½−θ)|r|} + e^{−(½+θ)|r|}, a sum of Cauchy kernels (G.3c).
- **Verdict.** The expectation "fails OS2 only, via the constituent filter" is correct **provided OS2a is imposed on GL(1) constituents**.
  - If OS2a were imposed only on the product (the GL(2) complementary series is unitarizable), F_θ would satisfy OS0, OS1, OS2a, OS3 and OS4 and fail only OS2b.
  - In that reading F_θ is the cleanest proof that OS2b is independent of the other axioms, and that positivity of the weights (F_θ's are 2Λ(n)cosh(θ log n) > 0) does not give OS2b (R3).
  - Either way, the failure sits in OS2.

**DH (PROVED).** DH(s) = ½(1−iκ) L(s,χ) + ½(1+iκ) L(s,χ̄), with χ mod 5, χ(2) = i, κ = (√(10−2√5) − 2)/(√5 − 1) ≈ 0.28408.
- Its coefficients are a(n) = Re χ(n) + κ Im χ(n).
- So a(2) = κ, a(3) = −κ (χ(3) = χ(2)³ = −i), and a(6) = χ(1) = 1. Hence a(2)a(3) = −κ² ≈ −0.0807 ≠ 1 = a(6).
- A Dirichlet series with a(1) = 1 has an Euler product iff a is multiplicative. So **DH has no Euler product, and OS3 fails.** (Checked numerically this session.)
- OS1 and OS2a hold sector-wise: χ and χ̄ are unitary, and DH's FE is by construction.
- OS2b fails: there are zeros in σ > 1 (Davenport–Heilbronn 1936; spec Cor. 5.3 via Lemma 5.2).
- OS4: L(s,χ) is entire, so there is no boundary term, and OS4 holds only vacuously. **DH also has no pole state (no R4).** That is an addition to the expectation, not a contradiction.

**DMV / P_X (PROVED, with a refinement).**
- The prime atoms of Δ_P beyond X sit at log p_j for the Beurling primes p_j, which are not the rational primes. Otherwise the system would equal ζ's past X.
- **Generic case (p_j ∉ ℚ).** No ℚ^*-sum realizes the atoms, so there is no lattice and **OS1 fails**.
- **Degenerate case (some p_j ∈ ℚ but not prime).** Then ℚ^* exists, but the local factors are not indexed by the places of ℚ, and OS3 fails instead.
- In both cases there is no Poisson summation, hence no FE (R2: "DMV, Euler product, no FE"). So the global part of OS2a, which is transported Poisson (G.4a), is undefined.
- **R1.** P_X agrees with ζ on every window L ≤ (log X)/2 (Goddess of PNT, `team/prime-number-theory/finite_windows_dmv.md`). So every axiom restricted to bounded windows holds for P_X: the axioms must be imposed for all L.
- The expectation "fail OS1 (no ℚ^* lattice)" is **verified** in the generic case. In the degenerate case the failure moves to OS3.

**Epstein Z(s) for x² + 5y² (h(−20) = 2) (PROVED).**
- Let r(n) be the number of representations, and a(n) = r(n)/2. Then r(1) = 2, r(2) = r(3) = 0, r(6) = 4 (6 = 1 + 5). So a(2)a(3) = 0 ≠ 2 = a(6): not multiplicative, hence no Euler product, and **OS3 fails**. (Checked numerically this session.)
- OS1 holds over K = ℚ(√−5). Z = (w/2)(ζ_K + L(s,χ_{−4}) L(s,χ_5)) is a sum of two C_K-sectors with unitary genus character, so OS2a holds too.
- OS4 holds: the pole of ζ_K gives MF ≥ 0 (G.3c).
- OS2b fails: off-line zeros (Davenport–Heilbronn 1936; spec quartet numerics, NUMERICAL).
- The expectation was "determine". Determined: **OS3**, agreeing with ARP G.4e ("breaks at G.1d/G.4b").

**H₋ₜ (t > 0) (PROVED for the natural realization; the intrinsic form is CONJECTURE).**
- H₋ₜ(z) = ∫ e^{−tu²} Φ(u) cos(zu) du, where Φ is the theta-side (E-side) kernel and u ∝ log|g|.
- Lift the multiplier to 𝔸: m_t(x) = e^{−tc·log²|x|}. It is ℚ^*-invariant (|qx| = |x|), so M_t E(φ⁰) = E(m_t φ⁰) formally. So OS1 holds formally, and OS2a holds (m_t is even under inv).
- **OS3 fails in this realization.**
  - log|x| = Σ_v a_v with a_v = log|x_v|, so log m_t = −tc(Σ_v a_v)².
  - Then ∂_{a_v}∂_{a_w} log m_t = −2tc ≠ 0 for v ≠ w. A function is a product of one-place factors iff all these mixed derivatives vanish.
  - So m_t φ⁰ is not a pure tensor, and is not even Bruhat–Schwartz: the deformation couples ∞ with every p.
  - This is the precise form of ARP G.4e's "γ_∞ is not Tate's".
- OS2b fails: H₋ₜ has non-real zeros for t > 0 (Rodgers–Tao 2020, Λ_dBN ≥ 0; with de Bruijn 1950 the set of t where all zeros are real is [Λ_dBN, ∞)).
- **CONJECTURE:** no other realization gives H₋ₜ an Euler product.
- OS4 is convention-dependent: the deformation acts on the entire Ξ, and G.3c's "pole kept" is a convention.

**Expectation summary.** No expectation was wrong. Three need qualifiers:
- (i) F_θ fails "OS2 only" only under the GL(1)-constituent convention for OS2a. Otherwise it fails OS2b only, which is still OS2.
- (ii) DMV/P_X fail OS1 in the generic case. If the Beurling primes are rational the failure moves to OS3.
- (iii) DH fails OS3 and also lacks the pole state (OS4 vacuous, no R4).

---

## Q.3 OS0–OS4 is RH-equivalent; OS2b is THE GAP; the proved OS2′

**Theorem Q.3 (PROVED, given the citations).** For ζ: (OS0 ∧ OS1 ∧ OS2a ∧ OS2b ∧ OS3 ∧ OS4) ⇔ RH.

Proof:
- OS0, OS1, OS2a, OS3 and OS4 are PROVED unconditionally (Q.1).
- OS2b, Q_L(f) ≥ 0 for all L and all f ∈ 𝒯_L, is exactly Weil's criterion, ⇔ RH (Weil 1952; Bombieri 2000).
- (⇐) The RH direction also gives the OS reconstruction. Q_L ≥ 0 makes Δ a positive-definite distribution on C. The GNS construction gives the physical Hilbert space with a unitary C-action: PROVED-LEAN `PositiveDefiniteGNS.gns_from_positive_definite` for an arbitrary group, algebraic level.
- The zero-side form of the equivalence is PROVED-LEAN `WeilPositivityCriterion.rh_iff_weil_pairedForm_nonneg`, which is listed in ARP §0.3 and was checked this session. Scope: it is the finite paired form over the zero set with the FE involution, **not** the explicit-formula window form. The explicit formula is not in GPPVerify. ∎

R-line: R1, the equivalence needs all L | R2, both halves enter, through OS1+OS2a (FE) and OS3 (Euler) | R3, OS2b is not a weight-sign statement | R4, OS4 | R5, Q.2 | R6, exact.

**Hence OS2b is THE GAP, identical to G.4d.** With Q_L = MF_L + Fl_L exactly (G.3d, ψ-fluctuation; R6):

> **GAP (G.4d = OS2b).** For every L > 0 and every real f ∈ C_c^∞((−L, L)):
>   Fl_L(f) ≥ −MF_L(f),
> where MF_L(f) = ∫∫ f(u) f(v) e^{−|u−v|/2} du dv and Fl_L(f) = W_∞(g) − ∫_0^∞ (g(u) + g(−u)) dE_Λ(u), with g = f ∗ f̃ and E_Λ = Σ_n Λ(n) n^{−1/2} δ_{log n} − e^{u/2} du. This must hold **uniformly in L**. Equivalent to RH.

**OS2′ (mean-field OS2; PROVED).** For every L and every f ∈ 𝒯_L:
- (Toeplitz) MF_L(f) = ∫ |f̂(t)|² dt / (t² + ¼) ≥ 0 (up to the 1/2π normalization of f̂). The Fourier density of e^{−|r|/2} is 1/(t² + ¼) > 0.
- (Half-line Hankel, i.e. genuine OS reflection) ∫∫_{u,v≥0} f(u) f(v) e^{−(u+v)/2} du dv = (∫_0^∞ f e^{−u/2})² ≥ 0. This is rank one, and the vector is the pole state e^{−u/2}.
- Citations: Astra's Ward identity, replayed by Grok (ARP G.3a). PROVED-LEAN for the Gram shape: `ArithmeticOSGram.rankOneGram_nonneg`. The specific identity is Lean ask G.6 #5, still open.
- With OS4 consumed (S103), **OS2′ is all of OS2 that the pole and the PNT mean can give.** It is fake-blind: P_X, DMV, F_θ, Epstein and H₋ₜ satisfy it, DH vacuously (G.3c). So it is RH-silent.

R-line: R1, uniform | R2, it uses neither the Euler product nor the FE beyond the pole | R3, ok | R4, **yes** | R5, fake-blind, as it must be for a PROVED statement | R6, exact.

---

## Q.4 Sign check: range vs cokernel (the E-transfer); verdict

Box copy: `/workspace/queen/work/arp/sign_check.md`. Sources read this session: Connes 1999, Selecta Math. 5, 29–106 (arXiv math/9811068), §III, §V Thm 3, §VII Thm 4, §VIII (16), Thm 5 with proof, Lemma 3. Meyer 2005 (Duke 127, 519–595; arXiv math/0311468) is used for the realization only, not for positivity.

**Connes' positivity (PROVED, Connes 1999 VIII (21)–(25)).**
- Let S_Λ be the window projection on L²(C) (|g| ∈ [Λ^{−1}, Λ]). Let B_{Λ,0} be the span of f ∈ 𝒮(𝔸)₀ with f and f̂ vanishing for |x| > Λ, and Q′_{Λ,0} = E Q_{Λ,0} E^{−1}.
- E(f)(g) = E(f̂)(1/g) (Θ transported). Hence E(B_{Λ,0}) ⊂ S_Λ, i.e. Q′_{Λ,0} ≤ S_Λ.
- Therefore Δ_Λ(h) := Tr((S_Λ − Q′_{Λ,0}) V(h)) is of positive type for every Λ, with Tr(S_Λ V(h)) = 2h(1) log′Λ exactly.
- Scope: the proof is written for function fields. For ℚ, Connes states the analogue with prolate B_Λ, which is only approximately space-and-band limited. That is stated by Connes and not re-checked here.

**Verdict.**
1. **Wrong sign for "additive positivity transfers through E".** The Poisson comb's positive-definiteness (G.4a) is positivity of the **range** Q′_{Λ,0} ≥ 0. The range enters Δ_Λ with a **minus** sign. So additive positivity yields only Weil ≤ (window term), the trivial side. It does not bound the Weil form from below.
2. **The positive object is the cokernel complement S_Λ − Q′_{Λ,0}.** Its positivity is free: a projection inequality that uses only the Θ-symmetry of the doubly localized space.
3. **So the GAP is an identification, not a positivity.**
   - Weil = lim_Λ Δ_Λ is Connes' global trace formula (16), ⇔ RH (Thm 5).
   - Exactly (Lemma 3, function fields; stated for ℚ): lim Δ_Λ(f) = Σ_ρ m_ρ ∫ ĝ dμ_ρ, where μ_ρ is the harmonic (Poisson-kernel) measure of ρ with respect to the critical line.
   - So the residual δ_Λ(f) := Δ_Λ(f ∗ f^*) − Q_L(f) tends to Σ_{ρ off-line} m_ρ [∫ ĝ dμ_ρ − ĝ(γ_ρ)].
   - Off-line zeros appear as Poisson-smeared resonances on the positive side and as point values on the Weil side. **Sign-correct GAP in Connes' language:** δ_Λ(f) → 0 for every f; one-sidedly, lim sup δ_Λ(f) ≤ lim Δ_Λ(f ∗ f^*).
4. **Consistency with G.4d.** "Fl_L ≥ −MF_L uniformly" and "δ_Λ → 0" are both RH-equivalent, hence consistent. No identification of (MF, Fl) with (Δ_Λ, δ_Λ) is claimed (CONJECTURE-free: none is made). **Two wording fixes are proposed for G.4d:**
   - (a) The positive part is Tr((S_Λ − Q′_{Λ,0}) V(h)V(h)^*), a trace over the complement of the range. It is not "‖P_Λ^⊥E(φ)‖²", which is the norm of a range vector, the wrong space.
   - (b) The 2 log Λ term is not part of δ_Λ: it cancels exactly against Tr(S_Λ V(h)).

R-line: R1, the GAP is the Λ → ∞ uniformity, i.e. Connes' own "uniformity" dead end (§0.2) | R2, E is the single operator carrying both halves (G.4b) | R3, ok | R4, B_{Λ,0} ⊂ 𝒮(𝔸)₀ removes the pole; Connes' D term restores it ("the nuance between B_Λ and B_{Λ,0}") | R5, DMV/P_X have no E; DH/Epstein have E only sector-wise; F_θ's Δ_Λ is not of positive type on the constituents | R6, the identities are exact; only the limit is open.

---

## Q.5 Where an OS2b failure shows up: quartet (1,1), Lemma 5.2, Turán

- **Quartet identity (PROVED, elementary given the explicit formula; spec §5, corrected S100).** An off-line pair {ρ₀, 1−ρ̄₀}, ρ₀ = ½ + β + iγ₀, contributes pair(f) = 2|⟨f, χ_L⟩|² − 2|⟨f, φ_L⟩|². Here χ_L = e^{−iγ₀x} cosh βx and φ_L = e^{−iγ₀x} sinh βx: **signature (1,1)**, the same shape as the pole's cosh/sinh pair.
  - Abstract finite form: PROVED-LEAN `QuartetNegativeIndex.quartet_negative_family`, `negativeIndex_le_one_forces_line`, `NoGo.QuartetSignature.Q_neg_plane`, `onAxis_psd`.
  - So every OS2b failure has a quartet-sized negative plane.
  - R-line: R1 n/a | R2 n/a | R3 ok | R4, it mirrors the pole's (1,1) | R5, it locates DH's and Epstein's negativity (NUMERICAL residual 1.3% and 0.5%) | R6, exact.
- **Lemma 5.2 (PROVED, qualitative; Landau).** η_D ≥ 2Θ_D − 1, with η_D = lim sup L^{−1} log⁺(−λ_min Q_{D,L}) and Θ_D = sup Re ρ.
  - Landau's theorem is PROVED-LEAN `Analysis.LandauDirichlet.landau`.
  - Consequence: an OS2b failure grows like e^{(2Θ−1)L}. Any proof of OS2b must be uniform in L (R1) and exact (R6), since ζ's margin leaves no room.
  - R-line: R1 yes | R5, DH has η ≥ 1.2407 (Cor. 5.3), F_θ has η ≥ 2θ | R6, an inequality, but it is used only to *detect* failure.
- **Turán Theorem A (rightmost attained target): PROVED in sketch** (`2026-10-22_turan_lemma52.md`). λ_min(Q_L) ≤ −(1−ε) e^{2βL}/(2β) on a relatively dense set of L. The only unwritten piece is explicit-formula validity for Lipschitz g, which is standard. **Theorem B (non-rightmost target): PARKED** (CONJECTURE, three gaps stated in the draft).

---

## Q.6 Cross-references (read-only; not re-derived)

- **Astra, thermal circle** (bridge `work/2026-10-26_astra_prime_os_thermal_circle.md`). The per-prime thermal-circle S_p is OS-positive, but it is RH-silent: per-prime positivity is an OS2′-type statement, and the global OS2b is not implied. This is consistent with G.1c (no place-by-place positivity for Weil's prime part). Related PROVED-LEAN: `CriticalRiggedGram.critical_gram_nonneg` (the ζ(1+t+u) OS kernel, untruncated, also RH-silent).
- **Astra, mirror-OS superselection** (`work/2026-10-26_astra_mirror_prime_vacuum_superselection.md`). Π_p R_p(σ) = 0 for σ ≠ ½ is PROVED (Kakutani-type; Lean ask G.6 #6). It is zero-blind: it selects the line for a σ-family of product states, not the zeros. It fails R2, since it uses no FE/Poisson, and R5, since F_θ- and DMV-type products give the same dichotomy. Abstract analogue: PROVED-LEAN `HolographicReflectionPositivity.reflectionForm_nonneg_iff` (a reflection form is positive iff the involution fixes every point).
- **Astra, doubled prime** (`work/2026-10-27_astra_doubled_prime_os_sewing_audit.md`). The raw two-opposite-prime OS form is indefinite exactly, so it is not OS-positive. Any compensator must use the full prime tail. This is the operator form of R1 and G.4d(ii)/T2.
- **Turbo DOM±.** These forms restate window positivity (OS2b) and are not re-derived here.

---

## Q.7 Proposed edits to ARP (for the lead; nothing pushed to main)
1. G.4d: the two wording fixes in Q.4 (a) and (b).
2. G.5 (T1): replace "it cancels against θ's −√x" by the S103b ψ-repair wording.
3. §0.3: confirm `rh_iff_weil_pairedForm_nonneg` is cited with its scope, i.e. finite paired form over the zero set, not the explicit-formula window form.

## Q.L Lean names used (all present in GPPVerify `5f54359`; no `sorry`/`axiom` in the cited files)

| Lean name | file |
|---|---|
| `primeSide_eq_truncation`, `primeSide_eq_zero_of_support_lt_log_two` | RiemannHypothesis/WeilSupportLadder |
| `ArithmeticOSReflection.osReflection`, `osReflection_eq_self_iff` | CelestialHolography |
| `ArithmeticOSGram.rankOneGram_nonneg` | CelestialHolography |
| `HaarSquareNotReflection.hankel_neg`, `GaussianHeatOS.gaussian_hankel_det_neg` | NoGo |
| `euler_factor_bridge` | RiemannHypothesis/PadicEulerFactorBridge |
| `RationalIdeleGroup` (type only) | RiemannHypothesis/IdeleGroup |
| `PositiveDefiniteGNS.gns_from_positive_definite` | RiemannHypothesis |
| `WeilPositivityCriterion.rh_iff_weil_pairedForm_nonneg` (scope as stated) | RiemannHypothesis |
| `quartet_negative_family`, `negativeIndex_le_one_forces_line` | RiemannHypothesis/QuartetNegativeIndex |
| `Q_neg_plane`, `onAxis_psd` | NoGo/QuartetSignature |
| `landau` | Analysis/LandauDirichlet |
| `critical_gram_nonneg` | RiemannHypothesis/CriticalRiggedGram |
| `reflectionForm_nonneg_iff` | RiemannHypothesis/HolographicReflectionPositivity |

Not cited, per instructions:
- `ArithmeticTemperedPrincipalSeries`;
- the `HaarPositivityWeil` `open_*` True-placeholders;
- `ArithmeticCoPoissonOrientationTemplate`.

## Q.8 Off-line residual theorem and the R6 identity
See [`ARP_queen_Q8.md`](ARP_queen_Q8.md) (2026-10-28). It also records corrections: (c1) the Gram step of Connes' Lemma 3; (c2) for ℚ, Q′ ≤ S_Λ is not exact with the prolate B_Λ, so §Q.4 item 2's "free positivity" is exact only for function fields, and the exact ℚ object is Δ^♮_Λ; (c3) the F_θ OS2a wording is reconciled with §P.3.

## References
- Weil 1952, "Sur les 'formules explicites' de la théorie des nombres premiers".
- Bombieri 2000, "Remarks on Weil's quadratic functional".
- Barner 1981.
- Tate 1950 (thesis).
- Connes 1999, Selecta Math. 5, 29–106, arXiv math/9811068.
- Meyer 2005, Duke Math. J. 127, 519–595, arXiv math/0311468.
- Connes–Consani–Moscovici, arXiv:2310.18423 (Ann. Funct. Anal. 2024).
- Connes–Consani 2021, Selecta 27:77.
- Yoshida 1992.
- Davenport–Heilbronn 1936.
- de Bruijn 1950.
- Rodgers–Tao 2020 (Forum Math. Pi 8, e6).
- Diamond–Montgomery–Vorhauer 2006 (Math. Ann. 334).
- Team notes: ARP §G; spec `2026-10-22_noncircular_positivity_spec.md` §5; `2026-10-22_turan_lemma52.md`; the Astra bridge notes above.
