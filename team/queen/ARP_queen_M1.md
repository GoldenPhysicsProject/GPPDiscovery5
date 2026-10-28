# M1 — The prolate leak (1 − S_Λ)Q′_Λ: a reduction (Queen of Mathematics; appendix to ARP_queen_Q8.md)

Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
Queen of Mathematics (RH sub-team), 2026-10-28. Against GPPDiscovery5 main `21067e4` (G.4d stated with Δ^♮_Λ; "the unsharpened
form only given M1") and team/queen `663cedf` (§Q.8 after the Golden Eye audit fixes).

**RH is OPEN. M1 proves nothing about RH and cannot.** Analytic only: no numerics were run for this note.
Labels: PROVED / PROVED-conditional (an implication whose hypotheses are named) / CONJECTURE / GAP; HEURISTIC marks orientation.
R-lines as in §Q (R1 uniform in L; R2 Euler product + FE; R3 not via weight signs; R4 soft mode/pole; R5 fakes; R6 exact identity).

## M1.0 What M1 is, and what it is for (aligned with main 21067e4)

- §Q.8.3 defines Δ^♮_Λ(f ∗ f^*) = Tr(A_Λ V(f)V(f)^*), A_Λ = S_Λ − S_ΛQ′_ΛS_Λ; it is positive type at every Λ (PROVED). Main's G.4d
  now states the GAP with Δ^♮_Λ. Connes' unsharpened Δ_Λ = Tr((S_Λ − Q′_Λ)V(h)) is not positive type over ℚ (§Q.8 (c2)).
- **The leak functional.** Lk_Λ(h) := Δ_Λ(h) − Δ^♮_Λ(h) = −Tr((Q′_Λ − S_ΛQ′_ΛS_Λ)V(h)), and
  Q′ − SQ′S = (1−S)Q′(1−S) + SQ′(1−S) + (1−S)Q′S, so Lk_Λ is controlled by (1 − S_Λ)Q′_Λ alone.
- **M1 (statement used by main 21067e4).** For each fixed real f ∈ C_c^∞(ℝ_u) (ζ-sector), Lk_Λ(f ∗ f^*) → 0 as Λ → ∞.
  Given M1, lim Δ_Λ and lim Δ^♮_Λ coincide (when either exists), which is exactly what licenses main's parenthetical
  "lim Δ_Λ = Weil ⇒ RH holds only given M1". M1 is also the input that lets the M2 argument treat S_ΛQ′_ΛS_Λ as asymptotically a
  projection (its defect from a projection is SQ′(1−S)Q′S, trace ≤ ℓ_Λ below).
- **M1 is not RH content.** It compares two operator traces built from the same E, S_Λ, Q′_Λ; no zero of ζ enters it. Even
  M1 + M2 + M3 (= Conj Q.8.C) leave RH ⇔ Res ≡ 0 (Thm Q.8.A, a restatement). R-line: R1 see M1.5 | R2 not engaged | R3 n/a |
  R4 the pole data sit in the leak region (Lemma Q.8.2) and are killed exactly by φ(0) = φ̂(0) = 0 | R5 fake-blind (any datum with an
  E has the same leak) | R6 n/a.

## M1.1 Parameters (Golden Eye split; Spectrum's convention: Λ and L independent)

- **Λ** = Sonin/prolate cutoff: S_Λ = 1_{|u| ≤ log Λ} on L²(ℝ_u) (ζ-sector, u = log|g|); B_Λ = span of the even prolates ψ_n,
  n ∈ I_Λ, of the time–band limiting pair (P_Λ, P̂_Λ) on L²(ℝ)^ev, bandwidth parameter c = 2πΛ² (CC 2106.01715 §3:
  ψ_{m,λ}(x) = PS_{2m,0}(2πλ², x/λ)); B_{Λ,0} = {φ ∈ B_Λ : φ(0) = φ̂(0) = 0}; Q′_Λ = orthogonal projection onto K_Λ := E(B_{Λ,0}).
  Placement (P): ψ_n truncated to [−Λ, Λ] (Prop Q.8.3: (P) ≡ (P̂) on real f).
- **L** = support half-width of the test function: f ∈ C_c^∞((−L, L)), h = f ∗ f^* supported in (−2L, 2L).
  **L is independent of Λ** (Spectrum's LAMBDA_L dictionary convention; the earlier guess L = log Λ is withdrawn). Conj Q.8.E's
  rule Λ = e^L is one diagonal choice, treated separately in M1.5.
- **Support past p = 2** means 2L > log 2, where the prime side of the Weil form Q_L first becomes nonzero
  (primeSide_eq_zero_of_support_lt_log_two). That is Celestial's track; M1.6 states the interface only.
- **CC 2021 (2006.13771)** is a third, separate setting: their Thm 1 / Thm 6.11 restrict supp g ⊂ [2^{−1/2}, 2^{1/2}] (a support
  condition, i.e. an L-condition with 2L ≤ log 2 for g ∗ g^*) **and** fix the cutoff at Λ = 1 (P_1, P̂_1). The two are unrelated
  parameters; their positivity at Λ = 1 says nothing about Λ → ∞ and vice versa.
- **Verdict on the split.** Every bound in this note is a bound in **Λ** at fixed or controlled L. M1 **pushes Λ only**: it is
  **RH-silent** (archimedean prolate data, plus the integers through E's ℤ-sum and |ζ(½+is)| on the line, M1.4). It does **not**
  push the support past log 2, and it is the same statement for every L, so it **cannot carry RH content**. R2 is not engaged.

## M1.2 Prior-art verdict (read this session; cite, don't duplicate)

No paper on the reading list bounds (1 − S_Λ)Q′_Λ, in any norm, for Λ > 1 in the E-transported setting. Specifically:
- **Connes 1999 (math/9811068) §VIII (31)–(33):** introduces B_Λ (n ≲ 4Λ², counting both parities), notes no nonzero f has f and f̂
  both supported in [−Λ, Λ], and states the ℚ analogue of Lemma 3 without proof. No leak bound.
- **CC 2106.01715 §3 (lines ~1395–1475 of the text):** for λ > 1 uses the Slepian–Pollak "approximate intersection" P_λ ∩_0 P̂_λ
  (≈ 2λ² very small angle eigenvalues) and, quoting, "act[s] as if E(φ_n)(u^{−1}) = (−1)^n E(φ_n)(u)". **That "as if" is exactly the
  leak**; it is used without a bound. Thm 1.1 (spectrum on (Σ_μ E(S_0^ev))^⊥ = critical-line zeros) is an exact statement on the
  circle C = ℝ_+^*/μ^ℤ and does not involve the leak.
- **CC 2021 Selecta (2006.13771):** Λ = 1 throughout. Prop 2.2 gives the **exact** trace identity Tr(ϑ(f) P P̂ P) = ∫ f(ρ^{−1})(δ − τ) d^*ρ
  with δ explicit (Si formula (49)); Lemma 4.1 δ(ρ) = Tr(ϑ(ρ^{−1})P̂_1P_1); trace class via [P, ϑ(f)] trace class for Schwartz f
  (Prop 1.5(iv)); Thm 3/4.7 Tr(ϑ(f)S) = W_∞(f) + ∫fε through the pair-of-projections decomposition. These are the technical model
  for M1 (commutators with a cutoff are trace class and boundary-localized), but there is no E and no Λ → ∞.
- **CCM 2310.18423:** semilocal prolate W_λ, Sonin space = negative spectrum up to a finite discrepancy, Thm 1 (Mellin form (7) with
  the local Euler factors), Thm 2 (θ_S isomorphism of Sonin spaces). The E-images of prolates give the extremely small eigenvalues of
  the Weil form. No leak bound; but (7) is the source of the conditioning issue in M1.4 (H-κ).
- **Connes–Moscovici 2112.05500:** W_sa self-adjoint, commutes with F (Thm 2.6); semiclassical count of negative eigenvalues
  (Prop 4.2); Thm 6.1. Spectral, not a leak bound. Points to Slepian [14–16] and Osipov–Rokhlin–Xiao [13].
- **Prolate eigenvalue literature (archimedean; to be supplied in our normalization by Spectrum, citations to be verified there):**
  Slepian–Pollak 1961 and Slepian 1965 (plateau λ_n ≈ 1 for n < 2c/π, asymptotics); Landau–Widom 1980, J. Math. Anal. Appl. 77
  (plunge count #{n : ε < λ_n < 1−ε} ~ (1/π²) log c · log((1−ε)/ε)); Osipov–Rokhlin–Xiao, *Prolate Spheroidal Wave Functions of Order
  Zero* (Springer 2013) and Osipov 2013 (nonasymptotic upper bounds); Bonami–Karoui 2017 and Bonami–Jaming–Karoui 2021 (nonasymptotic
  decay beyond 2c/π); Karnik–Romberg–Davenport 2021 (improved plunge bounds).
- **Spectrum's branch (origin/team/spectrum-measure 872aa60):** S(b)4/S(c)4 (no slack) and GAP2 (Bombieri window count, which this
  note cites as **PROVED under H_fin**, though it does not use it). No leak bound and no prolate tail bound yet; LAMBDA_L_dictionary.md
  has not landed (only the Λ/L independence is used here, per the lead).

## M1.3 Theorem M1.A (leak estimate in terms of one number ℓ_Λ; PROVED)

Define the **leak** ℓ_Λ := ‖(1 − S_Λ)Q′_Λ‖²_HS = Tr(Q′_Λ(1 − S_Λ)Q′_Λ) = Σ_j ‖(1 − S_Λ)e_j‖² for any orthonormal basis (e_j) of K_Λ
(= Σ sin²θ_k over the principal angles θ_k between K_Λ and ran S_Λ; basis-free).

**Theorem M1.A.** For every Λ > 1 and every f ∈ L¹ ∩ L²(ℝ_u), with h = f ∗ f^*:

  |Lk_Λ(h)| ≤ ‖f‖_1² ℓ_Λ + 2 ‖f‖_1 ‖f‖_2 (2 log Λ)^{1/2} ℓ_Λ^{1/2}.

*Proof.* Put V = V(f), so V(h) = VV^* and ‖V‖ ≤ ‖f‖_1. Q′ has finite rank, so all traces below are of finite-rank operators.
- Tr((1−S)Q′(1−S)VV^*) = ‖Q′(1−S)V‖²_HS ≤ ‖Q′(1−S)‖²_HS ‖V‖² ≤ ℓ_Λ ‖f‖_1² (and it is ≥ 0).
- Tr(SQ′(1−S)VV^*) = ⟨Q′(1−S)V, Q′SV⟩_HS, so |·| ≤ ‖Q′(1−S)V‖_HS ‖Q′SV‖_HS ≤ ℓ_Λ^{1/2}‖f‖_1 · ‖Q′SV‖_HS.
- ‖Q′SV‖²_HS = Tr(V^*SQ′SV) = Tr(SQ′S·VV^*) ≤ Tr(S VV^*) because SQ′S ≤ S and VV^* ≥ 0. And Tr(S_Λ V(h)) = h(0)·|{|u| ≤ log Λ}|
  = 2‖f‖_2² log Λ (kernel h(u − v) on the window).
  *(Grok audit: checked. SQ′S ≤ S because Q′ is an orthogonal projection, so SQ′S = (Q′S)^*(Q′S) ≤ S·S = S. V(h) = VV^* ≥ 0. Tr(S_Λ V(h)) is finite by Mercer, since the kernel is continuous on the compact window. Constants 2(2 log Λ)^{1/2} confirmed.)*
- The term (1−S)Q′S is the adjoint case and gives the same bound. Summing the three pieces gives the claim. ∎

**Corollary M1.B (sufficient condition for M1; PROVED-conditional).** If ℓ_Λ log Λ → 0, then M1 holds, with the explicit rate
|Lk_Λ(f∗f^*)| ≤ ‖f‖_1² ℓ_Λ + 2^{3/2}‖f‖_1‖f‖_2 (ℓ_Λ log Λ)^{1/2}. A weaker, f-localized sufficient condition is
ℓ_Λ(f) := ‖V(f)^*(1−S_Λ)Q′_Λ‖²_HS = o(1/log Λ) for each f (same proof).

R-line: R1 uniform in f through ‖f‖_1, ‖f‖_2 only (so in L, see M1.5) | R2 not engaged | R3 n/a (operator inequality SQ′S ≤ S, not a
weight sign) | R4 n/a | R5 fake-blind | R6 n/a.

Remark (the log Λ is sharp for this method). The cross term pays (Tr S_ΛV(h))^{1/2} ~ (log Λ)^{1/2}, the same 2‖f‖² log Λ that
appears in Conj Q.8.C; ℓ_Λ → 0 alone does not suffice for this argument.

## M1.4 Theorem M1.C (reduction of ℓ_Λ to prolate data; PROVED-conditional) and hypotheses for Spectrum

**Lemma M1.1 (exact leak formula; PROVED given Lemma Q.8.2 on B_{Λ,0}).** Normalize E(φ)(u) = e^{u/2} Σ_{n≠0} φ(n e^u) and
Fφ(ξ) = ∫ φ(x)e^{−2πixξ}dx. For φ ∈ B_{Λ,0} in placement (P):
- E(φ)(u) = 0 for u > log Λ (all |n e^u| > Λ), so the leak sits at u < −log Λ only.
- There, Lemma Q.8.2 with φ(0) = φ̂(0) = 0 (pole terms vanish) gives E(φ)(u) = E(Fφ)(−u), so
  ‖(1 − S_Λ)E(φ)‖² = ∫_Λ^∞ |Σ_{n≠0} Fφ(ny)|² dy.
- Fφ on |ξ| > Λ is the **out-of-band tail** of φ; for φ = ψ_n it equals μ_n ψ_n(ξ) (whole-line prolate), with
  ∫_{|ξ|>Λ}|Fψ_n|² = 1 − λ_n (λ_n = |μ_n|², ψ_n L²[−Λ,Λ]-normalized).
- Caveat (GAP-minor): Lemma Q.8.2 is proved on S_0^ev; extending it to truncated prolates (in L²(X)_0, not Schwartz) is by
  density; the sum Σ_n Fφ(ny) converges only conditionally (see below), so the extension must be done in L²(dy), not pointwise.

**Lemma M1.2 (jump/sawtooth split; PROVED).** For even φ ∈ C²[−Λ, Λ] extended by 0:
Fφ(ξ) = φ(Λ) sin(2πΛξ)/(πξ) + R_φ(ξ), with |R_φ(ξ)| ≤ (|φ′(Λ)|/2 + ‖φ″‖_{L¹}/4)/(π²ξ²) (two integrations by parts).
For y > Λ with Λy ∉ ℤ, Σ_{n≥1} sin(2πnΛy)/n = π(½ − {Λy}) (sawtooth), so the jump part of Σ_{n≠0}Fφ(ny) equals
(2φ(Λ)/y)(½ − {Λy}) and has modulus ≤ |φ(Λ)|/y. The triangle inequality alone would give Σ 1/n = ∞ here; the cancellation is
the only place the ℤ-structure of E enters the leak, and it enters only through the sawtooth.

**Archimedean inputs (RH-silent; to be filled by Spectrum in this normalization).** Let (ψ_n)_{n∈I_Λ} be the even prolates retained in
B_Λ, L²[−Λ,Λ]-orthonormal.
- **(S1) index set.** |I_Λ| and the split of I_Λ into plateau (λ_n ≥ 1 − ε) and plunge (ε < λ_n < 1 − ε), with the plunge count
  ≤ A log Λ · log(1/ε) + B, explicit A, B (Landau–Widom 1980; nonasymptotic: Osipov / ORX, Karnik–Romberg–Davenport).
- **(S2) endpoint sum.** ε_Λ := Σ_{n∈I_Λ} |ψ_n(Λ)|². (Classical route, to be verified by Spectrum: for unit-interval prolates,
  ∂_c λ_n = (2/c) λ_n ψ_n(1)², so endpoint values are controlled by the c-derivative of λ_n; they are tiny on the plateau.)
- **(S3) remainder tail.** With R_n := R_{ψ_n} from Lemma M1.2, ρ_Λ(X) := Σ_{n∈I_Λ} ∫_{|ξ|>X} |R_n(ξ)|² dξ, and the hypothesis
  ρ_Λ(X) ≤ ρ_Λ · (Λ/X)^{1+2σ} for all X ≥ Λ, some σ > 0. (Supporting data: 1 − λ_n and the Bonami–Karoui tail bounds.)

**Arithmetic input (critical-line only; RH-neutral, NOT archimedean; owner: Queen, currently GAP).**
- **(H-κ) conditioning of E on B_{Λ,0}.** κ_Λ := sup_{0≠φ∈B_{Λ,0}} ‖φ‖²_{L²(ℝ)} / ‖E(φ)‖²_{L²(du)}.
- By Mellin–Plancherel (PROVED on S_0^ev; Connes 1999 §III, cf. CCM (7) with all primes):
  ‖E(φ)‖² = 2 · [∫ |ζ(½+is)|² |Mφ(½+is)|² ds / ∫ |Mφ(½+is)|² ds] · ‖φ‖², where Mφ(z) = ∫_0^∞ φ(x)x^z d^*x.
- So 1/(2κ_Λ) is the minimum, over the Mellin profiles of B_{Λ,0}, of the |Mφ|²-average of |ζ(½+is)|². It is finite for each Λ
  (E is injective on the finite-dimensional B_{Λ,0}: a nonzero Mφ is not supported on the discrete zero set). It involves only ζ on
  the critical line, so no RH input or output; but it is arithmetic, and it is **not** something Spectrum's archimedean bounds supply.
  HEURISTIC: B's Mellin profiles resolve s to ~1/log Λ, comparable to the zero spacing at height ~Λ², so κ_Λ can be large; how
  large is the open point.
- *Grok audit (clarification; label unchanged, GAP).*
  - "RH-neutral" means (H-κ) involves only ζ on the line, and no implication to or from RH is known. It does not mean independent of RH.
  - κ_Λ is a **minimum over a family** of weights |Mφ|² (the whole of B_{Λ,0}), not a mean value against one fixed weight. Fixed-weight or mollified second-moment asymptotics (Balasubramanian–Conrey–Heath-Brown type) therefore do not bound it directly. They give averages, and the minimizing profile can adapt to ζ.
  - The minimizer can concentrate its Mellin mass near critical zeros. So a lower bound on 1/κ_Λ most likely needs zero-spacing or local-density information on the line, at resolution ~1/log Λ near height ~Λ². Zero-free-region information off the line does not help here.
  - *Follow-up (Queen, ARP_queen_Hkappa.md):* unconditionally κ_Λ ≤ Λ^B (two-constants + Kovrijkine), with no zero-spacing input; the local information needed is thickness of a superlevel set of |ζ| on the line. Λ^{o(1)} reduces to (Z-avg) there.

**Theorem M1.C (PROVED-conditional on (S2), (S3); κ_Λ always finite).**
  ℓ_Λ ≤ 2 κ_Λ ( ε_Λ / Λ + 4 ζ(1+σ)² ρ_Λ ).

*Proof.* Let (φ_j) be an E-orthonormal basis of B_{Λ,0}, so e_j = E(φ_j) is an orthonormal basis of K_Λ. By Lemma M1.1,
ℓ_Λ = Σ_j ∫_Λ^∞ |Σ_{n≠0} Fφ_j(ny)|² dy. By Lemma M1.2 and (a+b)² ≤ 2a² + 2b²:
- jump part ≤ 2 Σ_j |φ_j(Λ)|² ∫_Λ^∞ y^{−2} dy = (2/Λ) Σ_j |φ_j(Λ)|²;
- remainder: by Cauchy–Schwarz, (Σ_{m≥1}|R(my)|)² ≤ ζ(1+σ) Σ_m m^{1+σ}|R(my)|², and ∫_Λ^∞ Σ_m m^{1+σ}|R(my)|² dy
  = Σ_m m^σ ∫_{mΛ}^∞ |R|², so the remainder part ≤ 2·4ζ(1+σ) Σ_m m^σ Σ_j ∫_{mΛ}^∞|R_{φ_j}|².
- Each Σ_j q(φ_j) is the trace, in the E-inner product, of a positive form q on B_{Λ,0}. Since ⟨·,·⟩_{L²} ≤ κ_Λ⟨E·,E·⟩, it is
  ≤ κ_Λ × the L²-trace over B_{Λ,0}, which is ≤ the L²-trace over span{ψ_n : n ∈ I_Λ} = Σ_{n∈I_Λ} q(ψ_n) (B_{Λ,0} is a subspace;
  q ≥ 0). This gives Σ_j|φ_j(Λ)|² ≤ κ_Λ ε_Λ and Σ_j ∫_{mΛ}^∞|R_{φ_j}|² ≤ κ_Λ ρ_Λ m^{−1−2σ}.
- Σ_m m^σ m^{−1−2σ} = ζ(1+σ). ∎

**Corollary M1.D (M1 reduced; PROVED-conditional).** If κ_Λ (ε_Λ/Λ + ρ_Λ) log Λ → 0, then M1 holds (Theorem M1.A).

R-line for M1.C/D: R1 Λ-uniform constants, no L | R2 not engaged (E's ℤ-sum enters only via the sawtooth and ζ(1+σ); ζ on the
line enters only via κ_Λ) | R3 n/a | R4 the pole terms are removed exactly by B_{Λ,0}'s two conditions | R5 fake-blind | R6 n/a.

**Status of the reduction.** Archimedean side: (S1)–(S3) are standard-type prolate statements; expected (HEURISTIC, not used)
to make ε_Λ/Λ + ρ_Λ polynomially small in Λ, since the plateau states contribute exponentially small endpoint values and only
~log Λ plunge states contribute. The whole of M1 would then rest on **(H-κ): κ_Λ = o(Λ^a / log Λ)** for the a delivered by
Spectrum's bounds. (H-κ) is the genuine GAP, and it is where E's non-isometry (P.5–P.6, Meyer Lemma 5.3) bites. Labels:
(S1)–(S3) **inputs (Spectrum)**; (H-κ) **GAP (Queen)**; M1 itself **CONJECTURE**.

## M1.5 Uniformity in L (for Conj Q.8.E's diagonal Λ = e^L)

With f ∈ C_c^∞((−L, L)), ‖f‖_1² ≤ 2L‖f‖_2², so Theorem M1.A gives
|Lk_Λ(f∗f^*)| ≤ ‖f‖_2² (2L ℓ_Λ + 4 (L log Λ ℓ_Λ)^{1/2}) (PROVED).
- At fixed L (Λ, L independent): Corollary M1.B as stated.
- On Q.8.E's diagonal Λ = e^L, uniformly on the unit ball of 𝒯_L: needs ℓ_{e^L} = o(L^{−2}), i.e. ℓ_Λ = o(1/log²Λ). This is the
  only extra price of uniformity; it is still RH-silent.

## M1.6 Interface for Celestial (support past log 2)

- **What M1 needs from f:** only ‖f‖_1 and ‖f‖_2. Nothing in Theorems M1.A/M1.C changes when 2L > log 2. The leak ℓ_Λ is a
  property of (Λ, E, B_Λ) alone.
- **What M1 gives:** Δ_Λ(h) − Δ^♮_Λ(h) = O(‖f‖_1²ℓ_Λ + ‖f‖_1‖f‖_2(ℓ_Λ log Λ)^{1/2}) for every L, so on 2L > log 2 any identity
  proved for the unsharpened Δ_Λ transfers to the positive-type Δ^♮_Λ with this error, and conversely.
- **What M1 does not give:** any information about the prime side of Q_L. Matching lim Δ^♮_Λ to the Weil form when 2L > log 2
  is M2/M3 plus Celestial's W(k∗k^*) identity (Principal series); primes (R2) enter there, not in the leak.
- If Celestial couples Λ to L, use M1.5 (needs ℓ_Λ = o(1/log²Λ) along the coupling).

## M1.7 Summary table

| Item | Label |
|---|---|
| M1.A: |Lk_Λ| ≤ ‖f‖_1²ℓ_Λ + 2‖f‖_1‖f‖_2(2 log Λ)^{1/2}ℓ_Λ^{1/2} | PROVED |
| M1.B: ℓ_Λ log Λ → 0 ⇒ M1 | PROVED-conditional |
| Lemma M1.1 exact leak formula on B_{Λ,0} | PROVED on S_0^ev; extension to truncated prolates GAP-minor (density in L²(dy)) |
| Lemma M1.2 jump/sawtooth split | PROVED |
| M1.C: ℓ_Λ ≤ 2κ_Λ(ε_Λ/Λ + 4ζ(1+σ)²ρ_Λ) | PROVED-conditional on (S2), (S3) |
| (S1)–(S3) prolate inputs | inputs from Spectrum (RH-silent, archimedean) |
| (H-κ) conditioning of E on B_{Λ,0} | GAP (arithmetic, critical line only, RH-neutral) |
| M1 | CONJECTURE (RH-silent; pushes Λ only) |
| M1.5 diagonal Λ = e^L needs ℓ_Λ = o(log^{−2}Λ) | PROVED (reduction) |

Not cited anywhere in this note: ArithmeticTemperedPrincipalSeries, rh_iff_weil_pairedForm_nonneg, or any True-placeholder.
