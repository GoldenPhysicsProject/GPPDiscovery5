# GAP 2 — Bombieri's negative-eigenvalue count on the support window [−L, L]

Spectrum & Measure seat, 2026-10-28. **RH is OPEN.** Nothing here proves or assumes RH.
Labels: PROVED / CONJECTURE / GAP (plus NUMERICAL for project data, which is orientation only). Each step carries an R1–R5 line
(ARP §0.1: R1 uniform in L; R2 Euler product + FE jointly; R3 not via prime-weight signs; R4 soft mode / pole state; R5 where each
fake breaks, G_θ′ must pass). No numerics were run for this note.
Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
Parent: ARP_spectrum-measure.md S(b)3.6–3.8 (this file opens the S(b)3.7 GAP). Conventions are those of S.0 there.

## 0. Set-up and multiplicity conventions

- f ∈ V_L := C_c^∞([−L, L]), g = f ∗ f̃, f̂(z) = ∫ f e^{izu} du, ĝ(z) = f̂(z)·conj f̂(z̄).
- Zeros ρ = ½ + iγ_ρ, γ_ρ = (ρ − ½)/i. Write Γ_ℝ for the distinct real γ (on-line zeros) and m_ρ for multiplicity.
- Weil form (explicit formula; PROVED, Weil 1952; this normalization in Suzuki arXiv 2606.09096 eq. (3.1) and Alpöge–Furman
  arXiv 2608.13637 §1.2, both with the factor m):  Q_L(f) = Σ_ρ m_ρ ĝ(γ_ρ) = Σ_ρ m_ρ f̂(γ_ρ) conj f̂(conj γ_ρ).
- Off-line pairs. For a zero with Re ρ > ½ write γ_ρ = t − ib, b = Re ρ − ½ > 0. Its FE partner 1 − ρ̄ has γ = t + ib = conj γ_ρ.
  A *pair* is the set {ρ, 1 − ρ̄}; it contributes m_ρ·2 Re[f̂(γ_ρ) conj f̂(conj γ_ρ)].
- **Counts.** p := #{distinct ρ : Re ρ > ½} = number of distinct pairs = ½·#{distinct off-line zeros}.
  p_m := Σ_{Re ρ > ½} m_ρ = ½·#{off-line zeros counted with multiplicity}. p ≤ p_m, with equality iff every off-line zero is simple.
  Since ζ(σ) ≠ 0 on (0, 1) (classical), every off-line ρ is non-real, so the quadruple {ρ, ρ̄, 1−ρ, 1−ρ̄} has four distinct members
  and splits into the two distinct pairs of ρ and ρ̄; hence p is even (or infinite).
- κ_Q(L) := negative index of Q_L on V_L (max dim of a subspace on which Q_L < 0); κ_Q(∞) := sup_L κ_Q(L).
- **Relation to main's μ_ρ (G.4d, 7328e39).** Main's μ_ρ is the Poisson (harmonic) measure of ρ relative to the critical line.
  The off-line term used here is the exact point value m_ρ ĝ(γ_ρ), not ∫ĝ dμ_ρ. Exact split:
  Q_L = Σ_ρ m_ρ ∫ĝ dμ_ρ (≥ 0) + Σ_{off-line} m_ρ[ĝ(γ_ρ) − ∫ĝ dμ_ρ] (minus main's residual). μ_ρ = μ_{1−ρ̄}, so the residual pairs
  exactly like the pair form below. We use only the point-value form. The limit identity lim_Λ Δ_Λ = Σ m_ρ∫ĝ dμ_ρ is PROVED for
  function fields only (Connes 1999 §VIII Lemma 3), CONJECTURE for ℚ (G.4d), and **is not used here**.

## 1. What Bombieri states, and what the secondary sources say (citations verified 2026-10-28)

- **Bombieri 2000** (Rend. Lincei (9) 11, 183–233), abstract, verbatim: "if the Riemann Hypothesis is false but only with finitely
  many non-trivial zeros off the critical line we show that the number of negative eigenvalues is precisely one-half of the number
  of zeros failing to satisfy the Riemann Hypothesis, provided the truncation is big enough."
  Introduction (verbatim excerpts, bdim.eu text via search index): the functional's "Fourier transform … is easily interpreted as a
  quadratic form in infinitely many variables", "the eigenvalue 0 cannot occur for the finite dimensional approximations", and
  "the results obtained for the finite dimensional case carry over to the infinite dimensional case, provided we deal with negative
  eigenvalues bounded away from 0". The full PDF (bdim.eu) returned HTTP 500 again this session, so the theorem number, the exact
  truncation, and the multiplicity convention are **not verified from the source**.
- **Suzuki, arXiv 2606.09096 v3 (23 Sep 2026), "Weil's quadratic form via the screw function"** (read in full this session):
  - restates Yoshida 1992 as cited: [Yoshida, Lemma 2] positivity for small a; [Lemma 3] positivity on K_N(a) of finite codimension;
    [Theorem 2] RH ⇔ nondegeneracy on the completion of K(a); [Proposition 1] odd/even sectors;
  - restates Bombieri as cited: [Bombieri 2000, Problem 1] (H_0^1(−a,a)) and [Problem 2] (L²(𝓔)); [Theorem 3] minimizer exists;
    [Theorem 4] negative infimum attained in L²(−a,a); [Theorem 5] continuity of the even/odd infima (Suzuki: "details … not fully
    provided"); [Lemma 1] Bombieri's Lagrangian; [Lemma 3] the explicit formula. Bombieri 2003 (CPAM 56, 1151–1164) Thms 4.3, 4.4;
  - proves: Q_L = form of a self-adjoint A_L (Connes–Consani–Moscovici arXiv 2511.22755, Thm 3.6: lower bounded, discrete spectrum),
    A_L = Friedrichs extension of D*G_L D (Thm 1.1); C_c^∞(−L, L) is a form core (Cor. 1.2); λ_L = min spec A_L is continuous in L
    (Thm 1.3); for small L, λ_L > 0, simple, even eigenvector (Thm 1.4). The span of e^{iπnx/L} is a form core (Connes–Consani 2023,
    Prop. 2.3, as cited by Suzuki).
  - Suzuki does **not** restate the negative-eigenvalue count.
- **Alpöge–Furman, arXiv 2608.13637 (24 Aug 2026)** (accessible secondary restatement): "That the negative index of finite truncations
  of Weil's form equals the number of off-line zero pairs seen by the truncation was observed by Bombieri [Bom00] (Introduction)"
  (§7.1; also §1.3). Their Prop. 4.1 proves the inertia bound n_±(off-line part) ≤ #pairs via the pull-back of ⊕ m_ρ diag(1, −1),
  i.e. **per distinct pair, independent of m_ρ**. (They also note DH and Epstein as controls, §1.4.)
- **Yoshida 1992** (Adv. Stud. Pure Math. 21, 281–325): cited only through Suzuki 2206.03682 and 2606.09096 (numbers as above).

## 2. Statements

### G2.0 Window index = number of negative eigenvalues of A_L. PROVED (min–max; Suzuki Thm 1.1/Cor. 1.2, CCM Thm 3.6)
κ_Q(L) = #{negative eigenvalues of A_L, with multiplicity}. Reason: A_L is lower bounded with discrete spectrum; C_c^∞(−L, L) is a form
core; for a closed lower-bounded form, the max dimension of a negative subspace in a core equals the number of negative eigenvalues
(spectral subspace E(−∞,0) is finite-dim; approximate its basis in form norm from the core).
- **Truncation lemma (PROVED, Cauchy interlacing + the same approximation).** For any increasing sequence of finite-dim V_N ⊂ V_L whose
  union is a form core, ind_−(Q|V_N) is nondecreasing in N and equals κ_Q(L) for N ≥ N₀(L). So any "finite truncation" in this sense
  sees exactly κ_Q(L) once it is big enough. (Bombieri's "0 cannot occur" for finite truncations is not needed.)
R1 per window | R2 — | R3 — | R4 — | R5 fake-blind (form-theoretic).

### G2.1 Upper bound. PROVED (elementary; ARP S(b)3.6; Alpöge–Furman Prop. 4.1 for their compressions)
κ_Q(L) ≤ κ_Q(∞) ≤ p (distinct pairs). Each pair term m·2Re(A conj B) = (m/2)(|A+B|² − |A−B|²) has ≤ 1 negative square **whatever m is**;
on-line terms are ≥ 0. Hence a count "with multiplicity" (p_m) cannot be attained when some off-line zero is multiple.
R1 uniform in L | R2 FE used (pairing ρ ↔ 1−ρ̄) | R3 — | R4 — | R5 fake-blind; holds for every fake with an FE.

### G2.2 Exact window criterion. PROVED (elementary; written out here). Hypothesis H_fin: p < ∞.
Let z_1, …, z_{2p} be the distinct off-line points {γ_ρ, conj γ_ρ : Re ρ > ½}, ℓ(f) := (f̂(z_k))_k ∈ ℂ^{2p}, and
P_L(f) := Σ_{γ∈Γ_ℝ} m_γ |f̂(γ)|² ≥ 0 (absolutely convergent). Then Q_L(f) = P_L(f) + ℓ(f)* J ℓ(f) with J = ⊕_pairs m_ρ [[0,1],[1,0]]
(signature (p, p), negative eigenvalues −m_ρ). Define the **minimal real-zero energy**
  E_L(w) := inf{ P_L(f) : f ∈ V_L, ℓ(f) = w },  w ∈ ℂ^{2p}.
E_L is a positive semidefinite Hermitian form (squared quotient seminorm of P_L modulo ker ℓ; ℓ is onto because e^{iz_k u} are
linearly independent on [−L, L]). **Then κ_Q(L) = ind_−(J + E_L).**
Proof. (≤) If W ⊂ V_L is Q-negative, ℓ is injective on W (ℓf = 0 ⇒ Q(f) = P_L(f) ≥ 0), and for w = ℓf, (J+E_L)(w) ≤ P_L(f) + w*Jw =
Q(f) < 0; so ℓ(W) is a (J+E_L)-negative subspace of the same dimension. (≥) If U is (J+E_L)-negative of dim k, take a linear right
inverse R₀ of ℓ and, for a basis of U, elements of ker ℓ approximating the P_L-projection of R₀u onto (ker ℓ)^⊥; the linear map
R = R₀ − (those corrections) has ℓR = id_U and P_L(Rw) ≤ E_L(w) + ε|w|². Then Q(Rw) ≤ (J+E_L)(w) + ε|w|² < 0 on U for ε small
(compact unit sphere). ∎
- Monotonicity: E_L is nonincreasing in L (inf over a larger set), consistent with κ_Q nondecreasing.
- **Sufficient window condition:** λ_max(E_L) < min_{Re ρ>½} m_ρ ⇒ κ_Q(L) = p (Weyl: eigenvalues of J + E_L ≤ those of J + ‖E_L‖).
R1 per window; L* below is not effective | R2 FE (pairing) | R3 — (no prime weights used) | R4 E_L is the real-zero energy; at the
onset the negative direction is the P_L-softest admissible direction (cf. §3 DH) | R5 fake-blind under H_fin.

### G2.3 E_L ↓ 0. PROVED (written out here). Hypothesis H_fin.
E_∞(w) := lim_L E_L(w) = inf over all f ∈ C_c^∞(ℝ). Claim: E_∞ ≡ 0.
Proof. If E_∞(w) > 0 for some w, polarization λ(u) := E_∞(u, w) gives |λ(ℓf)|² ≤ E_∞(w)·E_∞(ℓf) ≤ E_∞(w)·P(f), so the nonzero
combination Λ := Σ_k ā_k ev_{z_k} = λ∘ℓ is P-bounded on C_c^∞(ℝ). Fix h ∈ C_c^∞ with ĥ(z_k) ≠ 0 for all k (generic; e.g. dilate) and
test on translates h_t := h(· − t): P(h_t) = P(h) for all t (|ĥ_t| = |ĥ| on ℝ), while Λ(h_t) = Σ_k ā_k ĥ(z_k) e^{iz_k t} is an
exponential polynomial with distinct non-real exponents. Bounded on ℝ forces every coefficient to vanish (take the largest
|Im z_k| group as t → ∓∞; after dividing by its growth, the remaining almost-periodic sum → 0, so its Bohr mean Σ|c_k|² = 0; induct).
So ā_k ĥ(z_k) = 0, i.e. Λ = 0, contradiction. ∎
Since E_L is a decreasing family of PSD forms on ℂ^{2p} with pointwise limit 0, ‖E_L‖ → 0.
R1 limit only, non-effective | R2 not used | R3 — | R4 — | R5 fake-blind under H_fin.

### G2.4 Main statement (Bombieri's count, window form). PROVED under H_fin (G2.1–G2.3); priority Bombieri 2000.
If ζ has only finitely many off-line zeros, then with p = #{distinct ρ : Re ρ > ½} (= #distinct pairs {ρ, 1−ρ̄}):
  (a) κ_Q(L) ≤ p for every L, and κ_Q(L) = ind_−(J + E_L) exactly;
  (b) κ_Q(∞) = ind_−(J) = p;
  (c) there is L*_p < ∞ with κ_Q(L) = p for all L ≥ L*_p; it suffices that λ_max(E_L) < min m_ρ;
  (d) multiplicity: κ_Q(∞) = ½·#(off-line zeros counted **without** multiplicity). It equals ½·#(with multiplicity) iff all off-line
      zeros are simple. Bombieri's "one-half of the number of zeros failing RH" is therefore correct as stated iff he counts distinct
      zeros or the off-line zeros are simple (his convention: GAP (citation), full text inaccessible).
  (e) Parity split (PROVED, same proofs with f_t = h(·−t) ± h(−·−t), which keep P bounded and give exponents ±z_k, all distinct and
      non-real): on even (resp. odd) f, ev_z = ±ev_{−z}, each quadruple {ρ, ρ̄, 1−ρ, 1−ρ̄} becomes one pair of weight 2m, so
      κ_even(∞) = κ_odd(∞) = p/2 = #quadruples. (Matches Yoshida's even/odd criterion, Prop. 1 via Suzuki; real zeros are excluded for ζ.)
- **Matching to Bombieri's truncation.** "Truncation big enough" = (support parameter t = L ≥ L*_p) and (finite-dim compression large
  enough, G2.0). With G2.0, Bombieri's statement for any truncation scheme that exhausts a form core of the window space is equivalent
  to (b)+(c). Whether his truncations are of that kind (likely the exponential/Fourier basis on [−t, t], which is a form core by
  Connes–Consani 2023 Prop. 2.3 via Suzuki) is **GAP (citation)** only; (a)–(e) do not depend on it.
R1 no effective L*_p; finite-window data cannot bound p | R2 FE only (pairing); no Euler product ⇒ no RH content | R3 — | R4 onset is
governed by E_L (soft-mode energy) | R5 see §3: any fake satisfying H_fin obeys the same count; fakes with p = ∞ are outside H_fin.

### G2.5 Transition points. PROVED as cited (Suzuki 2606.09096 Thm 1.3; Suzuki 2206.03682 Thm 6.1; Yoshida Thm 2 via Suzuki)
λ_L is continuous; if RH fails, Q is degenerate at the first L where λ_L crosses 0. In G2.2 language: L₁ := inf{L : λ_min(J + E_L) < 0}.
Monotone E_L and continuity give that each jump of κ_Q happens where J + E_L acquires a zero eigenvalue (CONJECTURE: this is the
same degeneracy point as Suzuki Thm 6.1; not checked).

### G2.6 Infinitely many off-line zeros (outside H_fin). Partly CONJECTURE.
- PROVED: κ_Q(L) < ∞ for each L for ζ (Yoshida Lemma 3 via Suzuki; = #neg. eigenvalues of A_L, G2.0).
- Note: G2.3 does not give a lower bound when p = ∞, because the remaining off-line terms are not ≥ 0 and cannot be absorbed into P.
- CONJECTURE: if p = ∞ then κ_Q(∞) = ∞ (expected; would follow from a version of G2.3 with an infinite J).

## 3. Fakes (R5)

| Fake | Off-line pairs | G2.4 applies? | Window index |
|---|---|---|---|
| DH (Davenport–Heilbronn, mod 5) | infinitely many (Davenport–Heilbronn 1936); project zero 0.808517 + 85.699348i (Queen's dh_sinh.py) | No (p = ∞) | finite per window (expected; GAP for DH: Yoshida-type finiteness not checked) |
| Epstein x² + 5y² | infinitely many (Davenport–Heilbronn 1936, J. LMS 11) | No | as DH |
| F_θ (project fake) | all zeros off line, infinitely many pairs | No | — |
| H₋ₜ (backward heat, t > 0) | non-real zeros exist (Rodgers–Tao 2020, Λ ≥ 0); finiteness unknown to us | if finite: yes | — |
| DMV Beurling / P_X with planted pairs | finitely many (planted) | Yes, given an explicit formula of the same shape (not checked) | = κ_ζ(L) for L below the agreement window (R1) |
| G_θ′ (shift ±iθ′) | 2p if ζ has p (both shifts) | Yes iff ζ's p finite | κ_G ≤ 2κ_ζ (pull-back of two ζ forms); **passes** |

**DH and its window onset L* = 1.7104.** Project data (not recomputed here):
- PROVED (computer-assisted: Arb-certified negative Rayleigh quotients at L = 1.7105, K300/K360; S75, GROK_RESEARCH_NOTES:992;
  not replayed here): κ_DH(L) ≥ 1 for L ≥ 1.7105 (κ is nondecreasing in L).
- NUMERICAL: DH positive at L = 1.71; onset L*_DH = 1.7104 ± 0.0001 (S75).
  The lower side is **not** certified, so "κ_DH(L) = 0 for L < 1.7104" is NUMERICAL.
- NUMERICAL (principal-series.md lines 123–161): the onset eigenvector (even sector) is **nodeless** (~80% mode 0); it is not the pair
  vector e^{−iγ₀x} sinh(β₀x) (2Lγ₀/π ≈ 93 nodes). S84: the off-line zero 0.8085 + 85.699i supplies the negative zero-side term
  (−1.234e−28 vs +1.19e−28 on-line).
- Reading through G2.2 (CONJECTURE, interpretation only — DH has p = ∞, so G2.2 does not literally apply): the onset is where the
  pair term outweighs the *minimal real-zero energy* E_L, and that minimum is attained on the softest direction (nodeless bump),
  not on a vector oscillating at γ₀. This matches the project's finding that DH "fails along a ζ-stiff direction" only after the
  soft mode is spent (S75/S84): L* is set by the soft mode, not by γ₀.
- Consequence for ζ (R1): the same mechanism means a hypothetical off-line ζ pair at large height could have its onset L*_p beyond
  any certified window. A certified window positivity (e.g. the 2.84 replay, 4cf1bfd) gives κ_ζ(L) = 0 there and nothing about p.
- What separates ζ from DH is not the count (G2.4 is fake-blind) but whether p = 0, i.e. RH itself.

## 4. Discrepancy flagged
- KREIN_LANGER_S97.md §A.1 says κ(∞) = #{ρ : Re ρ > ½} "(multiplicity)". By G2.1 the Weil form's negative index counts **distinct**
  zeros (a pole of i(ξ′/ξ)(½ − iz) in ℂ₊ is simple whatever m is, so the generalized-Nevanlinna count is also per distinct zero).
  Not edited here (one file per commit); suggested fix: replace "(multiplicity)" with "(distinct zeros)".

## 5. Remaining gap
1. **GAP (citation):** Bombieri 2000's exact truncation, theorem number, and multiplicity convention (full text inaccessible; bdim.eu 500).
   The mathematics of the window statement no longer depends on it (G2.4 PROVED under H_fin).
2. **GAP (effectivity, R1):** no bound on L*_p or on the first onset L₁ in terms of the zeros; E_L → 0 is non-effective. Bounding E_L
   from above needs quantitative interpolation at off-line points with small real-zero energy, which depends on unknown zero locations.
3. **CONJECTURE:** p = ∞ ⇒ κ_Q(∞) = ∞ (G2.6), relevant to DH, Epstein, F_θ.
4. No RH content: G2.4 is a counting identity, true for every L-function-like object with H_fin. RH ⇔ p = 0 ⇔ κ_Q ≡ 0 (Weil) is unchanged
   and OPEN.

## References (checked this session unless marked)
- E. Bombieri, Remarks on Weil's quadratic functional in the theory of prime numbers, I, Rend. Lincei (9) Mat. Appl. 11 (2000) 183–233
  (abstract and introduction only; full text not accessible). E. Bombieri, CPAM 56 (2003) 1151–1164 (as cited by Suzuki).
- M. Suzuki, Weil's quadratic form via the screw function, arXiv:2606.09096 v3 (2026). M. Suzuki, arXiv:2206.03682.
- L. Alpöge, R. Furman, arXiv:2608.13637 (2026), §§1.3, 1.4, Prop. 4.1, §7.1.
- H. Yoshida, On Hermitian forms attached to zeta functions, Adv. Stud. Pure Math. 21 (1992) (via Suzuki).
- A. Connes, C. Consani, H. Moscovici, arXiv:2511.22755, Thm 3.6 (via Suzuki). A. Connes, C. Consani (2023), Prop. 2.3 (via Suzuki).
- H. Davenport, H. Heilbronn, J. London Math. Soc. 11 (1936) (not re-checked). B. Rodgers, T. Tao, Forum Math. Pi 8 (2020) (not re-checked).
