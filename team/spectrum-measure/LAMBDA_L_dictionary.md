# Λ ↔ L dictionary: phase-space (Sonin/prolate) cutoff Λ vs support window L

Seat: Spectrum & Measure (team/spectrum-measure). Analytic only. **RH is OPEN**; nothing here proves or assumes it.
Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
Date: 2026-10-28. Companion files: GAP2_bombieri_window.md (GAP 2; cited as **PROVED under H_fin**), PROLATE_BOUNDS.md (prolate inputs).

Labels: **PROVED** (with a verified citation, or elementary and proved here), **PROVED (computer-assisted)**, **CONJECTURE**, **GAP**.
"verified" = statement quoted from the arXiv PDF read this session; "via restatement" = only a secondary source checked.
Track tags (lead's assignment):
- **[Λ-track]** (ours): Λ > 1 at fixed 2L < log 2. RH-silent. Feeds Principal series' projector dictionary and PROLATE_BOUNDS.md.
- **[L-track]** (Celestial's): the support extension L ↑ past ½ log 2 at fixed Λ. Described here only for the 2D picture.
- **[both]** / **[neither]** where applicable.

Ownership (no duplication): Celestial owns §C.13 (operator/trace side, the CC proof chain, 2 log Λ). Principal series owns the projector
dictionary GAP (Connes S_Λ ⊖ E(B_{Λ,0}) vs CC Sonin 1 − P₁∨P̂₁, ARP §P.8 at 9d520f8). Queen owns the leak (1 − S_Λ)Q′_Λ (M1) and
Conj Q.8.E (ARP_queen_Q8.md, merged to main in 1999c2a). This file only fixes **which parameter is which** and what follows elementarily.

**Withdrawn:** the earlier guess "L = log Λ" (spectrum-measure, previous turn) is **withdrawn**. It conflated two different cutoffs that
share a letter: the *support* cutoff λ of CC 2023 / CCM ([λ⁻¹, λ] ↔ [−L, L], so λ = e^L) and the *phase-space* cutoff Λ of CC 2021 /
Connes 1999. Below, (Λ, L) are **independent parameters**.

---

## 0. Conventions and verified source statements

- Ours: f̂(z) = ∫ f(u) e^{izu} du; g = f ∗ f̃, supp f ⊂ [−L, L] (L = half-width, binding convention), Q_L(f) = Σ_ρ m_ρ ĝ(γ_ρ),
  γ_ρ = (ρ − ½)/i. Multiplicative picture u = log x on ℝ*₊. ĝ(z) = f̂(z)·\overline{f̂(z̄)}.
- **CC 2021** (Connes–Consani, "Weil positivity and Trace formula the archimedean place", Selecta Math. 27:77, arXiv 2006.13771) — verified:
  - Abstract: "…compressed onto the orthogonal complement of the range of the cutoff projections associated to the cutoff in phase
    space, for Λ = 1."
  - Sonin space: "functions, which, together with their Fourier transform, vanish identically in the interval [−1, 1]".
  - Def 4.4: "S(α, β) := {ξ ∈ L²(ℝ)ev | ξ(q) = 0, ∀q, |q| ≤ α, (F_{e_ℝ}ξ)(p) = 0, ∀p, |p| ≤ β}", F_{e_ℝ}ξ(y) = ∫ξ(x)e^{−2πixy}dx (24).
  - (61) "(ϑ(λ)ξ)(v) := λ^{−1/2} ξ(λ^{−1} v)"; (62) P = 1 − P₁, P̂ = 1 − P̂₁.
  - Support restriction, verbatim: "support of the test function is contained in the interval (1/2, 2) ⊂ ℝ*₊ so that rational primes
    are not involved"; QW(g) := W∞(g ∗ g*) on "the vector space V of smooth functions with support in the interval [2^{−1/2}, 2^{1/2}]
    and whose Fourier transform vanishes at i/2".
  - Weil's inequality as restated by CC (Weil 1972 — **via restatement**): "For any smooth, positive definite function f with support
    in the interval (1/2, 2) and whose Fourier transform vanishes at ±i/2 one has: W∞(f) ≥ 0 where W∞ := −W_ℝ."
  - Theorem 1, verbatim: "Let g ∈ Cc∞(ℝ*₊) have support in the interval [2^{−1/2}, 2^{1/2}] and Fourier transform vanishing at i/2
    and 0. Then one has W∞(g ∗ g*) ≥ Tr(ϑ(g) S ϑ(g)*)." (4)
  - Theorem 6.11, verbatim: "Let g ∈ Cc∞(ℝ*₊) be a smooth function with support in the interval [2^{−1/2}, 2^{1/2}] and whose Fourier
    transform vanishes at −i/2. Let S be the orthogonal projection of L²(ℝ)ev onto the subspace of even functions which vanish as well
    as their Fourier transform in the interval [−1, 1]. Then W∞(g ∗ g*) ≥ Tr(ϑ(g) S ϑ(g)*) − c|ĝ(0)|², c = 4γ/log 2." (141)
  - **So CC 2021 has TWO independent restrictions (Golden Eye; both verified above):** (R-Λ) phase-space Sonin cutoff Λ = 1, i.e.
    S = S(1, 1); (R-L) test-function support [2^{−1/2}, 2^{1/2}], i.e. our L ≤ ½ log 2 (2L ≤ log 2). **The prime 2 is excluded by
    the support restriction (R-L), not by Λ.**
- **Connes 1999** (Selecta Math. 5, math/9811068) — verified: (21) "S_Λ = {ξ ∈ L²(C_k); ξ(g) = 0 ∀g, |g| ∉ [Λ⁻¹, Λ]}" (a *position*
  cutoff on the idèle class group); (31) "P_Λ = {ξ ∈ L²(ℝ); ξ(x) = 0 ∀|x| > Λ}", P̂_Λ = F P_Λ F⁻¹ with α(x) = e^{−2πix}; B_Λ = span of
  prolates ψ_n, n ≤ 4Λ²; (23) Q′_{Λ,0} ≤ S_Λ (function fields); (24) Δ_Λ(f) = Trace((S_Λ − Q′_{Λ,0})V(f)).
- **CC 2023** ("Spectral triples and ζ-cycles", Enseign. Math. 69 (2023) 93–148, arXiv 2106.01715) and **CCM** (arXiv 2511.22755) —
  verified: the form QW_λ lives on L²([λ⁻¹, λ], d*u) with primes "1 < n ≤ λ²" (CC 2023 (2.11), Prop 2.1). This λ is a **support**
  cutoff: λ = e^L.

## (a) Relation between Connes–Consani's Λ and our L

**(a1) PROVED (definitions above).** Four letters, two parameters:

| symbol | source | what it cuts | in our variables |
|---|---|---|---|
| λ (CC 2023, CCM; Suzuki's a = log λ) | support of the test function | [λ⁻¹, λ] | λ = e^L |
| support [2^{−1/2}, 2^{1/2}] (CC 2021) | support of the test function | fixed | L = ½ log 2 |
| Λ (CC 2021 Sonin S(Λ,Λ); Connes 1999 P_Λ, P̂_Λ, B_Λ) | phase space: vector *and* its Fourier transform | — | independent of L |
| Λ (Connes 1999 S_Λ, (21)) | position on C_ℚ (|g| ∈ [Λ⁻¹, Λ]) | — | independent of L |

Q_L involves only L. Λ enters only through the comparison operators (Sonin trace, Δ_Λ, Δ^♮_Λ).
**Caution (PROVED, elementary):** "Λ = 1" means different things. For ℚ, Connes' S_1 = {ξ supported on |g| = 1} = 0 (the set has Haar
measure 0), so S_1 is vacuous; CC's S(1, 1) is infinite-dimensional (it contains the Sonin space). The two "Λ = 1" are not the same object.

**(a2) Rescaling. PROVED (elementary; confirms Principal series §P.8 (1)).** With ϑ as in CC (61):
ϑ(λ) P_Λ ϑ(λ)⁻¹ = P_{λΛ}, ϑ(λ) P̂_Λ ϑ(λ)⁻¹ = P̂_{Λ/λ}. Proof: ϑ(λ)ξ is supported in |v| ≤ λΛ iff ξ in |x| ≤ Λ; and
F ϑ(λ) = ϑ(λ⁻¹) F (substitution in (24)). With λ = Λ⁻¹: (P_Λ, P̂_Λ) ≅ (P₁, P̂_{Λ²}). Likewise ϑ(λ) S(α, β) = S(λα, β/λ), so
S(Λ, Λ) = ϑ(Λ) S(1, Λ²): **only the product αβ (here Λ²) is a unitary invariant**; the prolate parameter is c = 2πΛ² (PROLATE_BOUNDS.md §1).
[Λ-track] | R1 n/a (no L) | R2 RH-silent (archimedean only) | R3 — | R4 — | R5 every fake has the same archimedean projections.

**(a3) Monotonicity in Λ. PROVED (elementary).** For Λ′ ≥ Λ, S(Λ′, Λ′) ⊂ S(Λ, Λ) (vanishing on a larger set). Hence for every g,
T_Λ(g) := Tr(ϑ(g) S(Λ, Λ) ϑ(g)*) is nonincreasing in Λ, and CC Theorem 1 / 6.11 remain true with S(Λ, Λ), Λ ≥ 1, in place of S
(weaker). Moreover T_Λ(g) → 0 as Λ → ∞ whenever T_1(g) < ∞ (it is, by Thm 1: T_1 ≤ W∞(g∗g*)): X_Λ := ϑ(g)S(Λ,Λ)ϑ(g)* is a decreasing
family of positive operators dominated by the trace-class X_1, and X_Λ → 0 strongly because ∩_Λ S(Λ, Λ) = {0}; dominated convergence
in Σ_n ⟨X_Λ e_n, e_n⟩. [Λ-track] RH-silent.

## (b) The two-sided Sonin space vs our one-sided window — precise relation

Celestial holds the dictionary at CONJECTURE until this is precise. Statement:

**(b1) They are different kinds of object; there is no inclusion. PROVED (definitions).** The window [−L, L] constrains the *test
function* f (a group-algebra element of ℝ*₊, acting by ϑ(g)). The Sonin conditions constrain *vectors* ξ ∈ L²(ℝ)ev. In the variable
u = log|x| the Sonin conditions are one-sided: ξ supported in u ≥ log α and F ξ supported in u ≥ log β ("two-sided" = both ξ and Fξ).
The window is a two-sided bound |u| ≤ L on f, and |u| ≤ 2L on g = f ∗ f̃. A statement "window ⊂ Sonin" or the reverse is not well-typed.

**(b2) Mapping relation. PROVED (elementary).** If supp g ⊂ [e^{−L}, e^{L}] then ϑ(g) := ∫ g(λ) ϑ(λ) d*λ maps
S(α, β) → S(α e^{−L}, β e^{−L}).
Proof: by (a2), ϑ(λ) S(α, β) = S(λα, β/λ) ⊂ S(α e^{−L}, β e^{−L}) for every λ ∈ [e^{−L}, e^{L}]; integrate (closed subspace).
Consequences:
- ϑ(g) S(Λ, Λ) ⊂ S(1, 1) whenever Λ ≥ e^{L}; the threshold curve is **Λ = e^{L}** (with ϑ(g*) the same).
- On CC's own support (L = ½ log 2) the threshold is Λ = √2.
- Queen's Conj Q.8.E also uses Λ = e^{L}, but for Connes' *position* cutoff S_Λ (21), not for S(Λ, Λ). Whether the two choices agree
  through the projector dictionary is part of Principal series' GAP (§P.8 (3)); not claimed here.
[Λ-track] | RH-silent | R5 fake-blind (archimedean).

**(b3) Pole/projector interaction. PROVED, as settled by Principal series §P.8 (2) (their result, cited).** Mellin of Ẽf is
f̂(0)/(s − 1) − f(0)/s; f(0) belongs to the P_Λ side, f̂(0) to the P̂_Λ side; the removed pair is G.3a's OS4 span{e^{±u/2}}. Consistent
with Queen's Lemma Q.8.2 (E(f)(g) = E(f̂)(g⁻¹) + f̂(0)|g|^{−1/2} − f(0)|g|^{1/2}). How Λ ↔ L interacts with it: the two pole
evaluations f(0), f̂(0) do not involve Λ, and the window L does not change which side each sits on; under (a2) the P̂-side is P̂_{Λ²}
after normalizing P_Λ to P₁. (Bookkeeping only; nothing beyond P.8 (2) is claimed.)
**GAP (Principal series' own, minor):** the projector dictionary Connes S_Λ ⊖ E(B_{Λ,0}) vs CC Sonin 1 − P₁ ∨ P̂₁. We add only (a2),
(b2): the CC side at general Λ is 1 − P_Λ ∨ P̂_Λ ≅ 1 − P₁ ∨ P̂_{Λ²}.

## (c) The lead's question, two-parameter form

**Question.** At Λ = 1 and 2L ≤ log 2 (CC's region), is κ_Q(L) = 0 there because no primes are seen, so that Q_L = arch + pole?

**(c1) No primes are seen. PROVED** (main ARP G.2b, PROVED-LEAN `primeSide_eq_zero_of_support_lt_log_two` for L < ½ log 2; at
L = ½ log 2 the convolution g is smooth and supported in [−log 2, log 2], so g(±log 2) = 0 as well). Hence
Q_L(f) = [pole pair] + [archimedean term W∞(g)] in CC's normalization (W∞ := −W_ℝ).
Normalization match (CC g ↔ our f via Δ^{1/2}, u = log x, t → −t): **PROVED (elementary)** at the level used here (pole points ±i/2,
the point 0, support), not re-derived line-by-line for every constant.
[L-track boundary, Λ-independent] | R1 bounded window only (P_X/DMV agree with ζ here) | R2 **absent**: no Euler product, so RH-silent.

**(c2) The pole is killed by one condition. PROVED (elementary).** ĝ(±i/2) = f̂(±i/2)·\overline{f̂(∓i/2)}, so the pole pair is
2 Re(f̂(i/2) \overline{f̂(−i/2)}): one negative square at most, and it vanishes on the codimension-1 class f̂(−i/2) = 0.

**(c3) What the literature gives on CC's region.**
- κ_Q(½ log 2) ≤ 2: **PROVED** from CC Theorem 1 (verified): on the codimension-2 class {f̂(−i/2) = 0, f̂(0) = 0} (CC's "vanishing at i/2 and 0" after t → −t),
  Q_L(f) = W∞(g ∗ g*) ≥ Tr(ϑ(g) S ϑ(g)*) ≥ 0. (Equivalently from Thm 6.11: on the codimension-1 class, Q_L ≥ −c|ĝ(0)|², one more square at most.)
- κ_Q(½ log 2) ≤ 1: from Weil's inequality — **via restatement** (Weil 1972 as quoted by CC).
- κ_Q(L) = 0 for all L ≤ ½ log 2: **PROVED (computer-assisted)** by the audited certificate at L = 71/50 = 1.42 (AUDIT_284.md,
  Q ≥ 2^{−4200}‖f‖²) plus monotonicity of κ_Q in L (windows nest). Small-L theorems: Suzuki 2606.09096 Thm 1.4 (verified, "for
  sufficiently small a > 0, the lowest eigenvalue λa is positive"), Yoshida Lemma 2 (via restatement) — both non-explicit in a.

**(c4) Answer, with the structural reason named honestly.**
- "No primes are seen" (c1) is why κ_Q on this region is a statement about **the archimedean place plus the pole alone**. That is a
  *reduction*, not the reason for the sign. PROVED.
- The analytic *source* of positivity that the literature offers is CC's Sonin-trace inequality W∞ ≥ Tr(ϑ(g)Sϑ(g)*) ≥ 0, valid off
  two directions (pole, ĝ(0)). PROVED (verified). It does **not** by itself give κ = 0 on the full window space; the remaining
  ≤ 2 directions (≤ 1 via Weil) are closed by the certificate, not by a structural argument. So: **κ = 0 there is PROVED, but the
  structural mechanism covers codimension 2 (CC) / 1 (Weil, via restatement); the full κ = 0 is computer-assisted.**
- In GAP 2 terms (G2.2, **PROVED under H_fin**): κ_Q(L) = ind_−(J + E_L), so κ = 0 ⇔ J + E_L ≥ 0. A restatement, not a mechanism.
- **Not** orthogonality: the evaluation map ℓ is surjective (G2.2, PROVED under H_fin), so off-line pairs are not invisible on small
  windows; their energy E_L is just too large (G2.3: E_L ↓ 0 only as L → ∞).
- **Not** vacuity: CC's Sonin space S(1,1) is infinite-dimensional. (Vacuity holds only for Connes 1999's S_1 for ℚ, (a1).)
- **Flag for the lead:** main ARP G.2b states "In that range Q_L(f) = pole + W_∞ ≥ 0" citing Yoshida 1992 and CC 2021. CC 2021 covers
  the pole-free class (Thm 1: codim 2; Weil via CC: codim 1). The full-space statement needs Yoshida Lemma 2 (via restatement, only
  "sufficiently small a"), Suzuki Thm 1.4 (same), or the certificate (AUDIT_284). Suggest G.2b cite AUDIT_284 / Mori for the full
  space at L = ½ log 2. Not edited (not our file).
- RH-silence: any object with ζ's archimedean factor and pole agrees with ζ on this region (P_X, DMV, F_θ, G_θ′ coincide with ζ for
  2L < log 2); DMV has off-line zeros yet κ = 0 here — consistent with G2.4 (κ = p only for L ≥ L*_p). R1 illustration.
R-line (c): R1 fails by design (bounded window) | R2 absent (no prime) | R3 — | R4 soft mode/pole handled by (c2) | R5 fakes below.

## (d) Where the first negative direction can enter (2D picture)

Coordinates: horizontal L (support), vertical Λ (phase-space comparison cutoff).

**(d1) Q_L and κ_Q(L) do not depend on Λ. PROVED (trivial: Λ does not occur in Q_L).** Λ changes only the comparison term T_Λ,
which is nonincreasing and → 0 (a3). **Hence no negative direction can enter by raising Λ at fixed L.** [Λ-track]
- Our track Λ > 1, 2L < log 2: κ = 0 throughout (c3); the content of the track is the *size* of the Sonin/prolate comparison
  (T_Λ, the defect Q_L − T_Λ, the projector dictionary), i.e. archimedean calibration. RH-silent. Inputs: PROLATE_BOUNDS.md.
- Hook for Queen: on this track A_Λ = S_Λ − S_Λ Q′_Λ S_Λ and Δ^♮_Λ = Tr(A_Λ V V*) ≥ 0 (Lemma Q.8.4, hers). For ℚ, Q′ ≤ S_Λ is not exact
  with the prolate cutoff (§Q.8); the leak (1 − S_Λ)Q′_Λ is Queen's M1 (not estimated here). Conj Q.8.E is the operator form of the GAP.
  PROLATE_BOUNDS.md supplies the prolate spectral inputs (plunge mass, tails, counts) in her normalization.

**(d2) Along L. PROVED under H_fin (GAP 2).** If RH fails with p distinct off-line pairs, κ_Q(L) jumps from 0 to p, first at
L₁ := inf{L : λ_min(J + E_L) < 0} (G2.2–G2.5), and κ = p for L ≥ L*_p (G2.4). **PROVED (computer-assisted):** L₁ > 1.42 (AUDIT_284).
No effective upper bound on L₁ (GAP, R1). [L-track — Celestial]
- Primes enter at 2L = log n: n = 2 at L = ½ log 2 ≈ 0.347 (inside the certified window 1.42); then 3 (≈ 0.549), 4 (≈ 0.693), 5 (≈ 0.805),
  …, up to e^{2.84} ≈ 17.1 (n ≤ 17) inside the certified window. So the certificate already covers the first ~7 primes (2…17) entering.
- **The real R2 test is the L-direction past ½ log 2 at Λ = 1** (where the Euler product first enters), not Λ > 1 at fixed small L
  (which never sees a prime). PROVED as a statement about which parameter carries arithmetic (d1, c1).
- Prime's heuristic (the first negative direction appears near e^{2L} ≈ X, not Λ²): consistent with (d1) — Λ cannot matter at all.
  The quantitative "≈ X" for P_X-type systems is **CONJECTURE/heuristic** (Prime's); the exact agreement P_X = ζ for 2L ≤ log p⁺(X) is
  main R1 (Prime's DMV theorem). Not developed here.

**(d3) Interaction (both parameters). CONJECTURE (not ours to settle).** For the CC-type identity at general (Λ, L) one expects
Tr(ϑ(f) S(Λ,Λ)) = W∞(f) + (correction built from P_Λ, P̂_Λ) for supp f in the window; by (a2) the correction at Λ is CC's at Λ = 1
with P̂_{Λ²}. Proving this belongs to Celestial §C.13 (trace side) and Principal series (projectors); the threshold Λ = e^L of (b2) is
where ϑ(g) stops pushing S(Λ,Λ) out of S(1,1).

## Fakes (R5) — where each breaks in this file

| fake | (a)/(b) [Λ-track] | (c) region 2L ≤ log 2 | (d2) L-track |
|---|---|---|---|
| P_X / DMV | same archimedean data: passes (blind) | identical to ζ (κ = 0) | differs from ζ from 2L > log p⁺(X); DMV fails RH (G.4e) |
| F_θ (positive weights) | blind | identical to ζ (no prime) | breaks at R3 once primes enter |
| G_θ′ (signed weights) | blind | identical | passes R3 (must pass) |
| DH (mod 5) | different gamma factor/conductor, no pole: different W∞ | not covered by CC; not checked (GAP) | off-line zeros, p = ∞ (outside H_fin) |
| Epstein | different archimedean factor | not covered by CC; not checked (GAP) | off-line zeros; no Euler product |
| H₋ₜ (de Bruijn–Newman side) | not a local explicit-formula deformation; not classified here | GAP | GAP |

## Status summary

| item | label | track |
|---|---|---|
| L = log Λ guess | **withdrawn** | — |
| (a1) parameter table; CC 2021 two restrictions (verbatim) | PROVED (verified) | both |
| (a2) rescaling (P_Λ, P̂_Λ) ≅ (P₁, P̂_{Λ²}); S(Λ,Λ) ≅ S(1,Λ²) | PROVED (elementary) | Λ |
| (a3) T_Λ nonincreasing, → 0 | PROVED (elementary + CC Thm 1) | Λ |
| (b2) ϑ(g)S(α,β) ⊂ S(αe^{−L}, βe^{−L}); threshold Λ = e^L | PROVED (elementary) | Λ |
| (b3) pole/projector | PROVED (Principal series P.8(2)); projector dictionary GAP (theirs) | Λ |
| (c1) Q_L = pole + arch for 2L ≤ log 2 | PROVED (Lean, main G.2b) | boundary |
| (c3) κ ≤ 2 (CC Thm 1); κ ≤ 1 (Weil via CC); κ = 0 | PROVED (verified) / via restatement / PROVED (computer-assisted) | — |
| (d1) κ independent of Λ | PROVED | Λ |
| (d2) onset along L at L₁ > 1.42; count p | PROVED under H_fin + computer-assisted; effectivity GAP | L |
| (d3) general (Λ, L) trace identity | CONJECTURE | both |
| RH | **OPEN** | — |

## References
- A. Connes, C. Consani, Weil positivity and Trace formula the archimedean place, Selecta Math. 27 (2021) 77, arXiv:2006.13771 — verified (Abstract, Intro, Thm 1, Def 4.4, (24), (61)–(62), Thm 6.11).
- A. Connes, Trace formula in noncommutative geometry and the zeros of the Riemann zeta function, Selecta Math. 5 (1999), math/9811068 — verified ((21)–(24), (31)–(32)).
- A. Connes, C. Consani, Spectral triples and ζ-cycles, Enseign. Math. 69 (2023) 93–148, arXiv:2106.01715 — verified (Prop 2.1, (2.11)).
- A. Connes, C. Consani, H. Moscovici, Zeta Spectral Triples, arXiv:2511.22755 — verified (Thm 3.6, Cor 3.7).
- M. Suzuki, arXiv:2606.09096 v3 — verified (Thm 1.4). H. Yoshida (1992) Lemma 2 — via restatement. A. Weil (1972) — via restatement (CC).
- Team: GAP2_bombieri_window.md (G2.0–G2.6, PROVED under H_fin); AUDIT_284.md; main derivation/ARP.md G.2b, R1–R5; team/principal-series ARP §P.8 (9d520f8); team/queen/ARP_queen_Q8.md §Q.8.
