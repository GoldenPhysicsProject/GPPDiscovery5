# AUDIT14: Golden Eye's audit of Rounds 14B–14F and 14E2

Auditor: Golden Eye (Golden Physics Project, RH team). Requested by Heart for Daniel Toupin. Date: 2026-11-18.
Scope: labels (PROVED / CHECKED / NUMERICAL / CONJECTURE / OPEN), overclaims, tautological controls, literature.
Standard used: **CHECKED** = code or identity check. **NUMERICAL** = evidence. A control only counts if it could have come out the other way.
I did not edit anyone else's files. My re-runs are tiny probes; their code is in the Appendix, so anyone can reproduce them from this directory.

**Verdict, one line each.**
- **14B:** the math holds as a restatement of Weil's criterion. The errata were applied accurately (§1).
- **14C:** the main no-go is correct but mostly by construction. Two items marked CHECKED are tautologies. "Rules out every semilocal route" is an overclaim (§2).
- **14D:** Faces 1+2 are correct; that is Tate's thesis. "E lands in the radical" is a correct theorem, but the numerical check doesn't test radical membership. "Positivity must live on the quotient ⇔ RH" is Weil's criterion with nothing added, and is partly misattributed to Connes 1999 (§3).
- **14F:** the algebra is correct and my falsifier was resolved correctly (unitarity is the separating axiom). Item 3 (additivity over products) is a tautology. Item 5 (ghosts are not spectrum of any Hecke-equivariant quotient) is vacuous as stated, though its algebraic restatement holds. Two ghost-column labels need fixing (§5).
- **14E:** the conclusion holds (no gap; a uniform gap is falsified at L = 8 by a proved Lehmer bound of 0.063, which I reproduced). But the Sylvester check is a one-sided restriction with low power, since a random subspace also keeps every ghost. The "floor" is a precision limit. The zeros in §5 are proof inputs, and `lehmer.py` doesn't run as committed. The Gaussian λ_min fit fails just outside its range: 80-digit data show double-exponential decay exp(−B e^L), B ≈ 13 (§6).
- **14E2:** reproduces (λ₁ and Q to 8 digits); no finite-L gap holds (Q/λ₁ 1.09–3.57). The Schur complement is the right quotient. But "S nearly orthogonal to u₁" is false: Q = λ₁/(1 − ‖P_S u₁‖²) with ‖P_S u₁‖² = 0.08–0.72. F_θ quotient ghost persistence is Haynsworth (a theorem), not evidence. The errata match §6 (§7). Round 15 is audited in AUDIT15.md.
- RH: OPEN. Nothing in 14B–14E2 moves it.

---

## 1. Round 14B: errata summary and check of how they were incorporated

The errata I raised on 14B (ea7490f), and how they appear in the block at the top of `round14b_euler/ROUND14B.md` (added in b581f8c), which 14C §0 summarises:

| # | my finding | incorporated? |
|---|---|---|
| framing | "No ghosts ⇔ RH" is Weil's criterion (1952) restated | yes, verbatim in spirit |
| a | ζ "0 ghosts, L ≤ 8" is a code check. The zeros in that band are known to be on the line, and the count depends on the threshold | yes |
| b | F_θ control is tautological: W_Fθ(g) = 2B_ζ(e^{−θu}g, e^{θu}g) is a cross term and indefinite by construction | yes. The appended remark "Euler + FE insufficient is standard (Ramanujan / Selberg class; Kaczorowski–Perelli)" is fair. The Kaczorowski–Perelli attribution is generic; no specific paper is pinned, so treat it as a pointer, not a citation |
| c | Spliced Beurling ghosts measure the borrowed-archimedean defect, not zeros | yes ("defect / not applicable") |
| d | A real off-line pair {β, 1−β} gives 1 ghost, not 2. The EPS L = 3 odd count is a threshold artifact | yes, and it adds "count with multiplicity" and "peak-frequency matching is heuristic". Accurate |
| e | Euclidean OS failure is general (Bernstein–Widder + spectrum condition) | yes, accurately stated |
| f | Tate "same time" is INTERPRETATION; at p the Fourier transform is the Weyl element | yes |
| — | Prolate alone cannot give exact ‖Tf‖² = W past log 2 (CC 2021) | carried into 14C §4, confirmed numerically there |

**Verdict on 14B errata:** incorporated accurately. No residual mislabel found.

---

## 2. Round 14C (b581f8c): prolate operator as T in ‖Tf‖² = W

**Verdict: mostly sound, labelled honestly.** The central conclusion holds: no T built from archimedean or finite-S data can give ‖Tg‖² = W past log 2. But it holds largely by construction. Two CHECKED items are tautologies, and §6's "rules out every finite-S / semilocal-prolate route" misdescribes the Connes–Consani semilocal program and overclaims.

### 2.1 Label table

| claim | 14C label | audit | note |
|---|---|---|---|
| CM 2021/22 and CCM 2024 statements (deficiency (4,4), Fourier-commuting W_sa, Sonin, Thm numbers) | CHECKED (read from arXiv) | **not independently verified by me** | Plausible and consistent with my reading of the abstracts. I did not re-check lemma and theorem numbers |
| CC 2021: Weil positivity at the archimedean place for support in [2^{−1/2}, 2^{1/2}] | fact | **verified** (arXiv 2006.13771, Introduction: test function support in (1/2, 2) for f = g ∗ g*, no rational prime involved) | Caveat: CC's statement (2) puts the conditions g̃(z) = 0 for z ∈ F ⊇ {0, 1} on the test functions, so the pole terms are removed. 14C's T_∞ = arch + pole on *all* window g is the same object once the poles are included, not literally CC's functional. Harmless, but say so |
| UV main term comes only from Γ_ℝ and the conductor | PROVED (elementary) | **earned** | Riemann–von Mangoldt main term θ(T)/π |
| T(ζ) = T(BEUR″), T(χ) = T(DH) | PROVED (construction), CHECKED (0.0 difference) | **PROVED earned; CHECKED is vacuous** | `pairtest.py` copies the arch data: BEUR″ borrows ζ's, and `chi5()` copies DH's `C, w`. A 0.0 difference was guaranteed |
| "Zero defect for L < log 2" | CHECKED | **tautology; relabel as "identity by support"** | The code masks atoms with log n < L, so for L < log 2 the prime matrix is empty by construction |
| "Defect is exactly the n = 2 atom at log 2 + 0.05" | CHECKED | **tautology** | `defect_is_prime_term = allclose(W − Ainf, −q)` where `W := Ainf − q`. That is an identity check of an assignment. The real content, that only n = 2 lies in (log 2, log 3), is arithmetic |
| T_∞ exists (arch + pole ≥ 0) up to ≈0.74, indefinite from 0.8 | NUMERICAL | **earned** | Reproduced from `pairtest.json`: λ_min(arch + pole) = +1.06e-3 at 0.743, −0.065 at 0.8 |
| Pair sign table, atom split | NUMERICAL | **earned**, with the label issue below | |
| Semilocal scan: W_S indefinite unless S ⊇ {p < e^L} | NUMERICAL | **earned as numbers, but uninformative** (see 2.2.3) | |
| "Rules out every finite-S / semilocal-prolate route" | PROVED (decomposition + CC 2021) | **overclaim** (see 2.2.4) | |

### 2.2 Issues

1. **Label clash: "L(χ₅)" is not the χ₅ of Round 13j.** In `pairtest.py`, `chi5()` uses χ(2) = i, the odd quartic character mod 5, and keeps Re χ(n). On real Galerkin vectors, W with coefficients Λ(n)·Re χ(n) is exactly W_{L(χ)} = W_{L(χ̄)}, so the model is correct. It is also the right partner for DH: same Γ_ℝ(s+1), q = 5. But 13j's χ₅ is the *even* quadratic character with Γ_ℝ(s). **Rename it to "L(χ), χ odd quartic mod 5"** so nobody mixes up the two.
2. **"As Tate vectors the two differ only at p = 5" is INTERPRETATION.** It is true at the level of adelic Schwartz data restricted to the archimedean ray, since both coefficient sequences are periodic functions on ℤ/5. It is not a statement about Tate integrals over C_ℚ, where averaging over Ẑ^× projects onto characters. The *finding* (the change shows up as non-Euler composite atoms, not at 5) is correct but expected: DH has no Euler product, so −DH′/DH has composite coefficients.
3. **The semilocal table is generic, not structural.** I added a ±ε copy of a single prime atom to ζ's full W, which is the effect of deleting or adding any atom. Probe `eps_atom.py` (Appendix):
   - At L = 2 (M = 37), W has **8** eigenvalues below 1e-8·max; at L = 3 (M = 51) it has **22**. These are the Connes–Consani "very small eigenvalues".
   - ±1e-9 × (atom at log 2) gives λ_min ≈ −1.2e-9 to −1.6e-9 **for either sign**, at L = 2 and at L = 3.
   - So W sits numerically on the boundary of the PSD cone, and *any* omitted atom of *either* sign makes it indefinite. "W_S ≥ 0 only when S is complete" is what you get for any perturbation. It carries no information about semilocal structure.
4. **Overclaim in §6, "rules out every finite-S / semilocal-prolate route".**
   - The CC program (CC 2021 Introduction) does *not* ask the truncated form arch + pole − prime_S to be ≥ 0 on a window that is too wide for S. It uses S = {∞, 2, …, p} for test functions supported in (p^{−1}, p): **S complete for the support**, which is exactly the case 14C dismisses as circular.
   - Their positivity is not T := W^{1/2}. It is an independent positive functional Tr(ϑ(g) S ϑ(g)*), with S the Sonin projection, which is positive by construction. The open part is controlling the difference from W, which CC show is compact (an infinitesimal) at the archimedean place.
   - So 14C rules out only "W_S ≥ 0 with incomplete S" (expected, item 3) and "T = W^{1/2}" (circular by definition). **The CC semilocal route is untouched.** Replace the sentence with: "rules out positivity of incomplete truncations; does not address the CC Sonin-compression route."
5. **The ζ/EPS "wrong direction" is real, but cheap:** EPS's archimedean constant includes log 20 versus log 1. A larger conductor makes arch more positive for every g. It is correctly reported; just don't read structure into it.

### 2.3 Falsifier for 14C's remaining claim
The only non-tautological claim left is "the positivity mechanism must act on the complete prime term at once." **Cheapest falsifier:** implement CC's Sonin-trace functional τ_S(g) = Tr(ϑ(g) S ϑ(g)*) at the archimedean place for L < log 2, and check that W − τ_S has the small compact defect CC prove. If it does, the next test is the {∞, 2} semilocal analogue at log 2 < L < log 3, which is the first case where a prime enters and S is complete for the support. That test directly addresses the route 14C claims to rule out.

---

## 3. Round 14D (7e89769): one global object, three faces

**Verdict: correct theorems, weak checks, and one framing overclaim.**
- Faces 1+2: PROVED, but this is Tate's thesis (1950), not a new object.
- 3a (range E_F ⊆ radical W_F): a correct theorem. It is the trivial direction of the explicit formula. The numerical "CHECK" does **not** test radical membership: it measures isotropy W(v, v) ≈ 0 with a low-power probe.
- 3b pencil: numbers reproduced, but they add nothing beyond λ_min(W). The "R indefinite past log 2 = CC 2021 boundary" reading is a coincidence of thresholds.
- "Positivity must live on the quotient, ⇔ RH": true, and it is exactly Weil's criterion. Quotienting by a subspace of the radical cannot change positivity. The attribution to Connes 1999 is partly wrong.

### 3.1 Label table

| claim | 14D label | audit | note |
|---|---|---|---|
| Faces 1+2: Tate Θ on C_ℚ; Poisson Θ(a^{−1}) = \|a\|Θ(a) + (\|a\|−1)f₀(0); Mellin gives Λ(s) | PROVED | **earned (classical)** | I re-derived the Poisson identity from Σ_{q∈ℚ} f(qa) = \|a\|^{−1} Σ f̂(q/a) with f̂₀ = f₀: correct. The ray restriction gives Jacobi; unfolding gives Γ_ℝ(s)ζ(s). Credit it as Riemann/Tate, and call it "Daniel's framing of Tate's thesis", not a new result |
| 3a: E_F H admissible and W_F((E_F H) ∗ Θg) = 0 | PROVED (generic) | **earned** | Mellin(E_F H) = F(s)Ĥ(s). It vanishes at every zero of F, with multiplicity, and Ĥ kills the poles; the explicit formula then gives 0. Admissibility: shifting the contour left gives rapid decay as x → 0, but g = e^{u/2}E_F H(e^u) has support (−∞, log b], so **it is not a window function**. See 3.2.4 |
| 3a numerics "CHECKED rel ≤ 1e-4" | CHECKED | **mislabelled: it tests isotropy, not radical membership, and has low power** | See 3.2.1–3.2.2 |
| BEUR″: "E_B is ~ in ζ's radical (6e-6)" | in the pair-test table | **not supported** | At this probe's power, a same-band *non*-radical vector also scores ~4e-6 (3.2.2). 14D's own caveat ("too small to resolve") is right; the table entry contradicts it |
| Plancherel normalization | CHECKED | **earned** | 3.0e-4 (L = 1) and 3.8e-5 (L = 4), reproduced from `pencil.json` |
| N_F = ‖E_F g‖² ≥ 0 for every F | PROVED | **earned, with a domain caveat** | N_F ≥ 0 is trivial (integral of \|F\|²\|ĝ\|²). The identity with ‖E_F g‖² holds only when ĝ vanishes at F's poles, otherwise E_F g ∉ L². The pencil uses all window g, so there N_F is a *defined* form, not a norm of E_F g. The script's docstring says this; the md doesn't |
| Pencil and R-sign tables | NUMERICAL | **earned as numbers** | Reproduced from `pencil.json`. ζ: μ_min tracks λ_min(W) (ratio ≈ 0.7–1.4) because λ_min(N) ∈ [0.25, 0.73]. So the pencil restates λ_min(W) → 0 (the CC small eigenvalues) |
| "R indefinite past log 2 (CC 2021 boundary), CHECKED consistent" | CHECKED | **overclaim; NUMERICAL coincidence** | R = N − W involves N, which CC say nothing about. DH's R is already indefinite **below** log 2 (R_min = −1.78 at log 2 − 0.05), and EPS's R stays ≥ 0 until L ≈ 2 (R_min = −0.0098 at L = 2). The sign change of R is model-specific. The CC-consistent fact is λ_min(W) > 0 for L < log 2, already in 14C |
| EPS "μ_min < 0 at L ≥ 3" | NUMERICAL | **borderline at L = 3** | μ_min = −1.7e-6 and λ_min(W) = −1.7e-5, at the noise / threshold level. 14B erratum (d) already calls the L = 3 EPS ghost a threshold artifact. Write "L > 3" |
| sign(μ_min) = sign(λ_min W) tautological | PROVED, tautological | **earned, and correctly flagged** | Good practice |
| 3c: W weighs \|ĝ\|² by the zero-counting measure θ′/π dt + dS_F | PROVED (elementary) | **informal; true as a positive measure only under RH** | Off-line zeros enter W through ĝ at complex points, not through a measure on the line. With S defined as arg F along the critical line, the off-line pair contributions to Re F′/F cancel. Label it INTERPRETATION. The conclusion of 3c (no injective T on range E) needs only 3a and is fine |
| Quotient H_W = admissible / closure(range E) with W; positivity ⇔ RH | OPEN (⇔ RH) | **⇔ RH is correct but is Weil's criterion; attribution and topology need fixing** | See 3.2.3–3.2.4 |

### 3.2 Issues

1. **Isotropy is not radical membership.** `radical.py::W_grid` computes only the *quadratic* value W(v, v). For ζ (W ≥ 0 in the band) isotropy implies radical, but for the fakes W is indefinite, and there W(v, v) = 0 does not imply W(v, ·) = 0. The theorem is fine; the "CHECK" doesn't test the claimed property. **Fix:** test the bilinear form W(v, φ) against probes φ that have spectral mass at the model's zeros (Appendix `radprobe.py` does this for ζ).
2. **The quadratic check has too little power to separate radical from non-radical.** H is a C^∞ bump on log x ∈ [0.5, 2.5], and its Mellin transform at the first ζ zero is only 1.4e-2 of its peak (\|Ĥ(½ + i·14.13)\|/max). W(v, v) is quadratic in that leakage, so *any* vector in this band has small W. My probe for ζ (`radprobe.py`) compares three vectors:

   | vector | W(v,v)/‖v‖² | bilinear W(v, φ₁₄)/(‖v‖‖φ‖), φ₁₄ = cos²-windowed cos(14.13u) |
   |---|---|---|
   | E_ζH (radical, by theorem) | −2.2e-6 | −4.7e-6 |
   | E with a(2) = 0.9 instead of 1 (same band, **not** in the radical) | **+4.3e-6** | **−1.05e-3** |
   | H itself (14D's control) | 6.1e-4 | −1.7e-2 |

   - The quadratic test cannot tell a radical vector from a same-band non-radical one (2.2e-6 vs 4.3e-6, both under 14D's 1e-4 tolerance). The bilinear probe separates them by 220×.
   - With a narrower bump on [0.5, 1.0], so that Ĥ is O(1) at height 14: quadratic 8.9e-5 vs 1.8e-3 (20×), bilinear 1.8e-5 vs 3.2e-2 (≈1800×).
   - The "H itself" control is a vector from a *different* band, so it cannot serve as a control for membership.
   - The same issue hits **F_θ**: in absolute terms, own W = 3.5e-8 and cross W_ζ = 6.2e-8. Since F_θ's zeros are ζ's ordinates shifted off the line, at heights ≥ 14, both are tiny for the same reason. EPS and DH separate only because they have low-lying zeros where Ĥ is large.
   - **Conclusion:** "CHECKED rel ≤ 1e-4" should read "consistent with the theorem at the grid floor (DU = 2e-3, ~1e-6); not a discriminating check."
   - Also, the denominator max(\|arch\|, \|pole\|, \|prime\|) rewards cancellation between large parts. Normalise by ‖v‖² or by the probe instead.
3. **"Positivity must live on the quotient ⇔ RH" adds nothing to Weil.** For e in the radical, W(x + e, x + e) = W(x, x). So W ≥ 0 on the quotient by closure(range E) ⇔ W ≥ 0 on the whole admissible space ⇔ RH. That is Weil's criterion; CC 2021 state it as their eq. (2) with the support refinement. "A positive form factors through its null quotient" is GNS, not a discovery. The faithful statement is: **"if W ≥ 0, its GNS space is a quotient of the admissible space by a subspace containing range E; whether that subspace equals closure(range E) is the spectral-realization question."**
4. **Attribution and topology.**
   - **What Connes 1999 does.** (Selecta Math. 5 (1999) 29–106; arXiv math/9811068.) The abstract: critical zeros appear as an *absorption spectrum* and non-critical zeros as *resonances*, in the L² (or Sobolev-weighted L²_δ) setting of the adèle class space. The paper "reduces the Riemann hypothesis to the validity of the trace formula." Theorem 5 (positive characteristic) proves that equivalence. Connes' cokernel carries the **L² (or L²_δ) inner product, not the W-norm**. 14D's "H_W = admissible / closure(range E) with the W-norm is Connes' spectral realization" conflates Weil's GNS with Connes' L² cokernel. They are related but not the same Hilbert structure.
   - **Topology.** W is unbounded relative to L² (its arch part grows like log\|t\| in frequency). So "W vanishes on range E" does **not** extend automatically to the *L²*-closure. The closure has to be taken in the form-domain norm, e.g. ‖g‖² + ∫\|ĝ\|² log(2 + \|t\|). Connes' handling of δ is exactly this kind of issue.
   - **Window issue.** E_F H has a half-infinite tail in u, so range E does not lie in any window space [−L/2, L/2]. I suspect range E ∩ C_c = {0}: compact support would force F(s)Ĥ(s) to have exponential type, against ζ's factorial growth to the left. That is plausible but **not proved here; not a citation.** Either way, the window-truncated range E is **not** in the radical, and that truncation error is precisely what CC handle with prolate functions in "Spectral triples and ζ-cycles" (arXiv 2106.01715; EMS Press lists L'Enseignement Math. 69 (2023), and the arXiv HTML template says J. Inst. Math. Jussieu, so check the venue before citing). There, Σ_μ E(S^ev_0) on a circle of length L fails to be dense exactly for "ζ-cycles", and the scaling spectrum on its complement gives critical zeros.
   - **So what is known vs new:**
     - Known: W ≥ 0 ⇔ RH (Weil 1952). Range E in the null space, and the cokernel spectral realization (Connes 1999). Archimedean positivity for support < log 2 (CC 2021). Small eigenvalues and the window version Σ_μ E (CC ζ-cycles 2021/23). Semilocal prolate / Sonin (CCM 2024, arXiv 2310.18423).
     - New in 14D: the explicit pair-test bookkeeping across the fakes, and the observation that 3a is "fake-blind". Both are correct and both are expected.
5. **3b's "fakes more coercive at small L"** is correctly reported and honest. Note that it reflects conductor (EPS q = 20, DH q = 5) and low-lying-zero structure, not arithmetic.

### 3.3 Falsifier for 14D's remaining claim
The claim worth testing is "closure(range E_F) = radical(W_F)" (equality, not just inclusion) in a stated topology. **Cheapest falsifier:**
- For ζ at L = 2–3, build V = span of window-*truncated* E_ζH_k for a family of H_k.
- Measure the bilinear residual sup_φ \|W(v, φ)\|/‖v‖‖φ‖ over v ∈ V, against the dimension of W's near-null space (8 at L = 2, 22 at L = 3, M = 37 and 51).
- If truncated range E cannot reproduce the near-null space to the residual floor, equality fails at the window level, and 14E's plan (next section) inherits the gap.

---

## 4. Forward look: biggest hole and cheapest falsifier for 14E and 14F

### 14E (λ_min of W on the quotient, uniform in L): biggest hole
**The quotient subspace V is not canonically defined at finite L, and W has no spectral gap to cut at.** Probe `spec.py` (Appendix), ζ, Dirichlet-sine Galerkin:

| L | M (fmax 45 → 90) | # eigenvalues < 1e-8·max | next four eigenvalues (fmax 45) |
|---|---|---|---|
| 2 | 37 → 66 | 8 → 8 | 1.2e-7, 1.2e-5, 4.2e-4, 1.3e-2 |
| 3 | 51 → 94 | 22 → 30 | 1.1e-7, 1.1e-6, 2.1e-5, 2.5e-4 |
| 4 | 66 → 123 | 38 → 59 | 1.5e-6, 2.5e-6, 2.2e-4, 3.9e-4 |
| 5 | 80 → 152 | 52 → 88 | 2.2e-7, 4.3e-7, 2.7e-5, 7.5e-5 |

- The spectrum decays geometrically *through* any threshold, with no plateau.
- The near-null count changes with basis size at fixed L.
- So "λ_min on the complement of V" depends on how V is chosen:
  - **V by eigenvalue threshold:** λ_min ≥ threshold is a tautology.
  - **V = window-truncated range E:** V is not in the radical (3.2.4), and the residual leaks into λ_min.
- Two more points:
  - Even a clean uniform gap for ζ at L ≤ 5–8 is a **code check**: in-band zeros are verified on the line.
  - Heuristically (Ingham/Landau type, **not proved here**), λ_min on the quotient is governed by the lower Riesz bound of {e^{iγu}} on [−L/2, L/2] near the critical height t ≈ 2πe^L, where zero spacing ≈ 2π/L. Close zero pairs (Lehmer pairs, e.g. near t ≈ 7005) should push it toward 0 around L ≈ 7. So a uniform-in-L gap would be a statement *stronger* than RH (it needs simple zeros with gap control) and may well be false.

**Cheapest falsifiers (seconds each):**
1. **Sylvester test on the fakes:** quotienting by a subspace of the radical cannot change the inertia. So for EPS (L = 4), F_θ (L = 2) and DH (L = 6) the negative index must be **identical** before and after quotienting by their own E_F. Any change measures V ⊄ radical, i.e. a pipeline error.
2. **Threshold and basis invariance for ζ:** report λ_min(quotient) for V chosen at thresholds 1e-6, 1e-8, 1e-10 and fmax 45 and 90. If it moves with either, the "gap" is an artifact.
3. **Lehmer probe:** at L ≈ 7, put g = (normalised) difference of windowed e^{iγu} at the Lehmer pair near 7005, project it off V, and compute the Rayleigh quotient. If it is far below the reported λ_min, uniformity is falsified.

### 14F (representation theory of the quotient; "single irreducible vs sum of two"): biggest hole
**The group must be named, and the criterion must not presuppose positivity.**
- The scaling group (ℝ₊* or C_ℚ) is abelian, so its irreducible unitary representations are 1-dimensional. ζ's quotient is already an infinite sum of characters \|·\|^{iγ}, one per zero, so "single irreducible" is meaningless for it unless a non-abelian group (which one?) acts.
- For the fakes W is indefinite, so the quotient is a Krein space and "unitary irreducible" is undefined. An off-line pair is a 2-dimensional neutral block of signature (1, 1), the (2,2) cross term of 14B.
- Structurally, W is **additive over products** (W_{FG} = W_F + W_G) and **not additive over sums**. The real dividing line among the fakes is product versus linear combination, plus unitarity, not irreducibility:
  - EPS = ζ_{ℚ(√−5)} + L(χ₋₄)L(χ₅), per `lineval.py`, i.e. genus theory with h(−20) = 2;
  - DH = a combination of L(χ) and L(χ̄);
  - F_θ = ζ(s−θ)ζ(s+θ), a *product* of non-unitary shifts.

**Cheapest falsifier:** run 14F's classifier on two controls already in the repo.
- **ζ_K for K = ℚ(√5)**, i.e. ζ·L(χ₅) with χ₅ even quadratic. 13j showed W_K = W_ζ + W_χ₅ exactly (2e-13) with no ghosts. It is reducible but ghost-free.
- **F_θ**: reducible (a product) and full of ghosts.

Any "irreducible ⇔ no ghosts" criterion misclassifies one of the two. If 14F's criterion separates them, it must be detecting unitarity (the θ shift), not reducibility. Then it should be compared with the Ramanujan / unitarity axiom, which is classical.

---

## 5. Round 14F (team/principal-series 104719a): representation theory of the quotient

Audited: `principal_series/ROUND14F.md`, `r14f_quotient.py`, `r14f_patch.py`, `r14f_collision.py`, `r14f_logderiv.py` and their JSON outputs. Re-runs are listed in 5.5.

**Verdict: a good round with honest bottom lines, but three overreaches.**
- The algebra is correct: the genus decomposition, DH's root number, the Hecke defects and the log-derivative atoms.
- My pre-audit falsifier was adopted and resolved correctly: ζ_{ℚ(√5)} is a ghost-free product, F_θ is a ghosted product, and the separating axiom is unitarity.
- The candidate criterion (an isobaric sum of unitary cuspidal automorphic representations) is the **classical** setting in which GRH is conjectured. The round says so ("selects, does not prove").
- **Overreach 1:** item 3, W_{FG} = W_F + W_G, is a true tautology presented with a CHECKED number.
- **Overreach 2:** item 5 ("ghosts are not spectrum of any Hecke-equivariant quotient") is vacuous as stated. In Connes' framework *no* off-line zero is ever spectrum, equivariant or not.
- **Overreach 3:** in the six-case table, the ζ / ζ_K ghost entries are labelled as evidence, and Beurling's "ghosts" ignore erratum (c).
- Two numerical details need sharpening (5.4).

### 5.1 (a) Item 3: "no cross term over products", W_{FG} = W_F + W_G
- **True, and tautological.** W_F is linear in −F′/F, in log Γ_F and in the pole divisor, and all three are additive under F ↦ FG. 14F labels it PROVED (correct) and "CHECKED to 1e-14". That check, max \|c_{ζL(χ₅)} − c_ζ − c_{L(χ₅)}\| = 9.8e-15, is a floating-point test of the identity log(FG) = log F + log G in the Dirichlet-convolution recursion. **Relabel it "identity; float check", not content.** The "non-additivity over sums" defect of 315 is likewise just the statement that log is not additive over sums.
- **No contradiction with 14B's polarisation identity.** The two identities are bilinearity in *different variables*:
  - **Additivity in F (14F):** W_{Fθ} = W_{ζ(·−θ)} + W_{ζ(·+θ)}.
  - **Twist as a multiplier on the test vector (14B):** shifting s by ±θ multiplies the autocorrelation A_g(t) by e^{±θt}, and e^{θt}A_g(t) = A_{e^{−θu}g, e^{θu}g}(t) is a *cross*-correlation of two different vectors. So W_{ζ(·∓θ)}(g) = B_ζ(e^{∓θu}g, e^{±θu}g), and summing gives **W_{Fθ}(g) = W_ζ(2cosh(θt)A_g) = 2B_ζ(e^{−θu}g, e^{θu}g).**
  - **I verified this term by term.** Prime side: c_{Fθ}(n) = 2cosh(θ log n)Λ(n), as in `models14b.fth`. Arch weight: w_{Fθ} = 2cosh(θt)w_ζ. Poles: cosh((½ ± θ)t). The archimedean constants close exactly: C_{Fθ} − 2C_ζ = ψ(¼+θ/2) + ψ(¼−θ/2) − 2ψ(¼) = −4∫₀^∞ w_ζ(t)(cosh θt − 1)dt, which matches to 1e-14 at θ = 0.1, 0.2, 0.3 (`twistid` in 5.5). So both statements hold exactly.
- **What the twist changes.** For a **unitary** twist \|·\|^{±iτ}, e^{iτt}A_g = A_{e^{iτu}g} is the autocorrelation of a *single* vector. So W_{ζ(·−iτ)}(g) = W_ζ(e^{iτu}g) ≥ 0 whenever W_ζ ≥ 0: positivity is transported. For a **non-unitary** twist \|·\|^{±θ}, the multiplier e^{±θu} is not an isometry. The form becomes a cross term between two different vectors, which is indefinite by construction. That is exactly 14B erratum (b).
- **The right reading of 14F's unitarity axiom is "the twist is a unitary multiplier, so it preserves the positive cone."** That is the one-line reason unitarity is the separating axiom between ζ_{ℚ(√5)} and F_θ.
- **Control run (new, `twist.py`).** F_{iτ} = ζ(s−iτ)ζ(s+iτ) with τ = 3 has the same product structure as F_θ and differs only in unitarity:

  | model | L = 2 (M = 37) | L = 4 (M = 66) |
  |---|---|---|
  | F_θ, θ = 0.2 | λ_min = −0.091, 7 ghosts | λ_min = −0.755, 20 ghosts |
  | F_{iτ}, τ = 3 | λ_min = +1.4e-10, 0 ghosts | λ_min = +6.1e-11, 0 ghosts |

  The pole term was built by quadrature and validated against ζ's closed form to 1e-13. The "0 ghosts" result is a code check, since the zeros are on the line, but it confirms the pipeline separates the pair exactly on unitarity. **Recommend adding F_{iτ} as 14F's seventh row.** It is the cleanest control for the unitarity column.

### 5.2 (b) Item 5: "ghosts are not spectrum of any ℋ-equivariant quotient" and Connes 1999
**What Connes 1999 proves** (Selecta Math. 5 (1999) 29–106, arXiv math/9811068; I re-read the arXiv text):
- **Thm 1 (§III).** In L²_δ(X) with **δ > 1**, the spectrum of D_χ on the cokernel H_χ is the set of imaginary parts of zeros of L(χ̃, s) **with real part ½**. The multiplicity is min(mult, largest n < (1+δ)/2). Connes himself notes that L²_δ "artificially eliminates the non-critical zeros" (§VIII).
- **δ = 0 (cutoff) version.** Critical zeros appear as an **absorption spectrum**, and non-critical zeros appear only as **resonances**, through their harmonic potential in the trace formula.
- **The global trace formula is not proved.** Connes shows its validity is *equivalent* to RH for all Grössencharakter L-functions (in positive characteristic: Thm 5). So nothing in Connes 1999 is unconditional beyond the semi-local trace formula and Thm 1.
- The setting is **GL₁ / Grössencharakters over a global field k**. There is no GL_n (n ≥ 2) cokernel in that paper.

**Audit of item 5:**
1. **Vacuous as stated.** Under Thm 1 (and in any setting with a unitary scaling action), an off-line zero is *never* spectrum of *any* cokernel: not EPS's, not DH's own non-equivariant cokernel, and not ζ's if ζ had one. So "the ghost is not spectrum of any ℋ-equivariant subquotient" holds equally for non-equivariant ones. Equivariance does no work in that sentence. **Overreach: relabel it from PROVED (sketch) to "true but vacuous".**
2. **The correct and non-vacuous content is algebraic, from Tate / Hecke theory, not from Connes.**
   - For an ℋ-isotypic (eigen) vector, the zeta integral is L(π_i, s) times local polynomial factors. So the zero set of any equivariant construction is contained in ∪ zeros(L(π_i)) ∪ (local-factor zeros).
   - EPS's and DH's off-line zeros are not in that union: the pieces are non-zero there (\|ζ_K\| = \|L_ψ\| = 0.933 at the EPS zero, \|L(χ)\| = \|L(χ̄)\| = 0.702 at the DH zero). 14F §2c's "same mechanism for ζ" (the 14A vector with the factor (1 − 2^{1−s})) is exactly this Tate gcd statement.
   - **Restate item 5 as:** "any off-line zero of an ℋ-equivariant zeta integral is an off-line zero of some L(π_i). The EPS / DH off-line zeros are not, so they come from non-equivariance (interference)." That is PROVED, given standard Tate / Jacquet–Langlands local theory.
3. **The isotypic split V = V_χ ⊕ V_χ̄ under ℤ₅^×** is correct (a compact group acting by distinct characters). The cokernel spectrum "Z_crit(L(χ)) ∪ Z_crit(L(χ̄))" follows from Thm 1 per character, and so contains critical zeros only, as Thm 1 says.
4. **"Equivalently: positivity of W_π on Connes' π-isotypic quotient is GRH for π"** is fine for GL₁ π. For GL_n with n ≥ 2 there is no "Connes quotient" in Connes 1999. For general π, cite Weil's criterion for L(π) (and its standard extensions), not Connes. **This is a label fix.**
5. **"Equal moduli at the zeros" is a tautology.** At any zero of αF₁ + βF₂, F₁/F₂ = −β/α, and here \|β/α\| = 1 (DH: β/α = ε(χ); EPS: α = β = 1). So \|F₁\| = \|F₂\| is forced, whether the zero is on or off the line. It should not be labelled NUMERICAL evidence. The non-trivial number is the *size* (0.70, 0.93), i.e. that the pieces are not near their own zeros.

### 5.3 (c) The six-case classification

| case | 14F verdict | criterion label PROVED? | ghost-column label | audit |
|---|---|---|---|---|
| ζ | passes | yes (GL₁ trivial) | "none found, NUMERICAL" | **relabel CHECKED (code check).** 14B erratum (a): zeros in band are known on the line |
| ζ_{ℚ(√5)} = ζ·L(χ₅) | passes, a product | yes. 1 ⊞ χ₅, χ₅ even quadratic | "none (13j), NUMERICAL" | same, CHECKED. Correct as my falsifier's positive control |
| Epstein x²+5y² | fails: not ℋ-equivariant | **yes**, and there is a simpler proof than the Hecke defect: a(6) = 2 ≠ a(2)a(3) = 0 (a = r_{Q₁}/2), so there is no Euler product, so it is not any L(Π) | quartet 0.933+15.67i, a verified zero | correct. Z_{Q₁} = ζ_K + L_ψ with constants exactly 1 (I re-checked the coefficient identity r_{Q₁} = b₁+b₂, r_{Q₂} = b₁−b₂ to n ≤ 300 with my own lattice count, 0 failures) |
| DH | fails: not equivariant under ℤ₅^× | yes. No Euler product (906 non-prime-power atoms ≤ 3000) | 0.8085+85.70i, verified | correct. β/α = ε(χ) correct, consistent with 14C's χ (odd quartic mod 5) |
| F_θ | fails: unitarity | **yes**. By Jacquet–Shalika the isobaric decomposition is unique, so L = L(\|·\|^θ ⊞ \|·\|^{−θ}) has no unitary-constituent representation | every zero off by θ | correct, with one label point: the column "unitary (Ramanujan at p)" conflates two axioms. For GL₁ they coincide (\|χ(p)\| = 1). For GL_n, unitary cuspidal π are not known to be tempered (Ramanujan is open for Maass forms). So c(p)/log p > degree detects non-temperedness, not non-unitarity in general. Fine for these six; rename the column |
| spliced Beurling | fails: not automorphic | yes. A generalised Dirichlet series over non-integers is not any L(Π) | "yes, proliferation (14B)" | **mislabel.** Erratum (c): there is no own explicit formula, so its form is a defect, not ghosts. The entry should be "N.A. / defect" |

- **My falsifier.** It predicted that any "irreducible ⇔ no ghosts" criterion misclassifies ζ_K or F_θ, and that a criterion separating them must be detecting unitarity. 14F confirms both and adopts unitarity explicitly. **The falsifier is resolved correctly.**
- **Consequence.** The criterion as a whole (automorphic + isobaric + unitary) is the classical domain of GRH (Langlands; Selberg-class conjectures; Kaczorowski–Perelli 1999 for degree ≤ 1, as 14F cites). 14F says so. The new content is expository and numerical: the interference picture, and the collision homotopy.
- **"PROVED for the classification" is earned for the predicate.** "Classifies correctly" against ghost-freeness is conditional on RH/GRH for ζ and ζ_K, and is NUMERICAL only for the three ghosted fakes. Beurling is not a ghost case at all.

### 5.4 Numerical details
- **Collision λ*.** 14F gives "≈ 0.14, bracket (0.06, 0.16), √-extrapolation 0.138". I counted on-line zeros of the real Λ_λ on t ∈ [14.5, 16.2] (step 0.005): the pair is at 15.22 / 15.255 at λ = 0.140, and there are no on-line zeros at λ = 0.145, 0.15 or 0.16. So **λ* ∈ (0.140, 0.145)**, collision height ≈ 15.24. The 0.138 extrapolation is slightly low. Tighten the bracket; the conclusion is unchanged. The mechanism (a real-coefficient, FE-symmetric family losing on-line zeros in colliding pairs) is classical.
- **Reproducibility trap.** `r14f_quotient.py` as committed tests `rQ == 2*(b1+b2)`. That is wrong, since r_{Q₁}(1) = 2 = b₁(1) + b₂(1), and it would report failures at almost every n. Only `r14f_patch.py` writes the correct key `r_Q1 = b1+b2` into the JSON, and `r14f_collision.py` then pops the bisection result. Running the scripts in a different order gives a different JSON. Fix the first script.
- **Verified as stated:** the Hecke defects of r_{Q₁}/2 (−4, −4, −4, −2, −4) and of b₁ (0); the DH translation ratios (0.284, −3.52); the FE check (9e-31); the s = 3 lattice sums.

### 5.5 Biggest hole and falsifier
- **Biggest hole (stated by 14F itself).** The criterion selects but does not prove: inside one unitary cuspidal L(π), W_π ≥ 0 is GRH. Item 5's restated content (5.2.2) excludes only *interference* ghosts.
- **Cheapest non-tautological test of the "selects" claim.** Add a pair inside the class that differs only in *unitarity*: F_{iτ} vs F_θ (done above: 0 vs 7–20 ghosts). Then add the hardest classical near-miss: a degree-2 Euler product with FE whose local parameters are unitary but which is **not** automorphic. The Selberg-class conjecture says none exists, so the criterion makes a falsifiable prediction only there. No cheap instance is known to me.
- **Probes (run from `round14c_prolate/`, ~50 s total):**
  - `twistid`: the digamma identity above, via scipy quad;
  - `twist.py`: F_{iτ} vs F_θ Galerkin with the quadrature pole term;
  - `f14.py`: the integer identities and the λ* scan.
  - The code is in my session scratch (`/tmp/ge/twist.py`, `/tmp/ge/f14.py`); key lines are reproduced below.

```python
# twist.py core: F_{i tau} model for the 14B Galerkin (pole term by quadrature, validated on zeta to 1e-13)
C  = 2*digamma(complex(0.25, tau/2)).real - 2*math.log(math.pi)
mo = dict(n=Z['n'], c=Z['c']*2*np.cos(tau*np.log(Z['n'])), C=C, w=lambda t: 2*np.cos(tau*t)*md.wz(t), poles=[])
kern = lambda t: 2*math.cosh(t/2)*math.cos(tau*t)       # poles at 1/2 +- i tau
# pole matrix P = 2 * int_0^L (A(t)+A(t)^T) kern(t) dt  (Gauss-Legendre, 400 nodes, wg.Aclosed)
# twist identity: digamma(1/4+th/2)+digamma(1/4-th/2)-2*digamma(1/4) == -4*int_0^inf wz(t)*(cosh(th t)-1) dt
```

## 6. Round 14E (main cbcf08a): W on the complement of the window-truncated range of E

Audited: `round14e_gap/ROUND14E.md`, `gap.py`, `run_gap.py`, `bilinear.py`, `fits.py`, `lehmer.py`, `lehmer_opt.py` and their JSON outputs. I ran everything from a detached worktree of origin/main, so nothing in it was edited. Re-runs are listed in 6.8.

**Verdict.** 14E took the pre-audit seriously, reports its vacuous cases honestly, and its main conclusion stands: no cut-independent positive gap is visible for ζ, and a uniform-in-L gap is falsified at L = 8 at the 0.06 level. But four points need correcting:
- (1) The Sylvester check is a restriction, which can only lower the negative index, and it has **almost no power**: a *random* subspace of the same dimension also preserves every ghost.
- (2) The "float floor" for L ≥ 1.5 is a precision limit, not a zero. In 80-digit arithmetic, λ_min(W_ζ) is positive and resolved down to 7.8e-30 at L = 2.
- (3) The Lehmer bound is correct (constant 12 and L³ verified, value reproduced). The zeros there are **proof inputs**, not just labels, and `lehmer.py::optimal` does not run as committed.
- (4) The Gaussian fit is **contradicted by high-precision data** just outside its range. The decay is double-exponential, ≈ exp(−B e^L) with B ≈ 11–13, consistent with Slepian/prolate decay at scale S = e^{L/2}. The extrapolation to L = 8 (≈ e^{−994}) is unsupported.

### 6.1 Errata blocks in 14C/14D (applied in cbcf08a): do they match this audit?
- **14C, items 1–6: match** §2 (tautological CHECKs, vacuous 0.0 difference, χ rename, semilocal claim withdrawn, Tate-vector interpretation). The withdrawal sentence in §6 is struck through and replaced correctly.
- **14D, items 1–6: match** §3.2: the easy direction; isotropy vs radical with the bilinear redo; BEUR″ entry withdrawn; quotient positivity = Weil plus GNS, with the L² vs W-norm and form-domain closure; the R-boundary coincidence; EPS "L > 3"; and Tate credit.
- **Three minor items from §3 are not carried over.** Add them for completeness:
  - (i) 3c's "zero-counting measure on the line" is a positive measure only under RH, so label it INTERPRETATION.
  - (ii) N_F = ‖E_F g‖² only when ĝ vanishes at F's poles.
  - (iii) Window-truncated range E is not in the radical, since E_F H has a half-infinite tail.
- I spot-checked the 14D erratum's new separations against `bilinear.json`. They match the table in 6.6.

### 6.2 (1) Sylvester "ghost-count invariance": which direction is proved?
- **What the code does.** `gap.py::analyse` computes Q = null_space(Sᵀ) and the eigenvalues of QᵀWQ. That is the **restriction (compression)** of W to the L²-complement of S, not a quotient form.
- **What is proved.**
  - (a) If S ⊆ radical exactly (WS = 0), then W = QQᵀWQQᵀ. Restricting to S^⊥ then equals the quotient W/S, and the inertia is exactly preserved. This is the classical statement 14E cites.
  - (b) For *any* S, Cauchy interlacing gives n₋(W) − dim S ≤ n₋(QᵀWQ) ≤ n₋(W).
  - So the pipeline can only *lose* ghosts. "Invariance" proves only the lower inequality is not active. BEUR″'s failure (10 → 8/7/5 and 40 → 6/3/3) is a loss, in the only possible direction. The ζ control (0 → 0) is automatic. EPS L = 4 and F_θ L = 2 are vacuous (dim S = 0), as 14E says.
- **Power: low. New control, my re-run.** At F_θ L = 4 (M = 66) I reproduced dim S = 14, 20 → 20 ghosts, λ_min −0.75496155371 → −0.75496155359, and bilinear residual 6.7e-6. Then:
  - **3 random 14-dimensional subspaces also give 20 → 20**;
  - only S = the span of the 14 lowest eigenvectors drops to 6, the interlacing floor.
  - So keeping the ghost count is the *generic* outcome of a small restriction. It is not evidence that S ⊆ radical; the evidence for that is the bilinear residual. The ghost count can only catch an S aligned with the negative eigenspace, which is what happened for BEUR″ (ζ's range overlaps BEUR″'s defect directions).
- **Relabel:** "inertia preserved: consistency check, one-sided (interlacing), low power. Radical membership: bilinear CHECK (6.6)." The exact claim 6.2 in 14E ("quotienting by range E_F cannot change the inertia") is PROVED for the *exact* quotient. The pipeline tests the restriction.

### 6.3 (2) ζ: cut and basis instability of λ_min
- **Confirmed and correctly reported.** The near-null count depends on the basis (L = 8: 97 vs 176), and λ_min(W|S^⊥) sits at the double-precision floor for L ≥ 1.5 for every ε and both bases.
- **Correction: the "floor" is a precision limit, not λ_min = 0.** I rebuilt ζ's Galerkin matrix in mpmath at 80 digits, with the same basis and formula. The archimedean integral uses Gauss–Legendre, 4 panels × 192 nodes; degree 7 and degree 8 agree to 8 digits at L = 2, M = 66. The prime and pole terms are closed form. W itself (no quotient) is positive with λ_min resolved far below 1e-15. One "floor" entry was really a value: at L = 1.5, 14E's 4.6e-15 / 3.2e-15 compare with true 6.4e-15 (M = 30) and 4.1e-15 (M = 73), so it sat right at the floor. The quotient construction S is built in double precision, so the double "floor" for λ_min(W|S^⊥) cannot tell 0 from 1e-30. **Rephrase: "undetermined at double precision".** 14E's own wording ("no cut-independent positive value") is fine.
- **Basis non-convergence grows with L.** At L = 2, λ_min = 1.28e-27 (M = 37, fmax 45), 1.01e-29 (M = 66), 7.8e-30 (M = 94). The fmax-45 values are Ritz upper bounds and are off by a factor of up to ~160 at L = 2. 14E's "basis-converged to ≤ 5% (11% at L = 1.25)" holds only for L ≤ 1.25.

### 6.4 (3) Lehmer-pair bound: derivation, constant, inputs
- **Derivation (checked).**
  - Assume V ⊂ L²[−L/2, L/2] with ê(γ₁) = ê(γ₂) = 0.
  - For h = g − e: ĥ(γ_i) = ĝ(γ_i), and |ĥ(γ₁) − ĥ(γ₂)| ≤ ∫|h(u)|·|e^{iγ₁u} − e^{iγ₂u}| du ≤ Δ∫|u||h| ≤ Δ‖h‖₂ (∫_{−L/2}^{L/2} u² du)^{1/2}.
  - Since ∫_{−L/2}^{L/2} u² du = **L³/12** (centring at the window midpoint is optimal), we get dist(g, V)² ≥ 12|Δĝ|²/(Δ²L³), hence **W(g)/dist(g,V)² ≤ W(g)Δ²L³/(12|Δĝ|²)**.
  - The constant 12 and L³ are correct. "Ideal" Δ²L³/24 uses W ≥ |ĝ(γ₁)|² + |ĝ(γ₂)|² ≥ |Δĝ|²/2, so q ≤ 2. That is correct.
  - The true radical lies inside {ê(γ₁) = ê(γ₂) = 0}, and the form-domain distance is ≥ the L² distance, so the bound transfers to the true quotient.
- **Caveats.**
  - (i) V must be *window-supported*. If e carries E's half-infinite tail, ∫|u||h| over ℝ is not controlled by √(L³/12)‖h‖. 14E's text omits this hypothesis; add it.
  - (ii) The bound is never better than the probe's own Rayleigh quotient: applying the same inequality to g gives bound ≥ W(g)/‖g‖². At L = 8 the probe has W/‖g‖² = 0.0094 and bound 0.063, so the bound's content is that dist(g, rad) ≥ 0.39‖g‖.
- **Reproduced.** Pair #6709/6710 at t = 7005.0629 / 7005.1006, Δ = 0.03770. Bound **0.1203 (L = 7), 0.0629 (L = 8)**, q = 0.338 and 0.964. These match 14E's 0.120 and 0.063. The data-side W exceeds the truncated zero side by 7.0% and 3.8%, as expected since the low zeros are omitted from the zero side; using the data value is conservative.
- **Zeros: inputs, not labels.** γ₁ and γ₂ (from `mp.zetazero`) are used to (a) build the probe and (b) make the bound valid at all: ê(γ_i) = 0 for e in the radical holds *because* the γ_i are zeros of ζ. That is legitimate, since an upper bound needs a witness and mpmath's zeros are reliable to the precision used. But 14E's header ("zeros are never proof inputs; only CHECK inputs in §5") is inaccurate for §5. **Relabel:** "zeros are proof inputs to the falsifier, verified numerically; W(g) itself comes from integer data."
- **Reproducibility bug.** `lehmer.py::optimal` uses an undefined `d` (NameError on import-and-call; I reproduced it). `lehmer_opt.json` therefore cannot come from the committed code. My fix is one line, `d = ph([γ₁])[:,0] − ph([γ₂])[:,0]`, and it reproduces the table.
- **Scope of the conclusion.** "No uniform gap ≳ 0.07 at L = 8" is correct for the window quotient (NUMERICAL, with a PROVED bound). It is weaker than the double-exponentially small λ_min(W) already visible at L = 2 (6.5): if the near-null directions left in S^⊥ are genuinely outside the radical, λ_min(quotient) is already ≤ 1e-29 at L = 2. Either way there is no gap. The Lehmer probe matters because it is a *proved* witness that doesn't depend on Galerkin artifacts.

### 6.5 (4) The λ_min fit: what was measured, and what the decay really is
- **What was actually measured above the floor.** 14E fits only L ∈ {0.25, 0.5, 0.75, 1.0, 1.25}: 5 points, a 3-parameter fit, 2 degrees of freedom. Every L ≥ 1.5 is at the double floor.
- **Code drift.** The quadratic fit is in `fits.json`, but `fits.py` as committed only computes the *linear* fit. Commit the code that produced it.
- **The fit's own formula gives e^{−994} at L = 8.** If anything downstream uses it, that is a pure extrapolation from L ≤ 1.25, about 400 orders below anything measured.
- **High-precision data (mine, 80 digits):**

| L | 0.25 | 0.5 | 0.75 | 1.0 | 1.25 | 1.5 | 1.75 | 2.0 | 2.25 | 2.5 |
|---|---|---|---|---|---|---|---|---|---|---|
| λ_min(W_ζ), best M | 0.329 | 3.60e-2 | 5.55e-4 | 1.06e-6 | 3.04e-10 | 4.12e-15 | 2.07e-21 | 7.82e-30 | ≤2.6e-40 | ≤3.4e-51 |
| M (fmax) | 16 (90) | 23 | 30 | 37 | 44 | 73 (135) | 84 (135) | 94 (135) | 73 (90) | 80 (90) |

  At L ≤ 1.25 these reproduce 14E's double values to all printed digits. Values at L ≥ 2.25 are Ritz upper bounds, not yet M-converged.
- **The Gaussian fails out of range.** The 14E quadratic fit predicts log λ = −32.1, −44.3, −58.5 at L = 1.5, 1.75, 2.0; measured: −33.1, −47.6, −67.0. At L = 2 it is off by a factor of ~5000, and the error grows with L.
- **Fits on the converged points L ≤ 2** (max log-residual):
  - geometric a + bL: 10.6 (rejected);
  - quadratic: 1.12;
  - a − Be^L: 2.1;
  - **a + bL − Be^L: 0.32, with B = 13.1.**
- **Slope test.** The local slope d log λ / d(e^L) is −10.4, −10.6, −11.3, −11.4, −11.9 over L = 0.75 to 2.0. It is nearly constant and slowly increasing in magnitude. A Gaussian in L would make this slope *decrease* in magnitude (from −10.3 at L ≈ 1.4 to −7.7 at L ≈ 2.1).
- **Conclusion: double-exponential decay, λ_min ≈ exp(−B e^L + O(L)), B ≈ 11–13.**
- **Theory (heuristic, stated as such).**
  - CC "Spectral triples and ζ-cycles" (arXiv 2106.01715) show the very small eigenvalues for support in [S^{−1}, S] have eigenvectors built from prolate functions at scale S. Here S = e^{L/2}.
  - Slepian's asymptotic for the prolate concentration defect is 1 − λ_n(c) ∝ c^{n+½} e^{−2c}, and the time–band product is c = 2πS² = 2πe^L for the e^{−2πixy} Fourier convention. That gives e^{−4πe^L}, i.e. **B → 4π ≈ 12.57**, consistent with the fitted 13.1 and the slope trend toward −12.
  - I have **not** derived λ_min(W) from Slepian's formula. The prefactor (b ≈ 7.9) is not explained, and I am not citing a theorem that gives this rate for the Weil form.
  - Bombieri studied Weil's functional on short supports ("Remarks on Weil's quadratic functional…", around 2000). I did not verify whether it states a decay rate, so I am not citing it for one.
- **Relabel 14E's "Gaussian decay in L (CC very small eigenvalues; qualitatively known)"** as: "super-exponential; on 0.25 ≤ L ≤ 2.0 the data fit exp(−B e^L + O(L)) with B ≈ 13, consistent with prolate decay at scale e^{L/2} (heuristic). The Gaussian fit is valid only on its 5-point range and fails at L ≥ 1.75."

### 6.6 (5) The bilinear radical check, redone
- **Construction (checked):** residual = ‖WS‖₂/‖W‖ = sup_{v∈S, φ∈V_M} |W(v,φ)|/(‖v‖‖φ‖‖W‖). The control is the same family pushed through F − 0.1·2^{−s}. That is exactly the test I asked for (3.2.1–2).
- **Numbers.** The JSON matches the table: ζ L = 4 and 6 at 6.4e3 and 1.2e4, EPS 1.3e2, F_θ 4.7e3, DH 7.4e2 at ε = 1e-8, and I reproduced F_θ L = 4's own residual of 6.7e-6. It drops to 47–1.9e3 at ε = 1e-6, as reported. **This is now a discriminating CHECK. Label earned.**
- **Minor points.**
  - (i) EPS's separation (1.3e2, with dim S = 9) is the weakest. The longer E tails at conductor 20 are the likely reason; 14E notes this.
  - (ii) The residual is measured inside V_M (the Galerkin span). Radical membership against probes *outside* the band, i.e. high-frequency φ near the window edge, is untested. That matters because the leftover near-null directions in S^⊥ are band-edge and window-edge (14E §4).
  - (iii) As 14E says, it proves the easy direction numerically and does not single out ζ.

### 6.7 Smaller items
- "Non-near-null count ≈ 2N(fmax) + 10–14 is Connes' semiclassical count (known)": the count depends on the threshold (1e-8‖W‖) and on leakage from the non-band-limited sine basis. Label it INTERPRETATION, not "known".
- The fakes' ghost-depth fits e^{(0.8–1.1)L} use 4 points each. Fine as NUMERICAL; don't extrapolate them.
- §7 "Next step" (Connes–Consani Sonin compression with S complete for the support) is the right next target, and it matches §2.3.

### 6.8 Falsifier and next check
- **Biggest hole.** At double precision, the quotient question for ζ (λ_min(W|S^⊥) at fixed L) cannot be settled. W's own spectrum below 1e-15 is a double-exponential ladder, and S is built in double precision.
- **Cheapest falsifier for "λ_min(quotient) = 0 at every fixed L".** Rebuild S in mpmath at L = 1.5–2. The line values come from Hurwitz sums, which are cheap at 40 digits for t ≤ 250. Then compute λ_min(W|S^⊥) at 80 digits for two cuts and two bases.
  - If λ_min(W|S^⊥) converges to a positive, cut-independent value well above λ_min(W) (about 1e-29 at L = 2), there *is* a finite-L gap on the quotient, and the Lehmer bound then caps it.
  - If it tracks λ_min(W) down the ladder, the uncaught near-null directions are genuine and the quotient has no gap at finite L.
- **Re-runs, ~20 core-minutes total:**
  - my high-precision Galerkin `mpw.py`, a line-by-line mpmath port of `weil_galerkin.operators` for ζ;
  - the Lehmer reproduction with the one-line fix;
  - the F_θ L = 4 Sylvester rerun plus the random-subspace control.
  - Scratch copies are in `/tmp/ge/` on the shared box, not committed. They are short; ask if you want them committed.

---

## 7. Round 14E2 (main 24e17c3): the true quotient at 80 digits, plus the 14E/14D errata

**Verdict.** It reproduces, and the conclusion stands: at L = 1.5–2 the quotient seminorm on S^⊥ tracks λ₁(W_ζ) (Q/λ₁ ∈ [1.09, 3.57] over the 41 non-empty ζ runs), so there is no finite-L gap. The Schur complement is the right object. But the **principal-angle diagnostic is the wrong quantity, and the "S nearly orthogonal to u₁" explanation is false**: S carries 8–72% of u₁'s energy. Q is close to λ₁ because S does not *contain* u₁, not because S misses it. The F_θ "7 → 7" is the restriction count; for the true quotient the count is fixed by a theorem (Haynsworth), so it isn't evidence either way. The errata match §6.

### 7.1 Is the Schur complement the right quotient object? Yes.
- For x ∈ S^⊥ (coordinates in an orthonormal basis Q of S^⊥), q([x]) = inf_{s∈S} W(x + s). When W_SS > 0 the infimum is attained at s = −W_SS^{−1}W_{SQ}x, and q = xᵀ(W_QQ − W_QS W_SS^{−1} W_SQ)x = xᵀ(W/W_SS)x. The quotient norm of the coset is ‖x‖, since x is its minimal-norm representative. So λ_min(W/W_SS) is the bottom of the quotient form for the quotient Hilbert norm. `mpgap.py` builds exactly this (`Sc = R − WQS·WSS⁻¹·WQSᵀ`, after checking that W_SS has no zero eigenvalue). Correct.
- Two caveats on what it means:
  - S is the **window-truncated** range of E. It is not in the radical (bilinear residual 4e-7 to 2e-3). So W/S is "W with the S-directions optimised away", not the 14D/Connes quotient by a radical. On the Galerkin space W_ζ > 0 has trivial radical, so a quotient by a true radical would be the identity here.
  - For W > 0 the bracket is exact and elementary: (W/W_SS)⁻¹ = (W⁻¹)_QQ. That gives λ_min(W/W_SS) ≥ λ₁(W). Testing the coset of u₁ gives λ_min(W/W_SS) ≤ λ₁/(1 − ‖P_S u₁‖²). The upper bracket λ_{k+1}(W) quoted in 14E2 also holds (Sc ≤ R, then interlacing).
- **The Gaussian H switch is legitimate.** E_F H is in the radical for Schwartz H whose Mellin transform vanishes at F's poles (the constraint rows Cn), and Gaussians are Schwartz. Compact support is not needed for that. The Gram matrices come from Parseval on the line, t ∈ [0, 260], step 0.1, 90 digits, with ζ(½+it) from mpmath. No zeros are used. Truncation at |t| = 260 costs ≲ e^{−400} for these σ = w/4 (σ ≥ 0.18, modulation k ≤ 99). The price of the switch is a longer tail, which is truncated by the window projection. That leak is exactly what the bilinear residual measures, and it is of the same size as in 14E.

### 7.2 Spot-check reruns (80/90 digits, scratch copy of `mpW.py` / `mpline.py` / `mpgap.py`, about 4 min)
- **W reproduces to all printed digits:** λ₁ = 6.394341e-15 (L = 1.5, M = 30) and 1.2843737e-27 (L = 2, M = 37). Every λ_min_quot in three gap files (ZETA 1.5/45/1e-13, 2.0/45/1e-13, 2.0/45/1e-40) matches the committed JSON to 8 digits. The 41-run range Q/λ₁ = [1.0922, 3.5666] is recomputed from the committed JSONs.
- **The missing diagnostic.** I added ‖P_S u₁‖² (u₁ = bottom eigenvector of W) to the run:

| L, τ | ε | dim S | ‖P_S u₁‖² | λ₁/(1 − ‖P_S u₁‖²) | Q = λ_min(W/W_SS) | sin(largest angle) |
|---|---|---|---|---|---|---|
| 1.5, 1e-13 | 1e-6 | 1 | 0.1417 | 7.4502e-15 | 7.4503e-15 | 0.926 |
| 1.5, 1e-13 | 1e-4 | 3 | 0.1741 | 7.742e-15 | 7.742e-15 | 0.999 |
| 2.0, 1e-13 | 1e-10 | 3 | 0.1626 | 1.5338e-27 | 1.5338e-27 | 0.987 |
| 2.0, 1e-13 | 1e-6 | 6 | 0.2317 | 1.6718e-27 | 1.6717e-27 | 0.997 |
| 2.0, 1e-13 | 1e-4 | 8 | 0.5512 | 2.862e-27 | 2.861e-27 | 0.996 |
| 2.0, 1e-40 | 1e-6 | 7 | 0.5776 | 3.041e-27 | 3.040e-27 | 0.999 |
| 2.0, 1e-40 | 1e-4 | 8 | 0.6281 | 3.453e-27 | 3.453e-27 | 0.995 |

- So **Q = λ₁/(1 − ‖P_S u₁‖²)** to 3–5 digits. That is the test-coset upper bound, saturated because λ₂/λ₁ ≈ 10³ (the correction is O(λ₁/λ₂)). Inverting it across all 41 runs gives ‖P_S u₁‖² ∈ [0.08, 0.72].
- **Why the angle diagnostic misleads.** "sin of the largest principal angle between S and the bottom-k eigenvectors ≈ 1" only says that *one* direction of S is nearly orthogonal to span(u₁…u_k). It says nothing about how much of u₁ lies in S. The same JSON has sin(largest angle vs bottom-2k) = 0.33 at L = 2: S sits largely *inside* the near-null cluster. The text "removing S removes none of the ladder" is therefore wrong. Correct statement: S overlaps u₁ substantially but never contains it. A finite-L gap of size G would need ‖P_S u₁‖² ≥ 1 − λ₁/G, i.e. S containing u₁ to relative accuracy λ₁/G, and the data are nowhere near that.
- The conclusion "no finite-L gap at L = 1.5–2" does not depend on the explanation; it is the measured Q/λ₁. **NUMERICAL, earned.** Basis caveat as in §6.3: the fmax-45 λ₁ is a Ritz value up to 127–160× too high at L = 2, and Q follows whichever λ₁ the basis gives.

### 7.3 F_θ control
- The restricted form R = QᵀWQ: 7 → 7 (7 → 6 at dim S = 3), consistent with interlacing, low power (§6.2).
- **The quotient count is a theorem, not a measurement.** With W_SS > 0, Haynsworth inertia additivity gives In(W) = In(W_SS) + In(W/W_SS). So n₋(W/W_SS) = n₋(W) **exactly, whatever S is**. The JSONs report n_neg_quot = 8 against n_neg_W = 7. That is not a contradiction; it is the threshold. W_Fθ has an 8th eigenvalue at −4.6e-7 (just under the −1e-7·‖W‖ cut) and a 9th at −3.2e-12, and the Schur complement moves them across the cut. Report the quotient's ghost persistence as PROVED (Haynsworth), and "more negative" (−0.7 to −2.6) as expected, since S is not in F_θ's radical. Neither is evidence about F_θ.

### 7.4 Errata vs §6
- 14E block items 1–5 match §6.2–6.7: the fit is withdrawn and e^{−994} dropped; the replacement is a + bL − Be^L with B = 13.13, and B → 4π is labelled heuristic; floor = precision; the Lehmer window-support hypothesis is added; zeros are proof inputs; the NameError is fixed and reruns bit-for-bit (0.1203 / 0.0629); Sylvester → "one-sided, low power", with the random-subspace control quoted; the semiclassical count → INTERPRETATION.
- The 14D items 7–9 (3c measure → INTERPRETATION, N_F domain, E tail) are the three minor items I flagged in §6.1, and they are worded correctly.
- **One wording to fix in 14E2 itself** (ROUND14E2.md line 11 and §6 item 3): replace "S is nearly orthogonal to the near-null eigenvectors" with "S contains 8–72% of u₁ but never all of it; Q = λ₁/(1 − ‖P_S u₁‖²)". Also relabel the F_θ quotient persistence as PROVED (Haynsworth).

---

## Appendix: reproduction probes (run from the stated directory; ~15 s total)
Run `radprobe.py` from `round14d_global/`, and `eps_atom.py` and `spec.py` from `round14c_prolate/`. They reuse `weil_galerkin.py`, `models14b.py` and `radical.py` unchanged.

```python
# radprobe.py (cwd round14d_global): quadratic vs bilinear radical test for zeta
import sys, os, math, numpy as np
sys.path[:0] = [os.getcwd(), '../round14b_euler', '..']
import radical as R, models14b as md
u = R.u
def bump(v, a, b):
    y = (2*v-a-b)/(b-a); o = np.zeros_like(y); m = np.abs(y) < 1; o[m] = np.exp(-1/(1-y[m]**2)); return o
def mkH(a, b, exps):  # Mellin of H vanishes at s in exps
    vv = np.linspace(a, b, 20001); bb = bump(vv, a, b); k = len(exps)
    A = np.array([[np.trapezoid(bb*vv**j*np.exp(s*vv), vv) for j in range(1, k+1)] for s in exps])
    c = np.linalg.solve(A, -np.array([np.trapezoid(bb*np.exp(s*vv), vv) for s in exps]))
    return lambda x: bump(np.log(x), a, b)*(1+sum(cj*np.log(x)**(j+1) for j, cj in enumerate(c)))
def Eu(H, n, a, lo_, hi_):
    x = np.exp(u); G = np.zeros_like(u)
    for k, xx in enumerate(x):
        lo, hi = np.searchsorted(n, math.exp(lo_)/xx), np.searchsorted(n, math.exp(hi_)/xx, side='right')
        if hi > lo: G[k] = np.sum(a[lo:hi]*H(n[lo:hi]*xx))
    return G*np.exp(u/2)
Z = md.zeta(); n, a1 = R.coeffs('ZETA'); Wq = lambda g: R.W_grid(Z, g)
Wbil = lambda g, h: 0.5*(Wq(g+h)['W']-Wq(g)['W']-Wq(h)['W'])
probe = np.where(np.abs(u) < 1, np.cos(math.pi*u/2)**2*np.cos(14.1347*u), 0.0)
for A, B in [(0.5, 2.5), (0.5, 1.0)]:
    H = mkH(A, B, [1.0]); ap = a1.copy(); ap[1] = 0.9      # a(2) = 0.9: same band, not in range E_zeta
    for lab, v in [('E_zeta H', Eu(H, n, a1, A, B)), ('E_pert H', Eu(H, n, ap, A, B)), ('H', np.exp(u/2)*H(np.exp(u)))]:
        w = Wq(v); print(A, B, lab, w['W']/w['norm2'], Wbil(v, probe)/math.sqrt(w['norm2']*Wq(probe)['norm2']))
```

```python
# eps_atom.py (cwd round14c_prolate): +-eps single-atom perturbation of zeta's full W
import sys, math, numpy as np; sys.path[:0] = ['../round14b_euler', '..']
import weil_galerkin as wg, models14b as md
Z = md.zeta()
for L in [2.0, 3.0]:
    M = int(math.ceil(45*L/math.pi))+8; a, p, q = wg.operators(Z, L, M); W = a+p-q
    P = wg.Aclosed(L, M, math.log(2)); P = P+P.T
    for e in [1e-3, 1e-6, 1e-9]: print(L, e, np.linalg.eigvalsh(W+e*P)[0], np.linalg.eigvalsh(W-e*P)[0])
```

```python
# spec.py (cwd round14c_prolate): low spectrum of zeta's W, basis dependence of the near-null count
import sys, math, numpy as np; sys.path[:0] = ['../round14b_euler', '..']
import weil_galerkin as wg, models14b as md
Z = md.zeta()
for L in [1., 2., 3., 4., 5.]:
    for f in [45., 90.]:
        M = int(math.ceil(f*L/math.pi))+8; a, p, q = wg.operators(Z, L, M); ev = np.linalg.eigvalsh(a+p-q)
        k = int(np.sum(ev < 1e-8*np.abs(ev).max())); print(L, f, M, k, ev[k:k+4])
```

— Golden Eye (auditor), for Daniel Toupin / Golden Physics Project.
