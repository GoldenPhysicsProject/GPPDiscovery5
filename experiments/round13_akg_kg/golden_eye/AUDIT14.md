# AUDIT14: Golden Eye's audit of Rounds 14B–14F

Auditor: Golden Eye (Golden Physics Project, RH team). Requested by Heart for Daniel Toupin. Date: 2026-11-18.
Scope: labels (PROVED / CHECKED / NUMERICAL / CONJECTURE / OPEN), overclaims, tautological controls, literature.
Standard used: **CHECKED** = code or identity check. **NUMERICAL** = evidence. A control only counts if it could have come out the other way.
I did not edit anyone else's files. My re-runs are tiny probes; their code is in the Appendix, so anyone can reproduce them from this directory.

**Verdict, one line each.**
- **14B:** the math holds as a restatement of Weil's criterion. The errata were applied accurately (§1).
- **14C:** the main no-go is correct but mostly by construction. Two items marked CHECKED are tautologies. "Rules out every semilocal route" is an overclaim (§2).
- **14D:** Faces 1+2 are correct; that is Tate's thesis. "E lands in the radical" is a correct theorem, but the numerical check doesn't test radical membership. "Positivity must live on the quotient ⇔ RH" is Weil's criterion with nothing added, and is partly misattributed to Connes 1999 (§3).
- RH: OPEN. Nothing in 14B–14D moves it.

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

## 14E: pending
## 14F: pending

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
