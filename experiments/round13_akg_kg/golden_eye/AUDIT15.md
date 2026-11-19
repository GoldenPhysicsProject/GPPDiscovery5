# AUDIT15: Golden Eye's audit of Round 15 (Hasse Gram for ℤ)

Auditor: Golden Eye (Golden Physics Project, RH team). Requested by Heart for Daniel Toupin. Date: 2026-11-18.
Object: `experiments/round13_akg_kg/round15_hasse/ROUND15.md` at main 4dbeb9d, with its scripts and JSON.
Standard (as in AUDIT14): **CHECKED** = code or identity check. **NUMERICAL** = evidence. A control only counts if it could have come out the other way.
I did not edit anyone else's files. My reruns used a scratch copy of the round's own scripts (§6).

## 0. Verdict

- **Part 1 (F_q side) is correct and standard.** I re-derived Γ′_j·Γ′_k = −q^j a_{k−j} and Γ′_j² = −2g q^j. T[j,k] = Σω_i^{k−j} ⪰ 0 is Hodge index, i.e. the Grothendieck / Mattuck–Tate route. The Toeplitz form equals 2·deg(m − nF) as an algebraic identity: the "verification" compares the form to m² − a·mn + q·n², not to an isogeny degree, so it is CHECKED-as-algebra, not a test. The synthetic genus-2 fake passing at 2×2 and failing at size 4 = 2g reproduces, and it is forced by a Vandermonde/Sylvester argument (§4).
- **The one-shift Hasse block is Cauchy–Schwarz for W.** W is translation-invariant, so W(τ_s g) = W(g). Then W ≥ 0 on span{g, τ_s g} ⇔ |W(g, τ_s g)| ≤ W(g). So:
  - it is **implied by positivity**;
  - every failure is an explicit vector g ± τ_s g with W < 0, i.e. a ghost at width ≤ l + s;
  - it is **strictly weaker** than W ≥ 0. EPS's W is indefinite from width 3.0, but the block passes until s = 4.40 (width 5.0).
  
  As a ζ-vs-fake separator it is a **tautological control**: a weaker ghost detector than 14B. ROUND15 itself says this in Part 3, which is to its credit; the Summary item 2 and §2.2 readout should say it too.
- **"Holding for all s ⇔ RH" (PROVED, sketch) is earned.** ⇒ is Cauchy–Schwarz under RH. ⇐ is the Laplace-transform argument, which is correct as sketched. So is the subexponential version.
- **The pole identity is right only for even g.** The general form is 2cosh(s/2)·ĝ(i/2)ĝ(−i/2), with ĝ(z) = ∫g e^{izu}, and that product is not |ĝ(i/2)|². CHECKED here: for sine₁ + sine₂, pole/2cosh = 0.24284 = ĝ(i/2)ĝ(−i/2), while ĝ(i/2)² = 0.28219. It is also an algebraic identity (c_±(τ_s g) = e^{±s/2}c_±(g)), so "measured at every s" is a code check.
- **ζ passing is a code check. F_{iτ} passing is expected.** Margins are at least 6 orders above the relevant numerical floor, so the passes and fails are resolved; they are not threshold artifacts. The l = 0.6 vs 0.3 dependence shows a finite scan is only a necessary condition.
- **Two literature/label fixes:**
  - The CC 2021 support is g ⊂ [2^{−1/2}, 2^{1/2}] (f = g*g* ⊂ (1/2, 2)). It also needs ĝ(i/2) = 0, so it does not cover the diagonal's pole part.
  - "Kronecker makes ρ → 1" should be relabelled HEURISTIC. F_{iτ} nearly saturates too (0.9963).
- RH: OPEN. Round 15 does not move it, and says so.

## 1. Claim 1: the one-shift block |W(g, τ_s g)| ≤ W(g)

### 1.1 What it is
- `hasse_blocks.py` tests W_l ± Bsym_s ⪰ 0 on the M-dimensional window space, with s ≥ l so that g and τ_s g are disjoint and ‖g ± τ_s g‖² = 2‖g‖². Since W(g ± τ_s g) = 2W(g) ± 2W(g, τ_s g), this is exactly:
  - "W ≥ 0 on g + τ_s g and on g − τ_s g for every window g", or equivalently
  - "the 2×2 Gram [[d, b], [b, d]] is PSD for every g".
- **Translation invariance is real.** W depends on g only through its autocorrelation, which τ_s does not change. The prime data cutoff is respected: s_max = log n_max − l − 0.05.
- **Implied by W ≥ 0.** If W ≥ 0 on functions of width l + s, Cauchy–Schwarz gives |W(g, τ_s g)| ≤ √(W(g)·W(τ_s g)) = W(g).
- **Tautological as a control.** A failure exhibits v = g ± τ_s g with W(v) < 0, an explicit ghost of width ≤ l + s. So a fake can fail only where its W is already indefinite.
- **Strictly weaker than W ≥ 0** (a test on a union of 2-dim subspaces, not even on all of V_l ⊕ τ_s V_l, since Bsym drops the antisymmetric part of B_s). The onset widths l + s versus 14B's ghost onsets:

| model | block onset s (l = 0.6) | width l + s | 14B ghost onset L | Gram-INT onset |
|---|---|---|---|---|
| FTH | 0.66 | 1.26 | 1.0 | 1.29 |
| BEUR″ | 3.60 | 4.20 | 4.5 (grid 4.0 → 4.5) | 4.10 |
| EPS | 4.40 | 5.00 | 3.0 | 3.08 |
| DH | 6.30 | 6.90 | 4.0 | 3.64 |

  At every block onset, W is already known to be indefinite at that width, as the logic requires: the 14B or Gram onset is ≤ l + s in every row. BEUR″'s 4.20 < 4.5 is consistent: 14B's grid only brackets the onset in (4.0, 4.5], and the Gram finds it at 4.10. FTH's Gram-INT onset (1.29) lies above 1.26 only because the INT family is a different subspace; 14B's onset is 1.0.

### 1.2 Reproduction (scratch copy, unchanged code)
- `hasse_blocks.py ZETA 0.6` reproduces exactly: 499 shifts, 0 violations, min μ = 2.1233e-5, max ρ = 0.999587.
- `hasse_blocks.py FTH 0.6` reproduces exactly: first violation s = 0.66, 377 violations, min μ = −4.4868, max ρ = 2.5384.

### 1.3 Numerical precision: margins vs floor
- **Floor.** B_0 agrees with 14B's W to ≤ 1.4e-13 (`check_core`). Panel doubling moves μ by < 7e-14 (`check_quad`). The independent grid check agrees to 1e-6–1e-7, which is that grid's own accuracy. The ghost threshold is −1e-7·max|W|.
- **Passes.** ζ: min μ = 2.1e-5 at l = 0.6 and 1.2e-2 at l = 0.3. F_{iτ}: 6.2e-4 and 7.5e-2. ζ's are 8–11 orders above the 1e-13 floor; F_{iτ}'s are ≥ 6 orders above its own 3.4e-10 quadrature floor. These passes are resolved, unlike 14B's ζ "λ_min at the float floor".
- **Fails.** At depth they are huge: EPS −16.9, FTH −4.49, DH −2.11, BEUR″ −0.114. At each onset, μ crosses zero between consecutive scan points: FTH +4.7e-3 → −3.2e-4, BEUR″ +6.4e-4 → −1.5e-4, EPS +3.5e-2 → −0.14, DH +1.8e-2 → −0.48. The smallest is |−1.5e-4|, still 9 orders above the floor. So the onset s is resolved to the scan step (0.02, or the prime spacing), and the sign is not in doubt.
- **The near-misses are at l = 0.3:** DH max ρ = 0.99929 (min μ = 2.9e-3), comfortably resolved. But "DH passes at l = 0.3" means only that this window cannot see DH's off-line zeros before s ≈ 8.25. It does not mean DH looks like ζ.

### 1.4 Window dependence (l = 0.6 vs 0.3)
- Halving l halves the frequency resolution of g (M = 17 → 13; 32 → 20 for DH). The diagonal "genus" grows: λ_min(W_l) = 0.0089 → 0.229 for ζ. A narrower g spreads |ĝ|² over more zeros, so an off-line pair's e^{δs} term must grow longer before it beats d(g).
- Result: at l = 0.3 only EPS (from s = 5.64) and FTH (from 1.78) fail; DH and BEUR″ pass to s ≈ 8.25.
- This is the expected behaviour of a necessary condition scanned over finite s. It is not a property of the fakes. ROUND15's reading ("the onset moves out as the genus grows") is a fair heuristic. Label the l-dependence NUMERICAL.
- At l = 0.3, BEUR″'s min μ and max ρ equal ζ's to 5 digits: its extremes are set by the shared small-n data, and the Beurling splice never beats them. As noted in AUDIT14 §5, BEUR″'s "ghosts" measure the archimedean/splice defect, not zeros.

### 1.5 ζ passing is a code check; F_{iτ} passing is expected
- **ζ.** The block needs W_ζ ≥ 0 on functions of width ≤ l + s ≤ 8.6 in the fmax-45 band. Their |ĝ|² weight sits on zeros with |γ| ≲ 50 plus a tail that decays polynomially, and those zeros are verified on the line to astronomically larger heights. So a ζ failure would mean a bug, and a pass cannot count as evidence. Same status as 14B's ζ ghost count (AUDIT14 §1), except that here the margin is resolved.
- **F_{iτ} = ζ(s−iτ)ζ(s+iτ).** Its zeros are ζ's, shifted by ±iτ along the line. Its form is W_ζ(e^{iτu}g) + W_ζ(e^{−iτu}g), the unitary twist (AUDIT14 §5.1). So W_{F_{iτ}} ≥ 0 whenever W_ζ ≥ 0 on complex g, and the pass follows from ζ's. That makes it a second code check: of the FTAU pole kernel 4cosh(t/2)cos(τt) and the twisted arch weight. It is not an independent control.
- **Saturation.** "ζ saturates, max ρ = 0.9996 ... Kronecker makes it approach 1" should be HEURISTIC:
  - Under RH, ρ(s) = sup_g |Σ_{γ>0} |ĝ(γ)|² cos(γs)| / Σ_{γ>0} |ĝ(γ)|². It comes near 1 whenever a window g can put its spectral weight on zeros with |cos γs| ≈ 1. That is limited by M and l (ρ_max drops to 0.956 at l = 0.3).
  - The unitary control does the same: 0.9963.
  - For "sup over s → 1", Dirichlet's simultaneous approximation is enough; Kronecker/LI is not needed.
  - The de Bruijn–Newman / "edge" reading is INTERPRETATION.

## 2. Claim 2: "holding for all s ⇔ RH" (PROVED, sketch): earned
- **⇒.** RH ⇒ W ≥ 0 ⇒ Cauchy–Schwarz (§1.1).
- **⇐.** With b_s(g) = Σ_ρ ĝ(γ_ρ)ĝ(−γ_ρ)e^{−iγ_ρ s}, an off-line zero ρ = ½ + δ + it has partner 1 − ρ̄ = ½ − δ + it, with γ = t + iδ and growth e^{δs}.
  - For compactly supported g, ĝ is entire and decays faster than any power in horizontal strips. So the sum converges absolutely, and the Laplace transform ∫_l^∞ b_s e^{−zs} ds = Σ c_ρ e^{…}/(z + iγ_ρ) is meromorphic with a pole at z = δ − it whenever c_ρ = ĝ(γ)ĝ(−γ) ≠ 0.
  - That holds for generic g in the window, since ĝ is entire and not identically zero. Distinct zeros give distinct poles, so there is no cancellation.
  - A bounded b_s makes the transform analytic on Re z > 0. Contradiction.
- **Subexponential version.** The same argument shows a bound C(g)e^{εs} for all ε > 0 also gives RH. This is the classical equivalence of RH with smoothed ψ(x) = x + O(x^{1/2+ε}).
- **Precision of the statement.** The ⇐ direction needs one g with ĝ(γ)ĝ(−γ) ≠ 0 at an off-line zero, and all s → ∞, not a finite scan. ROUND15 states both. "Finite checks prove nothing" is correct.

## 3. Claim 3: pole term = 2cosh(s/2)·|ĝ(i/2)|²
- **Derivation.** Write c_± = ∫g(u)e^{±u/2}du, so ĝ(−i/2) = c_+ and ĝ(i/2) = c_−. The pole part of W is P(g) = 2c_+c_−. Since c_±(τ_s g) = e^{±s/2}c_±(g), polarisation gives P(g, τ_s g) = ½[P(g + τ_s g) − 2P(g)] = **2cosh(s/2)·c_+c_− = 2cosh(s/2)·ĝ(i/2)ĝ(−i/2)**.
- For real g both factors are real, and their product is **not** |ĝ(i/2)|² unless c_+ = c_−, e.g. g even about the window centre. ROUND15's test function (first sine, which is even) hides this.
- **Rerun** (`hasse_core.Bs`, l = 0.6, s = 1, 3, 7):

| g | pole / 2cosh(s/2) | ĝ(i/2)ĝ(−i/2) | ĝ(i/2)² |
|---|---|---|---|
| sine₁ (even) | 0.488418389 | 0.488418387 | 0.488418387 |
| (sine₁ + sine₂)/√2 | 0.242837331 | 0.242837330 | 0.282188251 |

- **Labels.**
  - It is an identity of the pole kernel, so "CHECKED to 6 digits, constant at every s" is a code check, not a measurement.
  - The dictionary entry "ℤ image of (q + 1)/√q = √q + 1/√q" is fine as INTERPRETATION, with q = e^s.
  - Because c_+c_− can be negative for sign-changing g, the pole form is rank 2 and indefinite (signature (1, 1)). See §5 for why that matters to the CC citation.

## 4. Claim 4: the finite-field side
- **Intersection numbers (PROVED, re-derived).** With H = C×pt, V = pt×C, Γ_j·H = q^j, Γ_j·V = 1, Γ_j·Γ_k = q^j N_{k−j} (j < k) and Γ_j² = q^j(2 − 2g):
  - Γ′_j = Γ_j − H − q^jV is orthogonal to H and V;
  - Γ′_j·Γ′_k = q^j N_{k−j} − q^j − q^k = −q^j a_{k−j};
  - Γ′_j² = −2g q^j.
  
  Hodge index (the primitive part is negative definite) gives T[j,k] = Σω_i^{k−j} ⪰ 0 for every size. This is Grothendieck's 1958 proof and Mattuck–Tate's (1958, Castelnuovo–Severi from Riemann–Roch on C×C). Weil's own route was Castelnuovo positivity / positivity of the Rosati involution on the Jacobian. ROUND15's account matches.
- **"Toeplitz form = 2·deg(m − nF)".** For g = 1, (m, −n√q)·[[2, a/√q], [a/√q, 2]]·(m, −n√q)ᵀ = 2(m² − a·mn + q·n²) = 2·deg(m − nF), using deg(m − nF) = (m − nF)(m − nF̂). Correct. But `ff_hasse.py` checks it by comparing against m² − a·mn + q·n² for |m|, |n| ≤ 4. No isogeny degree is computed independently, so "verified" is an algebraic identity (CHECKED-as-algebra). The content is Hasse's theorem: |a| ≤ 2√q with a from point counts. The elliptic rows confirm it on examples.
- **Size 2g decides (PROVED, elementary).** Write T_n = V_n* C V_n, where V_n is the n×2g Vandermonde in the ω_i and C pairs ω with 1/ω̄ (the identity on unitary ω; a (1, 1) block on an off-circle pair). Once n ≥ 2g (distinct ω), V_n has full column rank, so by Sylvester In(T_n) = In(C) ⊕ zeros. Hence T_{2g} ⪰ 0 ⇔ all |ω_i| = 1. Repeated ω need a confluent Vandermonde; same conclusion. Ghosts can appear before size 2g but cannot be postponed past it.
- **Synthetic genus-2 fake (reproduced, `ff_hasse.py`).** {e^{±iφ}, r, 1/r} with (φ, r) = (2.2, 2), (2.6, 1.5), (1.0, 1.2) gives c₁ = 1.323, 0.453, 3.114, all below 2g = 4. So the 2×2 block [[4, c₁], [c₁, 4]] passes, and the first ghost is at K = 3 (size 4 = 2g) every time. The one-shift blocks first fail at k = 3, 4, 6, when |r^k + r^{−k} + 2cos kφ| > 4. As predicted.
  - Caveat: the fake has real but non-integer a_k, so it is not a curve. It tests the Toeplitz algebra, not geometry.
  - The commit message's "ghost at K = 2g − 1" means index K = 3, i.e. size 2g. Fine, but say "size 2g".
- **Real genus-2 curves** (p = 5, 7, 11, 13; brute-force N₁, N₂): |ω| = 1 to print precision, 0 ghosts to K = 10. That is Weil's theorem on examples. CHECKED.
- **Where the analogy holds or breaks.** Over F_q, positivity is inherited from Hodge index on a surface with an ample class, and the 2g-block decides because H¹ is finite-dimensional. Over ℤ, neither is available. ROUND15 §1.3 and Part 3 §4 say exactly this. Agreed.

## 5. Other items
- **CC 2021 citation (§1.2 table and Part 3 item 1).** Connes–Consani, Selecta 27 (2021), Thm 1 and eq. (5), and Weil/Yoshida positivity as they cite it, give W_∞(g*g*) ≥ 0 for supp g ⊂ [2^{−1/2}, 2^{1/2}], i.e. log-width log 2, matching l < log 2. They require **ĝ(i/2) = 0**; Thm 1 also requires ĝ(0) = 0, which (5) removes at the cost of c|ĝ(0)|². ROUND15 writes "supp g ⊂ (1/2, 2)", but (1/2, 2) is the support of f = g*g*.
  - Fix: on the subspace ĝ(i/2) = 0, Round 15's diagonal d = arch + pole reduces to W_∞, and CC applies.
  - On the full window, the pole part is a rank-2 indefinite form. So "the diagonal is known positive by CC 2021" is literally true only on that codimension-1 subspace. On the full window it is NUMERICAL here (λ_min(W_l) = 0.0089 > 0 at l = 0.6).
- **Full Hasse Gram (§2.3).** "Just W on a lattice of translates; reproduces 14B onsets; nothing new" is correct. So is "both are lower bounds on the inertia" (restriction, AUDIT14 §6.2).
  - ζ's λ_min at −1e-14 to −8e-14 is the float floor. Its 0 ghosts are a code check, threshold-limited as in 14B.
  - F_{iτ}'s −3.4e-10 is its quadrature floor. Fine.
- **Part 1's "add + mult" / Stepanov discussion.** Correctly labelled HEURISTIC / OPEN. Bombieri–Stepanov's use of x ↦ x^q being additive and multiplicative is accurately described. Bombieri's Sém. Bourbaki 430 (1973) is the right reference. "No published Stepanov argument on the arithmetic site" is labelled as a quick search, which is appropriate.
- **"N_k ≈ q^k + 1 with error bounded by the genus" (b_s remainder ±0.008 vs d = 0.0118).** This is the numerical face of the Cauchy–Schwarz bound in §1.1 for one fixed g. It is NUMERICAL, not an additional CHECK.

## 6. Biggest hole and cheapest falsifier for Round 16
- **Hole.** Round 16 proposes a Stepanov-type uniform bound on W(g, τ_s g) without zeros. By §2, any bound C(g)e^{o(s)} *is* RH. So the "finite toy" must locate the step where additivity enters, and show that step has a characteristic-1 analogue that acts on the prime-window sum in b_s. Otherwise it reproduces RH as a hypothesis.
- **Cheapest falsifier.** Run the Round 16 toy unchanged on the F_q "Epstein" Z₀ (probe B: functional equation, no Euler product, 1 ghost). A Stepanov-type count that does not use the multiplicativity of Frobenius on points would also bound Z₀'s a_k, and Z₀ violates the bound. If the toy's count passes Z₀, the toy is not using the coupling. Also run it on F_{iτ} vs F_θ: unitarity is the axiom that separated them in 14F.

## 7. Reproduction (scratch copy inside the worktree, not committed)
- Copy `round15_hasse/` to a sibling scratch folder so the `../round14b_euler` imports resolve; running from `/tmp` fails on `ModuleNotFoundError: weil_galerkin`. Then:
  - `python3 hasse_blocks.py ZETA 0.6; python3 hasse_blocks.py FTH 0.6` (about 2 min each; exact match);
  - `python3 ff_hasse.py` (needs `../clues/ff_epstein.json`; exact match).
- Pole-identity probe (cwd = the scratch copy):

```python
import numpy as np, math, hasse_core as hc
mo = hc.MODELS['ZETA'](); l = 0.6; M = hc.Mof(l, 45.); x = np.linspace(-l/2, l/2, 20001)
phi = lambda j: math.sqrt(2/l)*np.sin((j+1)*math.pi*(x+l/2)/l)
for tag, g in [('sine1', np.eye(M)[0]), ('sine1+sine2', np.r_[1, 1, np.zeros(M-2)]/math.sqrt(2))]:
    gx = sum(g[j]*phi(j) for j in range(M) if g[j]); cp = np.trapezoid(gx*np.exp(x/2), x); cm = np.trapezoid(gx*np.exp(-x/2), x)
    for s in [1.0, 3.0, 7.0]:
        D = hc.Bs(mo, l, M, s); P = 0.5*(D['pole']+D['pole'].T)
        print(tag, s, (g@P@g)/(2*math.cosh(s/2)), cp*cm, cm*cm)   # pole/2cosh, g^(i/2)g^(-i/2), g^(i/2)^2
```

RH: OPEN.
