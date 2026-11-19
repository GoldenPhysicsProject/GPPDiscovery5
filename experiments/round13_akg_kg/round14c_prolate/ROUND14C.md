# Round 14C — the prolate operator as "T" in ‖Tf‖² = W: pair test first

Framing credit: **No-Ghosts framing — Daniel Toupin.** The task is to write the ghost pairing ⟨f, f⟩_Θ = W(f ∗ Θf) from Round 14B as ‖Tf‖², with T built from the Connes–Moscovici prolate operator, and to say exactly where that fails.

Code: `pairtest.py`, which reuses `../round14b_euler/weil_galerkin.py` and `models14b.py`. Data: `pairtest.json`, `pairtest.log`. Zeros are never inputs.

> **Errata (Golden Eye, AUDIT14 on team/golden-eye 062e4e2; applied in Round 14E).**
> 1. "Zero defect for L < log 2" is an **identity by support**: the code drops every atom with log n ≥ L, so for L < log 2 the prime matrix is empty by construction. It is a code check, not evidence.
> 2. "Defect = the n = 2 atom at log 2 + 0.05" is an **identity check of an assignment** (W := A_∞ − q, then tested W − A_∞ = −q). The arithmetic content, that only n = 2 lies in (log 2, log 3), is elementary. Code check.
> 3. The T(ζ) = T(BEUR″) and T(χ) = T(DH) "0.0 difference" is vacuous as a CHECK: the code copies the arch data. PROVED by construction only.
> 4. **Rename:** the model called "L(χ₅)" below is **L(χ), χ the odd quartic character mod 5** (χ(2) = i; Γ_ℝ(s+1), q = 5), the right partner for DH. It is **not** 13j's χ₅ (even quadratic, Γ_ℝ(s)).
> 5. **Withdrawn:** "this rules out every finite-S / semilocal-prolate route". ζ's own W has 8 / 22 eigenvalues below 1e-8·max at L = 2 / 3, so a ±1e-9 change of *any* atom makes it indefinite; the semilocal table is uninformative. Replacement: "rules out positivity of incomplete truncations W_S and T = W^{1/2}; does **not** address the Connes–Consani Sonin-compression route with S complete for the support."
> 6. "As Tate vectors the two differ only at p = 5" is INTERPRETATION.

## 0. Honest framing of 14B (and its errata)
- **"No ghosts ⇔ RH" is Weil's criterion restated. It is classical.** 14B's ghost counts are finite-L NUMERICAL evidence. The ζ "no ghosts" row is a CHECKED code check.
- Golden Eye's audit is applied as an errata block at the top of `../round14b_euler/ROUND14B.md`:
  - (a) ζ row relabelled as a code check;
  - (b) F_θ's ghosts are built in, and "Euler + FE insufficient" is standard;
  - (c) Beurling moved to "defect / not applicable";
  - (d) ghosts counted with multiplicity;
  - (e) Euclidean OS is impossible for any L-function with zeros (Bernstein–Widder);
  - (f) the Tate "same time" reading is INTERPRETATION, and at p the Fourier transform is the Weyl element.

## 1. The papers: what is proved and what is numerical
**Connes–Moscovici, "Prolate spheroidal operator and zeta"** (arXiv 2112.05500; PNAS 2022, "The UV prolate spectrum matches the zeros of zeta"). The operator is W_λ = −∂(λ² − x²)∂ + (2πλx)² on L²(ℝ). It is purely **archimedean**: there is no prime anywhere in it.
- PROVED:
  - The deficiency indices of W_min are (4, 4) (Lemma 2.1).
  - There is a self-adjoint extension W_sa that **commutes with the Fourier transform** and has discrete spectrum, unbounded in both directions (Thm 2.6).
  - Eigenfunctions for negative eigenvalues lie in the **Sonin space** (Lemma 3.1, Cor. 3.2).
  - The Darboux / Riccati construction of an isospectral family of Dirac operators D (§5).
- SEMICLASSICAL, not a rigorous Weyl law: the count of negative eigenvalues, 2σ(E, λ) with σ ∼ (E/2π)(log(E/2π) − 1) + const (Prop. 4.2, stated as "the semiclassical approximation"). Theorem 6.1's counting statement N(E) ∼ (E/2π)(log(E/2π) − 1) + O(1) for √2·D rests on that proposition.
- NUMERICAL: Figs. 2–4 compare the spectrum with zeta zeros. This is a UV / counting-level similarity, not a match of individual zeros. The authors state explicitly that they "do not claim" every eigenvalue on Sonin space is negative; that is supported only by numerics.
- **Structural point (PROVED, elementary).** The main term (E/2π)(log(E/2π) − 1) of N(E) comes entirely from the Γ_ℝ factor. It is the same for every function with ζ's gamma factor and conductor, so an archimedean operator cannot carry any information that separates ζ from a fake with the same Γ factor.

**Connes–Consani–Moscovici, "Zeta zeros and prolate wave operators"** (arXiv 2310.18423; Ann. Funct. Anal. 2024).
- PROVED:
  - Thm 1, the semilocal Hardy–Titchmarsh transform for a finite set S of places. It uses the measure dm_S = ∏_{v∈S} |L_v(½ − is)|² ds, and F_S ∘ η_S = η_S ∘ F_ℝ.
  - Thm 2, the semilocal Sonin spaces are hilbertian-isomorphic for every finite S ∋ ∞ (stability in S).
  - In the archimedean case, V W_λ V* = −s² + 2πλ²(4N + 1) − ¼, where N is the orthogonal-polynomial grading.
  - Thm 5.2, the metaplectic discrepancy coefficients.
- NOT done there:
  - The semilocal prolate operator itself is announced for a forthcoming paper. The metaplectic description is "at the formal level".
  - The low-lying ("infrared") zeros from the radical of the Weil quadratic form, through E(f)(u) = u^{1/2} Σ_n f(nu), are a **numerical** spectral realization (Connes–Consani, ζ-cycles).
- Weil positivity is **proved only at the archimedean place**. That is Connes–Consani, Selecta Math. 27 (2021): support in [2^{−1/2}, 2^{1/2}], i.e. width L ≤ log 2, where no prime enters.

## 2. What data each candidate T uses, and what the fakes lack (stated up front)

| candidate T | data used | ζ | EPS x²+5y² | DH | F_θ | spliced Beurling BEUR″ |
|---|---|---|---|---|---|---|
| T_∞: ‖T_∞g‖² := (arch + pole)(g). This is the archimedean Weil form; prolate / Sonin is the CC-2021 proof of its positivity | Γ factor, conductor, poles | Γ_ℝ(s), q = 1 | **different**: Γ(s), q = 20 | Γ_ℝ(s+1), q = 5, **identical to L(χ mod 5)** | Γ_ℝ(s±θ) | **identical to ζ** (borrowed) |
| T_S: ‖T_Sg‖² := (arch + pole − prime_S)(g). Semilocal: Euler factors at p ∈ S | the above + L_p for p ∈ S | has them | has **no** Euler factors | has **no** Euler factors; its Tate vector differs from L(χ)'s only at 5 | factors (1 − p^{θ−s})^{-1}(1 − p^{−θ−s})^{-1} | **identical to ζ** for p < 55 |

## 3. The pair test (Principal Series's cheapest falsifier). Same g for both members, g = the fake's ghost eigenvector (NUMERICAL, `pairtest.json`)

| pair | L, M | T_∞ identical across the pair? | W(member 1) | W(fake) | ‖T_∞g‖² (1 / fake) | what changes W |
|---|---|---|---|---|---|---|
| ζ / BEUR″ | 4.5, 73 | **yes** (exactly) | +0.042 | **−0.587** | 1.528 / 1.528 | prime atoms ≥ 55 only: −0.629 |
| L(χ₅) / DH | 4, 168 | **yes** (exactly) | +2.773 | **−0.755** | 4.183 / 4.183 | prime term only: −3.528 |
| L(χ₅) / DH | 6, 247 | yes | +3.991 | **−3.675** | 4.222 / 4.222 | prime term only: −7.666 |
| ζ / EPS | 4, 66 | no (different Γ) | +2.230 | **−1.615** | 0.993 / **4.811** | arch +3.818 (*more* positive for the fake), prime −7.663 |
| ζ / EPS | 6, 94 | no | +0.216 | −8.643 | 1.064 / 4.949 | arch +3.885, prime −12.744 |
| ζ / F_0.2 | 2, 37 | no | +0.0004 | −0.091 | 1.382 / 2.764 | arch +1.381, prime −1.474 |

- **The pairs with identical archimedean data fail at once** (PROVED by construction, CHECKED in numbers). For ζ / BEUR″ and L(χ₅) / DH, T_∞ is the same operator, so ‖T_∞g‖² is the same number. W has opposite signs. So no T built from archimedean data can satisfy ‖Tf‖² = W for both members.
- **Where the archimedean data differ, they differ in the wrong direction.** For ζ / EPS, EPS's archimedean part is *larger* on its own ghost (4.81 vs 0.99). The ghost is created entirely by the prime term, as in 13e. The same holds for F_θ.
- **L(χ₅) / DH, as Principal Series asked.** As Tate vectors the two differ only at p = 5. I split the prime-term difference by atom type, on the first two ghosts at L = 4:
  - atoms with 5 | n: **0.0** (Λ_DH(5^k·m) = 0, the same as L(χ));
  - prime powers p^k, p ≠ 5: −2.08 / −2.06;
  - **composites prime to 5** (atoms an Euler product would forbid): −2.86 / −2.74;
  - L(χ₅) itself: −1.41.
  - At L = 6 the composite atoms dominate: −8.10 of the first ghost's −7.67 total.
  - **So the change "at the place 5" never shows up at 5.** It shows up as non-Euler coefficients at every other prime and composite. **We are back at the prime term of Weil's criterion.**

## 4. Golden Eye's falsifiers (NUMERICAL)
- **The unchanged construction on EPS at L = 4 and on F_θ at L = 2. W = ‖Tf‖² must fail, and it does.**
  - EPS: λ_min(W) = −1.615 (M = 66) and −1.625 (M = 100), with 4 ghosts. λ_min(arch + pole) = −3.36 / −3.38, so T_∞ does not even exist as a real operator there.
  - F_θ: λ_min(W) = −0.0909 (M = 60) and −0.0910 (M = 100), with 7 ghosts. λ_min(arch + pole) = −2.72 / −2.73.
- **ζ at L = log 2 + 0.05. The identity breaks just past the first prime atom.** λ_min(arch + pole) for ζ against L:

| L | 0.3 | 0.5 | 0.6 | 0.65 | log 2 − 0.01 | **log 2 + 0.05** | 0.8 | 1 | 2 | 4 |
|---|---|---|---|---|---|---|---|---|---|---|
| λ_min(arch + pole) | 0.225 | 0.036 | 0.0086 | 0.0033 | 0.0017 | **+0.0011** | −0.065 | −0.31 | −1.31 | −4.90 |
| W − ‖T_∞g‖² | 0 | 0 | 0 | 0 | 0 | **= −(atom n = 2)**, ‖·‖ = 0.070 | −(n = 2) | −(n = 2) | −(n ≤ 7) | −(n ≤ 53) |

  - For L < log 2 the identity W = ‖T_∞g‖² holds **exactly**, with zero defect.
  - At L = log 2 + 0.05 the defect is exactly the n = 2 prime atom. T_∞ still exists there, since arch + pole is still positive, but it is no longer W.
  - From L ≈ 0.8 on, arch + pole is itself indefinite, so T_∞ no longer exists.
  - **So the prolate / archimedean T is Connes–Consani 2021 (Selecta), not something new** (CHECKED against their support bound [2^{−1/2}, 2^{1/2}]).
- **Semilocal T_S for ζ (Euler factors at p ≤ P):**

| L | P = 2 | 3 | 5 | 7 | 11 | 23 | 53 |
|---|---|---|---|---|---|---|---|
| 2 | −1.03 (4 ghosts) | −0.68 (6) | −0.11 (2) | **0 (S ⊇ {p < e²})** | 0 | 0 | 0 |
| 3 | −2.53 (7) | −2.21 (12) | −1.81 (14) | −1.50 (10) | −0.99 (8) | **0 (complete)** | 0 |
| 4 | −5.06 (10) | −4.95 (14) | −4.58 (19) | −4.25 (20) | −3.83 (20) | −2.39 (12) | **0 (complete)** |

  - W_S ≥ 0 happens **only** when S contains every prime below e^L. At that point W_S *is* W, and T_S = W^{1/2}, which exists iff W ≥ 0. That is circular.
  - With S incomplete, the semilocal form is indefinite, so T_S does not exist. The defect is exactly the omitted primes' atoms (PROVED, by the definition of the split; NUMERICAL, the signs).

## 5. Verdicts
1. **Papers (§1).** The PROVED content is: a self-adjoint prolate operator commuting with Fourier; Sonin eigenvectors for its negative spectrum; the semilocal Hardy–Titchmarsh transform and Sonin-space stability; and archimedean Weil positivity (CC 2021). The UV match is **semiclassical counting plus numerical figures**, and the counting is purely archimedean. No prime enters the prolate operator. The semilocal prolate operator is not yet constructed (OPEN in the literature).
2. **Reproduce the UV spectrum: NOT RUN, by decision.** The pair test (§3) left no room for any archimedean T. The UV counting term is determined by Γ_ℝ alone and is shared by ζ, BEUR″, and any function with ζ's gamma factor (PROVED, elementary). So reproducing it cannot discriminate anything. This can be done later as a pure code check if wanted.
3. **‖Tf‖² = W.**
   - **Exact only for L < log 2** (T_∞, archimedean: Connes–Consani 2021).
   - For L > log 2 the **precise defect is the prime term** −2Σ_{p^k ∈ (log 2, L)} Λ(p^k) p^{−k/2} A(log p^k). At finite S it is exactly the omitted primes' atoms.
   - With all primes below e^L included, T = W^{1/2}, which is circular.
   - (PROVED: decomposition; CHECKED: zero defect for L < log 2; NUMERICAL: indefiniteness.)
4. **Fakes.**
   - The construction **does not distinguish** BEUR″ from ζ, or DH from L(χ₅): T is identical across each pair (PROVED).
   - It "distinguishes" EPS and F_θ only through a different Γ factor, and in the wrong direction: their archimedean part is more positive on their own ghosts.
   - Every ghost of every fake lives in the prime term. **The prolate construction fails exactly where the ghosts are, because it does not contain the place where they are.**

## 6. Best lead (exact claim) and next step
- **Exact claim (PROVED, an elementary decomposition plus CC 2021).** For every finite set of places S ∋ ∞ and every g supported in [−L/2, L/2]:
  W(g) = ‖T_S g‖²_formal − 2 Σ_{p ∉ S, p^k < e^L} Λ(p^k) p^{−k/2} A_g(log p^k)
  where ‖T_S g‖²_formal = W_S(g) is the semilocal Weil form.
  - Hypotheses: ζ's data, and A_g the autocorrelation of g.
  - **Fact:** W_S is positive semidefinite for S = {∞} and L ≤ log 2 (Connes–Consani 2021).
  - **Fact:** numerically, for ζ and L = 2, 3, 4, W_S is indefinite unless S ⊇ {p < e^L}.
  - So any T must contain *all* primes up to the window. The positivity mechanism has to act on the complete prime term at once, which is the Euler product as a whole, not place by place.
  - ~~This rules out every finite-S / semilocal-prolate route as a source.~~ **Withdrawn (errata 5):** it rules out positivity of incomplete truncations W_S and T = W^{1/2} only; the Connes–Consani Sonin-compression route (S complete for the support) is untouched.
- **Next step (OPEN).** Test a **global** T that uses all integers at once, built from the Connes–Consani map E(f)(u) = u^{1/2} Σ_n f(nu). It is the only integer-built map in their framework that contains every prime.
  - Question: is W(g) = ‖P_Sonin E(g)‖² − (explicit remainder) with a remainder that has a sign?
  - Falsifiers first: E for EPS (Σ over the lattice ℤ² with Q) and for DH (Σ a(n) f(nu)). The pair test applies again.

## Labels
- §1 literature statements: CHECKED (read from the arXiv texts: 2112.05500v1; 2310.18423v2).
- Γ-only UV main term: PROVED (elementary).
- Pair identities T(ζ) = T(BEUR″), T(χ₅) = T(DH): PROVED (construction), CHECKED (0.0 difference).
- Sign table, atom split, Golden Eye falsifiers, semilocal scan: NUMERICAL.
- Zero defect for L < log 2, and the defect being exactly the n = 2 atom at log 2 + 0.05: code checks (tautological by construction; errata 1–2).
- Global-E route: OPEN.
- RH: OPEN. Nothing here closes RH, ARP or H_ψ.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5). Audit: Turbo.
