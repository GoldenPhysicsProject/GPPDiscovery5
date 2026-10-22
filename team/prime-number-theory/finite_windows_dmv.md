# Finite windows cannot separate ζ from DMV-type Beurling systems

Goddess of Prime Number Theory, 2026-10-22. **RH is OPEN.** Nothing here bears on the truth of RH. The result
is a limitation on one kind of proof strategy.

Daniel Toupin / Golden Physics Project, with AI assistance from Grok.

Labels: PROVED (literature, with citation, or an elementary argument given in full here), CERTIFIED (Arb/interval),
NUMERICAL, CONJECTURE.

## 0. Correction to the brief
The DMV paper is **Math. Ann. 334 (2006) 1–36**, doi:10.1007/s00208-005-0638-2. It is not in Ramanujan J.
The authors also remark after their Theorem 1: *"By adjusting finitely many primes we could arrange to have κ = 1
in (3) without affecting the other assertions."* The gluing in §3 below is an explicit version of that remark.
The remark gives no proof, so §3 supplies one.

## 1. Definitions

**Beurling system.** This is a non-decreasing right-continuous prime measure Π_P on (1,∞). The discrete case puts
mass 1 on each g-prime q_j > 1. The integer measure is dN_P = exp(dΠ_P), with multiplicative convolution as in DMV §1.
The zeta function is ζ_P(s) = ∫_{1^-}^∞ x^{-s} dN_P(x) = exp ∫ x^{-s} dΠ_P(x) for σ > 1. Here Π_P is the *full* (Riemann-type) prime measure, so log ζ_P(s) = ∫ x^{-s} dΠ_P(x), and for a discrete system it puts mass 1/k at each q_j^k. The Chebyshev measure is dψ_P = log x · dΠ_P. Equivalently, −ζ_P'/ζ_P(s) = ∫ x^{-s} dψ_P(x). For the rational primes, dψ = Σ_n Λ(n) δ_n.

**Window form (the team's convention: s62/cosritz.py, Zhu arXiv:2608.24827).** Fix L > 0. For
f ∈ L²(−L, L), write g = f * f̃ with f̃(u) = conj f(−u). Then supp g ⊂ [−2L, 2L], and g is continuous and vanishes at ±2L.
Define

  Q_L^P(f) = ĝ(i/2) + ĝ(−i/2) − ∫_{1}^{∞} x^{-1/2} (g(log x) + g(−log x)) dψ_P(x) + W_∞(g),

where W_∞ is **ζ's** archimedean functional (the Γ_ℝ term, exactly as in cosritz.py's k1/k2 kernels and c0).
For P = rational primes, Q_L^ℙ = Q_L^ζ, the Weil quadratic form on the window. The Weil criterion is
PROVED (Weil 1952; Bombieri 2000): RH ⇔ Q_L^ζ ≥ 0 for every L.
For a general Beurling P, Q_L^P is **a definition, not a consequence of an explicit formula**. ζ_P has no
functional equation, so W_∞ is borrowed from ζ.

**Equal below the cutoff.** P agrees with the rational primes up to X when ψ_P and ψ coincide as measures
on [1, X]: ψ_P(x) = ψ(x) for all x ≤ X. With the convention above, Q_L^P depends on dψ_P only through
[1, e^{2L}] (g(±log x) = 0 for x ≥ e^{2L}). So:

> **Lemma 0 (PROVED, elementary).** If ψ_P = ψ on [1, X] and 2L ≤ log X, then Q_L^P(f) = Q_L^ζ(f) for every
> f ∈ L²(−L, L). In particular, the Gram matrix in any basis, its LDLᵀ pivots, and its λ_min are all identical.

## 2. Theorem

> **Theorem (PROVED, given DMV Theorem 1 and the elementary gluing in §3).** Fix 1/2 < θ < 1,
> a > (4/e)(1 − θ), and any X ≥ 2. There is a discrete Beurling system P_X such that
> 1. **(cutoff)** the g-primes of P_X in [1, X] are exactly the rational primes ≤ X. Hence ψ_{P_X} = ψ on [1, X],
>    and Q_L^{P_X} ≡ Q_L^ζ for every L ≤ (log X)/2;
> 2. **(integers)** N_{P_X}(x) = κ_X x + O_X(x^θ) with κ_X > 0;
> 3. **(continuation)** ζ_{P_X} is analytic on σ ≥ θ except for a simple pole at s = 1;
> 4. **(off-line zeros)** ζ_{P_X} has infinitely many zeros on σ = 1 − a/log t, t ≥ 2. Every zero with
>    t > e^{2a} has real part > 1/2, so "RH for ζ_{P_X}" fails;
> 5. **(PNT, de la Vallée Poussin quality)** ψ_{P_X}(x) = x + O(x exp(−2√(a log x))) and
>    ψ_{P_X}(x) − x = Ω_±(x exp(−2√(a log x))). This is DMV (iv)–(v) up to an O_X(log x) change.
>
> **Corollary (PROVED).** Let 𝒯 be any test whose input is Q_L^P for L in a bounded set [0, L₀] only. Examples are a
> positivity certificate at L = 1.7105, any finite list of Gram pivots, or λ_min(L) for L ≤ L₀. Then 𝒯 gives the
> same answer on ζ and on P_X with X = e^{2L₀}, and P_X has zeros off the critical line. So no bounded-window test can
> imply RH by itself. An argument that reaches RH must also use either (a) windows with L → ∞,
> uniformly, or (b) structure that P_X lacks, chiefly ζ's functional equation, i.e. the match between W_∞ and
> the prime side. Note that integer regularity N(x) − κx = O(x^{θ'}) with θ' < 1/2 does **not** rescue this. See
> integer_regularity.md: such systems with off-line zeros exist unconditionally.

## 3. Proof

Let P be the DMV system of Theorem 1 for the given (θ, a). It is discrete (DMV §§6–9, probabilistic
discretisation of a continuous Π_C). Let q_1 < q_2 < … be its g-primes. Define P_X = {rational primes p ≤ X} ∪ {q_j > X}.

*Step 1 (Euler factors).* For σ > 1,
 ζ_{P_X}(s) = ζ_P(s) · E_X(s), E_X(s) = ∏_{p ≤ X} (1 − p^{-s})^{-1} · ∏_{q_j ≤ X} (1 − q_j^{-s}).
Both products are finite. Each factor 1 − c^{-s} with c > 1 vanishes only when c^{-s} = 1, i.e. on σ = 0.
So E_X is holomorphic and zero-free on σ > 0 and extends ζ_{P_X} = ζ_P·E_X meromorphically to σ ≥ θ.
Then **3** holds, and the zeros of ζ_{P_X} in σ > 0 are exactly those of ζ_P, which gives **4** from DMV (ii)–(iii).

*Step 2 (integers).* dN_{P_X} = dN_P * dμ_X, where dμ_X is the finite multiplicative convolution of the geometric
measures Σ_k δ_{p^k} (p ≤ X) with δ_1 − δ_{q_j} (q_j ≤ X). Its Mellin transform E_X(s) converges absolutely for
σ > 0, so ∫ t^{-θ} |dμ_X|(t) < ∞. The standard convolution estimate (e.g. Diamond–Zhang, *Beurling Generalized
Numbers*, AMS Math. Surveys 213 (2016), Ch. 2–3) gives
 N_{P_X}(x) = ∫ N_P(x/t) dμ_X(t) = κ E_X(1) x + O(x^θ ∫ t^{-θ}|dμ_X|) = κ_X x + O_X(x^θ),
with κ_X = κ E_X(1) > 0. This is **2**.

*Step 3 (primes).* ψ_{P_X} − ψ_P = Σ_{p≤X, p^k≤x} log p − Σ_{q_j≤X, q_j^k≤x} log q_j = O_X(log x). This gives **5**
from DMV (iv)–(v).

*Step 4 (cutoff).* Below X, the only g-primes and g-prime powers of P_X are the rational primes p ≤ X and their powers
p^k ≤ X. Every q_j > X has all powers > X. So dψ_{P_X} = Σ_{n ≤ X} Λ(n) δ_n on [1, X], which is **1** by Lemma 0. ∎

Nothing in Steps 1–4 goes beyond finite products and DMV's stated theorem. That is why the label is PROVED.
If one prefers to cite only DMV's own remark, the claim is "DMV + elementary". Same label.

## 4. Gaps and caveats (honest list)
1. **Q^P is not canonical.** For P_X, positivity of Q_L^{P_X} for all L is **not** known to be equivalent to
   RH for ζ_{P_X}. With no functional equation, there is no Weil explicit formula pairing W_∞ with P_X's prime side.
   So the theorem says only that bounded windows cannot see the difference. It does **not** say that the windows of P_X fail
   to be positive for large L. **Open:** is Q_L^{P_X} ≥ 0 for all L? A natural CONJECTURE is "no, for L ≳ (log X)/2 + O(1)".
   We have not tested it. A numerical test is possible with cosritz-style assembly plus synthetic g-primes above X.
2. **DMV's zeros are near σ = 1, not near 1/2.** This makes the fake "far" from ζ in the zero picture, but invisible in the window.
   Fakes with zeros only slightly off the line under RH exist too (BDR, see integer_regularity.md). They give the
   same corollary with θ < 1/2.
3. **Endpoint convention.** If a team script includes n with log n = 2L exactly, use 2L < log X. For
   L = 1.7105, e^{2L} ≈ 30.6, so X = 31 suffices.
4. **Discreteness.** This is used only to call P_X "a number system". The proof works verbatim for measure systems.
5. **Scope.** The theorem is about tests that are functions of Q_L for L ≤ L₀. Tests that use ζ's analytic structure (zeros,
   functional equation, Γ-factor consistency, Li coefficients defined via ξ) are outside its scope.

## 5. Fakes battery relevance
| fake | separated from ζ by a bounded-window test? |
|---|---|
| **DMV (P_X glued)** | **No, by construction.** It is the content of the theorem. |
| DH (Davenport–Heilbronn) | Yes in principle. DH's coefficients b(n) differ from Λ(n) already at n = 2,3. Its window form differs at small L, and the team control shows DH is negative at L ≈ 1.71 (see precision/README.md). DH lacks an Euler product, so it is not a Beurling system. |
| H₋ₜ (team: zero-swap / de Bruijn–Newman-type deformation, GROK_HEURISTICS S72) | This fake is built on the zero side and keeps ζ's functional equation. Its window form differs from Q^ζ through the moved zeros, so bounded windows can separate it. It is **not** a Beurling system, and the theorem says nothing about it. |
| F_θ = ξ(s−θ)ξ(s+θ) | Euler product, but its prime weights are Λ(n)(n^θ+n^{−θ}). These differ from Λ at every n, so bounded windows separate it (the team has seen failure at modest L). It also differs from Λ below every cutoff, so it is **not** of P_X type. |
| G_θ' | I found no definition in GPP-bridge/GPPDiscovery5 under that name, so it is not classified here. Rule: it is inseparable by windows L ≤ L₀ iff its prime measure equals Λ on [1, e^{2L₀}] (Lemma 0). |
| P_X itself (new fake proposed for the battery) | It agrees with Λ exactly below X, has an Euler product, PNT with dlVP error, N-error x^θ, and off-line zeros, but no functional equation. Any mechanism that passes ζ at L ≤ (log X)/2 also passes P_X there. So a mechanism must use the functional equation or uniformity in L. |

The bottom line for Daniel's directive ("where do the integers enter"): at any finite window, the integers enter **only** through Λ on
[1, e^{2L}]. Agreement there is cheap to arrange in a Beurling system with off-line zeros. So the integers' decisive role must come
through behaviour uniform in L, or through the functional equation, which ties the integers to the Γ-factor.

## References
- H. G. Diamond, H. L. Montgomery, U. M. A. Vorhauer, *Beurling primes with large oscillation*, Math. Ann. 334 (2006) 1–36.
  https://doi.org/10.1007/s00208-005-0638-2 (open copy: https://deepblue.lib.umich.edu/items/80a3d22b-2022-435e-9bcb-360f281fb998)
- H. G. Diamond, W.-B. Zhang, *Beurling Generalized Numbers*, AMS Math. Surveys and Monographs 213 (2016).
- A. Weil, *Sur les "formules explicites" de la théorie des nombres premiers*, Comm. Sém. Math. Lund (1952).
- E. Bombieri, *Remarks on Weil's quadratic functional in the theory of prime numbers I*, Rend. Mat. Acc. Lincei 11 (2000) 183–233.
- F. Broucke, G. Debruyne, Sz. Gy. Révész, *Some examples of well-behaved Beurling number systems*, arXiv:2309.01567.
- Team convention: GPPDiscovery5 s62/cosritz.py (Grok S62); Zhu, arXiv:2608.24827.
