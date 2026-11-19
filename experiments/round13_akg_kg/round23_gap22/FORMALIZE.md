# Round 23: formalizable statements, a handoff for the verify / Lean track

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).

These are self-contained statements from ROUND23.md (and a few from R20/R22). They are ordered from easiest to hardest.
None of them asserts RH. Notation follows ROUND23.md §0: b_k = 2k + ½, and Λ is von Mangoldt.

## Tier 1: elementary real analysis (Mathlib-ready)

**F1 (convexity of the archimedean part, Lemma 5a).** For real L ≥ log 2: (e^L + 1)(1 − e^{−2L}) > 1. Hence
e^{L/2} + e^{−L/2} − e^{−L/2}/(1 − e^{−2L}) > 0.
*Lean hint:* set x = e^L ≥ 2. The inequality is (x + 1)(1 − 1/x²) > 1, i.e. x³ + x² − x − 1 > x², i.e. x³ − x − 1 > 0 for x ≥ 2.
Use `nlinarith`/`positivity`.

**F2 (geometric series for T'').** For L > 0: Σ_{k≥0} e^{−(2k+½)L} = e^{−L/2}/(1 − e^{−2L}).
Also the tail bound: Σ_{k≥K} e^{−b_k L}/b_k² ≤ e^{−b_K L}/(b_K²(1 − e^{−2L})).

**F3 (two-tangent lemma, Lemma 5c).** Let f: [a,b] → ℝ be convex and differentiable, with f'(a) < 0 < f'(b). Then for all
x ∈ [a,b], f(x) ≥ f(a) + f'(a)(L* − a), where L* = (f(b) − f(a) + f'(a)a − f'(b)b)/(f'(a) − f'(b)). The one-sided version
for piecewise-smooth f also holds. *Mathlib:* `ConvexOn` and the tangent-line inequality.

**F4 (piecewise-affine prime sum, Lemma 5b).** For q_j ≤ e^L < q_{j+1} (consecutive prime powers):
S(L) := Σ_{n≤e^L} Λ(n) n^{−½}(L − log n) = α_j L − β_j, with α_j = Σ_{n≤q_j} Λ(n) n^{−½} and β_j = Σ_{n≤q_j} Λ(n) n^{−½} log n.
S is continuous at each L = log q. (It is a finite-sum identity.)

**F5 (Laplace transforms of the pieces).** For σ > ½:
(a) ∫₀^∞ e^{−σL}·4(e^{L/2} − 2 + e^{−L/2}) dL = 2/(σ(σ² − ¼));
(b) ∫₀^∞ e^{−σL}(L − a)₊ dL = e^{−σa}/σ²  (a ≥ 0);
(c) 1/(σ + ½) + 1/(σ − ½) = 2σ/(σ² − ¼), which gives the 1/s + 1/(s−1) term of ξ'/ξ at s = ½ + σ.

**F6 (box Gram identity, Theorem 7 algebra).** For real x, y: 1 − cos x − cos y + cos(x−y) = Re[(1 − e^{ix})(1 − e^{−iy})].
Hence, for any real γ_1..γ_N and weights w_k ≥ 0, the matrix
G_{ij} = Σ_k w_k(1 − cos γ_k a_i − cos γ_k a_j + cos γ_k(a_i − a_j)) is positive semidefinite.
*Mathlib:* a sum of `Matrix.PosSemidef` rank-one Hermitian terms, plus `Complex.exp` trig.

**F7 (per-prime mean zero, R22).** For n ≥ 2: lim_{T→∞} (1/T)∫₀^T cos(γ log n) dγ = 0. So no single prime term can be a
passive (one-signed) element.

## Tier 2: classical analysis (needs some library building)

**F8 (Dirichlet ⇒ liminf, Theorem 6 core).** Let w_k ≥ 0 with Σ w_k < ∞, and let γ_k be real. Then
liminf_{L→∞} Σ_k w_k(1 − cos γ_k L) = 0.
*Needs:* Dirichlet's simultaneous approximation theorem (pigeonhole). Check whether Mathlib has the multi-dimensional version.

**F9 (sign of ζ on (0,1)).** For 0 < σ < 1: ζ(σ) < 0. *Route:* ζ(σ) = η(σ)/(1 − 2^{1−σ}), with η(σ) > 0 by alternating series.
(Mathlib has `riemannZeta`; this real-axis sign lemma may need building.)

**F10 (Pringsheim–Landau for Laplace transforms).** Let f ≥ 0 be locally integrable on [0,∞), with ∫ e^{−σL} f dL
converging for σ > σ_c and diverging for σ < σ_c, where σ_c is finite. Then F(σ) = ∫ e^{−σL} f(L) dL cannot be continued
analytically to a neighbourhood of σ = σ_c.
*Standard proof:* expand F in a Taylor series at a real point σ₀ > σ_c. Its coefficients are (−1)^m ∫ L^m e^{−σ₀L} f; by
positivity the series can be summed inside the integral (Tonelli).

## Tier 3: the theorems (Tier 1–2 plus the explicit formula)

**F11 (closed form, Theorems 1–2).** D(L) := 4(e^{L/2} − 2 + e^{−L/2}) + c₁L + C₂ − Σ_k e^{−b_k L}/b_k² − S(L), with
c₁ = (ψ(¼) − log π)/2 and C₂ = ψ'(¼)/4, satisfies ∫₀^∞ e^{−σL} D(L) dL = (ξ'/ξ)(½+σ)/σ² for σ > ½.
*Needs:* the digamma partial-fraction series (ψ(z) = −γ + Σ[1/(k+1) − 1/(k+z)]) and the Dirichlet series of ζ'/ζ for Re s > 1.
Then the identity extends by analytic continuation.

**F12 (Landau reduction, Theorem 3).** If D is bounded below on [0,∞), then ξ has no zeros with Re s > ½.
*Uses:* F9, F10, F11, the Euler product (no zeros on Re s > 1), and ξ(1) = ½. The **converse** (RH ⇒ D ≥ 0) needs the
Hadamard product. It can be left as a separate statement.

**F13 (strip version, Theorem 4, "⇐" direction).** If D(L) ≥ −C e^{εL} for some 0 ≤ ε < ½, then ξ ≠ 0 on Re s > ½ + ε.
(This is F12 applied to D + C e^{εL}.)

## Tier 4: certificate checking (optional)

**F14 (finite certificate ⇒ positivity).** For X₁ = 10⁸, the list of 5,762,859 gaps, with ball enclosures of
(D(a), D'(a⁺), D(b), D'(b⁻)), certifies min D > 0.0214 via F1–F4. A Lean check of the whole list is heavy. A small range
(e.g. X ≤ 100, 35 gaps) could be checked with `norm_num`/interval tactics, as a proof of concept. The certificate is
reproduced by `r23_interval.py <X1>`, in about 3 min for 10⁸.

## Suggested priority

F1, F3, F4, F6 (immediate) → F5, F2 → F9, F8 → F10 → F11–F13. F12/F13 are the conceptually important ones: "GAP-22 is
exactly RH", and "an exponential lower envelope D ≥ −Ce^{εL} gives a zero-free strip".
