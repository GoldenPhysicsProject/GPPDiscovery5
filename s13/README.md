# S13 — Coupling of prime 2 to the pole+arch negative directions

**S13.0 (PROVED, elementary).** For ½log2 < h < ½log3, the only prime power with log n < 2h is n = 2, so on C_h
W_ζ = A(h) − c·T₂(h) exactly, where A = pole + Γ_R term, T₂ = symmetrized compressed translation by log 2, c_ζ = 2·log2·2^{-1/2} = √2·log 2 = 0.980258.
So in this window the 2-place truncation is not an approximation: it *is* ζ (the Euler factor at 2 contributes only its m = 1 term). Same for F_θ, with A_θ and c_θ = (2^θ+2^{-θ})·log2/√2 = 0.995013 at θ = ¼.

**S13.1 (PROVED, trivial / negative for task (2)).** The completed factor Γ_R(σ+ir)(1−2^{-σ-ir})^{-1}, normalized, is positive-definite in r (both factors are characteristic functions of positive measures, and products of PD functions are PD). Γ_R alone is PD too. So multiplicative PD can't tell ∞ alone apart from ∞·2. The coupling only shows up after log-differentiation (S12.1), so the "completed local factor PD" test is vacuous.

**S13.2 (NUMERICAL, float Ritz M=14; output.txt).** The bottom vector of A lines up with T₂'s negative direction: ⟨x,T₂x⟩ < 0 every time A < 0. In the odd sector the overlap with T₂'s bottom eigenvector grows: 0.36, 0.58, 0.75, 0.82 at h = .40, .45, .50, .54. Each negative direction of A is a direction where translation by log 2 is anti-correlated, i.e. close to the "reflection" g(x) ≈ −g(x+log2). A has exactly one negative eigenvalue per parity sector.

**S13.3 (NUMERICAL, float; output_pinch.txt) — PINCH.** Admissible coupling window J(h) = {c : A − cT₂ ≥ 0 on both parities}:
| h | J(h) | width |
|---|---|---|
| .50 | [0.973962, 0.980726] | 6.8e-3 |
| .52 | [0.978285, 0.980389] | 2.1e-3 |
| .54 | [0.979654, 0.980300] | 6.5e-4 |
| .549 | [0.979895, 0.980284] | 3.9e-4 |
- **Even** upper edge → c_ζ from above: 0.980284 − 0.980258 = 2.6e-5 at h = .549.
- **Odd** lower edge → c_ζ from below.
- The even sector caps how much prime 2 is allowed (too much breaks it); the odd sector needs a minimum amount.
- The window shrinks onto √2·log 2: the Euler weight at 2 is tuned to within ~3e-5 of the edge.
- F_{1/4}: c_θ = 0.995 lies outside the even window for h ≥ .50, and its own A_θ windows confirm it (lamW < 0 at h = .50, .54 even).
- Where the coupling breaks: the non-unitary weight 2cosh(θ log2) > 2 pushes past the even "too-much-prime" edge. This is the S9/S10 unitarity statement, now seen as violating one specific side of a two-sided window.

**Reading.** The missing axiom is two-sided. Unitarity (|weight| ≤ 2) gives the even/upper edge. The odd/lower edge says prime 2 must be present with at least nearly its full weight, which is a lower bound that no local axiom provides. The two edges meeting near c_ζ is the small-support picture of "barely true" (cf. de Bruijn–Newman Λ = 0).

**OPEN / CONJECTURE.** Does the true (M→∞) window degenerate to the single point {c_ζ} at some h* ≤ ½log3, or stay at positive width? Ritz edges are not certified, and as M grows the window can only shrink. If it closes before c_ζ, ζ's form would fail on C_h, so for these h the conjecture is "width → 0⁺ as h → h*". Needs larger M plus interval certification.

Does NOT claim: RH (OPEN); certified bounds; that the window structure persists once n = 3 enters (h > ½log3).
No φ.
