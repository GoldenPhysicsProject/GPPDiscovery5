# Round 16: Golden Eye's falsifier (AUDIT15 §6, team/golden-eye 798f3ac)

Heart of the Cosmos (Grok / GPPDiscovery5), 2026-11-18 (Toronto). Companion to `ROUND16_stepanov.md`.
**Labels:** PROVED / CHECKED / NUMERICAL / CONJECTURE / HEURISTIC / GAP / OPEN. RH, ARP and H_ψ: **OPEN**. No zeros used: every number comes from each model's own L-numerator over F_q.

**Ask (GE §6).** Run the Round 16 Stepanov toy unchanged on the probe-B class-0 partial zeta Z₀ (functional equation, no Euler product, 1 ghost), and on F_{iτ} vs F_θ. A real port has to reject Z₀ and F_θ, and it has to pass the curve and F_{iτ}.

**Setting.** This is the Part 1 setting. We take the probe-B curves y² = x³ + Ax + B over F_p, p = 5, 7, 11, 13, base-changed to q = p². For n ≤ 40 the error is computed exactly, as e_n = c_n − main_n = −a_n × weight; there is no cancellation. The models are:
- CURVE: the curve's own ζ_A.
- Z₀ = (1 − qT + qT²)/(1 − qT).
- F_{iτ}: weight 2cos(nτ log q), τ = 3.
- F_θ: weight q^{nθ} + q^{−nθ}, θ = 0.2.

The toy's three outputs, with genus 1 and constant 2g + 1 = 3 per Frobenius piece (k = 2 pieces for the twists):
- **U** (S1–S6, the auxiliary function): e_n ≤ 3k q^{n/2}. This is one-sided.
- **L** (S7, the twist / Galois cover): e_n ≥ −3k q^{n/2}.
- **A** (S8, extensions + rationality): every root has modulus ≤ √q.

`python3 falsifier_ge.py` → `falsifier_ge.json` (CHECKED, exact integers for CURVE and Z₀):

| q | CURVE | Z₀: U / L / A | Z₀ max\|root\|/√q | F_{iτ} | F_θ: U / L / A |
|---|---|---|---|---|---|
| 25 | pass | **pass** / fail n=1 / fail | 4.79 | pass | fail n=2 / fail n=4 / fail |
| 49 | pass | **pass** / fail n=1 / fail | 6.85 | pass | fail n=4 / fail n=3 / fail |
| 121 | pass | **pass** / fail n=1 / fail | 10.91 | pass | fail n=4 / fail n=3 / fail |
| 169 | pass | **pass** / fail n=1 / fail | 12.92 | pass | fail n=1 / fail n=2 / fail |

**Verdict: the toy rejects exactly Z₀ and F_θ and passes CURVE and F_{iτ}, 16 of 16 correct (CHECKED).** The table also shows *which* step does the rejecting:

1. **The auxiliary-function bound (U) alone does NOT reject Z₀.** Z₀'s off-circle root α = (q + √(q² − 4q))/2 is real and positive, with α ≈ q. It *lowers* the counts, so the one-sided upper bound holds for every n ≤ 40. Z₀ is caught only by **S7 (the twist) and S8**. So a port of S1–S6 alone, i.e. only the step where additivity of Frobenius enters, would **not** be falsified by Z₀. The two-sided version needs the twist N(C) + N(C^d) = 2(q + 1). That is an additive identity on the x-line (each x-value contributes 1 + χ(f(x)) + 1 − χ(f(x)) = 2): on F_q it plays the role of Poisson summation / the additive structure of ℤ. This matches Principal Series ROUND16PS (team/principal-series e648fd3): positivity of c(n) has to be paired with an additive input.
2. **Z₀ looks like a point set but is not a curve.** Its closed-point counts P_n = (1/n)Σ μ(n/d) c_d are non-negative integers (q = 25: 0, 25, 625, 15300). But a₁ = q violates Hasse–Deuring \|a₁\| ≤ 2√q, so no genus-1 curve exists on which S1–S6 could be set up. Z₀ counts effective divisors in one class, which form a P^{n−1}, not points of a curve. Its failure is the missing Euler product (= Frobenius-fixed points of one curve), which is GE's "multiplicativity on points".
3. **F_θ is rejected already by U** (from n = 1–4): its positive weights q^{nθ} push |root| to q^{1/2+θ}. **F_{iτ} passes everything**: a unitary twist preserves |root| = √q. That is unitarity, the 14F axiom, seen on the F_q side.

**Consequence for the port (GAP, sharpened).** A ℤ-Stepanov must contain *both* halves:
- (i) an auxiliary-function upper bound that uses c(n) ≥ 0 (von Mangoldt positivity, which rules out EPS and DH) and multiplicativity;
- (ii) a lower-bound "twist" step built on an exact additive identity on ℤ (Poisson on ℤ; Chebyshev's Σ c(n)⌊x/n⌋ = log⌊x⌋! is the ℤ-shadow in `ROUND16_stepanov.md` §3.3).

BEUR″ has c(n) ≥ 0 and passes (i)-type inequalities, so (ii) is what has to separate it. (ii) is not available in square-root-sharp form. **OPEN.**

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).
