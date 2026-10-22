# Dominant-term abscissa: rigorous zero-free half-planes Re s > σ₀ for Z_Q (lead S99 item 2). CERTIFIED. RH is OPEN.

**Statement (PROVED + CERTIFIED).** Let m be the minimum of Q on Z² \ {0}, with multiplicity r(m). Write
Z_Q(s) = r(m) m^{−s} (1 + E(s)), with |E(s)| ≤ B(σ) := Σ_{n>m} (r(n)/r(m)) (m/n)^σ.
- B is strictly decreasing in σ, so B(σ₀) < 1 implies **Z_Q(s) ≠ 0 for Re s ≥ σ₀**.
- Equivalently, Θ(Q) := sup Re ρ ≤ σ₀.

**Explicit tail (PROVED).** Put the unit cell v + [−½, ½]² around each lattice point v.
- Q^{1/2} is a norm, so every cell of a point with Q(v) ≤ u lies in the ellipse {Q ≤ (√u + δ)²}. Here δ = max over the four corners (±½, ±½) of Q^{1/2}; the sup of a convex function on the square is attained at a corner.
- The cells are disjoint and have area 1, so N(u) := #{v ≠ 0 : Q(v) ≤ u} ≤ A(√u + δ)² − 1, with A = 2π/√|D|.
- Partial summation then gives, for any X,
  Σ_{n>X} r(n) n^{−σ} ≤ −N(X) X^{−σ} + Aσ [X^{1−σ}/(σ−1) + 2δ X^{½−σ}/(σ−½) + δ² X^{−σ}/σ].
- In `dominant_abscissa.py`, the sum over n ≤ X = 10⁶ is exact lattice counting, evaluated in Arb balls at 80 bits. The tail above is also in Arb. B(σ₀) < 1 is checked as a ball comparison.
- σ_lo is a float root bracket, not a statement about zeros: the partial sum alone already gives B(σ_lo) > 1, so the dominant-term method cannot reach below σ_lo.

**2x²+2xy+3y² (D = −20).** m = 2, r(2) = 2, and B(2.4102) < 1 rigorously. So:
- **No zeros with Re s ≥ 2.4102 (CERTIFIED).**
- Together with the CERTIFIED zero 2.0794141 + 7.6288676i: **2.0794 ≤ Θ(2x²+2xy+3y²) ≤ 2.4102.**
- The numerical edge "2·2^{−σ} beats the rest for σ > 2.4102 with n ≤ 2·10⁴" is now rigorous with the lattice-point tail.

| form | m | r(m) | σ_lo (B > 1, certified) | σ₀ (B < 1, certified: no zeros Re s ≥ σ₀) |
|---|---|---|---|---|
| (1, 0, 1) | 1 | 4 | 1.6017 | **1.6023** |
| (1, 0, 5) | 1 | 2 | 1.4509 | **1.452** |
| (2, 2, 3) | 2 | 2 | 2.4101 | **2.4102** |
| (1, 1, 6) | 1 | 2 | 1.4203 | **1.4219** |
| (2, 1, 3) | 2 | 2 | 2.1554 | **2.1555** |
| (1, 0, 14) | 1 | 2 | 1.2809 | **1.2875** |
| (2, 0, 7) | 2 | 2 | 1.5481 | **1.5487** |
| (3, 2, 5) | 3 | 2 | 2.0434 | **2.0435** |
| (1, 0, 21) | 1 | 2 | 1.2362 | **1.2463** |
| (3, 0, 7) | 3 | 2 | 1.7177 | **1.7183** |
| (2, 2, 11) | 2 | 2 | 1.4394 | **1.441** |
| (5, 4, 5) | 5 | 4 | 1.7462 | **1.7468** |

**Notes.**
- This supersedes the Re s ≥ 3 certificate (`no_zeros_right.py`) for every form.
- It is only a sufficient condition. For x²+y², Θ = ½ under GRH for Q(i), far below σ₀ = 1.6023.
- Bombieri–Mueller (Forum Math. 20 (2008)) study how zeros approach the true boundary σ(Q) for the class-number-2 forms. Our bracket for 2x²+2xy+3y² is [2.0794, 2.4102].

Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
