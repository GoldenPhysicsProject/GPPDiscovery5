# S14 — n=2 window: high-M Ritz + Arb certification of the prime-2 coupling window

Setup (as S13): on C_h, ½log2 < h < ½log3, W_ζ = A − c_ζ T₂ exactly, c_ζ = √2·log2. Basis is Legendre·(1−u²)⁴, j < M.

**S14.0 (PROVED, standard / EXACT construction).** The entries are now exact. F_jk(x) = h·G_jk(x/h) with G_jk an exact rational polynomial. The archimedean term is written in real space using Gauss's digamma integral:
arch(F) = (−γ−logπ)F(0) − F(0)log(1−e^{−4h}) + ∫₀^{2h} [(F(0)−F)/x]·e^{−x}·x/sinh x dx − 2∫₀^{2h} F·e^{−x/4} sinh(3x/4)/sinh x dx.
Moments are computed by Arb's acb.integral, which is rigorous. At M = 14, h = .5 this matches the float S13 matrices to 1.6e−15. (s14_exact.py)

**S14.1 (CERTIFIED, Arb).** For each M we take the Ritz null vector, round it to rationals, call that explicit g, and compute c* = ⟨g,Ag⟩/⟨g,T₂g⟩ as an Arb ball.
- Even sector (⟨g,T₂g⟩ > 0): A − cT₂ is not PSD on C_h for every c > c*.
- Odd sector (⟨g,T₂g⟩ < 0): A − cT₂ is not PSD for every c < c*.

At h = 0.549, M = 110:
- even: c* − c_ζ = 3.317137951e−6 ± 8e−16
- odd: c* − c_ζ = −8.797147555e−5 ± 5e−15

So the true admissible window satisfies J(0.549) ⊂ [c_ζ − 8.798e−5, c_ζ + 3.318e−6] (certified), with width ≤ 9.13e−5. At h = ½log3 (0.54930614, M = 80) it is ⊂ [c_ζ − 9.173e−5, c_ζ + 3.468e−6]. F_{1/4} (c = 0.99501) is excluded for certain.

**S14.2 (NUMERICAL, Ritz edges, mpmath dps 120).** Even gap at h = .549 for each M:
| M | even gap |
|---|---|
| 14 | 2.60e−5 |
| 20 | 1.41e−5 |
| 30 | 7.37e−6 |
| 40 | 5.25e−6 |
| 56 | 4.10e−6 |
| 80 | 3.54e−6 |
| 110 | 3.32e−6 |

The successive drops shrink by a factor of about 0.4 each step, so the gap is saturating, and extrapolation gives about 3.1–3.2e−6 > 0. The odd gap saturates near −8.5e−5. The edges never cross c_ζ, which is consistent with Weil positivity on this window. Also consistent with Connes–Consani's numerics ("Spectral triples and ζ-cycles", arXiv 2106.01715 §2.2–2.3), which show positivity for log2 ≤ L < log3 by numerics only, not by proof. The proved range is only L ≤ log2, i.e. h ≤ ½log2 (Connes–Consani, Selecta, arXiv 2006.13771).

**Conjecture (3), "edge → c_ζ", is NOT SUPPORTED.** The even edge converges to about c_ζ + 3e−6, not to c_ζ, even at h = ½log3. ζ's prime-2 weight sits just inside the admissible window, about 3e−6 below its upper edge, and is not the maximal admissible weight. The window also does not degenerate within the n = 2 range.

**Prior art.** Connes–Consani 2021 §2.3 already found the sharpness phenomenon by varying p near 2 at L = log3, where p is admissible only within an interval < 1e−3. What is new here is that the window is in the weight c, with explicit certified rigorous bounds (Arb), and is split by parity: even gives the cap, odd gives the floor.

Does NOT claim: RH (OPEN); positivity of W_ζ on C_h (the lower bounds on the true even edge and the upper bounds on the true odd edge are not certified); anything beyond n = 2. No φ.
