# Round 25: how the zeros line up — a visual attack

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).

Labels: PROVED, CHECKED, NUMERICAL, CONJECTURE, HEURISTIC, GAP, OPEN. Zeros are used
here only as plot data (illustration), never as proof inputs. **RH is OPEN; nothing
in these figures closes it.**

Controls: $F_\theta=\xi(s+\theta)\xi(s-\theta)$ with $\theta=\tfrac14$; Epstein
$Z_Q$ for $Q=x^2+5y^2$; Davenport–Heilbronn $f$ (period-5 coefficients). All four
completed functions are real on the critical line.

## Plain-words summary

Daniel asked: how do the zeros line up like they do, and can a picture show something
nobody has seen? The six figures below are the honest answer.

The lining-up is real, and it has a precise mechanism, but the mechanism is classical
and it is shared by the fakes. I did not find a new structural fact that would force
the crossings onto the line for $\zeta$ and only for $\zeta$.

1. **The imaginary part does not disappear only at the zeros. It disappears along
   the whole critical line.** That is just the functional equation plus real
   coefficients: $\xi(\tfrac12+it)$ is real for every $t$, not only at a sign change.
   A zero on the line is an ordinary zero of this real function — i.e. a place where
   the *real* part vanishes too. This is the right picture behind “they live where
   the imaginary part disappears,” but the “where” is the entire line, not the zeros
   themselves. (CHECKED, classical.) Shared by $F_\theta$, Epstein and DH: being
   real on the line does not force the zeros onto it. Visible in Figs. 3–5.

2. **Zeros of $\xi$ are the crossings of the two level-curve families
   $\mathrm{Re}\,\xi=0$ and $\mathrm{Im}\,\xi=0$.** For every model here, the
   critical line is itself a component of $\mathrm{Im}\,\xi=0$. The curves
   $\mathrm{Re}\,\xi=0$ hit that line at right angles at each simple zero
   (Cauchy–Riemann). Extra $\mathrm{Im}\,\xi=0$ branches meet the line at the
   *extrema of $\xi(\tfrac12+it)$*, not at the extrema of Hardy’s $Z(t)$. The two
   families of extrema differ because $\xi(\tfrac12+it)=A(t)Z(t)$ with a decaying
   archimedean factor $A(t)\asymp t^{9/4}e^{-\pi t/4}$. At low height the shift is
   large: the first $\xi$-extremum is at $t\approx 15.586$, while the first peak of
   $Z$ is at $t\approx 17.883$. (CHECKED, classical Cauchy–Riemann plus the known
   shape of $A(t)$. I have not seen the shift emphasised in the usual $Z$-pictures;
   it is not a new theorem.) Visible in Fig. 4.

3. **The fakes use exactly the same lining-up architecture, and still have off-line
   zeros.** An off-line zero is simply a red/blue crossing that happens *away* from
   $\sigma=\tfrac12$, on an extra $\mathrm{Im}=0$ branch. In this window:
   - $F_\theta$ has zeros at $\tfrac12\pm\tfrac14+i\gamma$, i.e. $\zeta$'s zeros
     rigidly shifted off the line;
   - Epstein’s first off-line zero sits at $0.933+15.668i$, between two on-line zeros;
   - DH’s first off-line zero sits at $0.809+85.699i$, between on-line zeros at
     $83.11$ and $87.65$, with a very clear extra blue branch meeting the red curve
     off the line.
   (NUMERICAL, illustration; zero locations from R13 / S84, used only as plot
   marks.) This is the most revealing picture. It does not suggest a new constraint
   that $\zeta$ obeys and the fakes break. Visible in Figs. 2 and 5.

4. **The minus signs are phases, not negative integers.** The $\eta$-walk (Fig. 1)
   and the Riemann–Siegel cosine sum (Fig. 6) show cancellation of rotating arrows
   $n^{-it}$ with *positive* prefactors $n^{-\sigma}$. Starting from only the
   positive integers and multiplying cannot, by itself, prevent that cancellation:
   $F_\theta$ is the counter-example (positive coefficients, Euler product, the
   same shadow symmetry, off-line zeros). At a zero of $\zeta$ the $\eta$-walk
   closes on the origin; at the same height off the line, or on the line between
   zeros, it misses. (NUMERICAL, illustration.)

5. **At low height the Riemann–Siegel remainder is not small.** At the first zero,
   the main sum is a single term $2\cos\vartheta(t)\approx-0.314$, and the remainder
   is $R(t)\approx+0.314$. The first zero is *not* “just $\cos\vartheta(t)=0$.”
   Even at $t=143$ (four terms) one still has $|R(t)|\approx0.24$, consistent with
   the classical size $t^{-1/4}$. (NUMERICAL, illustration; remainder theory
   classical: Riemann–Siegel, Gabcke, Arias-de-Reyna.)

6. **No new structural fact that looks like a proof.** Gram’s law is visible at
   low height in Fig. 3 and is a heuristic, already known to fail (first failure
   between $g_{126}$ and $g_{127}$; Hutchinson, Titchmarsh). The portraits match
   the classical level-curve picture (Spira; the $\zeta$-curve literature;
   Arias-de-Reyna’s high-precision $Z$). I did not see a pattern in the red/blue
   curves that would be fatal to the fakes and provable for $\zeta$. The dividing
   line remains the stronger arithmetic (coefficients exactly $1$, each integer
   counted once — Hamburger uniqueness), not “positive coefficients” and not
   “real on the line.”

RH is OPEN.

## The figures

All produced by `r25_viz.py` (matplotlib + mpmath `fp`, no zeros as proof inputs).

* `fig1_arrow_walk.png` — Partial sums of $\eta(s)=\sum(-1)^{n-1}n^{-s}$ in the
  complex plane, at the first two zeros, at Gram $g_0$, and at the same height off
  the line. The true value $\eta(s)$ is the gold star.
* `fig2_domain_coloring.png` — Wegert portrait of $\xi$ in the strip, versus the
  completed Epstein function. Colour-wheels on the line for $\zeta$; off-line
  colour-wheels (red rings) for $x^2+5y^2$. Conjugate symmetry across
  $\sigma=\tfrac12$ is visible in both.
* `fig3_hardy_Z_gram.png` — The real wave $Z(t)$, its zeros and extrema, Gram
  points $g_n$ with $\vartheta(g_n)=n\pi$, and $\vartheta(t)/\pi$.
* `fig4_level_curves_zeta.png` — The central picture: $\mathrm{Re}\,\xi=0$ (red)
  and $\mathrm{Im}\,\xi=0$ (blue) in the strip, with $Z(t)$ underneath. Gold
  zeros sit on red/blue crossings *on the line*. Violet marks are extrema of
  $\xi(\tfrac12+it)$, and the extra blue branches meet the line there — visibly
  shifted from the peaks of $Z$.
* `fig5_level_curves_fakes.png` — The same colour language for $\zeta$, $F_\theta$,
  Epstein and DH. This is the figure that shows the lining-up *mechanism* and its
  failure for the fakes, side by side.
* `fig6_cosine_cancellation.png` — Riemann–Siegel cosine contributions at a first
  zero, a sixth zero, a later zero ($t=143$), and at $g_0$. Every prefactor
  $n^{-1/2}$ is positive.

## Files

* `r25_viz.py` — generators (evaluators, grids cached locally under
  `/workspace/r25viz/data/`, not committed).
* The six PNGs listed above.
* `data/meta.json` — Gram points, $\xi$-extrema, and the plot marks used.

References: Gram (1903) · Hutchinson (1925) · Titchmarsh, *The Theory of the
Riemann Zeta-Function* · Spira, “The Riemann zeta function near the line
$\sigma=1$” and later level-curve work · Gabcke, dissertation on the
Riemann–Siegel formula · Arias-de-Reyna, “High precision computation of
Riemann’s zeta function by the Riemann–Siegel formula” · Edwards, *Riemann’s
Zeta Function* · Wegert, *Visual Complex Functions*.

RH: OPEN.
