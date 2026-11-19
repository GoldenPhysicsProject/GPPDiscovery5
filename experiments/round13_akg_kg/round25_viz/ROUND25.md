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
* `fig7_mirror_pairs.png` — Round 25b. Same four windows as Fig. 5, with
  $\mathrm{Re}\,\xi=0$ split by side of the line, plus a zoom of the first
  $F_\theta$ fold. See the Round 25b section below.
* `fig8_weil_vs_hook.png` — Round 25b. Fejér–Weil form $W_X(t)$ scanned across
  hook heights and ordinary crossings. See the Round 25b section below.

## Files

* `r25_viz.py` — generators (evaluators, grids cached locally under
  `/workspace/r25viz/data/`, not committed).
* The original six PNGs, plus `fig7_mirror_pairs.png` and `fig8_weil_vs_hook.png`.
* `r25b_viz.py` — Round 25b generators (Fig. 7 from the cached grids; Fig. 8 from the R22 explicit formula).
* `r22_models.py` — copy of the Round 22 L-data (von Mangoldt coefficients, archimedean factors, poles), used by Fig. 8.
* `data/meta.json` — Gram points, $\xi$-extrema, and the plot marks used. The two Fig. 5 violet marks are the $\xi$-extrema at $t\approx 15.5857$ and $22.0980$.
* `data/r25b_weil_scan.json` — numerical summary of the Fig. 8 scans.

References: Gram (1903) · Hutchinson (1925) · Titchmarsh, *The Theory of the
Riemann Zeta-Function* · Spira, “The Riemann zeta function near the line
$\sigma=1$” and later level-curve work · Gabcke, dissertation on the
Riemann–Siegel formula · Arias-de-Reyna, “High precision computation of
Riemann’s zeta function by the Riemann–Siegel formula” · Edwards, *Riemann’s
Zeta Function* · Wegert, *Visual Complex Functions*.

RH: OPEN.

## Round 25b — crayons, and a Weil-form test of the hook

Turbo asked for two things: name the two unlabelled violet marks on the $\zeta$
panel of Fig. 5, and redraw the level curves with each red curve and its mirror
in different shades. After Fig. 7 was in hand, she also asked whether a hook's
birth lines up with a sign change of the Weil form $W(g\ast\widetilde g)$.

### The violet marks (Fig. 5, $\zeta$ panel)

They are **extrema of $\xi(\tfrac12+it)$**, the points where an extra blue
branch of $\mathrm{Im}\,\xi=0$ meets the critical line (Cauchy–Riemann). In the
Fig. 5 window they sit at

$$
t = 15.585708589924,\qquad t = 22.097977280453.
$$

(CHECKED, recomputed from $\frac{d}{dt}\xi(\tfrac12+it)=0$.) They are *not*
Gram points ($g_0\approx 17.846$, $g_1\approx 23.170$) and they are *not* the
extrema of Hardy's $Z(t)$ (first $Z$-peak $\approx 17.883$). The shift is the
decaying archimedean factor in $\xi(\tfrac12+it)=A(t)Z(t)$, already noted in
the Round 25 summary, item 2.

They were plotted in Fig. 5 — `draw_level_panel` was called with
`extrema=xi_ext` — but they were left out of that figure's legend (Fig. 4's
legend had them). That is why they were hard to name from the picture alone.
Fig. 7 labels them.

### Fig. 7 — a curve meeting its mirror, or hooking away

Same four windows as Fig. 5, with $\mathrm{Re}\,\xi=0$ coloured coral on
$\sigma<\tfrac12$ and rose on $\sigma>\tfrac12$. A zoom of the first $F_\theta$
fold sits underneath.

The colour split is not decorative. For every model here the completed function
satisfies $\xi(1-\sigma+it)=\overline{\xi(\sigma+it)}$, so the two shades are
literal mirrors. This was checked on the plotting grids (NUMERICAL): the
maximum relative residual $|\xi(1-\sigma+it)-\overline{\xi(\sigma+it)}|/|\xi|$
is $2.3\times 10^{-12}$ ($\zeta$), $1.0\times 10^{-12}$ ($F_\theta$),
$1.1\times 10^{-12}$ (Epstein), $1.5\times 10^{-11}$ (DH), and
$5.5\times 10^{-12}$ on the zoom grid.

In $\zeta$ the two shades meet on the line, squarely, once at each simple zero.
In every fake they turn away from each other, and the tip of the hook catches
an off-line blue branch — the off-line zero. The $F_\theta$ zoom shows this
cleanest: a rose loop with its tip at $\tfrac34+i\gamma$ and its coral mirror
with its tip at $\tfrac14+i\gamma$, $\gamma\approx 14.135$. They do not meet on
the line. The other two members of the quartet, at $t=-\gamma$, are out of
view.

**Heuristic verdict.** Reading a hook as a broken pairing, and a square
crossing as a true inner product / reflection positivity, is **HEURISTIC**.
The drawing shows the geometry of an off-line zero — a red curve that fails to
meet its functional-equation mirror on the critical line — not the Weil form
itself. Both a hook and a failure of the Weil form are equivalent to an
off-line zero, but Fig. 7 does not exhibit the Weil pairing, and this note
does not derive one from the other. Fig. 8 tests one precise version of the
metaphor.

### Fig. 8 — does a hook's birth line up with a sign change of $W$?

Test function: a Fejér wave packet in $u=\log x$, Fourier support
$[-\log X,\log X]$, peaked at height $t$ on the critical line, main-lobe width
$2\pi/\log X$. The quantity plotted is the R22 Weil form $W_X(t)=W(g_t\ast\widetilde g_t)$,
evaluated from each model's own explicit formula (arithmetic side: $\Lambda_M(n)$,
archimedean factor, poles). No zeros are inputs. Off-line zeros appear only as
crimson plot marks, and as a NUMERICAL prime-side vs zero-side check for $\zeta$
and $F_\theta$ (200 zeta zeros: max abs residual $1.2\times 10^{-3}$).

Windows: the same ranges as Fig. 5. Widths $X=30,100,300,1000$ (lobe widths
$\approx 1.85,\,1.36,\,1.10,\,0.91$). $\zeta$ is the control.

**What we see (NUMERICAL).**

- $\zeta$: $W_X>0$ throughout. Valleys between on-line zeros stay positive
  (minima $\approx 0.0046$–$0.0084$ in this window). Peaks sit on the zeros
  and grow with $X$.
- The fakes do **not** go negative *at* the hook. At the hook height itself
  $W_X$ is a peak, and a large one: for $F_\theta$, $W_{100}(14.13)\approx 3.27$,
  $W_{1000}(14.13)\approx 5.62$; for Epstein, $W_{100}(15.668)\approx 4.21$,
  $W_{1000}\approx 8.91$; for DH, $W_{100}(85.699)\approx 3.59$,
  $W_{1000}\approx 6.38$.
- The sign change lives in the **shoulders**. Once the packet is narrow enough
  to resolve the hook, the valleys immediately beside it go negative, and the
  minimizer moves closer to the hook as $X$ grows:
  - $F_\theta$: already negative at $X=30$ (min $-0.025$ at $t=15.86$,
    distance $1.73$ from the first hook), deeper and closer at larger $X$
    (at $X=1000$, min $-0.518$ at $t=25.80$, distance $0.79$ from the nearest
    hook). Every $\zeta$-height in this window is a hook, so the negative
    valleys sit between hooks.
  - Epstein: still positive at $X=30$; negative from $X=100$
    (min $-0.374$ at $t=14.54$, distance $1.13$ from the first hook). At
    $X=1000$, min $-2.34$ at $t=14.98$, distance $0.69$. On-line zeros in the
    same window keep positive peaks.
  - DH: still positive at $X=30$; negative from $X=100$
    (min $-0.134$ at $t=84.48$, distance $1.22$ from the hook at $85.699$).
    At $X=1000$, min $-0.875$ at $t=84.94$, distance $0.76$. Both shoulders
    go negative; the on-line zeros on either side stay positive.

So the pictures line up with a characteristic *peak-plus-negative-shoulder*
signature at a hook, which ordinary on-line crossings do not have, and the
shoulders get sharper with narrower packets. They do **not** line up with
"min $W$ goes negative exactly at the hook." That is a numerical observation
for this Fejér family, not a theorem that a hook is a failure of reflection
positivity.

RH is OPEN.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).
