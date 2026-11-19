# Round 23: attacking GAP-22, the Fejér-network inequality

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).

Labels: PROVED, CHECKED, NUMERICAL, CONJECTURE, HEURISTIC, GAP, OPEN. Zeta zeros are never used as proof inputs. Below they appear
only in illustration and validation rows, and those rows are marked. RH is **not** closed here.

## Plain-words summary

* **GAP-22 reduces to one variable.** The inequality for all (X, γ) holds iff RH holds. It also holds iff its γ = 0 slice
  holds for all large X (Theorem 3, PROVED, Landau's classical argument). So the whole content of GAP-22 sits in one
  explicit real function D(L), with L = log X, built only from the integers: D(L) = 4(e^{L/2} − 2 + e^{−L/2}) + c₁L + C₂ − T(L) − Σ_{n≤e^L} Λ(n) n^{−½}(L − log n).
* **Best partial result (route b): PROVED by computer, in ball arithmetic.** D(L) > 0.0214 for every real X ∈ [2, 10⁸].
  For X ∈ [2, 10⁷] the bound is D(L) > 0.0272. This is the γ = 0 slice of GAP-22 for every X up to 10⁸, including non-integer X.
  The proof covers 5,762,859 gaps between consecutive prime powers, with no zeros used.
* **Uniform margin (route c): impossible, PROVED.** Unconditionally, liminf_{L→∞} D(L) ≤ 0. So no positive lower bound that is
  uniform in X can exist. The ~1/log X decay of the R22 Fejér margin ρ_X is only normalisation (ρ_X(0) = 2D/(πL)). Under
  RH, D has mean exactly 2c₀ = ξ''(½)/ξ(½) = 0.046210, and it returns arbitrarily close to 0 infinitely often (Dirichlet).
* **Positive identity (route a):** D is the Weil quadratic form evaluated on the box 1_{[0,L]}. The Gram matrix
  G(a,b) = D(a) + D(b) − D(|a−b|) is Weil's Gram matrix on box functions, and it is PSD on every finite set iff RH
  (CHECKED, standard argument). In tests, ζ's Gram matrices are PSD (min eigenvalue > 0.003 on 230×230). All three
  controls are indefinite, DH included, even though DH's scalar slice is still positive at 10⁵. An integer SOS identity
  for G is exactly Weil positivity on boxes: OPEN.
* **de Branges below ½ (route d):** pushing E_h = Ξ(z+ih) below h = ½ is equivalent to D(L) ≥ −C e^{εL} for some ε < ½
  (Theorem 4, PROVED). That in turn is a zero-free vertical strip, the quasi-RH. The Euler product gives exactly ε = ½ and
  nothing below it. OPEN.
* **Next sharp gap (GAP-23):** prove, from integers only, that D(L) ≥ −C e^{εL} for **one** ε < ½. Equivalently, prove that
  Weil's box-Gram matrix has a lower bound of the form G ≥ −C e^{ε·max} with ε < ½. The ε = 0 case is RH.

## 0. Notation

L = log X > 0. Λ is von Mangoldt. b_k = 2k + ½. c₁ = (ψ(¼) − log π)/2 = (−γ_E − π/2 − 3 log 2 − log π)/2 = −2.686091709613,
C₂ = ψ'(¼)/4 = (π² + 8G)/4 = 4.299332288627, where G is Catalan's constant.

* T(L) = Σ_{k≥0} e^{−b_k L}/b_k²
* S(L) = Σ_{n≤e^L} Λ(n) n^{−½}(L − log n)
* A(L) = 4(e^{L/2} − 2 + e^{−L/2}) + c₁L + C₂ − T(L)
* D(L) = A(L) − S(L)

The R22 Fejér network is ρ_X(γ₀) = Σ_ρ K(γ − γ₀), with K(r) = (L/2π) sinc²(rL/2). Its prime side is
(2/π) Σ Λ(n) n^{−½}(1 − log n/L) cos(γ₀ log n). It satisfies ρ_X(0) = 2D(L)/(πL).

## 1. The γ = 0 slice in closed form

**Theorem 1 (PROVED).** With D(L) := (πL/2) ρ_{e^L}(0) as in R22, D has the closed form above for all L > 0, and D(0⁺) = 0.

**Theorem 2 (PROVED), the Laplace identity.** For Re σ > ½,
∫₀^∞ e^{−σL} D(L) dL = (ξ'/ξ)(½+σ) / σ².
*Proof.* Write ξ'/ξ(s) = 1/s + 1/(s−1) − ½ log π + ½ ψ(s/2) + ζ'/ζ(s) at s = ½ + σ and divide by σ². The pieces transform
term by term:
* 1/(σ+½) + 1/(σ−½), divided by σ², equals 2/(σ(σ² − ¼)). That is the Laplace transform of 4(e^{L/2} − 2 + e^{−L/2}).
* −Σ Λ(n) n^{−½−σ}/σ² is the transform of −S, since ∫ e^{−σL}(L − a)₊ dL = e^{−σa}/σ².
* The ψ-series ½ψ(¼ + σ/2) = ½ψ(¼) − Σ_k [1/(k + ¼ + σ/2) − 1/(k + ¼)]/2 gives c₁L + C₂ − T(L), after partial fractions in σ.
  The 1/σ² and 1/σ coefficients produce c₁ and C₂.

Validation: the closed form agrees with the R22 explicit-formula code to about 1e−13 at L = 0.8, 1.5, 2.3, 4, 6, 9
(`r23_slice.json`, CHECKED). The zero side, using 400 zeros (illustration only), sits 0.00266 below at every L. That matches
the missing tail Σ_{γ>γ₄₀₀} 2/γ² ≈ (1/π)(log(γ₄₀₀/2π) + 1)/γ₄₀₀ ≈ 0.0028.

## 2. GAP-22 is equivalent to its γ = 0 slice: the Landau reduction

**Theorem 3 (PROVED, classical argument).** The following are equivalent:
(i) RH;
(ii) GAP-22: ρ_X(γ₀) ≥ 0 for all X ≥ 2 and γ₀ ≥ 0;
(iii) D(L) ≥ 0 for all L > 0;
(iv) there is L₀ with D(L) ≥ 0 for all L ≥ L₀;
(v) D is bounded below on [0, ∞).

*Proof.*
* (i)⇒(ii): under RH, ρ_X(γ₀) = Σ_ρ K(γ − γ₀) with K ≥ 0.
* (ii)⇒(iii)⇒(iv)⇒(v): trivial (D ≥ 0 on [0, L₀] is not needed for (v), because D is continuous).
* (v)⇒(i): let D ≥ −C. Then f = D + C ≥ 0, and its Laplace transform is F(σ) = (ξ'/ξ)(½+σ)/σ² + C/σ. By Landau's theorem
  (Pringsheim–Landau for Laplace transforms of non-negative functions), the real point σ_c of the abscissa of convergence
  is a singularity of F. On the real axis σ > 0, F is analytic, because ξ has no zeros on (½, ∞):
  - on (1, ∞) the Euler product rules them out;
  - ξ(1) = ½ ≠ 0;
  - on (½, 1), ζ(σ) = η(σ)/(1 − 2^{1−σ}) < 0, since η(σ) > 0 (alternating series with decreasing terms).

  Hence σ_c ≤ 0, the integral converges, and F is analytic on Re σ > 0. So ξ'/ξ(½ + σ) is analytic on Re σ > 0, and ξ has
  no zeros with Re s > ½. By the functional equation, RH follows. ∎

So the whole difficulty of GAP-22 is the single real function D, which comes from integers. (The control models need
their own real-axis check, ξ_M(σ) ≠ 0 on (½, ∞). They are not needed below because the controls fail numerically.)

**Theorem 4 (PROVED), the strip version, which is route (d)'s exact form.** Fix 0 ≤ ε < ½. Then
ξ(s) ≠ 0 for Re s > ½ + ε **iff** for every ε' > ε there is C with D(L) ≥ −C e^{ε'L} for all L ≥ 0.
*Proof.*
* (⇒) Use the explicit formula D(L) = Σ_ρ [e^{(ρ−½)L} − 1]/(ρ−½)², summed over all zeros (the inverse Laplace transform of Σ_ρ 1/((σ − (ρ−½))σ²); the L/(ρ−½) terms cancel in ρ ↔ 1−ρ pairs). The sum converges absolutely,
  since Σ 1/|ρ|² < ∞. Each term is O(e^{εL}).
* (⇐) Apply Landau to D + C e^{ε'L}, whose extra transform C/(σ − ε') is analytic for σ > ε'. Then
  (ξ'/ξ)(½ + σ) is analytic for Re σ > ε'. ∎

The unconditional exponent is ε = ½: all β < 1 gives D = o(e^{L/2}). That is exactly the R22 de Branges frontier h = ½,
the R20 Euler-product region u > ¼, and the R21 β > 1. By R22 Theorem D, E_h is Hermite–Biehler iff h > sup|β − ½|. So
"push de Branges from ½ toward 0" **is** "improve the exponent ε in D ≥ −C e^{εL}".

## 3. Route (b): computer-assisted proof for X ≤ 10⁸ (PROVED, γ = 0 slice)

**Lemma 5a (convexity, PROVED by hand).** A''(L) = e^{L/2} + e^{−L/2} − Σ_k e^{−b_k L} = e^{L/2} + e^{−L/2} − e^{−L/2}/(1 − e^{−2L}).
This is > 0 for L ≥ log 2. Indeed, A'' > 0 ⇔ (e^L + 1)(1 − e^{−2L}) > 1, and the left side is increasing, with value 9/4 at log 2.

**Lemma 5b (PROVED).** On each gap [log q_j, log q_{j+1}] between consecutive prime powers, S(L) = α_j L − β_j, where
α_j = Σ_{n≤q_j} Λ(n) n^{−½} and β_j = Σ_{n≤q_j} Λ(n) n^{−½} log n. S is continuous. Hence D is convex on each gap. The
right slope drops by Λ(q)/√q at each prime power.

**Lemma 5c (two-tangent bound, PROVED).** Let f be convex on [a, b] with f'(a⁺) < 0 < f'(b⁻). Then
min f ≥ f(a) + f'(a⁺)(L* − a), where L* is where the two endpoint tangents cross.

At L = log q everything is algebraic in q: e^{L/2} = √q, and T(log q) = Σ q^{−b_k}/b_k² with a rigorous geometric tail.

**Theorem 5 (PROVED, computer-assisted; `r23_interval.py`, python-flint `arb` balls, 160 bits).** D(L) > 0 for all
L ∈ [log 2, log 10⁸]. Equivalently, for every real X ∈ [2, 10⁸], the γ = 0 slice of GAP-22 holds:
Γ-term + pole-term ≥ (2/π) Σ_{n≤X} Λ(n) n^{−½}(1 − log n/log X).

Certified lower bounds:

| range of X | certified min D | where |
|---|---|---|
| [2, 10] | 0.03261 | |
| [2, 10²] | 0.02979 | |
| [2, 10³] | 0.02802 | |
| [2, 10⁷] | 0.02723 | gap ending at q = 3109 (float min 0.027521 at X ≈ 3100) |
| [2, 10⁸] | 0.02150 | gap ending at q = 34,186,499 |

The run checked 5,762,859 gaps, 136 of them subdivided (the wide gaps at small q). It took 196 s. The inputs are the
integer sieve, plus γ_E, π, log 2 and Catalan's G as `arb` constants. JSON: `r23_interval_X1e05/7/8.json`.

Scope, said honestly:
* This is the γ = 0 slice only. By Theorem 3 the slice carries all of RH, but only **as L → ∞**. A finite-range check
  proves nothing asymptotic.
* For γ ≠ 0 at finite X there is no proof here. R22's numerical γ-sweep to X = 10⁶ (`r22_fejer_ZETA_ext.json`) remains
  NUMERICAL.
* The trust base is the arb library and the sieve.

## 4. Route (c): the margin, and why it decays

The NUMERICAL scan (step 0.0005 in L, X ∈ [2, 10⁷], `r23_slice.json`) gives the minimum of D in each decade:

| decade | min D |
|---|---|
| [1, 10) | 0.03262 |
| [10, 10²) | 0.02979 |
| [10², 10³) | 0.02802 |
| [10³, 10⁴) | 0.02752 |
| [10⁴, 10⁵) | 0.02965 |
| [10⁵, 10⁶) | 0.02806 |
| [10⁶, 10⁷) | 0.02783 |

The mean is 0.046150. Under RH it should be 2c₀ = ξ''(½)/ξ(½) = 0.046210: D = Σ_ρ (1 − cos γL)/γ² averages to Σ 1/γ².

**Theorem 6 (PROVED).** Unconditionally, liminf_{L→∞} D(L) ≤ 0. Under RH, liminf_{L→∞} D(L) = 0.
*Proof.* If RH fails, D is unbounded below by Theorem 3. Under RH, D = Σ (1 − cos γL)/γ² (all terms ≥ 0, summable).
Pick N so that the tail beyond the first N ordinates is < ε/2. Dirichlet's simultaneous approximation theorem gives
arbitrarily large L with every γ_k L within δ of 2πℤ for k ≤ N, so the head is < ε/2. ∎

So **no uniform-in-X positive margin exists**, and the observed minima of about 0.027–0.033 must eventually be undercut. In
ρ-normalisation the margin is ρ_X(0) = 2D/(πL). The ~1/log X decay seen in R22 is this factor, not a weakening.
The physicist reading: D is a sum of independent oscillators (1 − cos γL)/γ², each ≥ 0. Positivity of each oscillator is
the content; the minimum floor is not. HEURISTIC: Dirichlet's L for N oscillators is of size ~δ^{−N}, so
the liminf 0 is approached only at astronomically large X; the finite-range minima ≈ 0.027 are typical values, not near the limit.

## 5. Route (a): a positive identity, the Weil box-Gram matrix

The box f_L = 1_{[0,L]} has f̂_L(γ) = (e^{iγL} − 1)/(iγ) and |f̂_L(γ)|² = 2(1 − cos γL)/γ². So D(L) = ½ Σ_ρ |f̂_L(γ)|² is the
**Weil quadratic form of the box**. Its polarisation is G(a,b) = D(a) + D(b) − D(|a−b|) = Σ_ρ Re f̂_a(γ) conj(f̂_b(γ))
(with the trig identity 1 − cos x − cos y + cos(x−y) = Re[(1 − e^{ix})(1 − e^{−iy})]).

**Theorem 7 (CHECKED, standard).** G is PSD on every finite set of boxes **iff** RH.
* (⇒) Under RH, G is a sum of Hermitian rank-one terms.
* (⇐) G PSD means D is of negative-definite type with D(0) = 0. By Schoenberg/Lévy–Khintchine, D(L) = ∫(1 − cos tL)/t² dμ(t)
  + cL², with μ ≥ 0. So D has at most polynomial growth, hence D ≥ 0, and RH follows by Theorem 3.

NUMERICAL test (`r23_gram.py`, grids L_i = i·h up to log 10⁵):

| model | h = 0.25 | h = 0.1 | h = 0.05 | first indefinite block (h = 0.05) |
|---|---|---|---|---|
| ζ | min eig +0.01161 | +0.00436 | +0.00322 | none (230×230 PSD) |
| FTH | −1.208 | −3.348 | −6.874 | L = 2.75 |
| EPS | −5.788 | −14.91 | −30.49 | L = 3.40 |
| DH | +0.01240 | −0.01092 | −0.08859 | L = 7.7 (X ≈ 2200) |

The Gram test detects DH by X ≈ 2200, while DH's scalar slice D_DH is still positive at 10⁵ (§7). Cross-terms see the
off-line zeros before the diagonal does.

OPEN: an **integer sum-of-squares** certificate G = Σ (integer data)(integer data)* would be Weil positivity on step
functions, which is RH. The finite-grid PSD check is NUMERICAL, and the controls show it is a genuine discriminator.

## 6. Route (d): de Branges below ½ (OPEN, precisely located)

By Theorem 4 together with R22 Theorem D: E_h is HB for some h < ½ ⇔ zero-free strip Re s > 1 − δ for some δ > 0
(δ = ½ − h) ⇔ D ≥ −C e^{(½−δ)L}. The Euler product gives absolute convergence only on Re s > 1, so it reaches exactly h = ½.
The known zero-free regions (de la Vallée Poussin, c/log t; Vinogradov–Korobov, c/(log t)^{2/3}(log log t)^{1/3}) have
width → 0 as t → ∞, so they give **no** h < ½. No partial step below ½ was obtained. This is the classical quasi-RH
problem, OPEN.

## 7. Controls: each candidate must fail for F_θ, Epstein, DH

| test | ζ | FTH (θ = ¼) | EPS (x² + 5y²) | DH |
|---|---|---|---|---|
| scalar slice D_M ≥ 0 (NUMERICAL, X ≤ 10⁵) | passes (PROVED to 10⁸) | fails at X ≈ 2943, min −0.160 | fails at X ≈ 4533, min −0.110 | positive to 10⁵ (min 0.0726) |
| box-Gram PSD (NUMERICAL) | passes | fails, L = 2.75 | fails, L = 3.4 | fails, L = 7.7 |
| full γ-sweep ρ_X(γ₀) (R22) | passes to 10⁶ | fails X = 30 | fails X = 100 | fails X = 100 |
| de Branges h-scan (R22) | no violation to h = 0.05 | violation at h ∈ (0.20, 0.25] | (0.40, 0.45] | (0.30, 0.35] |

DH's slice must eventually go negative. By Theorem 3, it carries off-line zeros. Illustration only, using DH's known zero
0.8085 + 85.70i: its term 2e^{0.3085L}/85.7² overtakes the 0.07 floor near L ≈ 18, X ~ 10⁸. That is HEURISTIC and order of
magnitude. The γ = 0 slice is therefore a *weak finite-X* detector for zeros with large γ, while the Gram matrix and the
γ-sweep are strong ones. All four candidate statements fail for all three controls, given enough X.

## 8. Ledger

* PROVED:
  - Theorem 1 (closed form) and Theorem 2 (Laplace identity);
  - Theorem 3 (GAP-22 ⇔ slice ⇔ RH);
  - Theorem 4 (strip ⇔ exponential lower envelope);
  - Lemmas 5a–c;
  - Theorem 5 (computer-assisted: D > 0.0214 on X ∈ [2, 10⁸]);
  - Theorem 6 (no uniform positive margin; liminf ≤ 0).
* CHECKED: Theorem 7 (box-Gram PSD ⇔ RH, standard Schoenberg argument); closed form = R22 code to 1e−13.
* NUMERICAL: decade minima, mean ≈ 2c₀, Gram eigenvalues, the controls' failure points.
* HEURISTIC: the DH crossover near X ~ 10⁸; the LI picture of when liminf is approached.
* OPEN / **GAP-23**: from integers alone, prove D(L) ≥ −C e^{εL} for one ε < ½. Equivalently: a uniform zero-free strip;
  de Branges at some h < ½; Weil box-Gram positivity up to an e^{εL} defect. ε = 0 is RH.

## Files

* `r23_slice.py`, `r23_slice.json`, `slice.log`, `r23_D_zeta.npy`: closed form, cross-checks, ζ scan to 10⁷, controls' scalar slice.
* `r23_interval.py`, `r23_interval_X1e05.json`, `r23_interval_X1e07.json`, `r23_interval_X1e08.json`, `interval_*.log`: Theorem 5.
* `r23_gram.py`, `r23_gram.json`, `gram.log`: Theorem 7 test.
* `r22_models.py`, `r22_fejer.py`: copied from round22_synthesis, unchanged.
* `FORMALIZE.md`: formalizable lemmas, as a handoff for the verify/Lean track.
