# Round 13f — de Branges / Lagarias test of the "projection curvature" R = ∂∂̄ log K(z,z)

Inspired by the Principia coherence functional R = ∇∇ log ∫K(x,y)² dμ(y). For a reproducing kernel, ∫K(z,y)K(y,z)dμ = K(z,z), so this is the Bergman/Kähler potential ∂∂̄ log K(z,z). Code: `debranges.py` (ζ, Epstein, spliced Beurling), `beur_zero_probe.py`. Data: `debranges_*.json`. Zeros were not used. They are test data only, and none were inputs here.

## 1. The space and the exact criterion
- **E = ξ(1/2 − iz) itself is degenerate.** The functional equation and real coefficients give E^# = E, so K ≡ 0.
- **The correct family** (de Branges 1986 for h = 1/2; Lagarias, *Zero spacing distributions for differenced L-functions*, Acta Arith. 120 (2005), §2 and §6, arXiv math/0601653) is
  - E_h(z) = ξ(1/2 + h − iz), h > 0, with s = 1/2 − iz, so Im z = σ − 1/2;
  - E_h^#(z) = ξ(1 + h − s) = ξ(s − h) (functional equation);
  - K_h(w,z) = (E_h(z)\overline{E_h(w)} − E_h^#(z)\overline{E_h^#(w)}) / (2πi(\bar w − z)).
- **Lagarias, Lemma 2.1:** E_h is Hermite–Biehler (HB) unconditionally for h ≥ 1/2. That input is ζ ≠ 0 on Re s ≥ 1, an Euler-product fact. Under RH, E_h is HB for every h > 0.
- **Converse.** If E_h is HB for all h ∈ (0, 1/2), then E_h has no zeros in ℂ⁺, so every zero has β ≤ 1/2 + h for every h, and RH follows. **So "E_h HB for all small h > 0" is RH-equivalent (CHECKED, classical).**
- **Precise positivity statements on ℂ⁺** (for E with no real zeros and |Θ| → 0 as Im z → ∞, which holds for E_h). Let Θ = E^#/E.
  - (a) HB ⇔ K(z,z) = (|E|² − |E^#|²)/(4π Im z) > 0 on all of ℂ⁺ ⇔ |Θ| < 1 on ℂ⁺ ⇔ K is a positive-definite kernel on ℂ⁺. For E_h this is D_h(s) := log|ξ(s+h)| − log|ξ(s−h)| > 0 on Re s > 1/2. (de Branges, *Hilbert Spaces of Entire Functions*, Thm 1–2. CHECKED, classical.)
  - (b) Exact identity on {|Θ| < 1}: **∂∂̄ log K(z,z) = 1/(4y²) − |Θ′|²/(1 − |Θ|²)².** The log|E|² term is harmonic and drops out. (PROVED, a one-line computation.)
  - (c) **HB ⇔ ∂∂̄ log K ≥ 0 on {z ∈ ℂ⁺ : K(z,z) > 0}.** (PROVED, short.)
    - ⇒: under HB, Θ maps ℂ⁺ into the disc, and Schwarz–Pick gives |Θ′|/(1−|Θ|²) ≤ 1/(2y).
    - ⇐: if HB fails, the set {|Θ| ≥ 1} ∩ ℂ⁺ is non-empty, while |Θ| < 1 near i∞. So the level curve |Θ| = 1 meets ℂ⁺. Approaching it from inside, (1 − |Θ|²) → 0 while Θ′ ≠ 0 except at isolated points, so the curvature → −∞.
    - Caveat: this needs E and E^# to have no common zero in ℂ⁺. For E_h, a common zero would require two zeros with the same ordinate exactly 2h apart.
- **Verdict on the Principia curvature: RIGHT as an equivalence, WRONG as a source.**
  - "Positive by construction" holds only once K is known to be a reproducing kernel, and that is HB, i.e. the RH-equivalent statement itself.
  - The curvature condition is Schwarz–Pick for Θ. It adds nothing beyond K(z,z) > 0 and is circular as a route.

## 2. Primes, pole, and arch only? Global?
- **Primes, pole, arch only: no.**
  - log Θ = log ξ(s−h) − log ξ(s+h) splits into three parts:
    - arch: log Γ_ℝ(s−h) − log Γ_ℝ(s+h);
    - pole: log[(s−h)(s−h−1)/((s+h)(s+h−1))];
    - primes: Σ Λ(n)/log n · (n^{−(s−h)} − n^{−(s+h)}).
  - The prime series converges only for σ > 1 + h. There the HB inequality is unconditionally true: by Hadamard, D_h = Σ_ρ log(|s+h−ρ|/|s−h−ρ|), and every term is positive when σ > β.
  - In the critical region 1/2 < σ ≤ 1 + h, the prime expression needs analytic continuation, which is the zero sum. (CHECKED.)
- **Global: yes, in the falsifier sense.** The condition at one point s involves ξ at s ± h, i.e. all primes, not a finite window. It is exactly a zero-location statement.

## 3. Runs on the three models
Grid: σ ∈ {0.51 … 1.2}, t ∈ [0.5, 60] (step 0.1 for D_h, 0.25 for curvature), h ∈ {0.05, 0.1, 0.2, 0.3, 0.45, 0.5}. SP = 2y|Θ′|/(1−|Θ|²) is the Schwarz–Pick ratio; curvature ≥ 0 ⇔ SP ≤ 1.
- **ζ.** D_h > 0 at every grid point, for every h. The minimum is at t = 0.5 and scales like 2·(σ−1/2)·∂_σ log|ξ|. SP_max = 0.99999 (→ 1 as y → 0, near t = 0.5). K(z,z) > 0 everywhere, and the curvature ≥ 0. (NUMERICAL.)
- **Epstein Q = x²+5y²** (its own Γ(s), conductor 20, and functional equation).
  - D_h < 0 at h = 0.05, 0.1, 0.2, 0.3, in t-clusters at 15.4–16, 30, 36.4, 44, 46.8, 52.7, 58.5. These are exactly the off-line zeros from Round 13e, and failure sets in iff β − 1/2 > h.
    - h = 0.1, σ = 0.7: D = −2.26 at t = 58.5.
    - h = 0.3, σ = 0.6: D = −1.32 at t = 30.0.
    - At σ = 0.8 and h = 0.3, D_h is positive again (largest β − 1/2 = 0.438).
  - h = 0.45 and 0.5: D_h > 0 on the whole grid, t ≤ 60.
    - Literature: Davenport–Heilbronn 1936 give Epstein zeros of class-number > 1 forms in σ > 1. So HB must fail even at h = 1/2 at larger t: the unconditional h ≥ 1/2 case fails without the Euler product. Not reproduced numerically: there is no σ > 1 zero for t ≤ 60.5 (argument-principle count in 13e).
  - Curvature: SP > 1 (negative curvature) next to every failure cluster, e.g. 85.8 at h = 0.1, σ = 0.8, t = 44; 14.0 at h = 0.3, σ = 0.6, t = 44. Away from them SP ≤ 1. This matches §1(c). (NUMERICAL.)
- **Spliced Beurling.** There is no functional equation, so E_h^# = ξ_B(1 + h − s) ≠ ξ_B(s − h). Its gamma factor is borrowed from ζ.
  - Validation: the Re s > 0 continuation used is accurate. log|ζ_B| = −26.3 and −28.1 exactly at 0.4 + 28i and 0.6 + 28i, with no other dips in t ∈ [26, 30].
  - D_h < 0 at **36–41% of all t**, for every h tested (0.05 → 0.3: 236 → 71 of 596 grid points at σ = 0.51). The failures are spread over the whole range, not just near the quartet; the deepest is at t ≈ 28.
  - SP up to 2.6e4.
  - The kernel cannot even be evaluated on all of ℂ⁺, because ξ_B(1+h−s) is needed at Re → −∞.
  - So without the functional equation, E_h is not a de Branges function at all, and the failure is not localized at zeros. (NUMERICAL + structural.)
- **Separation: yes.** The criterion separates all three where the finite-width terms could not separate ζ from the spliced Beurling model. But it does so because it is the zero-location statement (Hadamard form above), not because it brings in a new mechanism.
  - Epstein fails at its zeros.
  - Spliced Beurling fails structurally, through the missing functional equation.
  - ζ passes, which is only verified on the grid. Proving it is RH.

## 4. Literature obstructions
- **de Branges' positivity conditions** (1986, 1992, and the 1990s preprints) asked for more than HB: an operator positivity, Re⟨F, TF⟩ > 0 for a shift T on H(E), with E(z) = ξ(1 − iz) (h = 1/2). That would have forced RH.
- **Conrey–Li**, "A note on some positivity conditions related to zeta and L-functions", IMRN 2000(18) 929–940, arXiv math/9812166, showed it fails for ζ:
  - −Re ξ′(ρ)ξ(1+ρ) = −5.39·10⁻⁶⁹ < 0 at the 34th zero, ρ = 1/2 + 111.0295i;
  - for the space F(W), W = 1/ξ(1−iz): Re ξ(1+282i)/ξ(2+282i) = −1.3·10⁻⁴ < 0.
  - The analogous condition fails for L(s, χ₄).
- **So:**
  - the HB / curvature level is true under RH but equivalent to it;
  - the stronger de Branges positivity level, which would have been a route, is false for ζ itself.
  - The Principia curvature is the first level.

## 5. Labels
- §1(a): CHECKED (classical).
- §1(b), (c): PROVED (short).
- "E_h HB ∀h ∈ (0, 1/2) ⇔ RH": CHECKED (Lagarias + Hurwitz/zero location).
- §2: CHECKED.
- §3: NUMERICAL.
- Davenport–Heilbronn at h = 1/2: LITERATURE, not reproduced.
- RH: OPEN.

**Judgement:** the curvature route is right as a reformulation and wrong as a source of positivity. It is RH-equivalent, circular as "positive by construction", and the next, stronger de Branges step is refuted (Conrey–Li). It does fit the 13e sentence. The h ≥ 1/2 half needs the Euler product (Epstein fails it via Davenport–Heilbronn), and the construction itself needs the functional equation (Beurling fails it everywhere).

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5). Audit: Turbo.
