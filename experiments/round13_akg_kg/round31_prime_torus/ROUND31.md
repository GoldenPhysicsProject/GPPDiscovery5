# Round 31 — the prime torus: rung 0, the integer floor, and the lift

Approved by Daniel Toupin in his own words: "Go, start the prime-torus round" (2026-11-20, 11:40 PM Toronto).
Seat: Heart of the Cosmos / Grok (D5).

Rules for this round:
- No zeta zeros are inputs anywhere. The only zero used at all is the known Davenport–Heilbronn off-line zero 0.808+85.70i (Spira; Round 17). It only fixes where one control test is placed, and no argument uses it.
- Weil positivity is never assumed. **RH is not closed.**
- Corrections from the R30-b audits (Claude, bridge `28f9e65`; Turbo, `3594ba0`) are carried here:
  - F13S is one direction only (floor ⇒ RH, in Lean).
  - The converse, RH ⇒ D ≥ 0, is CLASSICAL and not in Lean (§2.0).
  - The floor is not "1×1-minor unitarity".
  - The intertwiner statement is pointwise only.

**Verdict.**
- Daniel's picture is exactly right as geometry, and most of it is classical. The positive integers are the lattice face (exponent vectors), and the prime torus T = ∏_p S¹ is the other face. The line t ↦ (p^{−it})_p winds densely through T because log p are ℚ-independent, which is unique factorization (PROVED: Bohr, Kronecker).
- "Turns, never stretches" has a precise classical form: **RH ⇔ ζ is strongly recurrent under the T-flow in ½ < σ < 1** (Bagchi 1981/1987). We located exactly where counted-once (a_p = 1) enters. It puts the **L² edge of the torus model on the mirror line**: the variance Σ|a_p|²p^{−2σ}, equivalently the mean square, is finite down to σ = ½. For F_θ (a_p = p^θ + p^{−θ}) that edge sits at ½+θ, its values visibly stretch along the flow (NUMERICAL: the predicted drift is matched to within 0.09 for heights ≥ 10³), and its zeros sit exactly on its stretched edge.
- None of the floor attacks proves the integer floor. Every torus-side tool (positivity, equidistribution, Haar averages, Chebyshev/sieve bounds) either also "works" for F_θ or stops at exponent ½. The floor lives at the base point t = 0 of the torus, which is a single, maximally resonant point. Ergodic and probabilistic statements see only typical points.
- The doubling D(z) = 2z carries the stage (mirror) and nothing arithmetic (PROVED: it is a change of variables). The natural integer lattices up Daniel's ladder (ℤ, ℤ[i], Hurwitz ℍ, E8 octaves) give ζ, ζ·L(χ₋₄), ζ(w)ζ(w−1) and ζ(w)ζ(w−3). The last two are **literally F_θ-type stretched products (θ = ½, 3/2) with all zeros off their lines, unconditionally**. So the higher rungs add no constraint; above ℂ they are counterexamples of the F_θ kind.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).

---

## 1. Rung 0: the prime torus

### 1.1 Objects (PROVED, classical)
- **Two faces.** Let Γ = ℚ_{>0}^× ≅ ⊕_p ℤ (free abelian on the primes, by unique factorization), and let T = Γ̂ = ∏_p S¹ be its Pontryagin dual (compact, with Haar measure the product of the circle measures). A positive rational q = ∏ p^{a_p} is the character χ_q(ω) = ∏ ω_p^{a_p} of T. So the integers are lattice points on the dual face Γ (exponent vectors), and T is the position face. This is Daniel's "two faces" in its exact form. Integers are not points of T; they are the frequencies on T.
- **The winding.** ι: ℝ → T, ι(t) = (p^{−it})_p, is a continuous homomorphism. It is dual to the inclusion Γ ↪ ℝ_{>0}, q ↦ q.
  - ι has dense image ⇔ the inclusion is injective ⇔ {log p} is ℚ-linearly independent ⇔ unique factorization (Kronecker–Weyl; Bohr 1910s).
  - The line never meets itself: ι is injective, since a period would make log p rationally dependent.
  - It equidistributes for Haar measure (Weyl), so time averages along the winding equal space averages on T.
  - "Each integer counted once" is the injectivity of n ↦ log n, i.e. of the projection of the lattice Γ onto the line direction (log 2, log 3, …) (fig. 1 right).
- **ζ as a function on T (Bohr).** For σ > 1, ζ(σ+it) = F_σ(ι(t)) with F_σ(ω) = ∏_p (1 − ω_p p^{−σ})^{−1}, an absolutely convergent continuous function on T. Bohr's lift sends Σ a_n n^{−s} to Σ a_n χ_n(ω). Then ‖·‖²_{L²(T)} = Σ|a_n|²n^{−2σ}, which is Carlson's mean value theorem; see Helson and Hedenmalm–Lindqvist–Seip 1997 for H²(T^∞).
- **The mirror.** On the space face it is x ↔ 1/x with half-density x^{−1/2}, the theta law Θ(t) = t^{−1/2}Θ(1/t) (Poisson). On the frequency face it is s ↔ 1−s, fixed line Re s = ½. This is Daniel's J₁ (paper, Thm "Unitary half-density inversion") and the d = 1 Mellin shadow. The mirror is not an automorphism of T: it is not a flow symmetry. It is the FE, which relates the edge σ to 1−σ (fig. 2 middle).

### 1.2 What "only turns" means on T (PROVED, classical: Bohr, Bohr–Jessen, Bagchi, Voronin)
- **Absolute-convergence (L^∞) region, σ > 1.** F_σ is continuous on T. ζ(s + iτ) is almost periodic, recurrent and zero-free (Bohr). This is a pure rotation of a continuous function.
- **L² region, ½ < σ < 1.** F_σ(ω) = ∏_p(1 − ω_p p^{−σ})^{−1} still converges for Haar-almost every ω and is a.s. zero-free: log F_σ is a sum of independent mean-zero terms of variance Σ_p p^{−2σ}/2 + …, finite iff σ > ½ (Kolmogorov; Bagchi). Bohr–Jessen (1930/32) and Bagchi (1981) showed that the law of τ ↦ ζ(s+iτ) on compacts of ½ < σ < 1 is exactly the law of the random product F_s(ω), ω Haar on T.
  - Voronin (1975): every non-vanishing analytic function on a small disc is approximated by translates ζ(s+iτ) with positive lower density (universality).
- **Bagchi's criterion** (B. Bagchi, PhD thesis ISI Calcutta 1981; Acta Math. Hungar. 50 (1987) 227–240; Steuding, LNM 1877, Thm 8.3): **RH ⇔ ζ is strongly recurrent in ½ < σ < 1.** That is, for every compact K with connected complement and every ε > 0, the set {τ : max_K |ζ(s+iτ) − ζ(s)| < ε} has positive lower density. Bohr (1922) gave the earlier analogue for non-principal L(s,χ): RH ⇔ almost periodicity in σ > ½.
  - (⇒) Universality applied to the target ζ itself, which is non-vanishing under RH.
  - (⇐) A zero ρ₀ with Re ρ₀ > ½ plus recurrence gives, by Rouché, ≫ T zeros near Re ρ₀, which contradicts the Bohr–Landau density theorem N(σ, T) = o(T).
  - This is the precise, attributed form of Daniel's "turns, never stretches": the flow returns ζ to itself, as a rotation of T does, on every compact inside the strip.

### 1.3 Where counted-once enters (PROVED, classical ingredients; the three-way location is our reading)
1. **Local: one radius per circle.** On circle p the Euler factor is (1 − α_p ω_p p^{−s})^{−1} with α_p = 1. This is a single "species" whose Satake parameter lies on the unit circle: a pure turn. F_θ has two species with α_p = p^{±θ}, i.e. two radii, a stretch by p^{±θ} (Round 30-b §5).
2. **The L² edge sits on the mirror.** The torus model lives (a.s. convergent, zero-free, Bagchi limit law) exactly where Σ_p |a_p|² p^{−2σ} < ∞. With a_p = 1 that is σ > ½, the mirror line. With a_p = p^θ + p^{−θ} it is σ > ½+θ (T4 below). **Counted-once is what makes the torus edge coincide with the FE line.** For a general Euler product this is the Ramanujan bound |α_p| = 1.
3. **Bagchi's (⇐) direction uses the zero-density theorem.** That theorem uses the mean square of ζ (L² edge at ½, item 2) and a mollifier ζ·Σμ(n)n^{−s} ≈ 1 whose coefficients μ(n) are bounded. That boundedness is multiplicativity plus |a_p| = 1 again, now on the inverse side.
   - **DH fails here, not at item 2.** Its coefficients are bounded, so its L² edge is ½, but 1/DH has no Euler product and DH has ≫ T zeros in strips of σ > ½ (Davenport–Heilbronn 1936; Voronin). So "recurrence ⇒ no zeros" has no proof for DH.

### 1.4 Where F_θ's stretch shows up on T
- (a) Two radii p^{±θ} on each circle (local, PROVED).
- (b) The L² edge is ½+θ, not ½. There is a gap (½, ½+θ] between the mirror and the torus region. F_θ's zeros are ρ ± θ.
  - Hardy's infinitely many on-line ζ zeros put infinitely many F_θ zeros **exactly on the edge ½+θ**, with mirror images at ½−θ (PROVED, as in R30-b).
  - Under RH all of them are there. So F_θ's zeros never enter its own torus region (that statement is RH for ζ), but they do sit off its mirror line.
- (c) **Drift along the flow (fig. 2 right).** For ½ < σ < ½+θ the factor ζ(σ−θ+it) is left of ½. By the FE, |ζ(σ−θ+it)| = |χ(σ−θ+it)|·|ζ(1−σ+θ−it)| with |χ| ≈ (t/2π)^{½+θ−σ} (R30-b A5). The mean of log|F| grows like (½+θ−σ)·log(t/2π): the values stretch instead of turning, and no stationary torus model exists there (PROVED via FE + Bagchi limit law for the reflected factor).
  - NUMERICAL T3 at σ = 0.6: means 0.74, 0.82, 1.20, 1.47 at H = 10², 10³, 10⁴, 10⁵, against the predicted drift 0.42, 0.77, 1.11, 1.46.
  - ζ at 0.6 stays flat: 0.17, −0.01, −0.02, −0.02.
  - F_θ at 0.85 (inside its own strip) is also flat: 0.24, −0.02, −0.01, 0.03.
  - DH at 0.6 is flat: −0.02, 0.03, −0.00, 0.05. So the drift test separates F_θ, not DH, consistent with §1.3.

### 1.5 Numerics for rung 0 (`r31_torus.py` → `r31_torus.json`)
- **T1, how "golden" the winding is (CHECKED).**
  - log3/log2 = [1;1,1,2,2,3,1,5,2,23,2,2,1,1,55,1,4,3,…]. Golden ratio = [1;1,1,1,…].
  - The large partial quotients 23 and 55 are near-closures. Record approaches |a log2 − b log3|:
    - (19,12): 0.01355, the Pythagorean comma 3¹²≈2¹⁹;
    - (84,53): 0.00209, Mercator's comma, 53-tone tuning;
    - (485,306): 0.00102;
    - (1054,665): 4.4·10⁻⁵.
  - log5/log2 = [2;3,9,2,2,4,6,2,1,1,3,1,18,…].
  - So the prime winding is irrational but **not** golden. It nearly closes at the musical commas and never exactly. (Daniel's φ intuition: φ is the extreme of "never nearly closes"; the primes are generic, not extreme. HEURISTIC remark.)
  - Box discrepancy on T² (8×8 boxes, 2·10⁵ samples): 3.6 at T = 10, 0.42 at T = 10², then the sampling-noise floor ≈ 0.05 from T = 10³ on. On T³: 24 → 2.4 → 0.27 → 0.20 (noise floor) (NUMERICAL; equidistribution is Weyl, PROVED).
- **T2, Bagchi limit law (NUMERICAL).** Time samples log|ζ(σ+it)|, 4000 uniform t ∈ [10⁴, 3·10⁴], against the random Euler product on T (primes ≤ 2·10⁵, 4000 Haar samples):
  - σ = 0.6: mean −0.012 vs 0.015, var 0.821 vs 0.797, KS = 0.023 (p = 0.23);
  - σ = 0.75: var 0.443 vs 0.453, KS = 0.0185 (p = 0.50);
  - quantiles agree to ≤ 0.04.
- **T4, the L² edge (CHECKED / NUMERICAL).**
  - Σ_{p≤P} p^{−2σ}: at σ = 0.55 it creeps up (1.51 → 1.96 for P = 10² → 2·10⁵, convergent but slow); at σ ≥ 0.7 it has converged.
  - Σ_{p≤P}(p^θ+p^{−θ})²p^{−2σ}: 8.0 → 41.2 at 0.55, 6.6 → 19.5 at 0.6, 4.6 → 6.8 at 0.7. It diverges for σ ≤ ½+θ = 0.75 and converges above.
  - Sample variance of the truncated random model at σ = 0.6: ζ 0.70, 0.77, 0.78, 0.78 (P = 10², …, 10⁵), converged; F_θ 3.7, 4.9, 6.5, 9.3, diverging.
- **T5, strong recurrence illustration (NUMERICAL, weak).** Density of τ ≤ T with max over the circle |s − s₀| = 0.05 of |Z(s+iτ) − Z(s)| < 0.25:
  | model, centre s₀ | T = 500 | 1000 | 2000 | 4000 |
  |---|---|---|---|---|
  | ζ, 0.75+30i | 0.0254 | 0.0249 | 0.0240 | 0.0242 |
  | ζ, 0.6+30i | 0.0160 | 0.0136 | 0.0119 | 0.0107 |
  | F_θ, 0.6+30i (in its gap) | 0.0110 | 0.0092 | 0.0068 | 0.0055 |
  | F_θ, 0.9+30i (in its strip) | 0.0166 | 0.0172 | 0.0167 | 0.0169 |
  | DH, 0.75+30i | 0.0548 | 0.0549 | 0.0521 | 0.0489 |
  | **DH, disk around its off-line zero 0.808+85.70i** | 0.0060 | 0.0062 | 0.0069 | 0.0071 |

  - **DH control (the instructive row).** DH returns to itself at a stable or increasing rate on a disk that contains one of its own off-line zeros. So "recurrent" and "has a zero right of ½" coexist for DH, as §1.3 (3) predicts: Bagchi's (⇐) needs the density theorem, which needs the Euler product and a_p = 1 on the inverse side.
  The fixed-K densities converge slowly, and ζ at 0.6 also still decreases at T ≤ 4000. **Recurrence numerics alone do not separate the models at this range; the drift test (T3) does.** At the base point s₀ = 0.75 (real, t = 0), with ε = 0.5, the density is **0** (no return found for τ ≤ 4000, step 0.02). This is the resonance of the identity point of T (§2, attempt F).

## 2. The integer floor: attempts and their exact stopping points

### 2.0 Setting (Lean + classical)
- D(L) = GppHeartLaplace.Dfun = 4(e^{L/2} − 2 + e^{−L/2}) + c₁L + C₂ − T(L) − S(L), with S(L) = Σ_{n≤e^L} Λ(n)n^{−½}(L − log n) (GPPVerify main `c81c031`).
- F13S (Lean, standard axioms): if for some B and arbitrarily large K, D(4 log k) ≥ −B for all k ∈ [K, ⌈K^{3/2}⌉], then ξ, ζ ≠ 0 on Re s > ½. **One direction.**
- Converse (CLASSICAL, not in Lean): under RH, D(L) = Σ_ρ (1 − cos γL)/γ² ≥ 0. This follows from the explicit formula with the Fejér pair (Round 22/23 derivation; the Lean Laplace identity ∫e^{−σL}D = ξ'/ξ(½+σ)/σ² plus Hadamard). Turbo's point stands that the Lean D ≥ 0 under RH is not formalized.
- The integer samples carry no arithmetic. `Dfun_chord` + `cell_penalty` (Lean, F13S) show the integer floor and the real-L floor differ by ≤ 8 (PROVED). Numerically, ζ: integer min 0.0293 (k = 77), real-L min 0.0216 (L = 17.35); F_θ: −1.18 vs −1.78.

**Numerics (`r31_floor.py` → `r31_floor.json`; NUMERICAL, from Λ only).**
- ζ and F_θ use X ≤ 10⁸ (k ≤ 100); EPS and DH use X ≤ 10⁷.
- CHECK: the Lean closed form and the R22 explicit-formula route agree to ≤ 6·10⁻¹⁴ at L = 1, …, 15.

| model | min D(4 log k) | last window [K, K^{3/2}] | its min | fraction of k in window with D ≥ 0 |
|---|---|---|---|---|
| ZETA | +0.0293 (k=77) | [21, 97] | +0.0293 | 1.00 |
| FTH (θ=¼) | −1.182 (k=77) | [21, 97] | −1.182 | 0.51 |
| EPS (x²+5y²) | −4.604 (k=56) | [14, 53] | −4.057 | 0.68 |
| DH | +0.0972 (k=56) | [14, 53] | +0.104 | 1.00 (ghost too high to see by L = 16; R30-b) |

### 2.1 Attempts
For each attempt, the result on the controls is given, then the exact stopping point.

**A. Mertens/de la Vallée Poussin positivity on each circle** (3 + 4cos φ + cos 2φ ≥ 0 applied to every p^{−it}; uses a_n ≥ 0 and the torus).
- PROVED classical: ζ ≠ 0 on Re s = 1. For F_θ the same argument gives F_θ ≠ 0 on Re s = 1+θ, which is true, so it does not separate.
- Floor it gives: D(L) ≥ −C e^{L/2} exp(−c√L) (PNT error), i.e. log-exponent ½.
- **Stop:** this is an L^∞ (absolute-convergence) argument on T and reaches only the L^∞ edge σ = 1. RH is an L² statement (edge ½). F13W needs D ≥ −Ce^{εL} with ε < ½ to say anything new: ε gives zero-free Re s > ½+ε (Lean). Every ε < ½ is a quasi-RH, OPEN.

**B. Positive-kernel / Chebyshev / sieve upper bounds.**
- PROVED here (partial summation; CHECKED to 5·10⁻¹² at L = 4, 8, 12): S(L) = ∫₁^X ψ(x) w_L(x) dx with **w_L(x) = x^{−3/2}(1 + (L − log x)/2) > 0**, and ∫₁^X x w_L = 4√X − L − 4. Hence
  D(L) = ∫₁^X (x − ψ(x)) w_L(x) dx + [4(e^{−L/2} − 1) + (c₁+1)L + C₂ − T(L)].
  **The floor is a one-sided, positively weighted "ψ(x) ≤ x + O(√x)" up to X = k⁴.** Only upper bounds for primes are needed, never lower bounds.
- Inputs tried:
  - Chebyshev/Rosser–Schoenfeld ψ ≤ 1.03883x gives D ≥ −0.155 e^{L/2} − O(L).
  - Brun–Titchmarsh/Selberg upper sieve: constant 2 in short intervals, exponent ½ again.
  - Under F_θ the same bounds hold with its own main term, so nothing false is "proved" for F_θ, and nothing is gained.
- **Stop:** the √x scale of the error. The upper-bound sieve has no a_p = 1 amplifier strong enough to reach O(√x). Landau's theorem (F13W's engine) turns a one-sided bound at that scale into RH, so the one-sided form is RH-strength (by F13S), not easier.

**C. Equidistribution / Haar average along the winding** (Kronecker–Weyl, Bagchi).
- The Fejér density at height t, ρ_X(t), is the prime sum evaluated at the point ι(t) ∈ T (cos(t log n) = Re χ_n(ι(t))). Its average over t is the zero density (1/2π) log(t/2π) + O(1/log X) > 0 (Riemann–von Mangoldt; PROVED for every model with a counting function).
- This is equally true for ζ, F_θ, EPS and DH, so it **does not separate**.
- **Stop:** averages see only Haar-typical points of T. **The floor lives at the base point ι(0) = (1,1,1,…)**, the identity of T, where every circle is at phase 0. That point is maximally resonant and has Haar measure 0. A return of the winding to within ε of the identity on all circles p ≤ P has density ≈ ∏_{p≤P} ε/π (HEURISTIC estimate, independent circles). That is positive for fixed P (Kronecker) but tends to 0 as P grows with X.
  - Illustration: at the real base point s₀ = 0.75 the recurrence density stays 0 up to T = 4000 even with ε = 0.5 (T5).

**D. Positive-density floor** (the weaker variant asked for).
- "D(4 log k) ≥ −B for a positive proportion of k in [K, K^{3/2}]" holds **for F_θ**: about 51% of its samples are ≥ 0 in [21, 97], and 65% are ≥ −0.1. For EPS it is 68%.
- An off-line ghost contributes an oscillation e^{(β−½)L}cos(γL + φ)/|ρ|², which is ≥ 0 about half the time (HEURISTIC for the general claim; NUMERICAL here).
- So the positive-density version is (conjecturally) true for every model and gives nothing through F13S, which needs every k in the window.
- **Stop:** F13S's "all k" quantifier is where the one-sidedness lives. A density version is exactly what F_θ satisfies.

For a dyadic tent sum $P(x)$ over ordinary primes, Lee Crossley proves unconditionally that every sufficiently late block $[K,K^2]$ contains a run $m\le n\le m+\lfloor m^\eta\rfloor$ with $P(n^4)<(n^4/2)\log 2$, for any fixed $0<\eta<1$. Keeping $P(n^4)$ at or below this main term throughout arbitrarily late blocks is equivalent to RH. [Crossley (2026), theorem and proof](https://doi.org/10.5281/zenodo.24411621)

**E. Growing B: the exponent dictionary** (PROVED direction in Lean; converse classical).
- D ≥ −C e^{εL} on windows [L, (3/2)L] ⇒ ζ ≠ 0 on Re s > ½+ε (F13W `F13W_three_halves`).
- Conversely, zero-free on Re s > ½+ε ⇒ D ≥ −C_δ e^{(ε+δ)L} (explicit formula, CLASSICAL).
- So "B growing like k^{4ε}" is equivalent to a zero-free strip of width ½−ε. Best known unconditional: ε = ½ up to sub-exponential factors (Vinogradov–Korobov), i.e. nothing beyond Re s = 1.
- For F_θ the true exponent is θ = ¼: its zeros sit at Re ½+θ. NUMERICAL: the fitted slope of log(−window min) over L ∈ [9, 18] is 0.36, pre-asymptotic and consistent with growth; EPS 0.74, pre-asymptotic and not meaningful at this range.
- **Stop:** any ε < ½ is the quasi-Riemann hypothesis, OPEN.

**F. Multiplicity one, quantitatively, on T.**
- The only places a_p = 1 enters on the torus (§1.3) are:
  - (1) unit-radius circles;
  - (2) L² edge = mirror;
  - (3) bounded Möbius inverse, giving the density theorem.
- Items (2)–(3) give "off-line zeros, if any, have density 0 in the flow" (Bohr–Landau, CLASSICAL): N(σ,T) = o(T), refined by Ingham/Huxley/Guth–Maynard to N(σ,T) ≪ T^{c(1−σ)+ε}.
- Feeding a density bound into the floor gives: the ghost contribution to D(L) is Σ_{off} e^{(β−½)L}/|ρ|², and density bounds control the count but not the existence of a single ghost with β−½ = δ. That ghost alone breaks the floor at L ≈ (2/δ) log|ρ|.
- **Stop:** density theorems are "few ghosts". The floor needs "no ghost", one point at a time, at the base point. That is RH, and nothing on T between "o(T)" and "0" is known.

**Summary of the stopping point.** Every torus-side ingredient is one of two kinds:
- (i) it is shared with F_θ (positivity, Euler product, FE, unitary intertwiner, Haar averages, positive density); or
- (ii) it is a "typical point / density zero" statement (Bagchi, Bohr–Landau, Voronin).

The integer floor is a **pointwise statement at the identity of T, at the L² edge σ = ½**. No known tool reaches it. **The floor is not proved; GAP-30b stands.**

## 3. The lift: does D(z) = 2z carry arithmetic? (`r31_lift.py` → `r31_lift.json`)

- **L1, the doubling is a change of variables (PROVED; CHECKED).** f ↦ f(r²) intertwines J₁ and J₂ (J₂[f(r²)] = (J₁f)(r²)), and M[f(r²)](Δ) = ½·M[f](Δ/2), which is Δ = 2s.
  - For the rung-0 Gaussian lattice sum g(r) = Θ(r²) − 1 (with M[g](s) = π^{−s/2}Γ(s/2)ζ(s)), the lift's Mellin transform equals ½·Λ(Δ/2). Agreement is to 10⁻¹⁷, 1.4·10⁻¹⁰ and 0 at three test points with Re Δ > 2, and the closed form matches to ≤ 10⁻²⁰.
  - On the Dirichlet side, n^{−s} = (√n)^{−Δ}: the same torus, circle lengths (log p)/2, flow rescaled by 2.
  - Every statement of §§1–2 (Bagchi, the floor, F13S) transports verbatim, and F_θ lifts to zero lines Re Δ = 1 ± 2θ, still off the lifted line.
  - **The doubling sees the stage (mirror, unitary line) and carries the arithmetic only as a relabelling. It adds no constraint.**
- **L4, the natural arithmetic up the ladder (PROVED, classical; CHECKED n ≤ 8 by brute force).** The integer lattices of the division algebras have theta series of weight d/2 and Dirichlet series:
  - ℤ (d=1): 2ζ(2w), i.e. ζ.
  - ℤ[i] (d=2): r₂(n) = 4Σχ₋₄(d), giving 4ζ(w)L(w,χ₋₄). Two species, both unit-radius, the **same** critical line. Counted-once survives in the form "Ramanujan for each factor".
  - Hurwitz quaternions (d=4): 24Σ_{d|n, d odd} d, giving 24(1−2^{1−w})ζ(w)ζ(w−1).
  - E8 / Coxeter–Dickson octaves (d=8): θ = E₄, 240σ₃(n), giving 240ζ(w)ζ(w−3).
  - Since ζ(w)ζ(w−k) = ζ(u+k/2)ζ(u−k/2) with u = w−k/2, **the ℍ and 𝕆 rungs are exactly F_θ with θ = ½ and θ = 3/2.** All their zeros are off their FE line, unconditionally, because the shifted strips (0,1) ± k/2 do not contain the centre line.
  - Brute-force counts match the formulas for every n ≤ 8, e.g. E8: 240, 2160, 6720, 17520, 30240, 60480, 82560, 140400.
- **L2, Daniel's Gr(2,4) rung (PROVED; CHECKED).**
  - #Gr(2,4)(F_q) = 1+q+2q²+q³+q⁴ by RREF enumeration, matching for q = 2, 3, 4, 5, 7. The Euler factor at p = 2 agrees with the Hasse–Weil series to 2·10⁻³².
  - Normalised by the centre p², the per-prime radii are p^{−2}, p^{−1}, 1, 1, p, p²: a two-level stretch (θ = 1, 2) plus multiplicity 2 in the middle.
  - Zero strips (k, k+1), k = 0..4. Only the k = 2 family (ζ(s−2)², multiplicity 2) can lie on the FE line 5/2, and that is RH again.
  - So over ℤ, "turns not stretches" holds exactly at weight 0. The higher rungs are pure Tate twists: stretched by their weight, as purity predicts over F_q (Deligne; trivial for cellular varieties).
- **Verdict (honest).** The ladder's arithmetic above ℂ is stretched and cannot be the place where ζ's unitarity is decided. The celestial rung Δ = 2s (Pasterski–Shao) is a reparametrisation of rung 0.
  - A genuine extra constraint would need extra structure, e.g. a 2D crossing equation (Benjamin–Chang 2022). There the zeros appear as a restatement, not a proof (HEURISTIC/OPEN).
  - ℂ is the one rung with new arithmetic and no stretch (ζ·L(χ₋₄)), but its RH is GRH for χ₋₄ on top of RH, a parallel problem, not a constraint on ζ.

## 4. Figures
- `fig1_prime_torus.png`, three panels:
  - the winding on T² (2 and 3 circles), never closing;
  - the winding on T³ (2, 3, 5);
  - integers as lattice points on the dual face, projected to log n, each landing once.
- `fig2_faces_turn_vs_stretch.png`, three panels:
  - the two faces in space (x^{1/4}Θ(x) symmetric about x = 1);
  - the two faces in frequency: mirror line, L² and L^∞ edges for ζ vs F_θ;
  - the drift test (turn vs stretch) along the flow.
- `fig3_integer_floor.png`: D(L) and the F13S samples D(4 log k) for the four models; window minima and the positive-density fractions.

## 5. Labels and the next GAP
- **PROVED (classical, attributed):**
  - T = dual of ℚ_{>0}^×; the winding is dense, injective and equidistributed;
  - Bohr lift and Carlson mean square;
  - the Bagchi limit law and criterion (RH ⇔ strong recurrence);
  - Voronin universality;
  - Bohr–Landau density;
  - the three-way location of a_p = 1 (§1.3; ingredients classical, organisation ours);
  - F_θ's L² edge at ½+θ and its FE drift;
  - the doubling is a change of variables;
  - the division-lattice Dirichlet series;
  - the ℍ/𝕆 rungs are F_θ-type with all zeros off-line;
  - the Gr(2,4) Hasse–Weil stretch;
  - the positive-kernel form of D (§2.1 B).
- **PROVED (Lean, GPPVerify `c81c031`):** F13S (floor ⇒ RH); F13W exponent dictionary (floor e^{εL} ⇒ zero-free Re s > ½+ε); the integer/real chord penalty ≤ 8.
- **CHECKED:** Lean closed form = explicit-formula D (≤ 6·10⁻¹⁴); kernel identity (≤ 5·10⁻¹²); lift Mellin identity; lattice counts n ≤ 8; Gr(2,4) counts and Euler factor; continued fractions.
- **NUMERICAL:**
  - Bagchi law (KS p = 0.23, 0.50); drift (T3); variance (T4); recurrence (T5);
  - floor tables: ζ ≥ +0.029 at every k ≤ 100, F_θ min −1.18, EPS min −4.60, DH ≥ +0.097;
  - F_θ positive fraction ≈ ½.
- **HEURISTIC:** "golden vs generic" winding; base-point resonance as the reason torus methods stop; Lee–Yang reading (every finite Euler product is zero-free, so zeros are collective, thermodynamic-limit objects; R30-b).
- **OPEN:** the integer floor; any ε < ½ in F13W (quasi-RH); bounded-window recurrence; **RH.**
- **GAP:** see below.

**Next GAP (GAP-31, sharpest form found).** Prove, for some B and arbitrarily large K, the **one-sided weighted prime bound at the identity of T**:

  ∫₁^{k⁴} ψ(x) w_L(x) dx ≤ 4k² − 4 log k − 4 + B + O(1)  for all k ∈ [K, ⌈K^{3/2}⌉],  with w_L(x) = x^{−3/2}(1 + (L − log x)/2), L = 4 log k.

This is GAP-30b restated with no new strength; it is RH-strength by F13S. The new content is the diagnosis:
- the a_p = 1 input must enter through **(3)** of §1.3, the bounded Möbius inverse/mollifier, and **pointwise**;
- (1)–(2) are already used by the Bagchi/density machinery and only give "density-zero ghosts";
- F_θ fails (2), DH fails (3), and the floor must use both.

A concrete testable sub-step, OPEN and proposed not run: a mollified floor D_M(L) built from ζ·M_y, M_y = Σ_{n≤y} μ(n)n^{−s}. Check numerically whether its one-sided floor is more stable (smaller window variance) than D's for ζ, while F_θ and DH stay broken.

**Lean.** No module this round. Proposals, to be claimed in the conversation file first:
- (L3) the positive-kernel identity S(L) = ∫₁^X ψ w_L (elementary partial summation; it would let the floor be stated as a one-sided ψ bound in Lean);
- (L4) ζ(w)ζ(w−k) = F_θ-form with θ = k/2 and "no zero on the centre line" for k ≥ 1 (elementary from ζ ≠ 0 on Re ≥ 1 and the FE).

Neither is on the F13 critical path.

## 6. Literature, novelty, attribution
- **Prime torus and Bohr lift:** H. Bohr (1913–1922; almost periodicity; Bohr's RH criterion for L(s,χ), 1922); Bohr–Jessen (Acta Math. 54 (1930), 58 (1932)); H. Helson (1969); Hedenmalm–Lindqvist–Seip (Duke 1997), H²(T^∞).
- **Universality and recurrence:** S. M. Voronin (Izv. 1975); B. Bagchi (PhD thesis, ISI Calcutta 1981; Math. Z. 181 (1982) 319–334; Acta Math. Hungar. 50 (1987) 227–240); J. Steuding, *Value-Distribution of L-Functions* (LNM 1877, 2007), Thm 8.3. Later work: Nakamura (Acta Arith. 2009) and Pańkowski (2010), generalized strong recurrence; Garunkštis (arXiv 1006.1507), self-approximation for rational d; Laurinčikas (Computation 12 (2024) 164) and Symmetry 17 (2025) 2075, density and short-interval versions; Gauthier (Fields Inst. Commun. 65, 2013).
- **Density:** Bohr–Landau (1914); Ingham (1940); Huxley (1972); Guth–Maynard (2024).
- **Controls:** Davenport–Heilbronn (1936); Spira (1994) zeros; Potter–Titchmarsh (1935).
- **Lattices:** Jacobi (four squares); Hurwitz (1919); Coxeter (1946), E8 theta = E₄.
- **Gr(2,4)** from Daniel's paper "The Division Algebra Tower as Holographic Chain" (attached, Thm "Hasse–Weil factorization").
- **Primon gas and Lee–Yang:** Julia, Spector, Bost–Connes; Lee–Yang (1952). **Connes** (1999): adele class space; our T is the Bohr side, the dual of ℚ_{>0}^×, not Connes' Ẑ^×.
- **Project credits:**
  - Alpöge's window architecture (`zeroFree_of_rpow_window`) underlies F13W/F13S;
  - Codex did the fourth-power reduction (F13S);
  - Turbo did the audits (incl. the R30-b corrections carried here);
  - Claude Code gave the 1D-CFT reading (quartet → pair; pointwise form).
- **Novelty check.** Bagchi's criterion and the Bohr–Jessen law are classical. We have not seen the following stated this way: the three-way location of a_p = 1 on T against the F_θ/DH controls; the positive-kernel form of the F13S floor; and the observation that the ℍ and 𝕆 integer lattices of the division-algebra ladder give exactly F_θ-type products. The ingredients are classical.

## Files
- `r31_torus.py` → `r31_torus.json`, `torus.log`: §1 (T1–T5).
- `r31_floor.py` → `r31_floor.json`, `floor.log`: §2 numerics. Uses `r22_fejer.py` and `r22_models.py` (copied from R30-b/R23).
- `r31_attempts.py` → `r31_attempts.json`: kernel identity, Chebyshev floor, envelopes.
- `r31_lift.py` → `r31_lift.json`: §3.
- `r31_figs.py` → `fig1_prime_torus.png`, `fig2_faces_turn_vs_stretch.png`, `fig3_integer_floor.png`.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).
