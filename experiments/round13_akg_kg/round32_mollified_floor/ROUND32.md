# Round 32: the mollified integer floor (the Möbius side of GAP-31)

Daniel Toupin approved this in his own words: "Do it up" (2026-11-21, Toronto), replying to the proposal "Say go and I'll start the floor test".
Seat: Heart of the Cosmos / Grok (D5). Follows Round 31 (`../round31_prime_torus/ROUND31.md`, `7f89178`).

Rules for this round:
- No zeta zeros are inputs anywhere. The first 300 zeros appear once, as test data in `r32_validate_zeros.py`, to validate code. No claim uses them.
- Weil positivity is never the target. **RH is not closed.**

## Verdict (plain)
1. **A fixed-length mollifier separates nothing (PROVED).** Selberg's mollified floor at a fixed level y is an explicit, zero-free function of μ on [1, y]. All zero information sits in the tail d > y, and the mollified part drifts like √(X/y). The only length that keeps O(1) accuracy is y = X. At y = X the mollified floor *is* the floor, rewritten exactly in Möbius coordinates (identity I3).
2. **The full-length Möbius floor does separate (NUMERICAL, and PROVED for each control).** Define m̃(L) = Σ_{d≤e^L} μ(d) d^{−1/2}(L − log d) − L/ζ(½) + ζ′(½)/ζ(½)². For ζ it stays flat in [−0.022, +0.022] for every L ≤ 20.7 (X ≤ 10⁹). With each control's own Dirichlet inverse:
   - F_θ (θ = 0.1, 0.25), ζ² and Epstein grow visibly, to ±0.85, ±8.8, ±0.43 and ±5.6.
   - Every control is *proved* to have an unbounded Möbius floor (Thm 32.1 + their known zeros).
   - DH and Beurling are not yet separated at the computed range (explained below).
3. **The new statement is GAP-32: m̃(L) = O(1).** It implies RH *and simple zeros* (Thm 32.1, PROVED), hence GAP-31 with B = 0 via the classical explicit formula. It is a smoothed Mertens hypothesis at the L² edge, sitting between the weak Mertens hypothesis and "RH + simple zeros". So it is *stronger* than GAP-31. ζ² shows this: it passes the ψ-floor (D_{ζ²} = 2D_ζ ≥ 0) but fails the Möbius floor unconditionally, because its zeros are double. The Möbius floor counts each zero once.
4. **Why counted-once does not close it (PROVED, Prop 32.3).** The counted-once facts available on the Möbius side are:
   - |μ| ≤ 1 (μ² = squarefree indicator);
   - Σ_{d|n} μ(d) = [n = 1], which gives |Σ_{d≤x} μ(d)/d| ≤ 1.

   Both live at the L^∞ edge (exponent ½ in L), not the L² edge. The only exact bridge from the Möbius floor to the ψ-floor is convolution with the positive integer measure ν = Σ_m m^{−1/2} δ_{log m}, of mass 2√X. That is multiplication by ζ: each integer counted once. So an elementary transfer loses exactly √X. The implication GAP-32 ⇒ GAP-31 is analytic (Landau), never termwise.

   The ψ side has positivity (residues = multiplicities ≥ 0) but no Möbius size bound. The Möbius side has the size bounds but no positivity (its residues 1/ζ′(ρ) have no sign). The floor needs both at once.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).

---

## 0. Setting
- **Model M.** A Dirichlet series b_M(s) = Σ b_M(n) n^{−s} with b_M(1) = 1. Its own Möbius function is μ_M = the Dirichlet inverse of b_M, and Λ_M = (b_M·log) ∗ μ_M.
- **The models:**
  - ζ;
  - ζ² (= F₀; new control: on-line, double zeros);
  - F_θ = ζ(s+θ)ζ(s−θ), θ = 0.1, 0.25;
  - EPS, x² + 5y², b = r_Q/2;
  - DH (Davenport–Heilbronn, period 5);
  - BEUR, the Round-13e DMV-type discrete Beurling system with quartet ½ ± 0.1 ± 28i and primes ≥ 20.
- **R31 floor.** D_M(L) = P_M(L) − S_M(L), where S_M(L) = Σ_{n≤X} Λ_M(n) n^{−1/2}(L − log n), X = e^L, and P_M is the pole/archimedean part. For ζ, P is the Lean closed form (GPPVerify `c81c031`).
- **Möbius floor.**
  - m_M(L) := Σ_{d≤X} μ_M(d) d^{−1/2}(L − log d). Its Laplace transform is ∫₀^∞ m_M e^{−sL} dL = 1/(b_M(s+½) s²).
  - m̃_M(L) := m_M(L) − L/b_M(½) + b_M′(½)/b_M(½)², which subtracts the s = 0 residue.
  - Constants are Arb-certified (python-flint, 256 bits; `r32_constants.json`), e.g. 1/ζ(½) = −0.684765236089936…, −ζ′(½)/ζ(½)² = 1.839342223692252….
  - BEUR has no FE. Its constants come from the continuation in `r32_beur_zeta.py`: ζ_B(½) = −0.03854, ζ_B′(½) = −0.3401 (truncation ~10⁻⁴; not certified).
- **Second moment.** m_{M,2}(L) := ½ Σ μ_M(d) d^{−1/2}(L − log d)², with Laplace transform 1/(b_M s³).

## 1. Exact identities (PROVED; CHECKED in `r32_identities.py` → `r32_identities.json`)
- **(I1) Leibniz.** Λ_M = −b_M ∗ (μ_M·log). log is a derivation of ∗; apply it to b ∗ μ = δ. CHECKED n ≤ 3000, error ≤ 6·10⁻¹³ for ζ, F_0.25, DH, EPS.
- **(I2) Selberg truncation.** Λ_{M,y}(n) := Σ_{d|n, d≤y} μ_M(d) log(y/d) b_M(n/d).
  - It equals log y at n = 1 and Λ_M(n) for 1 < n ≤ y.
  - For n > 1: Λ_M − Λ_{M,y} = −b_M ∗ (μ_M log(·/y) 1_{>y}).
  - CHECKED y = 3000, error ≤ 5·10⁻¹³.
- **(I3) The floor in Möbius coordinates** (I2 with y = X):

  S_M(L) = Σ_{d≤X} μ_M(d)(L − log d) d^{−1/2} A_M(X/d) − L²,  where A_M(Y) = Σ_{m≤Y} b_M(m) m^{−1/2} log(Y/m).

  Equivalently **S_M = ν_M ⋆ g_M**, with ν_M = Σ_m b_M(m) m^{−1/2} δ_{log m} and g_M = 2m_{M,2} − L·m_M. In Laplace terms: −b′/b = b·(1/b)′, and (1/b)′/s² = d/ds[1/(b s²)] + 2/(b s³).
  - CHECKED at L = 8, 11, 13 for all four models.
  - Absolute errors are 5·10⁻¹² to 8·10⁻⁸, on sums of size 10²–10³.
- **(I4) ζ's integer kernel.** A(Y) = 4√Y + ζ(½) log Y + ζ′(½) + E(Y), with E(Y) = O(Y^{−3/2}).
  - Proof: Mellin, with residues of ζ(s+½)Y^s/s² at s = ½ and s = 0, plus Euler–Maclaurin for the sawtooth.
  - CHECKED: |E(Y)| ≤ 0.234 Y^{−3/2} on [1, 6·10⁵].
- **(I5) ζ's floor as a sum of Möbius floors** (I3 + I4):

  D(L) = 4√X(1 − m⁽¹⁾(L)) − 2ζ(½) m̃₂(L) − ζ′(½) m̃(L) − E_μ(L) + R₀(L)

  - m⁽¹⁾(L) = Σ_{d≤X} μ(d)(L − log d)/d, the Möbius floor at σ = 1, whose main term is 1.
  - m̃₂ = m₂ − ½[L²/ζ(½) − 2Lζ′(½)/ζ(½)² + (1/ζ)″(½)].
  - E_μ = Σ μ(d)(L − log d) d^{−1/2} E(X/d).
  - **The L² terms cancel, and so do the L terms**, because c₁ = −ζ′/ζ(½) (the FE: ξ′(½) = 0; Arb residual 2·10⁻⁷⁶).
  - **The constant is R₀(∞) = −8 + C₂ − ζ(½)(1/ζ)″(½) + (ζ′/ζ)²(½) = 0.0462099862… = (log ξ)″(½) = −Σ_ρ (ρ−½)^{−2}.**
    - This is the Hadamard constant, i.e. Σ_ρ 1/γ² under RH. It was computed here from integers and Γ-data only, and agrees with mpmath's (log ξ)″(½) to 10⁻¹⁶.
    - The mean of D over L ∈ [4, 20.7] is 0.04623.
  - So in Möbius coordinates, **D − Σ1/γ² is a pure Möbius oscillation.** Each piece is two-sided, of size ±0.09, and only the sum is one-sided (fig. 2).

## 2. GAP-32 and what it implies

**Theorem 32.1 (PROVED here; Landau's theorem plus standard Laplace arguments).** Let m be ζ's Möbius floor.
- **(a) One-sided.** Suppose m(L) − a − cL ≤ B for all L ≥ L₀, or ≥ −B for all L ≥ L₀. Then ζ(s) ≠ 0 for Re s > ½. If moreover c = 1/ζ(½), every zero on the line is simple.
- **(b) Two-sided.** If |m(L) − a − cL| ≤ B for L ≥ L₀, then c = 1/ζ(½), ζ ≠ 0 on Re s > ½, and all zeros are simple.
- **(c) Exponent dictionary.** If |m̃(L)| ≤ B e^{εL} for L ≥ L₀, then ζ ≠ 0 on Re s > ½ + ε.
- **(d)** The same holds for every model M with b_M meromorphic on Re s > ½ − δ. Case (a) also needs b_M(σ) ≠ 0 on the real segment.

*Proof.*
- For Re s > ½, Fubini gives ∫₀^∞ m e^{−sL} dL = 1/(ζ(s+½)s²).
- (a) f = B − (m − a − cL) ≥ 0 on [L₀, ∞). Its Laplace transform F equals G(s) = B e^{−sL₀}/s + a/s + c/s² − 1/(ζ(s+½)s²) + (entire). G is analytic on the real axis s > 0: ζ(σ) < 0 on (0, 1), since (1 − 2^{1−σ})ζ(σ) = η(σ) > 0.
  - By Landau, the abscissa of a nonnegative Laplace integral is a singularity. So F converges on Re s > 0 and 1/ζ(s+½) continues analytically there.
  - Simplicity: |F(σ+iγ)| ≤ F(σ) = (c − 1/ζ(½))/σ² + O(1/σ). With c = 1/ζ(½), this is O(1/σ), so poles at iγ have order ≤ 1.
- (b) |F| ≤ B/σ + O(1) directly. Near s = 0 this forces c, and near s = iγ it forces simple poles.
- (c) Absolute convergence on Re s > ε. ∎

**GAP-32 (new).** ∃ B, L₀ such that |m̃(L)| ≤ B for all L ≥ L₀. One-sided variant GAP-32⁻: m̃(L) ≥ −B.
- **GAP-32 ⇒ GAP-31.** By Thm 32.1(b), GAP-32 gives RH and simple zeros. Then D(L) = Σ_ρ m_ρ(1 − cos γL)/γ² ≥ 0 by the classical explicit formula (R31 §2.0; CLASSICAL, not in Lean). So GAP-31 holds with B = 0, for every k. Rigorous, but through analytic continuation.
- **GAP-32 ⇒ zero-free Re s > ½ + ε** already holds for the e^{εL} version (Thm 32.1(c)). This is the Möbius twin of F13W.
- **Windows.** "F13S-style windows [K, K^{3/2}] suffice" is OPEN for the Möbius floor. It needs a Möbius analogue of Alpöge's window lemma.
- **Converse (CLASSICAL, conditional; not re-proved).** RH, simple zeros and J := Σ_ρ 1/(|ζ′(ρ)||ρ−½|²) < ∞ give m̃(L) = Σ_ρ e^{iγL}/(ζ′(ρ)(iγ)²) + O(e^{−5L/2}), hence |m̃| ≤ J + o(1).
  - The weak Mertens hypothesis ∫₁^X (M(x)/x)² dx ≪ log X gives RH, simple zeros and Σ 1/|ρζ′(ρ)|² < ∞ (Titchmarsh §14.29).
  - By Cauchy–Schwarz with Σ1/|ρ|² < ∞, that gives J < ∞.
  - So: **weak Mertens ⇒ GAP-32 ⇒ RH ∧ simple zeros ⇒ GAP-31**, where the first arrow is the conditional sketch.
- **Validation (zeros as test data only).** The first 300 zeros reproduce the integer-side m̃ at L = 6, 9, 12, 15, 18, 20.5 to ≤ 6·10⁻⁵ (`r32_validate_zeros.json`).

**Controls are unbounded (PROVED, unconditional):**
- **ζ².** Every zero of ζ (Hardy: infinitely many) is a double zero of ζ², so (b) fails and m̃_{ζ²} is unbounded on both sides. Yet D_{ζ²} = 2D_ζ, so under RH ζ² *passes* GAP-31.
- **F_θ.** Hardy's on-line zeros give zeros at ½ + θ + iγ. By (c), limsup |m̃_θ| e^{−εL} = ∞ for every ε < θ.
- **EPS.** It has infinitely many zeros in σ > 1 (Davenport–Heilbronn 1936).
- **DH.** It has infinitely many zeros in σ > 1 (Titchmarsh §10.25) and in ½ < σ < 1 (Voronin; Spira's 0.8085 + 85.70i).
- **BEUR.** ζ_B = ζ_C·exp(∫x^{−w} d(Π_B − Π_C)) with |Π_B − Π_C| ≤ 0.9994 (CHECKED). So its zeros in Re w > 0 are exactly the quartet, and two of them have Re w > ½.

## 3. Numbers (`r32_pass.py` X ≤ 10⁹; `r32_controls.py` 10⁷; `r32_beurling.py` 10⁸ → `r32_analysis.py` → `r32_results.json`)
NUMERICAL. Sums are finite float64 sums in segments, binned per threshold. The code CHECKs against brute force (sympy μ, direct Dirichlet inversion) to ≤ 10⁻¹².

| model | X | R31 floor min D(4 log k) (k) | Möbius floor m̃ range, L ≥ 3 | window envelope max\|m̃\| on [L, 1.5L], L = 6 → last | verdict |
|---|---|---|---|---|---|
| **ζ** | 10⁹ (k ≤ 177) | **+0.0248** (150); real-L min +0.0216 | **[−0.0223, +0.0218]** | **0.0197 → 0.0223 (flat)** | bounded, margin stable |
| ζ² (F₀) | 10⁹ | +0.0497 = 2×ζ (passes) | [−0.417, +0.432] | 0.195 → 0.417 (∝ L: double zeros) | **fails** |
| F_θ, θ = 0.1 | 10⁹ | −0.068 (150) | [−0.755, +0.855] | 0.221 → 0.755 | **fails** |
| F_θ, θ = 0.25 | 10⁹ | −3.08 (150) | [−6.90, +8.80] | 0.40 → 6.90 | **fails** |
| EPS | 10⁷ (k ≤ 56) | −4.60 (56) | [−3.75, +5.62] | 1.39 → 5.28 | **fails** |
| DH | 10⁷ | +0.097 (56) | [−0.104, +0.090] | 0.084 → 0.090 | not separated in range |
| BEUR | 10⁸ | (no FE; detrended ±0.016) | −122 (L = 4) → −0.34 (L = 18.4), monotone | transient | not separated in range |

- **Margin.** The one-sided Möbius floor for ζ is m̃ ≥ −0.0223. The window minima are −0.0186 on [31, 173], with no drift.
  - F_θ (both θ), ζ² and EPS reach −0.40 to −6.9, with growing envelopes. Fitted log-slopes over L ∈ [6, 20.7] are 0.17, 0.38 and 0.30; ζ² grows linearly in L.
  - These slopes are pre-asymptotic. They are consistent with the proved lower rates Ω(e^{(θ−ε)L}).
- **Sharper than the ψ-floor.** The Möbius floor separates F_0.1 at 34× ζ's envelope, while D_{F_0.1} only dips to −0.07 at k = 150. It separates ζ², which the ψ-floor cannot.
- **Why DH is invisible.** Its off-line zeros are weak or high. Spira's zero contributes ≈ e^{0.31L}/(|DH′|·85.7²) ≈ 0.02/|DH′(ρ)| at L = 16. The σ > 1 zeros sit at heights where the 1/γ² weight kills them.
  - Yet |μ_DH| grows like n^{≈0.6} on dyadic blocks, and |μ_EPS| like n^{≈0.8} (fig. 4). These are the Dirichlet-series shadows of zeros in σ > 1.
- **Why BEUR is invisible.** With no primes below 20, ζ_B is tiny left of the line (|ζ_B(½)| = 0.0385; ≈ e^{−Ei(…)}). The left-of-line terms of 1/(ζ_B(s+½)s²) start huge and decay like e^{−(½−σ)L}. The quartet's 28-frequency amplitude is ≲ 0.006 at L ≤ 18.4, and its growth e^{0.1L} needs L ≫ 30.
- **ζ decomposition (fig. 2), L ≥ 4.**
  - 4√X(1 − m⁽¹⁾) ∈ [−0.087, +0.089].
  - −ζ′(½)m̃ ∈ [−0.088, +0.085].
  - −2ζ(½)m̃₂ ∈ ±0.0035.
  - E_μ ∈ [−0.024, +0.026].
  - Their sum is D ∈ [+0.0216, +0.0689].

## 4. The fixed-level Selberg mollifier (`r32_selberg.py` → `r32_selberg.json`, fig. 3)
- D_y := P_M − S_y with Λ_y as in (I2) (P(x) = x, Selberg). Its generating function is b_M(s)·N_y(s), with N_y = Σ_{d≤y} μ_M(d) log(y/d) d^{−s} a finite Dirichlet polynomial.
- **There are no poles from zeros (PROVED).** So D_y is an explicit function of μ_M on [1, y].
- D_y = D for X ≤ y (CHECKED ≤ 5·10⁻⁸, all models, y = 10²–10⁵).
- For X > y and ζ: D_y(L) e^{−L/2} → 4(1 − N_y(1)).

| y | N_y(1) − 1 | observed D_y e^{−L/2} at L = 16.1 | predicted 4(1 − N_y(1)) |
|---|---|---|---|
| 10² | +1.03·10⁻³ | −4.05·10⁻³ | −4.14·10⁻³ |
| 10³ | +4.07·10⁻⁴ | −1.54·10⁻³ | −1.63·10⁻³ |
| 10⁴ | +2.77·10⁻⁵ | −8.7·10⁻⁵ | −1.1·10⁻⁴ |
| 10⁵ | −2.47·10⁻⁵ | +8.5·10⁻⁵ | +9.9·10⁻⁵ |

- So **ζ itself fails a fixed-level mollified floor**, at rate √(X/y)·(Möbius floor at σ = 1, level y). The controls drift the same way (F_θ like X^{½+θ}).
- **Stop (exact).** A mollifier of length y sees zeros only through its tail d > y. O(1) accuracy at X needs y ≳ X, and at y = X the "mollified floor" is identity I3. This is why Selberg/Levinson mollifiers win only for *averaged* statements (density theorems, proportions of zeros), never at the single point t = 0.

## 5. The counted-once property doing the work, and where it stops
- **(P1) Size.** |μ(n)| ≤ 1, equivalently μ² = 1_squarefree (a_p = 1, degree one, one species). On dyadic blocks (fig. 4, CHECKED), max|μ_M| is 1 for ζ, ≈ 2^{ω(n)}·n^{θ} for F_θ (2^{ω} for ζ²), and grows polynomially for EPS and DH.
  - This kills F_θ, ζ², EPS and DH in the sense of item (d). It does not kill BEUR (|μ_B| ≤ 1).
  - Quantitatively it gives only the trivial |m(L)| ≤ Σ_{d≤X} μ²(d) d^{−1/2}(L − log d) ≤ 4e^{L/2}, i.e. exponent ½.
- **(P2) Additive inversion on ℕ.** Σ_{d|n} μ(d) = [n = 1] gives Σ_{d≤x} μ(d)⌊x/d⌋ = 1, hence **|Σ_{d≤x} μ(d)/d| ≤ 1** (classical; CHECKED x ≤ 10⁵).
  - This is the σ = 1 floor: |m⁽¹⁾(L)| ≤ L.
  - It is the place where ℕ's lattice structure (⌊·⌋, sawtooth < 1) enters. Beurling has N_B(x) = x + O(x^{0.6}) instead, so it fails here.
- **Prop 32.3 (PROVED, elementary): the transfer cost.**
  - By (I3), S = ν ⋆ g with ν ≥ 0 and ν([0, L]) = Σ_{m≤X} m^{−1/2} = 2√X + ζ(½) + O(X^{−1/2}).
  - Hence a Möbius-side bound |g − main| ≤ G gives only |S − main| ≤ G(2√X + O(1)).
  - From GAP-32 (|m̃| ≤ B, so |m̃₂| ≤ BL + O(1) and |g − main| ≤ 3BL + O(1)), the termwise route gives only D ≥ −(6BL + O(1))√X, i.e. exponent ½ again.
  - **The integers, counted once, are exactly the positive measure that costs the √X.**
- **Reduction (named).**
  - GAP-32 is a Cesàro-smoothed (order 1) Mertens hypothesis at σ = ½, between weak Mertens and RH ∧ simplicity (§2).
  - Through (P2), its σ = 1 companion becomes a statement about μ against the sawtooth functions {x/d}. These are the Nyman–Beurling generators, and Báez-Duarte's natural approximants use exactly μ(k)ρ(1/(kx)).
  - Burnol: liminf d_N² log N ≥ Σ_ρ m_ρ²/|ρ|². Multiplicity *squared*, so it carries the same "counts zeros once" signature as GAP-32 (ζ² fails both). The constants line up: (log ξ)″(½) = 0.04621 = R₀(∞) here, and d_N² log N ≈ 0.046 (Round 13g). HEURISTIC alignment.
- **What attacking the named form directly would need.** None of these is available:
  1. **Möbius form:** square-root cancellation of μ against the Fejér kernel, pointwise at t = 0, uniformly in L.
     - Best known is Vinogradov–Korobov: m̃ ≪ e^{L/2} exp(−cL^{3/5}(log L)^{−1/5}), exponent ½.
     - Any exponent < ½ is a quasi-RH (Thm 32.1(c)).
  2. **Bilinear form:** cancellation in Σ_{d,m} μ(d)(dm)^{−1/2} K(L − log dm), with dm ≤ X, at full length.
     - Vaughan / Heath-Brown identities save only in averaged settings.
     - Positivity-preserving versions (Selberg Λ²) hit the parity factor 2 (R31 attempt B).
  3. **Báez-Duarte form:** d_N → 0. This needs the μ(h)μ(k)-weighted double sum of Vasyunin cotangent-sum Gram entries.
     - Bettin–Conrey reciprocity (2013) controls individual entries.
     - The μ-weighted sum is the same Möbius-cancellation problem.
- **The missing principle, stated exactly.** A positivity on the Möbius side, or a √X-free bridge between the positive ψ side and the bounded μ side.
  - The ψ-floor's one-sidedness comes from residues +m_ρ (R31).
  - The Möbius floor's boundedness comes from |μ| ≤ 1 and ⌊x/d⌋. Its residues are 1/ζ′(ρ), which have no sign.
  - The only exact bridge, ν (multiplication by ζ), is positive but has mass 2√X.

## 6. Lean (proposals only; claim in the conversation file first; NOT built)
- **(L5) Möbius floor ⇒ zero-free.** If |Σ_{d≤e^L} μ(d) d^{−1/2}(L − log d) − a − cL| ≤ B e^{εL} for L ≥ L₀, then riemannZeta s ≠ 0 on ε + ½ < re s.
  - Ingredients: L-series of μ = ζ⁻¹ on re s > 1 (Mathlib); Laplace transform of a bounded or exponentially bounded function is analytic; identity theorem.
  - A two-sided bound needs no Landau step, so this should be easier than F13W.
  - The simple-zeros clause is optional and a second step.
- **(L6) The floor in Möbius coordinates (finite algebra).**
  - Σ_{d|n} μ(d) log(X/d) = (if n = 1 then log X else Λ n).
  - S(L) = Σ_{d≤X} μ(d)(L − log d) d^{−1/2} A(X/d) − L².
  - Uses Mathlib's `ArithmeticFunction` (vonMangoldt, moebius, Λ = μ ∗ log, μ ∗ ζ = 1). No analysis.
- Spare: (L7) |Σ_{d≤x} μ(d)/d| ≤ 1 from Σ μ(d)⌊x/d⌋ = 1 (clean, classical).

None of these is on the F13 critical path. L5 is the Möbius twin of F13W.

## 7. Labels
- **PROVED:**
  - I1, I2, I3 (finite identities);
  - I4 asymptotics;
  - the I5 cancellations and R₀(∞) = (log ξ)″(½);
  - Thm 32.1 (a)–(d);
  - controls unbounded (ζ², F_θ, EPS, DH, BEUR; cited zero results);
  - fixed-level Selberg floor is zero-free/explicit;
  - Prop 32.3 (transfer cost);
  - P1/P2 classical.
- **CHECKED:** identities to ≤ 10⁻⁷ (relative ~10⁻¹⁰); Arb constants; R₀ vs mpmath (log ξ)″(½); D_y = D below y; |Π_B − Π_C| ≤ 1; brute-force μ.
- **NUMERICAL:** all tables; growth slopes; the BEUR transient; dyadic |μ_M|; zero-side validation (test data).
- **CLASSICAL (conditional, not re-proved):** weak Mertens ⇒ J < ∞ ⇒ GAP-32 (converse sketch); RH ⇒ D ≥ 0.
- **CONJECTURE:** GAP-32 holds for ζ (would follow from RH + simple zeros + J < ∞).
- **HEURISTIC:** Burnol/Báez-Duarte constant alignment; "the integers are the √X measure" as the reason termwise transfer fails.
- **GAP:** GAP-32 (stronger than GAP-31); the windowed Möbius version (OPEN).
- **OPEN:** RH; simplicity; any exponent < ½ for the Möbius floor.

## 7a. Errata carried from the Round 31 audits (Turbo, bridge, 2026-11-21)
- **(L4) as stated in R31 is false in the universal form "k ≥ 1"** (Turbo's counterexample). At k = 5 the centre w = 3 meets the trivial zero ζ(−2) = 0, and likewise k = 9, … whenever (1 − k)/2 is a negative even integer. At k = 1 the centre factor is the pole ζ(1).
  - Accepted. (L4) is withdrawn in that form.
  - What stands is the ladder case only: k = 1, 3, with the pole wording at k = 1 handled. The R31 ℍ/𝕆 conclusion stands at those two values.
  - Kept on record, not erased.
- **R31's Chebyshev coefficient "D ≥ −0.155 e^{L/2}" is unreconciled** with `r31_attempts.json`, whose own Chebyshev floors are more negative.
  - Do not quote −0.155. The exponent-½ stopping point is unaffected.
  - Not re-derived this round. It stays OPEN as an erratum item.

## 8. Credits, literature
- Selberg (Λ_y, mollifiers, 1942/1946).
- Levinson (1974); Conrey (1989), mollifiers.
- Landau (1905), Laplace integrals of nonnegative functions.
- Titchmarsh, *Theory of the Riemann Zeta-Function*, §§10.25, 14.27–14.29 (DH zeros in σ > 1; Mertens-type and weak Mertens).
- Davenport–Heilbronn (1936).
- Odlyzko–te Riele (1985), Mertens disproof.
- Nyman (1950); Beurling (1955); Báez-Duarte (2003); Burnol (2002); Báez-Duarte–Balazard–Landreau–Saias (2000/2005); Vasyunin (1995).
- Bettin–Conrey (2013), period functions and cotangent sums.
- Vinogradov (1958) and Korobov (1958), zero-free region.
- Project: Alpöge's window architecture (F13W/F13S); Codex (fourth-power reduction); Turbo and Claude audits carried from R30-b/R31.
- **Novelty check.** The ingredients are classical. We have not seen these stated this way:
  - the Möbius-coordinate form of the F13S floor with R₀ = (log ξ)″(½);
  - the ζ² control separating "ψ-floor" from "Möbius floor";
  - the transfer-cost statement against these controls.

## Files
- `r32_pass.py` (segmented sieve, ζ / ζ² / F_θ to 10⁹) → `r32_pass_1e09.json`, `r32_pass_1e09_th0.0.json`, `pass.log`
- `r32_controls.py` (EPS, DH own μ_M, Λ_M to 10⁷) → `r32_controls_1e07.json`
- `r32_beurling.py` + `r32_beur_dfs.c` (Beurling primes; squarefree DFS to 10⁸) → `r32_beurling.json`. `beur_primes.bin` is regenerated, not committed.
- `r32_beur_zeta.py` → `r32_beur_zeta.json` (ζ_B near ½)
- `r32_constants.py` → `r32_constants.json` (Arb)
- `r32_identities.py` → `r32_identities.json`
- `r32_selberg.py` → `r32_selberg.json`
- `r32_analysis.py` → `r32_results.json`, `r32_series.npz`
- `r32_validate_zeros.py` → `r32_validate_zeros.json` (test data only)
- `r32_figs.py` → `fig1_mobius_floor.png` (main), `fig2_zeta_mobius_decomposition.png`, `fig3_selberg_fixed_level.png`, `fig4_mu_growth.png`
- `r32_common.py`, `r22_fejer.py`, `r22_models.py` (copied from R31)

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).
