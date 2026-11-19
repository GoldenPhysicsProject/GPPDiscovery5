# Round 24: synthesis. The Haar-positivity paper, everything we know, and new mechanisms killed fast

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).

Labels: PROVED, CHECKED, NUMERICAL, CONJECTURE, HEURISTIC, GAP, OPEN. Zeros are never proof inputs; they appear only in rows
marked *illustration*. **RH is OPEN; nothing here closes it.** Controls: F_θ (FTH; θ = ¼ unless stated), Epstein x² + 5y² (EPS),
Davenport–Heilbronn (DH). Notation: Ψ(t) = D(t) is the R23 γ = 0 slice (see §0).

## Plain-words summary

1. **The paper.**
   - Its abstract positivity theorems (Haar squares are positive type; GNS; Haar projection) are correct and standard.
   - Its Weil instance is a **restatement**, not a proof step. The "shadow-positive datum" (Def. 7.1) *assumes* q(A,A) ≥ 0. For
     ω = Weil's distribution that assumption *is* RH (by the paper's own Theorem 3.2). So Theorem 7.2(i) and Prop. 3.3 hold only
     conditionally on RH, and the second "⇔" in (8.1) is a tautology.
   - The forward-direction proof sketch of Theorem 3.2 contains a **false claim**: that the local terms are "always ≥ 0 by
     Tate's thesis". They are not. The archimedean term and every prime term are negative on the very tests where Weil
     positivity holds (§1.2, NUMERICAL table).
   - The paper gives no candidate Ω. I supply and test the natural one built from the integers (C6).
2. **A correction to our own record.** R23's D(L) is exactly −g(t) = Ψ(t), Masatoshi Suzuki's *screw function* of ζ (J. London
   Math. Soc. 2023, arXiv:2206.03682). The following R23 results are **Suzuki's theorems**, not new:
   - the closed form and Laplace formula (his Thm 1.1);
   - "RH ⇔ Ψ ≥ 0 ⇔ Ψ bounded" (Thms 1.6–1.7);
   - "RH ⇔ the box-Gram is PSD" (Thm 1.2, Kreĭn screw function);
   - the strip version (Thm 11.1).

   Also, "RH ⇔ e^{−Ψ} is an infinitely divisible characteristic function" is Nakamura–Suzuki (arXiv:2306.08317). What
   remains ours from R23:
   - the **ball-arithmetic certificate Ψ > 0.0214 on [log 2, log 10⁸]**, which extends Suzuki's Thm 4.1 (t₀ > log 2) to
     t₀ > 18.42;
   - liminf Ψ ≤ 0;
   - the fake controls.

   ROUND23.md now carries an attribution erratum.
3. **Eight candidate mechanisms** were stated exactly, tested against ζ and the fakes, and searched in the literature.
   - Every one is either RH-equivalent, circular, known, or broken.
   - The deepest output is a **locality obstruction** (C7). Every mechanism that sees only *local* data (Satake parameters,
     Hagedorn point, Λ ≥ 0, Sym^k/Deligne amplification, Euler-factor positivity) can catch F_θ. Such mechanisms are blind to an
     off-line zero of ζ, because ζ's local data are trivial (Satake parameter 1 at every p).
   - Every mechanism that sees *global* data (Weil/Suzuki/Gram/Cramér) is RH-equivalent and offers no extra leverage.
4. **Shortlist** (§5). Neither route proves RH. Both have a real chance of producing a *new unconditional theorem* that is fatal to
   the fakes.
   - (1) **A density version of the screw function.** A frequency-localised Ψ_{T,H}, whose negative mass is controlled by
     N(σ, T) and hence by Carlson/Ingham density theorems. It is F_θ- and DH-fatal and provable. Next lemma: L5.1.
   - (2) **Turning finite certificates into explicit zero-free boxes** (a quantitative Landau/Turán step for Ψ). Next lemma: L5.2.

## 0. The common object (synthesis of R13–R23 in one line)

Every positivity in R13–R23 is a face of **one positive measure μ = Σ_ρ δ_{γ_ρ}** (real iff RH), read through different test families:

| face | where | test family | what it is in terms of μ | status |
|---|---|---|---|---|
| Stieltjes ladder / S-fraction a_k > 0 | R20 | moments at u = 0 | Σ γ^{−2k−2} (Hamburger/Stieltjes) | RH-equiv. (all k) |
| Pick at s_{k,N} | R20 | Cauchy kernels | Σ 1/((u+γ²)(v+γ²)) | RH-equiv. |
| Fejér network ρ_X(γ₀) | R22 | sinc² bumps | μ ∗ Fejér | RH-equiv. |
| de Branges E_h, h = ½ frontier | R22 | Hermite–Biehler | sup \|β − ½\| < h | h < ½ = quasi-RH |
| Ψ(t) = D(L), box-Gram G | R23 = Suzuki | boxes 1_{[0,L]} | ∫(1−cos tγ)/γ² dμ | RH-equiv. (Suzuki) |
| Bost–Connes Hagedorn β = 1 | R21 | Euler product region | abscissa of Σ n^{−β} | local only (C4) |
| Tate theta / BPY law | R17 | Mellin of Φ | ξ = ∫ x^{−Δ} f_ζ | normalisability = RH (R17) |
| graviton Δ = 2 / "s IS Δ" | R17B | conformal primaries | (ψ_Δ, f_ζ) | restatement (R17) |
| Hamburger uniqueness | R19 | FE + Dirichlet series | ζ is the only member of its FE class | PROVED, classical |
| paired cancellations / entangled bulks | R21 | A–A pairs | prime terms negative in each re-centred moment | no-go per piece (R22) |

Pole + arch alone are not coercive (S12, R22); per-prime pieces are mean-zero (R22). So positivity is *only* global. That is the
locality obstruction, made precise in C7.

## 1. Paper audit: "Haar Positivity: From Weil to Wightman"

**Version note.** The brief names the primary version `cbf92c1f….tex` and an earlier `02d1d5b2….tex`.
- Neither .tex file exists on the box. The attachments folder holds other papers only.
- Two extracted text copies of the paper are in the bot's asset store: a complete 723-line copy, and a copy of lines 80–676 with
  the first 79 and last 48 lines elided.
- On their 597 overlapping lines they are **byte-identical** (`diff` empty), so I found no substantive difference, including
  in the Weil section.
- If the two .tex files really differ, the difference is outside lines 80–676 (the preamble/title block or the bibliography), or
  in a file I could not access. The audit below uses the complete copy.

### 1.1 Status of each statement

| item | paper's claim | verdict |
|---|---|---|
| Thm 2.2, Haar squares Ω^∨∗Ω are positive type | proved | **PROVED** (standard; correct proof) |
| Thm 2.3, GNS from positive type | proved | **PROVED** (standard, Dixmier 13.4.5) |
| Thm 3.2, Weil's criterion | cited, plus sketch | **Statement correct (Weil 1952 / Bombieri). Sketch wrong:** see 1.2 |
| Prop 3.3, Weil positivity ⇒ Hilbert space | proved | **PROVED but conditional:** its hypothesis *is* RH |
| Wightman, OS reconstruction (positivity parts) | cited | correct citations; restatements |
| Thm 6.1, Haar projection; Prop 6.2; Peter–Weyl | proved/cited | **PROVED** / standard |
| Thm 7.2(i), "Weil positivity instantiates the universal construction" | "proved" (Prop 8.1 iii) | **Restatement, conditional on RH.** Def. 7.1 *requires* ω(Θ(A)A) ≥ 0; for ω = D_ℚ that is RH |
| (8.1), RH ⇔ D(Ω∗Ω^∨) ≥ 0 ⇔ ‖Ω‖²_D ≥ 0 | equivalences | first ⇔ = Weil; **second ⇔ is a tautology** (definition of ‖·‖_D) |
| §8.3, "mass gap M > 0 is unconditional by Seiler cluster expansion and RG invariance" | asserted | **Unsupported, and a red flag.** The continuum Yang–Mills mass gap is the open Clay problem; cluster expansions give lattice/strong-coupling results only |

**Where RH enters, precisely.** In the Weil instance the Haar square is the *test* P = Ω∗Ω^∨, which is automatically of positive
type (Thm 2.2). What would be RH is that the *distribution* D_ℚ is of positive type. By Bochner–Schwartz, that means D_ℚ's
Fourier transform is a positive measure, which is μ = Σ_ρ δ_{γ_ρ} with all γ real.

A literal Haar-square factorisation D_ℚ = Ω^∨∗Ω with Ω a function on ℝ (or C_ℚ) is **impossible even under RH**. The FT of
D_ℚ is purely atomic, and no L²/L¹ function has |Ω̂|² equal to a sum of point masses. The factorisation exists only on the
**Bohr compactification** bℝ (a compact group, which fits the paper's compact-Haar theme), with Ω in Besicovitch B². That is
candidate C6.

Under RH, the paper's GNS construction applied to D_ℚ yields a unitary translation group whose generator has spectrum {γ_ρ}.
That is the Hilbert–Pólya operator, so the construction is *equivalent* to RH (folklore; cf. Connes 1999), not a route to it.

### 1.2 The false step in the forward-direction sketch (NUMERICAL, `r24_audit.json`)

The paper writes D(Ω∗Ω^∨) = Σ_ρ |Ω̂(ρ)|² + (local terms, always ≥ 0 by Tate's thesis). In fact:
* D(P) *is* the spectral sum. The local terms are its other expression, not an extra summand.
* The local terms are **not** non-negative. Take the box tests (Suzuki/R23), where D = pole + arch + prime:

| L | pole | arch | prime | total D = Ψ |
|---|---|---|---|---|
| 0.69 | 0.481 | −0.418 | 0 | 0.063 |
| 2 | 4.345 | −2.545 | −1.746 | 0.053 |
| 5 | 41.06 | −9.46 | −31.55 | 0.045 |
| 12 | 1605.7 | −27.9 | −1577.8 | 0.029 |

The archimedean and prime local terms are negative; only the pole term is positive. Positivity of the total is a cancellation
to 1 part in 5·10⁴ at L = 12. This matches S12 (pole + arch not coercive past h ≈ 0.40) and the R22 per-prime mean-zero lemma.

The converse sketch is also mis-stated. Off the line the spectral term is Ω̂(ρ)·conj Ω̂(1 − ρ̄), not |Ω̂(ρ)|² + |Ω̂(1 − ρ̄)|². The
correct converse is Weil's, by choosing Ω so that this cross term is negative and dominant.

**Verdict.** As a *conceptual* unification (one GNS mechanism, four languages) the paper is sound. As a contribution to RH it adds
no inequality. It relocates RH into the hypothesis of Def. 7.1. Fix the Thm 3.2 sketch, and either mark the Yang–Mills claim
conditional or remove it.

## 2. Literature synthesis: why proofs fail (one line each)

* **Weil, function fields:** positivity = Hodge index on C×C, which needs an ample class on a surface. ℤ has no such surface
  (R15: Hasse Gram ⇔ RH).
* **Deligne:** the tensor-power trick amplifies |α| ↦ |α|^k while the a-priori Rankin–Selberg bound keeps a *fixed* +1 (pole edge).
  The zeros of ζ are not Satake parameters, so there is nothing to raise to the k-th power (C1).
* **Connes trace formula / semilocal:**
  - archimedean Weil positivity is proved for support ⊂ [2^{−½}, 2^{½}] (Connes–Consani, Selecta 2021);
  - the team has Arb certificates to full support 2.6 (s36);
  - the global step needs positivity uniform in the set of places. OPEN.
* **Connes–Consani prolate/Sonin (2021–25):** Sonin-space traces give positive functionals whose difference from W is controlled
  only on small windows.
* **Bombieri / Yoshida / Suzuki:** positivity on small supports is unconditional; extension is RH.
* **Li/Keiper; Bombieri–Lagarias:** λ_n ≥ 0 ⇔ RH. λ_n is the Weil form on Laguerre-type tests, another face of μ.
* **Lagarias:** de Branges spaces and Hilbert spaces of entire functions.
* **de Branges' attempts:** the positivity condition he needed is false for ζ (Conrey–Li 2000).
* **Hilbert–Pólya / Berry–Keating xp:** semiclassical counting works; the self-adjoint realisation is equivalent to RH.
  Yakaboylu (arXiv:2211.01899, 2309.00405, 2408.15135) states that self-adjointness of his similarity-transformed H_BK "still needs
  rigorous treatment".
* **Lee–Yang / Newman:** the de Bruijn–Newman constant satisfies Λ ≥ 0 (Rodgers–Tao 2020) and Λ ≤ 0.2 (Platt–Trudgian 2021;
  Polymath 15: 0.22). RH ⇔ Λ = 0, so ζ sits **exactly at criticality**.
* **Jensen polynomials (GORZ 2019):** hyperbolic for fixed d and large n, for a wide class. It **does not separate** ζ from EPS
  at d ≤ 8 (R17).
* **Nyman–Beurling / Báez-Duarte:** RH ⇔ 1 ∈ the closed span of {α/t} − α{1/t}. This is the classical additive×multiplicative
  coupling.
* **Goldbach average (Granville 2007; Bhowmik–Schlage-Puchta; Bhowmik–Ruzsa):** Σ_{n≤x} G(n) = x²/2 + O(x^{3/2+ε}) ⇔ RH. The
  k-prime version is also RH-equivalent.
* **Kurokawa absolute tensor product ζ^{⊗r}** has zeros at ρ₁ + … + ρ_r. Koyama–Kurokawa give an Euler product for ζ^{⊗2}, valid
  only in Re s > 2, so the region scales with r.
* **Suzuki 2023; Nakamura–Suzuki 2023:** screw function Ψ; e^{−Ψ} infinitely divisible ⇔ RH.

## 3. Candidates: stated, attacked, tested, searched

Numerical files: `r24_tests.py` with `r24_{audit,exponent,cells,levy,cramer}.json` and logs. ζ uses the integers only; the fakes
use the R22 explicit formula from their own Dirichlet data.

### C1. Deligne's tensor-power trick in the integer setting (brief 3a)
**Statement.** Find a family L_k with zeros at ½ + k(ρ − ½), each with an Euler product absolutely convergent on a *k-independent*
half-plane. Then GAP-23 at any ε < ½ amplifies: k(β − ½) < C for all k forces β = ½.

**Attack.**
* (i) **GL(1) rigidity (PROVED, elementary).** Every Schur functor applied to ζ's Satake parameter (α_p = 1) is trivial, so every
  functorial lift of ζ is a product of ζ's (and shifts ζ(s + iτ) for Eisenstein tensors). Zeros are *translated*, never
  dilated, so the amplification factor is exactly 1.
* (ii) **Dilation scales the a-priori region too (PROVED).** ξ(½ + (s−½)/k) has the desired zeros, but its Dirichlet series
  converges only for Re s > ½ + k/2.
* (iii) **Kurokawa tensor power / Goldbach.** The k-prime additive convolution produces the zeros ρ₁ + … + ρ_k (Fujii;
  Granville Thm 1D). With Abel weights, Σ_n (Λ−1)^{∗k}(n) e^{−n/X} = F(e^{−1/X})^k exactly, where F(z) = Σ(Λ(n) − 1)z^n. The
  additive tensor power is literally the k-th power, the a-priori bound X^k scales with it, and no fixed "+1" survives.
  Granville's RH ⇔ x^{k−1/2+ε} is the honest residue (known).
* (iv) It **does** amplify F_θ. Sym^k of Satake p^{±θ} is p^{±kθ}, so the defect is local.

**Tests.**
* (iv) confirmed (NUMERICAL, `r24_exponent.json`): the oscillation envelope of Ψ_FTH grows like e^{0.241 L} for θ = 0.25 and
  e^{0.092 L} for θ = 0.125, against ζ's 0.002 (flat). The fit matches Suzuki Thm 11.1.
* F_{0.125} stays *positive* to X = e^{11.5}. Its defect is visible only after θ ↦ 2θ (a Sym²-type step).
* Schur powers G^{∘k} of the box-Gram (`r24_levy.json`), min eig/max eig for k = 1..4:

| model | k = 1 | k = 2 | k = 3 | k = 4 |
|---|---|---|---|---|
| ζ | +2.9e−4 | +2.6e−3 | +3.6e−3 | +3.3e−3 |
| FTH | −0.23 | −0.21 | −0.27 | −0.27 |
| EPS | −0.073 | −0.061 | −0.074 | −0.088 |
| DH | −1.8e−3 | −9.3e−4 | **+4.9e−4** | **+2.1e−4** |

Schur ("tensor") powers *hide* DH at k ≥ 3. Positivity of high tensor powers is strictly weaker than positivity.

**Literature.** Deligne (Weil I); Kurokawa 1992; Koyama–Kurokawa (ζ^{⊗2} Euler product in Re s > 2); Granville 2007.

**Fate: dead (structural).** A GL(1) tensor power cannot amplify a global defect. The trick catches exactly the local fakes.

### C2. Castelnuovo/Hodge index on ℕ×ℕ: the additive-cell Weil Gram (brief 3b)
**Statement.** Let I_n = [log n, log(n+1)] be the additive unit cells in multiplicative coordinates. The Weil Gram is
M_N(n,m) = Ψ(log((m+1)/n)) + Ψ(log(m/(n+1))) − Ψ(log(m/n)) − Ψ(log((m+1)/(n+1))), for n, m ≤ N. This is the box polarisation of
Suzuki's kernel; it couples + (cells n → n+1) with × (log). Hodge index analogue: M_N ⪰ 0 for all N.

**Attack.**
* RH ⇒ M_N ⪰ 0 (Suzuki Thm 1.2 restricted).
* Conversely, M_N ⪰ 0 for all N ⇒ Ψ(log N) ≥ 0 for all N (sum of all cells). If RH failed, Ψ would be unbounded below
  (Suzuki 1.6). Ψ' = O(e^{L/2}) and the grid spacing is e^{−L}, so Ψ(log N) < 0 for some N.
* So the family is **RH-equivalent** (CHECKED, sketch). The "ample class" that makes Hodge index work on C×C has no counterpart:
  the pole term is the only positive local piece, and it is rank 2 (R15, S12).

**Tests (NUMERICAL, `r24_cells.json`, correlation-normalised min eigenvalue):**
- ζ is PSD to N = 600 (min 2.0e−6, decaying smoothly);
- FTH is indefinite from N = 30;
- EPS is indefinite from N = 81;
- DH stays PSD to N = 600, as with the R23 slice. DH's off-line zeros sit at γ ≈ 85.7 and need longer windows.

**Literature.** Nyman–Beurling (step functions); Suzuki's kernel on L²(−a,a); Yoshida. I found no explicit statement of the
cell Gram, but it is a restriction of a known RH-equivalent kernel.

**Fate: dead as a mechanism** (reformulation with no new input). It is kept as a clean "+/× coupling" form of Weil positivity.

### C3. Lévy process with exponent Ψ (brief 3c)
**Statement.** Under RH, e^{−tΨ(L)} = E e^{iL Z_t}, where Z is the compound-Poisson process with jumps ±γ_ρ at rate 1/γ_ρ². Its
total intensity is Σ 1/γ² = 2c₀ = 0.0462, which equals the mean of Ψ. This is dual to Khinchine's zeta distribution: for σ > 1,
t ↦ ζ(σ+it)/ζ(σ) is infinitely divisible with **prime** jumps (log n), in the frequency variable. At σ = ½, Ψ would make the
**zeros** the jumps, in the position variable.

**Attack.**
* RH ⇔ e^{−Ψ} is an infinitely divisible characteristic function. That is **Nakamura–Suzuki 2023, exactly** (their g_ζ = −Ψ).
  KNOWN.
* New lemma (PROVED): **no factor-wise Lévy structure.** Ψ = A − S with A = pole + arch. A conditionally negative definite
  function is O(L²), but A ~ 4e^{L/2} and S ~ 4e^{L/2}. So neither factor e^{−tA} nor e^{tS} is positive definite. The
  "Euler-product factorisation" e^{−tΨ} = e^{−tA} Π_n e^{tΛ(n)n^{−½}(L − log n)₊} has no positive-definite factors. This matches
  Queen's spec: the Lévy measure Σ Λ n^{−½} δ_{log n} is not a Lévy measure on ℝ.

**Tests (NUMERICAL, min eigenvalue of [e^{−tΨ(|L_i − L_j|)}], h = 0.05, L ≤ log 10⁵):**
- ζ is positive for t ∈ {0.3, 1, 3, 10, 30};
- FTH and EPS are negative for all t;
- DH is negative for t ≤ 3 and **positive for t ≥ 10**. Large t is blind.

**Fate: dead (known, RH-equivalent).** The prime ↔ zero Lévy duality is a nice dictionary, but it gives no lever.

### C4. de Bruijn–Newman heat flow as physical time; graviton/Hagedorn fixes Λ = 0 (brief 3d)
**Statement.** The BC/Hagedorn transition at β = 1, or the massless graviton Δ = 2, forces the heat-flow "temperature" Λ to 0.

**Attack.**
* The Hagedorn point is the abscissa of absolute convergence of Σ a_n n^{−β}. It is 1 for ζ **whatever the zeros are**: it is
  fixed by the pole and by Satake temperedness (PROVED).
* For F_θ it is 1 + θ, so it sees F_θ's local defect. DH and EPS have no BC system.
* Newman's Λ lives on the critical line, the Hagedorn point on σ = 1. These are different thresholds, h = 0 versus h = ½ in R22
  language.
* Rodgers–Tao (Λ ≥ 0) says ζ is *exactly* critical. Any thermodynamic argument must prove marginality with zero slack, and no
  perturbative mechanism can do that.
* H_t has no Euler product for t ≠ 0, so the integer structure is lost at once.

**Tests (R17, `r17_heat.json`):** heat-flow collision times give Λ_EPS ≥ 0.079, Λ_DH ≥ 0.044, Λ_FTH ≥ 0.020; ζ has none; Laguerre
passes only for ζ (and FTAU).

**Fate: dead.** Hagedorn sees only local data. The graviton/Δ = 2 identity is a restatement (R17B).

### C5. Bootstrap: GAP-23 at ε plus multiplicativity ⇒ smaller ε (brief 3e)
**Statement.** A map ε ↦ ε' < ε on the exponent of Ψ ≥ −C e^{εL}, using the Euler product and the FE.

**Attack.**
* **Fixed-point obstruction (PROVED).** For every θ ∈ (0, ½), F_θ = ζ(s−θ)ζ(s+θ) has an Euler product, a functional equation
  (its own Γ-factor), Λ_θ ≥ 0, order 1 and finite conductor, and sup(β − ½) = θ exactly. So no bootstrap that uses only those
  axioms can exist. Numerically the exponent is recovered: 0.241 for θ = ¼.
* Suzuki's monotonicity (Ψ_ω ≥ 0 on [0, t₀) ⇒ Ψ_{ω+η} ≥ 0) propagates positivity **upward** in ω only, which is the wrong direction.
* What F_θ lacks is **temperedness** (|α_p| = 1; a_n ≪ n^ε). With temperedness added we are in the Selberg class, where no fake
  is known but GRH is open.
* Hamburger (R19): with ζ's exact FE the class is {cζ}. So an "axiomatic" bootstrap must use the specific Γ_ℝ and conductor 1.
  Every fake changes exactly those.

**Fate: dead without temperedness. With temperedness it is the Selberg-class GRH (OPEN).**

### C6. The paper's Ω, built from integers: a Bohr-compact Haar square (brief 1 and 3f)
**Statement.** E(x) := e^{−x/2}(ψ(e^x) − e^x). Under RH, E = −Σ_ρ e^{iγx}/ρ + O(e^{−x/2}) is Besicovitch almost-periodic. Its
Haar (Bohr-mean) autocorrelation is
  M_x[E(x)E(x+y)] = Σ_ρ cos(γy)/|ρ|² (= Σ_ρ e^{iγy}/|ρ|² over the full symmetric zero set),
the Weil distribution smoothed by 1/|ρ|². So the paper's P = Ω^∨∗Ω **holds literally on bℝ with Ω = E, an integer object**.
Positivity is free, because it is an autocorrelation. The content moves to *size*: RH ⇔ the Bohr mean M[E²] is finite. That is
Cramér 1922 (⇒) plus Landau (⇐); Brent–Platt–Trudgian 2020 note that I(X)/X² is unbounded if RH fails.

**Attack.** Known equivalence, so it is RH-equivalent. GAP-23 becomes ∫₀^X E² ≪ e^{2εX} with ε < ½: quasi-RH again.

**Tests (NUMERICAL, `r24_cramer.json`; ψ from a sieve to 10⁸, fakes to 10⁶).**
- ζ: the fit M_X[E²] = a + C/X (X ≥ 9) gives **a = 0.04594** from primes alone. The exact Σ_ρ 1/|ρ|² = 2 + γ_E − log 4π =
  0.04619 (illustration).
- The autocorrelation matches the 400-zero side (illustration): 0.0073 vs 0.0065 at y = 3; 0.0163 vs 0.0171 at y = 12.
- Fakes: M[E²] grows. FTH: 2.06 → 5.87 (x = 7 → 13). EPS: 2.05 → 78.8. **DH: 0.166 → 0.397**, rising from x ≈ 9 (DH's zero
  0.8085 + 85.70i: e^{0.31x}/86).
- This is the only finite-range integer test this round that **detects DH by X ~ 10⁵–10⁶**, where the Ψ-slice and the cell
  Gram do not.

**Fate: dead as a proof (RH-equivalent, classical).** Alive as the correct reading of the paper's Haar square, and as the best
cheap fake-detector.

### C7. Locality obstruction (synthesis; brief 3f)
**Statement (HEURISTIC meta-lemma, with PROVED instances).** Sort mechanisms by the data they consume.

| mechanism | data | catches F_θ? | catches DH/EPS? | sees an off-line ζ zero? |
|---|---|---|---|---|
| Sym^k / Deligne (C1) | Satake | yes (kθ) | n/a (no Euler) | **no** (ζ's Satake = 1) |
| Hagedorn / BC (C4) | abscissa | yes (1 + θ) | n/a | **no** |
| Λ ≥ 0 / Markov (Queen spec) | local signs | **no** (F_θ is Markov) | partly | no |
| Euler-factor PD (s12 X1) | local unitarity | yes | n/a | no; S12 obstruction |
| pole + arch coercivity | archimedean | — | — | fails past h ≈ 0.40 (S12) |
| Weil/Suzuki/Gram/Cramér (C2, C3, C6) | global | yes | yes | yes, but **RH-equivalent** |

ζ's local data are as tame as possible. So any proof must convert **local tameness** (Euler product + temperedness + Γ_ℝ +
conductor 1, which by Hamburger pins ζ down) into **global positivity**, with something stronger than Landau. No mechanism in the
literature or in R13–R24 does that. Over F_q, the converter is the Hodge index on a surface; over ℤ it is missing (R15/R16/R19).

**Fate:** not a proof route. It is the map of where one must be.

### C8. Schur/tensor positivity as a weaker but provable target (brief 3f)
**Statement.** Prove G^{∘k} ⪰ 0 for some k ≥ 3, which is strictly weaker than RH.

**Attack.** DH passes at k = 3, 4 on our window (C1 table), so the statement is not fake-fatal at finite range. It also gives no
zero information: G^{∘k} ⪰ 0 corresponds to positivity of μ^{∗k}, which allows off-line zeros. **Dead** (weaker than RH and not
discriminating).

## 4. Fate table

| # | candidate | RH-equivalent / circular / broken | fakes fail? | known? | fate |
|---|---|---|---|---|---|
| C1 | Deligne tensor power | broken: amplification factor 1 on GL(1) (PROVED) | catches only F_θ | Deligne; Kurokawa; Granville | dead |
| C2 | ℕ×ℕ cell Gram (Hodge index) | RH-equivalent (CHECKED) | FTH N = 30, EPS N = 81, DH not by 600 | restriction of Suzuki | dead (kept as a form) |
| C3 | Lévy, exponent Ψ | RH-equivalent | FTH, EPS all t; DH only t ≤ 3 | Nakamura–Suzuki 2023 | dead |
| C4 | dBN + Hagedorn | broken: Hagedorn is local | F_θ only | Newman, Rodgers–Tao | dead |
| C5 | ε bootstrap | broken by the F_θ fixed-point family (PROVED) | — | Selberg class | dead / OPEN with temperedness |
| C6 | paper's Ω = Chebyshev error on bℝ | RH-equivalent | **all three, DH by 10⁶** | Cramér, Landau | dead as proof; best detector |
| C7 | locality obstruction | meta | — | new synthesis | map |
| C8 | Schur positivity | weaker than RH | DH passes k ≥ 3 | — | dead |

## 5. Ranked shortlist: the routes with a real chance (of a *new theorem*, not of RH)

**Route 1: a density theorem for the screw function (F_θ- and DH-fatal, provable).** Frequency-localise Ψ by applying the
explicit formula to Fejér/Beurling–Selberg windows at height T:
  Ψ_{T,H}(t) := Σ_ρ φ((γ − T)/H)·(e^{(ρ−½)t} − 1)/(ρ−½)², computed from the integers by the R22 machinery.
Off-line zeros in [T, T+H] make Ψ_{T,H} grow. The *number* of such zeros is controlled unconditionally for ζ by zero-density
theorems: Selberg N(σ, T) ≪ T^{1−(σ−½)/4} log T, Ingham N(σ, T) ≪ T^{3(1−σ)/(2−σ)} log⁵ T, so the count is o(T) for every σ > ½. In contrast, DH has ≫ T zeros in Re s > 1
(Saias–Weingartner, Acta Arith. 2009, for periodic-coefficient Dirichlet series without Euler product; cited as in Queen's spec, re-verify before use) and F_θ has ≍ T log T. This is Queen's CONJECTURE N, recast on Suzuki's Ψ.

*Next lemma (L5.1, Ψ mean-value dictionary).* For δ > 0, H ≥ T^{c} and t ∈ [L, 2L],
  ∫_L^{2L} |Ψ_{T,H}(t)|² e^{−2δt} dt = Σ_{ρ: γ∈[T,T+H], β>½+δ} w_ρ(L) + O_δ(H log T · L),
with explicit weights w_ρ ≍ e^{2(β−½−δ)L}/|ρ−½|⁴. It follows from the Montgomery–Vaughan Hilbert inequality, using
γ-spacing ≫ 1/log T on average. Combined with Ingham's density theorem, it gives an unconditional, ζ-specific bound on the
off-line mass seen by Ψ, which fails for DH and F_θ.

**Route 2: from finite certificates to explicit zero-free boxes (quantitative Landau/Turán for Ψ).** R23 certified Ψ > 0.0214
on [log 2, 18.42]. Landau's theorem is qualitative.

*Next lemma (L5.2).* If ξ(ρ₀) = 0 with β₀ = ½ + δ, then for every T ≥ T₀(δ, γ₀) there is t ∈ [T, T + K(γ₀)] with
Ψ(t) ≤ −c·e^{δt}/|ρ₀ − ½|² + O(Σ_{ρ≠ρ₀} …). The constants come from Turán's second main theorem on power sums, with
K ≍ log(γ₀)/δ. Feeding in the certified positivity on [0, t₀] then excludes zeros in an explicit box {β > ½ + δ, |γ| ≤ H(δ, t₀)}.

This is honest but modest: H(δ, 18.4) will be far below the Platt–Trudgian verification height. Its value is that it is
**formalizable end to end**: certificate (R23) → Turán → box, which suits the verify/Lean track. It also tells us how far
the certificate must be pushed for the box to become competitive.

**What would actually be needed for RH:** a converter from local tameness to global positivity (C7), i.e. an arithmetic
surface with an ample class. Neither route provides it; both are the strongest honest next steps I can see.

## Files
* `r24_tests.py`, with parts `audit`, `exponent`, `cells`, `levy`, `cramer`. Outputs `r24_*.json` and `*.log`.
* `r22_models.py`, `r22_fejer.py`, `r23_slice.py`: copied unchanged from rounds 22/23.

References: Suzuki, J. London Math. Soc. 108 (2023) (arXiv:2206.03682) · Nakamura–Suzuki, arXiv:2306.08317 · Connes–Consani,
Selecta Math. 27 (2021) · Granville, Funct. Approx. 37 (2007) · Bhowmik–Ruzsa, arXiv:1711.06442 · Kurokawa (1992); Koyama–Kurokawa,
CNTP 3 (2009) · Brent–Platt–Trudgian, arXiv:2008.06140 · Cramér, Math. Z. (1922) · Rodgers–Tao (2020) · Platt–Trudgian (2021) ·
Griffin–Ono–Rolen–Zagier, PNAS (2019) · Conrey–Li, IMRN (2000) · Yakaboylu, arXiv:2211.01899, 2309.00405, 2408.15135 · Weil (1952) ·
Deligne, Publ. IHES 43 (1974).
