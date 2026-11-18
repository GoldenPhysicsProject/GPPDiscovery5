# Round 14E — the Weil form on the complement of the window-truncated range of E

Framing credit: **Daniel Toupin** (No-Ghosts, 14B; one global object, 14D). Audit standard and the reframing used here: **Golden Eye** (AUDIT14, team/golden-eye 062e4e2).
Labels: PROVED / CHECKED (code or identity check) / NUMERICAL (evidence) / INTERPRETATION / OPEN. A control counts only if it could have come out the other way. Zeros are never proof inputs; they appear only as CHECK inputs in §5.

**RH: OPEN. Nothing here moves it.**

## 0. Verdict first
- **No cut-independent positive gap survives for ζ.** For every L ≥ 1.5, λ_min of W_ζ on the complement of the data-built truncated range S_L(ε) sits at the float floor (|λ| ≤ 6e-14, random sign), for every cut ε ∈ {1e-10, 1e-8, 1e-6} and both bases (fmax 45, 90). The only positive values are at L ≤ 1.25, where S is empty and the number is just λ_min(W) (the Connes–Consani archimedean regime).
- Golden Eye's warning is confirmed: W_ζ's small eigenvalues decay geometrically through any threshold, and the near-null count is basis-dependent (L = 8: 97 at fmax 45 vs 176 at fmax 90).
- **A uniform-in-L gap is falsified at the level the data can see.** At L = 8 a window probe built at the Lehmer pair (t ≈ 7005.06, gap 0.0377) has quotient Rayleigh quotient ≤ 0.063 (elementary bound, §5), against ≥ 2 for the band-limited bulk at fmax 45 and ≤ 33 for a median-gap pair at the same height. Any "gap" read off a band-limited Galerkin is a band artifact.
- **What does survive (cut- and basis-stable):**
  1. Sylvester inertia invariance for the fakes under quotienting by their *own* E (a pipeline check that can fail, and does fail for a borrowed range): EPS, F_θ, DH keep exactly their ghost counts.
  2. ζ's λ_min(W) at small L, converged in M: λ_min ≈ exp(−15.9 L² + 3.1 L − 0.9) on 0.25 ≤ L ≤ 1.25 (Gaussian in L, i.e. the CC "very small eigenvalues"; qualitatively known).
  3. The number of non-near-null modes, M − #{λ < 1e-8‖W‖}, is ~L-independent at fixed band: 26–29 (fmax 45), 62–64 (fmax 90), against 2N(F) = 16 and 50. That is Connes' semiclassical counting (known), not new.
- The radical check of 14D is redone with a **bilinear** probe and a same-band non-radical control: separation 1.3e2 to 1.2e4. That is now a discriminating CHECK (§2). Errata blocks for 14C and 14D were added to the tops of those files.

## 1. Construction (data only)
- Window g on u = log x ∈ [−L/2, L/2]; M Dirichlet sines, M = ⌈fmax·L/π⌉ + 8; W_F from 14B's `weil_galerkin.py` (own Λ_F(n), own Γ-factor, own poles). The Galerkin matrix is the **exact** restriction of W to the span of those sines.
- E_F H(u) = e^{u/2} Σ_n a_F(n) H(n e^u), with (E_F H)^(t) = F(½+it) Ĥ(t). F on the line is computed from Hurwitz sums (`clineval.py` → `cline_*.npy`, complex, t ≤ 250).
- Family H: bump((v−c)/w)·cos(k(v−c)+φ), with widths w ∈ {1.5, 1.0, (0.6)}, centres stepping w/2 across the window, k ≤ 1.1–1.2·fmax. Pole constraints Ĥ(½+a) = 0 at F's poles.
- Window coefficients a_i = ⟨φ_i, E_F H⟩ by Parseval from F on the line. Re of the integrand is even in t, so the trapezoid on [0, T] is spectrally accurate (CHECKED: t-steps 0.01/0.02/0.04 agree to 3 digits).
- **S_L(ε)** = span of P_window E_F H over the family directions with ‖P_M P_window E_F H‖² ≥ (1−ε)‖E_F H‖². This is a generalised eigenproblem (prolate-type concentration). Then W is compressed onto S^⊥ (`gap.py`).

## 2. Radical check, redone bilinear (fixes 14D 3a) — CHECKED
`bilinear.py`, `bilinear.json`. Quadratic residual = max|eig(SᵀWS)|/‖W‖. Bilinear residual = ‖WS‖₂/‖W‖ = sup over v ∈ S and φ ∈ V_M of |W(v,φ)|/(‖v‖‖φ‖‖W‖). Control: the same family pushed through F' = F − 0.1·2^{−s} (a(2) shifted by −0.1: same band, not in the radical).

| model | L | ε | dim S | own bilinear | control bilinear | separation |
|---|---|---|---|---|---|---|
| ζ | 4 | 1e-8 | 40 | 1.1e-5 | 6.8e-2 | 6.4e3 |
| ζ | 6 | 1e-8 | 80 | 9.6e-6 | 1.2e-1 | 1.2e4 |
| EPS | 6 | 1e-8 | 9 | 9.2e-5 | 1.2e-2 | 1.3e2 |
| F_θ | 4 | 1e-8 | 14 | 6.7e-6 | 3.2e-2 | 4.7e3 |
| DH | 4 | 1e-8 | 14 | 2.8e-5 | 2.1e-2 | 7.4e2 |

At ε = 1e-6 the separations drop to 47–1.9e3, because looser concentration lets more truncation in. **Range E_F (window-truncated) is in W_F's radical to 1e-5 bilinear, for every model with an explicit formula, fakes included.** That is the easy direction of the explicit formula (PROVED). It does not single out ζ.

## 3. (i) Pair test = Sylvester inertia test (Golden Eye falsifier 1) — CHECKED
`run_gap.py pair` → `pair.json`/`pair.log`; `fits.py` → `fits.json`. Ghost = eigenvalue < −1e-7‖W‖ (counted with multiplicity, as in 14B erratum d). Quotienting by a subspace of the radical cannot change the inertia, so any change would expose S ⊄ radical.

| model (own E) | L | M | dim S (ε = 1e-10 / 1e-8 / 1e-6) | ghosts W → W on S^⊥ | λ_min W → S^⊥ |
|---|---|---|---|---|---|
| EPS | 4 | 85 | 0 / 0 / 0 (**vacuous**) | 4 → 4 | −1.624 → −1.624 |
| EPS | 6 | 123 | 7 / 9 / 10 | 10 → 10 | −8.644 → −8.644 |
| F_θ | 2 | 37 | 0 / 0 / 0 (**vacuous**) | 7 → 7 | −0.0908 → −0.0908 |
| F_θ | 4 | 66 | 12 / 14 / 16 | 20 → 20 | −0.755 → −0.755 |
| DH | 4 | 168 | 12 / 14 / 15 | 2 → 2 | −0.7554 → −0.7554 |
| DH | 6 | 247 | 80 / 83 / 86 | 4 → 4 | −3.675 → −3.675 |
| BEUR″ (with ζ's E: see note) | 4.5 | 94 | 48 / 51 / 53 | 10 → 8 / 7 / 5 | −0.617 → −0.36 |
| BEUR″ (with ζ's E) | 6 | 123 | 77 / 80 / 83 | 40 → 6 / 3 / 3 | −2.41 → −1.1 to −1.2 |
| ζ (control) | 4, 6 | 85, 123 | 38–83 | 0 → 0 | floor → floor |

- **EPS, F_θ, DH keep their ghosts exactly** under their own quotient, for every cut. This is non-vacuous where dim S > 0: EPS L=6, F_θ L=4, DH L=4 and 6. At EPS L=4 and F_θ L=2 (Golden Eye's suggested points), no family vector reaches the concentration cut, so the test is vacuous there; I report that rather than hide it. That EPS has fewer concentrated vectors fits its conductor 20: longer E tails.
- **BEUR″: the borrowed-arch defect survives** (3–8 ghosts remain, λ_min ≈ −0.36 to −1.2). The inertia *does* change, and the test catches it: ζ's range is not in BEUR″'s radical (bilinear residual 0.10–0.41). BEUR″ has no admissible own E. Its own Σ_B H(nx) carries x^{−ρ} terms from the Re ρ = 0.6 quartet, so g grows like e^{0.1|u|} as u → −∞ unless Ĥ vanishes at the quartet (zeros as input). So "own quotient" is undefined: DEFECT / N.A., as in 14B/14D.
- Why this is a real control and not a tautology: it fails visibly (BEUR″) when S is not in the radical.

## 4. (ii) ζ: λ_min on S^⊥ against cut and basis (Golden Eye falsifier 2) — NUMERICAL
`run_gap.py zeta` → `zscan.json`/`zscan.log`; small L in `fits.json`.

| L | λ_min, fmax 45 (any ε) | λ_min, fmax 90 (any ε) | near-null count 45 / 90 | non-near-null modes M − count, 45 / 90 |
|---|---|---|---|---|
| 0.25 | 0.333 | 0.329 | 0 / 0 | 12 / 16 |
| 0.5 | 3.72e-2 | 3.60e-2 | 0 / 0 | 16 / 23 |
| 0.75 | 5.80e-4 | 5.55e-4 | 0 / 0 | 19 / 30 |
| 1.0 | 1.08e-6 | 1.06e-6 | 0 / 0 | 23 / 37 |
| 1.25 | 3.41e-10 | 3.04e-10 | 0 / 0 | 26 / 44 |
| 1.5 | 4.6e-15 (floor) | 3.2e-15 (floor) | – | – |
| 2 | floor (−1.9e-14) | floor | 8 / 8 | 29 / 58 |
| 3 | floor | floor | 22 / 30 | 29 / 64 |
| 4 | floor | floor | 38 / 59 | 28 / 64 |
| 5 | floor | floor | 52 / 88 | 28 / 64 |
| 6 | floor | floor | 67 / 117 | 27 / 63 |
| 7 | floor | floor | 82 / 147 | 27 / 62 |
| 8 | floor | floor | 97 / 176 | 26 / 62 |

- For L ≥ 2 the data-built S catches all but 7–9 of the near-null directions. The rest (window-edge and band-edge directions; the eigenvectors at L = 6 peak at frequency ≈ fmax) stay in S^⊥ at the floor, for every ε and both bases. **So λ_min(W|S^⊥) carries no cut-independent positive number for L ≥ 1.5.** The second eigenvalue on S^⊥ is also floor-to-1e-11 and ε-dependent.
- W_ζ's spectrum below the bulk is a geometric ladder with no plateau, e.g. L = 2: 6.9e-12, 1e-9, 1.2e-7, 1.2e-5, 4.2e-4, 1.3e-2, 0.33. Loosening ε eats it rung by rung (L = 6: ε = 1e-4 gives λ_min 1.5e-11, ε = 1e-2 gives 3.8e-8). There is no cut to make.
- **Fits.**
  - ζ, 0.25 ≤ L ≤ 1.25: log λ_min = −15.85 L² + 3.03 L − 0.86 (fmax 45, max residual 0.07) and −15.98 L² + 3.15 L − 0.90 (fmax 90, 0.03). A pure exponential fails (residual 2.0). Gaussian decay in L reaches the 1e-15 floor at L ≈ 1.5, consistent with that row. Basis-converged to ≤ 5% (11% at L = 1.25).
  - The fakes on their own S^⊥: |λ_min| = exp(a + bL) with b = 1.10 (EPS, L 3.5–6), 0.89 (F_θ, L 2–5), 0.79 (DH, L 4–6). Max log-residual 0.35 / 0.26 / 0.07. Ghosts deepen roughly like e^{L}; ζ shows none at any L.
- The non-near-null mode count is ~L-independent at fixed band (26–29 and 62–64). It sits about 10–14 above 2N(fmax) (16 and 50; zeros used only to compare). That is Connes' semiclassical count of the quotient (known; NUMERICAL here).

## 5. (iii) Lehmer-pair probe (Golden Eye falsifier 3) — PROVED bound + CHECKED explicit formula + NUMERICAL
`lehmer.py`, `lehmer_opt.py` → `lehmer.json`, `lehmer_opt.json`.

**Bound (PROVED, elementary).** Let V be any subspace whose elements e have ê(γ₁) = ê(γ₂) = 0 at two zero ordinates γ₁ ≠ γ₂ (true for range E_ζ, where Mellin = ζ·Ĥ). For window g:
|(g−e)^(γ₁) − (g−e)^(γ₂)| ≤ |γ₂−γ₁| ∫|u||g−e| ≤ |γ₂−γ₁| √(L³/12) ‖g−e‖, hence
**W(g)/dist(g,V)² ≤ W(g)·(γ₂−γ₁)²·L³ / (12·|ĝ(γ₁)−ĝ(γ₂)|²).**
If the closure is taken in the form-domain norm (Golden Eye), dist only grows, so the bound still holds.

**Explicit formula at height 7005 (CHECKED).** For fixed probes g = cos²-window·(cos γ₁u − cos γ₂u), W(g) computed from integer data (Λ(n), Γ_ℝ, pole; FFT autocorrelation at du ≈ 3.6e-6) agrees with the zero side Σ±|ĝ(γ)|² to 4 digits at L = 5, 7, 8. Examples: Lehmer pair at L = 7, 6.330e-3 vs 6.328e-3; the first two zeros at L = 8, 16.01 vs 16.01.

**Best probe in the window band** (maximise |ĝ(γ₁)−ĝ(γ₂)|²/W over sines within ±15 of the pair; zeros #6560–6860 as CHECK inputs; W(g*) recomputed from data, matching the zero side to ≤ 7%):

| pair | gap | L = 3 | L = 5 | L = 7 | L = 8 |
|---|---|---|---|---|---|
| Lehmer #6709/6710 (t ≈ 7005.06) | 0.0377 | ≤ 9.8 | ≤ 5.9 | **≤ 0.120** | **≤ 0.063** |
| median gap #6695/6696 | 0.884 | ≤ 14.3 | ≤ 19.5 | ≤ 22.8 | ≤ 33.3 |
| first two zeros #1/2 | 6.887 | ≤ 107 | ≤ 494 | ≤ 1356 | ≤ 2024 |

- At L ≥ 7 the window resolves the Lehmer pair (q = |Δĝ|²/W → 0.96 at L = 8, against the ideal maximum of 2). The quotient Rayleigh quotient is then ≤ 0.063, two orders below a median pair and far below the band-limited bulk. So **no uniform-in-L gap of size ≳ 0.07 exists on the full quotient at L = 8**, and any band-limited value is an artifact.
- At L ≤ 5 the window cannot separate the pair (q ≈ 0), so the bound is weak there. The ideal bound Δ²L³/24 → 0 along pairs with Δ → 0. Whether, at *fixed* L, q stays bounded below along such pairs as t → ∞ (the zero density grows like log t) is **OPEN**. So "λ_min = 0 on the quotient at every L" is a HEURISTIC (Ingham/Landau type, Golden Eye), not proved here.
- A uniform gap would also need simple zeros with gap control. That is stronger than RH (Golden Eye) and, by the table, false at L = 8 at the 0.07 level.

## 6. Exact claims
1. **PROVED (elementary):** the dist-bound of §5. For every L, every V ⊆ {e : ê(γ₁) = ê(γ₂) = 0} and every window g: W(g)/dist(g,V)² ≤ W(g)(γ₂−γ₁)²L³/(12|ĝ(γ₁)−ĝ(γ₂)|²).
2. **PROVED (classical; Sylvester + explicit formula):** quotienting by range E_F cannot change the inertia of W_F, for any F with an explicit formula. CHECKED in the pipeline: invariance for EPS, F_θ, DH wherever dim S > 0, with a visible failure for the borrowed (non-radical) range on BEUR″.
3. **CHECKED:** the window-truncated range of E_F is in W_F's radical to ~1e-5 bilinear, separated 1.3e2–1.2e4 from a same-band non-radical control (ζ, EPS, F_θ, DH).
4. **NUMERICAL:**
   - ζ: λ_min(W|S^⊥) has no cut- or basis-independent positive value for L ≥ 1.5.
   - ζ: for L ≤ 1.25, λ_min ≈ exp(−15.9L² + 3.1L − 0.9), basis-converged.
   - The fakes' ghosts deepen like e^{(0.8–1.1)L}.
   - At L = 8 the Lehmer pair forces the quotient Rayleigh quotient below 0.063.
5. **OPEN:** RH; whether the quotient has λ_min = 0 at every fixed L.

## 7. What this means for Daniel's line
Removing the radical removes nothing that matters for positivity (Weil: W(x+e) = W(x)). It also exposes no gap to prove: the quotient's lowest modes are set by how close the zeros are to each other (Lehmer pairs), not by arithmetic that separates ζ from the fakes. The fakes are separated only by their ghosts. That is Weil's criterion again, now shown stable under the quotient.

**Next step (suggested):** the only route Golden Eye flags as untouched is the Connes–Consani **Sonin-compression functional** τ_S(g) = Tr(ϑ(g) S ϑ(g)*), with S complete for the support. Concretely:
- First reproduce the compact archimedean defect W − τ_S for L < log 2.
- Then build the {∞, 2} semilocal analogue on log 2 < L < log 3, the first window where a prime enters with S complete.

That tests positivity as a *different* positive functional plus a controlled defect, instead of a gap.

## Files
`gap.py` (construction + analysis), `clineval.py` (F on the line, complex), `run_gap.py` (pair / zeta scans), `bilinear.py`, `fits.py`, `lehmer.py`, `lehmer_opt.py`; data `cline_*.npy`, `pair.json`, `zscan.json`, `bilinear.json`, `fits.json`, `lehmer.json`, `lehmer_opt.json`, logs.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5). Audit: Turbo. Round-14 audit standard: Golden Eye.
