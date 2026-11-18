# Round 14E2: the 14E quotient question at 80 digits (Golden Eye's falsifier, AUDIT14 §6.8)

Framing credit: **Daniel Toupin** (No-Ghosts, 14B; one global object, 14D). The falsifier and audit standard are **Golden Eye**'s (AUDIT14 §6, team/golden-eye eae661d). Her scratch `mpw.py` is ported here as `mpW.py`.
Labels: PROVED / CHECKED (code or identity check) / NUMERICAL (evidence) / CONJECTURE / GAP / OPEN. Proof inputs are integer data only: Dirichlet coefficients (through F(½+it) from mpmath's ζ), Γ-factors, poles. **No zeros enter 14E2.** (The §5 Lehmer falsifier in 14E does use zeros as proof inputs; see the 14E errata.)

**RH: OPEN. Nothing here moves it.**

## 0. Verdict first
- **The quotient tracks λ_min(W) down the ladder. There is no finite-L gap (NUMERICAL, 80 digits).** At L = 1.5, 1.75 and 2.0, for 2 bases (fmax 45, 90), 4 concentration cuts (ε = 1e-10 to 1e-4) and 2 family-span cuts (τ = 1e-13, 1e-40), the true quotient seminorm on S^⊥ has λ_min between **1.09× and 3.6× λ₁(W_ζ)**. At L = 2, fmax 90 that is 1.2e-29 to 3.6e-29, against λ₁(W) = 1.0e-29.
- The 14E pipeline number, the **restriction** λ_min(QᵀWQ), is not a gap but a **cut**. It moves over 10 orders with ε at fixed L (L = 2, fmax 90: 9.2e-19 at ε = 1e-10, up to 1.1e-9 at ε = 1e-4). It always sits below the interlacing cap λ_{dim S+1}(W). It rises only because S is *not* in the radical: S's vectors carry W-energy 1e-8 to 1e-5, and the L²-complement pays for that leak.
- **Why it comes out this way (diagnostic).** The window-truncated range S is nearly orthogonal to W_ζ's near-null eigenvectors: the sine of the largest principal angle between S and the bottom-dim S eigenvectors is 0.92–1.00 in every run. Removing S therefore removes none of the ladder 1e-29, 1e-26, 1e-23, … (L = 2).
- **Control (F_θ, θ = 0.2, L = 2).** The ghosts persist: 7 → 7 (or 7 → 6 at dim S = 3, above the interlacing floor 4). Restricted λ_min stays at −0.0908 to 5 digits, and the inf-over-coset form gets *more* negative (−0.7 to −2.6) because S is not in F_θ's radical either.
- **Reproduced (CHECKED):** all of Golden Eye's 80-digit λ_min(W_ζ) values (8 L values, same digits), her fit exp(a + bL − Be^L) with B = 13.13, and the Lehmer bound 0.1203 (L = 7) / 0.0629 (L = 8), bit-for-bit.

## 1. Construction (what was rebuilt, and the one deliberate change)
- **W at 80 digits** (`mpW.py`): line-by-line mpmath port of `weil_galerkin.operators` (Golden Eye's `mpw.py`, generalised to F_θ). Dirichlet-sine basis, M = ⌈fmax·L/π⌉ + 8. The arch integral uses Gauss–Legendre, 4 panels × 192 nodes (Golden Eye: deg 7 vs deg 8 agree to 8 digits at L = 2, M = 66). Prime and pole terms are closed form.
  - CHECKED: the λ_min values equal hers to all 8 printed digits (e.g. 6.394341e-15 at L = 1.5, M = 30; 1.2843737e-27 at L = 2, M = 37; 1.0135947e-29 at M = 66; 7.824487e-30 at M = 94).
  - F_θ's λ_min −0.0907676 (M = 37) matches 14E and 14F (−0.091).
- **F on the line at 90 digits** (`mpline.py`): ζ(½+it) and ζ(s−0.2)ζ(s+0.2), t = 0 to 260 in steps of 0.1.
- **S** (`mpgap.py`) is built exactly as in `gap.py::radical_space`, all in 90-digit arithmetic (python-flint arb for the t-integrals; mpmath for every eigenproblem):
  - same family layout (centres, half-widths w, frequencies k ≤ 1.1·fmax, phases; widths (1.5, 1.0); `run_gap.py zeta` settings, also used for F_θ);
  - same pole constraints Ĥ(½ + a) = 0;
  - same Parseval route a_i = ⟨φ_i, E_F H⟩ = (1/π) Re ∫₀^∞ φ̂_i · conj(F Ĥ) dt;
  - same Gram whitening, generalised concentration eigenproblem, cut ρ ≥ 1 − ε, and SVD cut 1e-8.
- **Deliberate change: Gaussian profile instead of the bump.** h(v) = exp(−(v−c)²/(2σ²))·cos(k(v−c)+φ), with σ = w/4 (the bump's s.d. is 0.40w).
  - Reason: the bump's transform decays only like exp(−c√t), so any t-truncated Parseval integral stalls near 1e-12. That is exactly why 14E's S was a double-precision object.
  - The Gaussian's transform is closed form, e^{sc}·√(2π)σ·½[e^{iφ}e^{σ²(s+ik)²/2} + c.c.] with s = ½+it. It decays like e^{−σ²t²/2}.
  - Double-precision cross-check of the change (`gap.py` with the Gaussian swapped in): dim S at L = 2 is 3–7 (bump: 2–5), and the restricted λ_min is again at the double floor. Same qualitative object.
- **Precision checks (CHECKED):**
  - `hcheck.py`: trapezoid h = 0.1 vs h = 0.2 on the pole-constrained combinations differ by ≤ 2.5e-91 (ζ) and ≤ 2.0e-90 (F_θ), against O(1) entries.
  - *Raw* family members differ by 1e-7 / 1e-4. Reason: F·Ĥ_j keeps F's pole at t = −i/2 unless Ĥ_j(1) = 0, which gives an e^{−π/h} trapezoid error. The pole constraint B removes it exactly, and every quantity is computed after B.
  - arb ball radii on the Grams are ≤ 4e-100.
- **Convention check (CHECKED, same-band non-radical control as 14E `bilinear.py`).** S was rebuilt from F′ = ζ − 0.1·2^{−s} (L = 2, fmax 45, τ = 1e-13). Its bilinear residual against W_ζ is 2.0e-2 to 3.0e-2, against 4.1e-7 / 1.2e-5 / 6.0e-5 / 1.2e-3 for ζ's own S at ε = 1e-10 / 1e-8 / 1e-6 / 1e-4. That is a separation of 4.9e4 / 1.9e3 / 3.6e2 / 24. The true-quotient λ_min for the control is again ≈ λ₁(W) (1.55e-27 to 1.82e-27), as §2 predicts for *any* S that misses u₁.

## 2. Two quotients, and what is forced (PROVED, elementary)
Let S ⊂ V_M with dim S = k, Q an orthonormal basis of S^⊥, and λ₁ ≤ λ₂ ≤ … the eigenvalues of W on V_M, with eigenvectors u_j.
- **Restriction (the 14E pipeline):** R = QᵀWQ. Cauchy interlacing gives λ_j(W) ≤ λ_j(R) ≤ λ_{j+k}(W).
- **True quotient seminorm:** q([x]) = inf_{s∈S} W(x+s), normalised by ‖x‖² for x ∈ S^⊥ (the L² distance of the coset to 0). For W > 0 this is the Schur complement W/W_SS on S^⊥. If S ⊆ radical, then q = R = the exact quotient W/S, which is 14E's Sylvester setting.
- **Forced brackets for W > 0:**
  - λ₁ ≤ λ_min(q), since W(x+s) ≥ λ₁‖x+s‖² ≥ λ₁‖x‖².
  - λ_min(q) ≤ λ₁ / (1 − ‖P_S u₁‖²): take x = u₁ − P_S u₁, whose coset contains u₁.
  - λ_min(q) ≤ λ_{k+1}(W) (interlacing applied to W⁻¹).
- **Consequence.** A finite-L gap on the quotient, meaning λ_min(q) ≫ λ₁, *requires* S to contain u₁ (and, for a gap above λ_j, u₁…u_{j−1}) to high accuracy. The falsifier therefore reduces to one geometric question: **does the window-truncated range of E contain W's bottom eigenvectors?**
- **Finite-L remark (NUMERICAL).** W_ζ is positive definite on V_M at every L tested (80 digits), so its radical *in the window Galerkin space is {0}*. No nonzero subspace of V_M can be quotiented exactly. At fixed L the honest "quotient gap" in V_M is λ₁(W_L) itself, which falls like exp(−13 e^L) (§4).

## 3. The λ_min table (ζ, 80 digits; `round14e2.json`, `summary14e2.md`)
Cell = dim S: **R** = restriction λ_min(QᵀWQ) / **Q** = true-quotient λ_min (Schur). "dim 0" means no family direction reaches the cut, so both equal λ₁(W).

τ = 1e-13 (14E's whitening cut):

| L | fmax (M) | λ₁(W) | ε = 1e-10 | ε = 1e-8 | ε = 1e-6 | ε = 1e-4 |
|---|---|---|---|---|---|---|
| 1.5 | 45 (30) | 6.39e-15 | dim 0 | dim 0 | 1: R 4.40e-12 / Q 7.45e-15 | 3: R 1.55e-09 / Q 7.74e-15 |
| 1.5 | 90 (51) | 4.67e-15 | dim 0 | dim 0 | 1: R 3.49e-12 / Q 5.44e-15 | 3: R 1.58e-09 / Q 5.64e-15 |
| 1.75 | 45 (34) | 5.22e-21 | 1: R 3.48e-18 / Q 5.91e-21 | 2: R 6.23e-16 / Q 6.04e-21 | 3: R 2.35e-13 / Q 6.08e-21 | 5: R 7.82e-11 / Q 7.57e-21 |
| 1.75 | 90 (59) | 2.24e-21 | 1: R 1.97e-18 / Q 2.52e-21 | 3: R 1.78e-13 / Q 2.61e-21 | 4: R 2.16e-11 / Q 2.71e-21 | 6: R 8.79e-08 / Q 3.73e-21 |
| 2.0 | 45 (37) | 1.28e-27 | 3: R 2.03e-19 / Q 1.53e-27 | 5: R 4.31e-15 / Q 1.57e-27 | 6: R 8.08e-14 / Q 1.67e-27 | 8: R 1.71e-09 / Q 2.86e-27 |
| 2.0 | 90 (66) | 1.01e-29 | 4: R 9.16e-19 / Q 1.21e-29 | 6: R 2.45e-14 / Q 1.31e-29 | 7: R 2.40e-12 / Q 1.35e-29 | 9: R 1.10e-09 / Q 2.33e-29 |

τ = 1e-40 (richer family span; robustness):

| L | fmax (M) | λ₁(W) | ε = 1e-10 | ε = 1e-8 | ε = 1e-6 | ε = 1e-4 |
|---|---|---|---|---|---|---|
| 1.5 | 45 (30) | 6.39e-15 | dim 0 | dim 0 | 2: R 3.09e-10 / Q 7.74e-15 | 3: R 3.48e-10 / Q 8.11e-15 |
| 1.5 | 90 (51) | 4.67e-15 | dim 0 | 1: R 1.19e-12 / Q 5.10e-15 | 2: R 2.53e-10 / Q 5.69e-15 | 3: R 4.15e-10 / Q 6.00e-15 |
| 1.75 | 45 (34) | 5.22e-21 | 2: R 6.39e-16 / Q 6.03e-21 | 3: R 7.59e-14 / Q 6.03e-21 | 4: R 3.31e-11 / Q 8.13e-21 | 6: R 2.84e-07 / Q 1.07e-20 |
| 1.75 | 90 (59) | 2.24e-21 | 2: R 4.26e-16 / Q 2.58e-21 | 3: R 1.92e-14 / Q 2.61e-21 | 5: R 8.98e-11 / Q 3.85e-21 | 6: R 2.59e-07 / Q 4.99e-21 |
| 2.0 | 45 (37) | 1.28e-27 | 4: R 1.06e-17 / Q 1.57e-27 | 6: R 9.98e-14 / Q 2.32e-27 | 7: R 3.73e-12 / Q 3.04e-27 | 8: R 2.13e-09 / Q 3.45e-27 |
| 2.0 | 90 (66) | 1.01e-29 | 5: R 3.93e-16 / Q 1.26e-29 | 7: R 1.34e-12 / Q 2.44e-29 | 8: R 1.09e-09 / Q 2.81e-29 | 9: R 1.35e-07 / Q 3.62e-29 |

- **Across all 41 non-empty ζ runs** (Q and angle ranges also include the 4 F′-control runs): Q/λ₁ ∈ [1.09, 3.57]; R/λ₁ ∈ [2.6e2, 1.3e22]; R/λ_{k+1}(W) ∈ [2e-4, 0.87] (cap always respected); sin(largest angle between S and u₁…u_k) ∈ [0.917, 1.000]; bilinear ‖WS‖/‖W‖ ∈ [2.3e-7, 2.0e-3].
- **Basis stability.** Q follows λ₁ of whichever basis is used: the fmax-45 Ritz value is 127× too high at L = 2, and Q inherits that. Q/λ₁ itself stays in the same 1.1–3.6 band for both bases. R is not stable in either cut or basis.
- **Cut stability.** Q moves by ≤ 3.6× over ε = 1e-10 to 1e-4 and both τ. R moves by up to 1e10.

## 4. Reproduction of Golden Eye's 80-digit data and the replacement fit (`fits80.py` → `fits80.json`) — NUMERICAL
| L | 0.25 | 0.5 | 0.75 | 1.0 | 1.25 | 1.5 | 1.75 | 2.0 |
|---|---|---|---|---|---|---|---|---|
| λ_min(W_ζ), best M | 0.3290 | 3.602e-2 | 5.552e-4 | 1.064e-6 | 3.044e-10 | 4.121e-15 | 2.065e-21 | 7.824e-30 |
| M | 16 | 23 | 30 | 37 | 44 | 73 | 84 | 94 |

- Identical to AUDIT14 §6.5 to the printed digits.
- Fits on these 8 points (max log-residual):
  - a + bL: 10.6, rejected;
  - quadratic: 1.12;
  - a − Be^L: 2.1;
  - **a + bL − Be^L: 0.32, with a = 14.09, b = 7.90, B = 13.13.**
- Local slope d log λ/d(e^L): −10.4, −10.6, −11.3, −11.4, −11.9 over L = 0.75 to 2.0.
- 14E's Gaussian, refit on its 5 points and extrapolated, predicts −32.1 / −44.3 / −58.5 at L = 1.5 / 1.75 / 2.0. Measured: −33.1 / −47.6 / −67.0.
- B → 4π ≈ 12.57 (Slepian prolate concentration at scale e^{L/2}, Connes–Consani ζ-cycles) is **heuristic**, not derived.

## 5. Exact claims
1. **PROVED (elementary).** For W > 0 on V_M and any subspace S: λ₁(W) ≤ λ_min(W/S-quotient seminorm on S^⊥) ≤ min(λ₁/(1 − ‖P_S u₁‖²), λ_{dim S+1}(W)). The restriction obeys λ₁ ≤ λ_min(QᵀWQ) ≤ λ_{dim S+1}(W). Hypotheses: finite-dimensional Galerkin space and W positive definite there.
2. **CHECKED.** 80-digit W_ζ reproduces Golden Eye's λ_min table. The trapezoid for the pole-constrained Grams is exact to 1e-90. The same-band non-radical control separates by 24 to 4.9e4. Lehmer reproduces bit-for-bit.
3. **NUMERICAL (80 digits; L = 1.5–2; 2 bases; 4 ε; 2 τ).** For ζ the true-quotient λ_min on S^⊥ is within a factor 3.6 of λ₁(W), which is 1e-29 at L = 2. The window-truncated range of E (Gaussian family) is nearly orthogonal to the near-null eigenvectors (sin angle ≥ 0.917). **No finite-L gap at these L.**
4. **NUMERICAL.** F_θ control: ghosts persist under the restriction (7 → 7, or 7 → 6 at dim S = 3). λ_min(R) = −0.0908 is unchanged to 5 digits.
5. **NUMERICAL.** λ_min(W_ζ) ≈ exp(14.1 + 7.9L − 13.1e^L) on 0.25 ≤ L ≤ 2. The B = 4π link is heuristic.
6. **OPEN.** RH. Whether some *other* subspace of the window space, e.g. prolate/Sonin-built rather than truncated E, captures W's ladder: that is the Connes–Consani route, not tested here. L > 2 at 80 digits (cost grows like M³; L = 3 needs M ≈ 94–135). Probes outside the Galerkin band.

## 6. What this means for Daniel's line
The "gap on the quotient" idea is now cleanly dead at finite L.
- Window-truncated E does not touch the near-null ladder of W_ζ. Since W_ζ is strictly positive on the window space, there is no radical there to quotient by. The quotient's floor is λ₁(W_L) ≈ e^{−13e^L}.
- Physically (INTERPRETATION): the ghost-free (unitary) sector of ζ is positive but *not gapped*. Its softest window modes, plausibly prolate-type edge modes (14E §4 saw them peak at the window and band edges), have energy that is double-exponentially small in L. The fakes differ by having negative modes (ghosts), not by having a gap.
- So positivity has to be proved as positivity (Weil), not as a gap. The next step stays Golden Eye's §6.7 pick: the Connes–Consani **Sonin-compression functional** with S complete for the support. In that framework these prolate edge modes are exactly what the compression is built to control.

## Files
- `mpW.py`: 80-digit W, ζ and F_θ; logs in `mp80/W.log`, `mp80/W_fit.log`. The matrix pickles `mp80/W_*.pkl` are regenerable and not committed.
- `mpline.py`: 90-digit F on the line (`mp80/line_*.pkl`, regenerable, not committed).
- `mpgap.py`: S, restriction, quotient, angles. Per-run `mp80/gap_*.json` and `.log`; `CONTROL=1` gives the F′ control.
- `run_mp.sh`: driver.
- `hcheck.py`: grid check, `mp80/hcheck.log`.
- `summarize14e2.py` → `round14e2.json`, `summary14e2.md` (full table, incl. F_θ and control).
- `fits80.py` → `fits80.json`.
- `lehmer.py`: `optimal` fixed, docstring relabelled; `lehmer_opt.log` is the rerun, and `lehmer_opt.json` is unchanged bit-for-bit.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5). Audit: Golden Eye.
