# λ_min(Q_L) for the Epstein window Weil form, compared with 2Θ − 1 (lead S99 item 1). NUMERICAL. RH is OPEN.

**Dedupe.** The probe-vector η fit is Queen's, not ours. See team/queen 7bdce47: η̂ = 3.20 against 2Θ − 1 ≥ 3.159 on the 2x²+2xy+3y² zero, and our x²+5y² probe numbers reproduced. The bound η_D ≥ 2Θ_D − 1 is her Lemma 5.2 (636005e). We do not redo either here. This file covers only the matrix minimum.

**Setup** (`galerkin.py`).
- The form is Queen's Lévy–Khintchine window form, built from the Dirichlet coefficients of −Z_Q′/Z_Q (`gen_coeffs.py`, which includes the rational-frequency atoms of the non-principal forms). It has the Epstein archimedean term and the pole term at s = 0, 1.
- Basis: orthonormal Legendre polynomials on **[−L, L] (half-width L)**, N functions.
- The kernel is real and symmetric, so λ_min over complex f equals λ_min over real f.
- Validation (`validate_galerkin.py`): the Galerkin Rayleigh quotient of the projected probe vector, for ρ₀ = 2.0794141 + 7.6288676i at L = 2, is −168.30488 from N = 32 on. The exact probe value is −168.3.

**Converged N.** We count an L as converged when the relative change of λ_min between the two largest N is < 1%.

| form (D) | L (half-width) | λ_min at largest N | change vs previous N | Θ used (largest-Re certified zero) |
|---|---|---|---|---|
| 2x²+2xy+3y² (−20) | 0.5 | −0.073082 (N=64) | 2e−5 | 2.0794141 + 7.6288676i |
| | 1.0 | −5.0936 (64) | 5e−5 | same |
| | 1.5 | −42.944 (64) | 2e−4 | same |
| | 2.0 | −190.15 (64) | 1e−4 | same |
| | 2.5 | −1060.4 (96) | 3e−3 | same |
| | 3.0 | −4857.7 (96) | 3e−4 | same |
| x²+5y² (−20) | 1.0 | +1.0532 (64) | 6e−5 | 0.9376669 + 29.983395i (t ≤ 100) |
| | 2.0 | −1.6379 (128) | 1e−3 | same |
| | 3.0 | −8.6439 (128) | 2e−4 | same |
| | 4.0, 5.0 | not converged at N = 64 (−28.3, −67.0, still moving 5–10%) | — | — |
| x²+y² (−4, h = 1, control) | 1.0 | +5.588e−6 (64) | 9e−4 | none off the line |
| | 2–5 | \|λ_min\| ≤ 8e−13 (float64 floor) | — | — |

**Fits.** Slope of log\|λ_min\| against L.

| form | ρ₀ for Θ | L range | fitted slope | 2Θ − 1 |
|---|---|---|---|---|
| 2x²+2xy+3y² | 2.0794141 + 7.6288676i | [1.5, 3] | **3.18** (local slopes 2.98, 3.44, 3.04) | 3.159 |
| 2x²+2xy+3y² | same | [2, 3] | 3.24 | 3.159 |
| x²+5y² | 0.9376669 + 29.983395i | [2, 3] | 1.66 | 0.875 (certified t ≤ 100 only) |

**Reading** (NUMERICAL):
1. **2x²+2xy+3y².** λ_min grows like e^{ηL} with η ≈ 3.18 ≥ 2Θ − 1 = 3.159, using Θ ≥ 2.0794 (CERTIFIED zero).
   - The local slopes oscillate around 2Θ − 1. That fits beating with the next zeros, 1.9939 + 23.02i and 1.9644 + 38.35i.
   - The fit is consistent with η = 2Θ − 1, with Θ the true sup, and with Lemma 5.2.
   - Θ is bracketed: 2.0794 ≤ Θ ≤ 2.4102 (`../epstein/dominant_abscissa.py`, CERTIFIED).
2. **x²+5y².** The slope 1.66 is well above 2·0.9377 − 1. This is no contradiction. Lemma 5.2 is a lower bound, and the true Θ of x²+5y² is > 1, because DH zeros with Re > 1 exist but have t > 100. So 2Θ_true − 1 > 1 is not determined here.
3. **Control x²+y².** λ_min ≥ −1e−12 (float floor) for L ≤ 5. This is consistent with Weil positivity for ζ_{Q(i)}, as expected.
4. **Lead's question** (does Θ ≈ 2.08 show up as η ≥ 2Θ − 1?): yes, numerically, for the matrix minimum at converged N.

Files: `galerkin.py`, `validate_galerkin.py`, `gal_*.log`, `galerkin_*.json`, `hiN/` (N = 80–160 reruns).
Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
