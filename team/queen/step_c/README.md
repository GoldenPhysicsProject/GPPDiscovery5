# Step (c): converged λ₂(R_e), λ₁(R_e), λ₁(R_o) (and λ₂(R_o)) for ζ, L = 0.5–2.0 (half-width), N = 24/32/40/48(/56/64)

Queen of Mathematics, 2026-10-22. Daniel Toupin / Golden Physics Project, with AI assistance from Grok. **RH is OPEN.** Everything here is **HIGH-PRECISION NUMERICAL** (mpmath, not interval).

## Setup
- **Form.** The S87 builder `s87/hp_forms.py` (ZETA), with exact polynomial moments and Gauss–Legendre kernel moments. R = Q − pole = A − P, per sector. Q_e = R_e + ppᵀ and Q_o = R_o − qqᵀ (s89 conventions).
- **Basis.** b_j(s) = √((2j+1)/L) P_j(2s/L − 1), with s = log x ∈ [0, L], extended by parity. Metric L²(ds) on [0, L]. **Half-width L**, primes n ≤ e^{2L}.
- **Grid.** L = 0.5, 0.6, …, 2.0. N = 24, 32, 40, 48 at every L. N = 56, 64 at L ∈ {0.5, …, 1.4, 2.0}.
- **Precision.** dps 120. The check at dps 170 is identical to ≥ 29 significant digits at (L, N) = (1.2, 48), (2.0, 48) and (1.0, 64). So working precision is not the limitation; **basis truncation is.**
- **Monotonicity (PROVED, Cauchy interlacing for nested subspaces).** The Legendre spaces are nested, so every Ritz value is non-increasing in N and is an **upper bound** for the corresponding operator eigenvalue (min–max). The tables confirm monotone decrease throughout.
- **Lower bound (PROVED, rank-one interlacing).** Since R_e = Q_e − ppᵀ, λ₂(R_e) ≥ λ₁(Q_e). Likewise λ₁(R_o) ≥ λ₁(Q_o). So wherever Q_{e/o,L} ≥ 0 (under RH, or at a window with a certificate), 0 ≤ λ₂(R_e) ≤ (Ritz value at any N).
- **Convergence criterion.** Relative change < 2% between the last two N (56 → 64 where available). Values are quoted as converged only under that criterion. The tolerance is stated per row in `tables_c.md`.

## Converged values (NUMERICAL; N = 64, rel. change N56 → 64 in brackets)
| L | λ₁(R_e) | λ₂(R_e) | λ₁(R_o) | λ₂(R_o) |
|---|---|---|---|---|
| 0.5 | −1.98746 (2e−6) | 1.1904e−2 (1.3e−4) | 1.97628e−2 (1e−6) | 0.38235 (1e−4) |
| 0.6 | −2.41691 (5e−6) | 6.620e−5 (4.5e−4) | 3.066e−3 (3e−4) | 6.6323e−2 (1.5e−5) |
| 0.7 | −2.86719 (6e−6) | 5.521e−8 (7e−4) | 6.526e−6 (1.0e−3) | 3.307e−2 (5e−4) |
| 0.8 | −3.31533 (6e−6) | 5.348e−12 (2.1e−3) | 9.826e−10 (2.2e−3) | 4.676e−5 (1.5e−3) |
| 0.9 | −3.82309 (9e−6) | 4.63e−17 (7.5e−3) | 1.534e−14 (6.5e−3) | (see tables) |
| 1.0 | −4.30576 (6e−6) | **1.42e−23 (1.8e−2)** | 8.17e−21 (7.3e−3) | 3.29e−15 (1.5e−2) |
| 1.1 | −4.83600 (1.5e−5) | NOT converged (≤ 1.01e−31; 15%) | NOT converged (≤ 9.04e−29; 6.7%) | – |
| 1.2–2.0 | converged, |Δ| < 4e−5 (rows in tables_c.md) | **NOT converged**: drops 10²–10⁷× per ΔN = 8 | **NOT converged** | NOT converged |

- **λ₁(R_e)** (the S98 nodeless soft mode, negative) is converged to ≤ 4·10⁻⁵ relative at every L in 0.5–2.0. At L = 1.5–1.9, only N ≤ 48 was run; the N40 → 48 change there is ≤ 3·10⁻⁵.
- **λ₂(R_e)** is converged (2%) **for L ≤ 1.0 only**.
  - For L ≥ 1.1, the Ritz values fall geometrically to super-geometrically in N, e.g. at L = 2.0: 1.5e−52 (N48), 9.6e−61 (N56), 1.2e−67 (N64).
  - There they are only **upper bounds**. With the interlacing lower bound, 0 ≤ λ₂ ≤ that value wherever Q_{e,L} ≥ 0.
  - Do **not** use the N = 24 values at L ≥ 1: they are 10⁴–10³⁹ too large.
- **Eigenvector convergence is much better than eigenvalue convergence.** 1 − |cos(v₂^{(N)}, v₂^{(N−8)})| is ≤ 2·10⁻¹⁰ for L ≤ 1.1, 1.6e−5 to 1.0e−4 at L = 1.2–1.4, 1.3e−3 at L = 2.0 (N56 → 64), and 5e−4 to 7e−3 at L = 1.5–1.9 (N40 → 48). This is in `lambda2_Re_eigenvectors.json`.
- **Super-exponential law (CONJECTURE, a fit to the converged points L ≤ 1.0).** log₁₀ λ₂(R_e) has slope d/dμ ≈ −3.8 to −4.9 in μ = e^{2L} (log₁₀ λ₂ = −1.92, −4.18, −7.26, −11.27, −16.33, −22.85 at μ = 2.72 … 7.39). That is consistent with Connes–Consani prolate-type smallness, λ ≈ exp(−c·μ). Extrapolating predicts ≈ 10^{−31.2} at L = 1.1, consistent with the N = 64 upper bound 1.0e−31, which is still moving 15%. It predicts λ₂ ≪ 10⁻⁶⁰ at L ≥ 1.5. Converging there would need N ≫ 64, i.e. a prolate or cosine basis as in Goddess of Prime Number Theory's K = 266–400 runs.

## For H99a / H100a (lead)
`lambda2_Re_eigenvectors.json` contains, per L: N, dps, the λ₂(R_e) Ritz value, a converged flag, |cos| against the N−8 vector, and the 40-digit coefficients in the basis above.
- Coefficient norm is L²(ds) on [0, L]. Divide by √2 for unit norm on [−L, L].
- The sign is arbitrary.
- For L ≥ 1.1, these are Ritz vectors of unconverged eigenvalues. They are stable as vectors (table above), but they are not certified eigenvectors.

Not duplicated, per coordination: DH node-count ladder (Principal series) and the 2.84 replay (Spectrum & Measure).
