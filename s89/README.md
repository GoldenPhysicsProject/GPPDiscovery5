# Grok S89 (2026-10-22): explicit in-band / out-of-band split of the Weil window form; H₋ₜ band test; S88 odd-sector sign correction. RH OPEN.

- `calib.py`: zero-side vs geometric-side calibration of the S87 builder (float64, 100k zeros, spherical-Bessel Legendre transforms).
  CHECKED: Q_e = A − P + Σ ppᵀ and **Q_o = A − P − Σ qqᵀ** (half zero sum Σ_{γ>0} u uᵀ matches to ~1e−3 = zero-tail size). S88's odd sector used +qqᵀ (bug).
- `band_split.py`: Q in the eigenbasis of the arch-only form A; G = arch-negative span, + = complement.
  Q ⪰ 0 ⇔ Q₊₊ ≻ 0 (out of band) and S_G = Q_GG − Q_G₊ Q₊₊⁻¹ Q₊G ⪰ 0 (in band). Also ρ_out = λ_max(A₊₊^{-1/2} P₊₊ A₊₊^{-1/2}).
  Forms ZETA, DH, FT(θ=0.1) at N = 24 (`band_main.json/log`); H_t (t = −1, −2, −4) = Q_ζ + Δ_t at N = 32 (`band_ht_32.json`, `band_ht32.log`),
  Δ_t swaps the first M = 50 ζ zeros for the first 50 zeros of H_t (zero lists from s64/hneg, real + non-real).
  Validated: H_{−2}, L 0.8 gives −2.80e−4 (even) / −1.81e−5 (odd) vs S64's independent cosine-basis Arb engine −2.85e−4 / −1.67e−5.
- `robust.py` / `robust.json`: onset points with M = 40 and N = 24 (signs, band, nodes unchanged).
- `table.md`: summary table. All HIGH-PRECISION NUMERICAL (mpmath dps 90), not interval.

Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
