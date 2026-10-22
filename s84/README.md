# Grok S84 (2026-10-22): near-null vectors of the ζ window Weil form, zero-side channels, and DH's failure direction
RH is OPEN. Note: GPP-bridge `work/grok/2026-10-22_s84_nearnull_zero_channels.md`. All shape and spectral numbers are NUMERICAL (mpmath). The RQ values come from Arb (`dhritz.py`, an unchanged copy of s75) and are CERTIFIED upper bounds.
- `selftest.py`: quadrature check of the closed forms for f̂.
- `spec_decomp.py`: ζ zero-side split of Q(f), with the explicit formula checked against the Arb RQ. Outputs `spec_*.json`.
- `band_mass.py`: L² mass of f̂ per frequency band, mode profile, and nodal count at high precision. Outputs `band_*.json`.
- `profile.py`: f̂ profile. Outputs `profile_results.json`. Its γ₁ entry uses a 20-digit γ₁ and is limited by that precision; use `spec_*.json` instead.
- `overlap.py`: coefficient overlaps and cosh(u/2) overlap. Outputs `overlap_results.txt`.
- `dh_channels.py`: ζ zero side of the DH witnesses. Outputs `dh_channels_results.json`.
- `dh_zeros.py` and `dh_zeros_range.py`: DH zeros, both on-line and off-line. `dh_spec.py` and `dh_spec_ext.py`: DH zero-side split.
- Ritz runs:
  - `hp/` (Arb prec 480, dps 130): ζ at L 71/50 with K 96 and 128, and at 137/80 with K 96 and 129.
  - `hp2/` (Arb prec 720, dps 210): ζ at 137/80 with K 129 and 160.
  - `dh/`: DH at 71/50 with K 96.
  - `*_dps60.json`: the first-pass dps-60 vectors. They are noise-floor limited, so do not use them.
Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
