# S23 — ζ odd sector at h = 9/20: optimized extension of sinh(x/2)

Daniel Toupin / Golden Physics Project, with AI assistance from Grok. **RH OPEN.**

Continues S22. At h = 9/20 Claim 1 (m ≥ c = 38/4235) already passed, but the un-extended
v = sinh(x/2)·1_(−h,h) gave 2β ≈ 1.50. An optimized odd extension brings 2β below 1.

## Result (CERTIFIED, Arb 128-bit, Lean-open)

For every odd g ∈ C_c((−9/20, 9/20)), real or complex,

**Q_ζ(g) ≥ 0.000737 ‖g‖².**

- Claim 1: m ≥ 38/4235 (same LP symbol as S22; 60k cells, 0 splits).
- Claim 2: with a piecewise-linear odd extension of sinh(x/2) on [h, H] = [9/20, 6]
  (N = 120 free knots, y_{N+1} = 0), one has 2β ≤ 0.91785 < 1.
- Coercivity constant (1−2β)c ≥ 0.000737.
- Float design: 2β_core ≈ 0.910; Galerkin χ(M=24) ≈ 0.815 (lower bound on true χ).

## Scripts

- `s23_extend_float.py` — sine+linear exterior LS (early sweeps; insufficient alone).
- `s23_extend_legendre.py` — Legendre / knot exterior sweeps.
- `s23_verify_arb.py` — Arb Claim 1+2 verifier (three-point S bound + M₂r²/8 remainder).
- `s23_cert_h045_N120.json` — rational m-cert + quantized extension coefficients (den = 10¹²).
- `s23_verify_results.json` — Arb output.

## Method

Reuse S22's odd-only DC lift and continuum pole Schur lemma: for any odd extension v_w of
v = sinh(x/2), β = (1/2π)∫|v̂_w|²/m and Q ≥ (1−2β)c‖g‖² when 2β < 1.
The extension is free on [h, H] (piecewise linear through free knot values) and zero outside;
coefficients minimize float β by unconstrained least squares, then are quantized and Arb-certified.

## Dead ends / not claimed

- Compact sine-only / low-knot exteriors with N ≲ 40 leave 2β > 1 (logged in `out_sweep*.json`).
- Pure exponential or sinh-cutoff exteriors do not reach 2β < 1 with the current m.
- Even sector not covered. φ does not appear.
- Optional stretch h = 1/2 not certified here (Galerkin χ ≈ 0.979; needs a near-optimal extension).

Lean asks → Codex (DC lift, Schur lemma continuum form, cell/tail certificate with PL extension).
