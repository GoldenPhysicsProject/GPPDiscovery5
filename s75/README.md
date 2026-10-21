# Grok S75 (2026-10-21): DH onset refinement, ζ-vs-DH separating test function, Mori-route cost at the DH window

RH OPEN. GRH OPEN. φ absent (DH's κ comes from q = 5 only).
Note: GPP-bridge `work/grok/2026-10-21_s75_dh_onset_and_mori_cost.md`.

- `dhritz.py` is an unchanged copy of `../s63/dhritz.py` (= s62 cosritz with FORM=DH). Arb Rayleigh quotient of a rational
  trig-basis Ritz vector, so a negative RQ is a rigorous negative witness of the DH Weil form (upper bound on inf Q/||f||^2).
- `sepritz.py` = the same plus env `VEC=<cosvec json>` to evaluate a FIXED rational vector under FORM (ZETA or DH).
- DH even onset scan: `dh_L*_K*.jsonl` (+ `cosvec_DH_even_*.json` witnesses).
  - CERTIFIED negative: L = 3421/2000 (1.7105) K300 [-2.06397269789e-31 ± 3.5e-43], K360 [-2.77559113650e-31 ± 2.7e-43];
    L = 1711/1000 K240 -1.3648e-30; 3423/2000 K240 -2.5358e-30; 214/125 K240 -3.5284e-30.
  - NUMERICAL (positive Ritz upper bounds, decreasing in K): L = 6841/4000 (1.71025) K300 +5.474e-31; 1.7105 K240 +7.09e-33;
    L = 171/100: K129 5.455, K180 2.550, K240 1.599, K300 1.358, K360 1.277 (e-30) -> extrapolated limit ≈ +1.24e-30.
  - So L*(DH) ∈ (1.70, 1.7105] (upper CERTIFIED, lower NUMERICAL), linear interpolation L* ≈ 1.7104 (support ≈ 3.4208 = log 30.6).
- Separation (`zeta_on_DHvec_L137_80_K129.jsonl`, control `dh_on_DHvec_L137_80_K129.jsonl`): the S63 DH witness f* at L = 137/80
  gives Q_ζ(f*)/||f*||² = [1.50064301068282e-8 ± 4.3e-23] (CERTIFIED > 0) while Q_DH(f*)/||f*||² = [-2.06863376668450e-30 ± 3.0e-45]
  (CERTIFIED < 0, reproduces S63). ζ's own even Ritz at 137/80, K129: RQ = [9.40325369652e-63 ± 4.5e-75] (CERTIFIED upper bound on
  ζ's even λ_min there).
- `mori_cost_dh.py` / `.json`: NUMERICAL resource extrapolation for the Mori route (Astra bound × cure factor 2.12–2.45).
  At L = 137/80: Ω ≈ 6250–7230, N ≈ 17,900–20,700, ≥ 23–26 GB per worker (linear-in-N extrapolation), dump ≥ 13–18 GB,
  exact LDLᵀ ≥ 39–61 h at 640 bits, assembly ≥ 210–1790 worker-hours. Not launched.
Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
