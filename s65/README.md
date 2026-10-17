# Grok S65 (2026-10-17): Perron-weight Schur bound for the prime operator; DH Schur bound; H_{-t} small-|t| edge

RH OPEN. Finite windows only. Note: GPP-bridge/work/grok/2026-10-17_s65_perron_schur_dh_backward_heat_edge.md

- `schur/perron_schur.py`: a weighted Schur test for C_p on L²[−L, L] with an optimised step weight (shifted nonlinear power iteration). It gives a rigorous upper bound: window maxima are widened by one cell, and the coefficients are rounded up from Arb. Modes: `zeta`, `dh` (|DH coefficients|, computed in Arb), `dhpos`, `dhneg`.
- `schur/arb_recheck.py`: recomputes the final sup in Arb, for rational L only. Results are in `arb_recheck.jsonl` (L = 7/5 and 71/50); the weights are in `w_*.npy`.
- `schur/cost_table.py`: a cost model, our one-stroke lift vs Mori/Liu bounded comparison with s = Schur bound. ESTIMATE only (calibrated on Mori VIII: it predicts Ω 716 / 810 modes vs actual 736 / 832).
- `hneg/`: `tc_track.py` tracks the γ₃₄/γ₃₅ real pair of H_t as t decreases; run_edge*.sh do the scans, the complex Newton and the K = 64 window tests at t = −0.3875 … −0.45. Summaries: `tc_summary_S65.json`, `window_summary_S65.json`.
