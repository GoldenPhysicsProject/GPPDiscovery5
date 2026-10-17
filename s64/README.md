# Grok S64 (2026-10-17): support-2.8 replay package; lift vs Mori/Liu Schur weight; DH design; H_{-t} window test

RH OPEN. φ absent. Note: GPP-bridge `work/grok/2026-10-17_s64_replay_lift_schur_backward_heat.md`.

- `replay/` (H64a): a standalone Level-1 PD replay (`replay_pd.sh`) plus the certificate JSONs, lift JSONs, PD checkers and `MANIFEST.sha256`. Levels 2 and 3 are in `replay/README.md`.
  - The self-test passed for both parities (`replay_out_*/REPLAY_RESULT.json`). It uses the same code as the original runs, so it is not independent.
- `lift_cmp/` (H64a ii):
  - `cp_norm.py` computes λ_max(C_p|L²[−L, L]) by Galerkin, and the Mori-VIII weighted Schur bound with a quadratic weight optimised. It reproduces Mori's R_w = 4.0915.
  - `dh_symbol.py` computes the size of the DH symbol at L = 1.70 for the H64b cost estimate.
- `hneg/` (H64c): the de Bruijn–Newman backward heat flow H_t, t < 0.
  - `ht_scan.py` (real zeros and count) and `ht_complex.py` / `ht_complex_grid.py` (non-real zeros, matched to the count deficit).
  - `ht_window2.py`: an Arb zero-side window Ritz test. It is validated at t = 0 against the geometric side: L = 1, K = 24 gives 2.3790252e−29 vs cosritz 2.37902545e−29.
  - Results are in `window2.jsonl` / `window_summary_S64.json`.

Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
