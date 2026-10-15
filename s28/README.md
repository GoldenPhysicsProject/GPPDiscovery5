# Grok S28 (2026-10-15) — RH OPEN, GRH OPEN
- `s28_onestroke_fast.py` — fast rigorous lifted one-stroke (Zhu 2608.24827 + S27 lift): parallel Arb assembly, blocked Schur-complement ball PD check.
- `s28_tails_arb.py` — Arb re-derivation of the final error budget (Codex 6204288080: exact ratio tails, rational rho = 347/100). Output `tails_arb_audit.json`.
- `cert_L119_T700_ng60.json` / `fast_L119_T700_ng60.log` — **CERTIFIED (Arb, Lean-open)**: zeta Weil positivity on [-1.19, 1.19] (support 2.38): even >= 6.8819e-48, odd >= 3.9149e-44.
- `cert_L119_T700_pp.json` / `fast_L119_T700_pp.log` — same at 40 nodes/panel: PD but NOT a certificate (quadrature radius >> lambda). Logged as a dead end.
- `cert_L10_T300_fast.json`, `cert_L08_T80_lift_fast.json` — reproductions of the S27 L=1 and L=0.8 certificates with the fast pipeline.
- `s28_pmax.py`, `pmax.json`, `pmax_conv.log` — lambda_max of the prime-shift operator P_L (floor of every invisible lift; NUMERICAL).
- `s28_torus_lp.py`, `s28_torus_bb.py`, `torus_*.json`, `tor_*.log` — cross-prime torus LP (NUMERICAL), float-B&B bound for the {3,5} toy. tor_L119_D3.log is the converged L=1.19 D=3 history (3.6019/3.6023); the D3b run was stopped.
- `s28_pairs_L2.py` — L=2 pair-group torus LP (not run to completion).
- `lift_L1.25.json`, `lift_L1.3.json` — per-prime lift designs for the next windows (A = 4.583, 5.313).
