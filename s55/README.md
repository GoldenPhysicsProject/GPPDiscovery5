# Grok S55 (2026-10-17). RH OPEN. GRH OPEN. φ absent.
See GPP-bridge work/grok/2026-10-17_s55_buffered_lift_frontier.md for results and labels.

- `s55_bundle_lp.py`: S36 bundle LP with a buffer. Characters need |k·log p| ≥ 2L + EPS (env EPS). NUMERICAL design.
- `s55_lp_lower.py`: unrestricted LP on a finite torus point set. This gives a lower bound on the best lift constant for a given character set. NUMERICAL.
- `lift/`: buffered L = 1.38 bundles. {2,3,5} D9: 3.0595. {7,11,13} D8: 1.8941. A ≈ 4.954, NOT B&B-certified. `lower_p235_D8.jsonl`: LP lower bound 2.883.
- `run13505/`: L = 2701/2000, T = 800, K = 881, ng = 92, with the certified S37 lift. Even and odd assembly (`../s53/s53_onestroke_par.py`), then PD (`pd.sh`: s36_cholres + s37_ldlt_exact + a negative control). The 97 MB Atil dumps stay on the box (sha256 in `atil_sha256.txt`).
- `s55_lam.py`: FIX for lam_estimate. It returns the smallest-modulus μ plus the Sylvester inertia (negative-pivot count) and flags INDEFINITE. `lamfix/` has the validation runs.
- `s55_knockout.py`, `h55b/`, `s55_ko_fit.py`: knockout margins δ_p for p ≤ 19, and fits (`h55b/ko_fit.json`).
- `s14_exact.py`: copy (unchanged).
