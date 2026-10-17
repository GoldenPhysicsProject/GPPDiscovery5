# Grok S58 (2026-10-17). RH OPEN. GRH OPEN. φ absent.
See GPP-bridge work/grok/2026-10-17_s58_lift_past_log16_even_knockout.md for results and labels.

- `run140/`: support 2.8 (L = 7/5). `chain_even.sh` (detached) = `lift/lift.sh` (buffered bundle LP warm-started from S55, then two
  B&B implementations per part, merge + group check) -> T from A_cert with beta_even = 0.15 -> even assembly (`run.sh`, s53_onestroke_par)
  -> `pd.sh` (s36_cholres + s37_ldlt_exact + negative control) -> s55_lam inertia. See `run140/STATUS.md`.
  Lift CERTIFIED A <= 5.189 (`run140/torus_L1.40_buf_D9D8_cert.json`). Even assembly RUNNING at session end.
- `ko/`: even-sector knockout margins delta_p^even (EDGEK = 0 free-edge Legendre, M = 40/56/80) and Ritz lambda_min(Q_zeta, even) at
  L = (1/2) log p (`s58_ko.py` = s57_ko.py with paths fixed). NUMERICAL.
- `wall/`: `s58_pmax.py` recomputes the S28 lift floor lambda_max(P_L) (all primes, and restricted to the prime groups {2,3,5}, {7,11,13})
  to measure how far the separable frontier lift sits above the floors. NUMERICAL (float Ritz lower bounds of the floors).
