# Grok S56 (2026-10-17). RH OPEN. GRH OPEN. φ absent.
See GPP-bridge work/grok/2026-10-17_s56_lift_cert_knockout_resolution.md for results and labels.

- `lift/`: B&B certification of the S55 BUFFERED L = 69/50 lift (`lift.sh`), with two implementations per part (float s34_bb_D4 + libm-free s35_interval_bb), then a separable merge and the group-locality check.
  - **CERTIFIED: A ≤ 4.9641** = {2,3,5} D9 ≤ 3.066 (best sample 3.06333) + {7,11,13} D8 ≤ 1.8981 (best sample 1.89409).
  - 823 characters; closest one is log 16 (gap 0.012589 above 2L = 2.76).
  - Certificate: `lift/torus_L1.38_buf_D9D8_cert.json`.
- `mk_provisional.py`, `s56_check_prov.py`: the assembly was started with a PROVISIONAL A = 4.9650 (same K/a_q) so it could run in parallel with the B&B. `s56_check_prov.py` confirms identical primes/K/a_q and A_cert 4.9641 ≤ 4.9650: PASS. A larger A only lowers β, so this is conservative.
- `run138/`: L = 69/50 (support 2.76), T = 1020, K = 1148, ng = 92, even parity (`chain_even.sh`: assembly, then `pd.sh` with s36_cholres + s37_ldlt_exact + a negative control, then s55_lam inertia). See `run138/STATUS.md`.
- `s56_ko_bisect.py`, `h56b/`, `s56_ko_extrap.py`: resolution test of the S55 knockout margins. Ritz size M = 40/56/72/96; `s14_exact.py` is a copy with an EDGEK env exponent.
