# Grok S36 — ζ Weil positivity on support 2.6 (L = 13/10), both parities; Rump-style PD step; {2,3,5} D=8 lift

RH OPEN. GRH OPEN. φ absent.

## H36a: CERTIFIED (Arb, Lean-open, single implementation)
**ζ Weil positivity on support 2.6 (L = 13/10) for all complex f ∈ L²[−1.3, 1.3]: Q_ζ(f) ≥ 3.6340e−62 ‖f‖².**
- Even ≥ 3.633998e−62 (binding); odd ≥ 1.4910e−58.
- Parameters: T=700, K=740, w=1/2, b=2/5, **ng=78**, assembly prec 384.
- Lift `../s34/torus_L1.3_D4_D20_cert.json`. A_cert ≤ 4.6565 = {2,3,5,7} D4 ≤3.4548 + {11,13} D20 ≤1.2017 (both S35 interval-confirmed). `s35_check_merge.py` OK, 546 chars group-local.
- β = 0.0552747. M_S ≤ 8.680e6, M0 ≤ 34.192. Wall 1086 s on 5 workers.

| parity | λ_est | PD certificate λ_min(Ã) ≥ | ‖E‖ (residual) | rE | eps_B (Arb audit) | lamD | certified_lower |
|---|---|---|---|---|---|---|---|
| even | 3.6376e−62 | 3.6340e−62 (l0 = 0.999 λ_est) | 3.8e−115 | 1.742e−68 | ≤1.14e−110 | 0.05527 | **3.633998e−62** |
| odd  | 1.4921e−58 | 1.4910e−58 | 1.3e−114 | 1.743e−68 | ≤3.51e−111 | 0.05527 | **1.4910e−58** |

### What changed vs. the failed S33 T=1000 run
1. **rE is quadrature-dominated, not precision-dominated.** In every past run rE = K·Eq exactly:
   - S33 L=1.28: 750 × 2.032e−65 = 1.524e−62;
   - S33 T=1000: 1090 × 3.557e−65 = 3.877e−62.
   Here M_S jumps to 8.7e6, because the 546-term lift contains cosh(0.4·|K·log p|). With ng=70 that would give rE ≈ 1e−59 > λ.
   Fix: ng=78 (factor ρ^16 ≈ 4.8e8). This gives rE = 1.74e−68, so λ/rE ≈ 2e6. `--prec 448` would not have helped, so prec stays 384.
2. **Ball blocked Schur fails because of radius growth (cond ≈ 1e62), even at 1536 bits.** All 10 tries per parity returned False, including l0 = 0.5 λ_est.
   It is replaced by `s36_cholres.py`, which uses approximate Cholesky plus a rigorous residual (Rump-style):
   - M = mid(Ã_s) − l0 I ≈ RᵀR, computed in point arithmetic at 640 bits;
   - then, in Arb balls, E = (Ã_s − l0 I) − RᵀR and ‖E‖₂ ≤ max row-sum |E|;
   - since RᵀR ⪰ 0, λ_min ≥ l0 − ‖E‖ − δ_str. Here δ_str ≤ 7.4e−112 bounds the decimal dump vs the 384-bit midpoints.
   - It has no elimination radius growth and takes **18 s per parity**.
   - Negative control: l0 = 1.02 λ_est and 1.1 λ_est both FAIL (approx Cholesky gives a negative pivot), so λ_min(even) ∈ (3.634e−62, 3.71e−62). Self-test on Hilbert(60) (λ_min 7.82e−90) passes at 7.02e−90 and fails at 5e−81.
3. Tails re-audited in Arb (`../s35/s35_tail_audit.py` → `tail_audit_L13_T700.json`): agrees in both parities.

Artifacts:
- `cert_L13_T700_K740_FINAL.json`: combined certificate (from `s36_finalize.py`).
- `cert_L13_T700_K740_ng78_p384s1536.{json,log}`: raw run, with ball-Schur nulls.
- `atil_L13_T700_K740_{even,odd}_cholres_p640.json` and `cholres_*.log`: PD step outputs.
- The dumps `atil_L13_T700_K740_{even,odd}.json` (68 MB each) stay on the box and are not committed; sha256 values are in `atil_sha256.txt`.

## H36b: {2,3,5} D=8 regrouped torus at L=1.3 — CERTIFIED envelope (float B&B + independent interval B&B)
- `s36_bundle_lp.py` is a proximal/trust-region cutting-plane LP. It fixes the S34/S35 exchange-LP oscillation (2.87–3.05).
  It started from the S35 LP and reached a grid-max of 2.7493. The grid NG=96 is too coarse: B&B found the true sup is 2.7770.
- **CERTIFIED: sup_{T³}(B−m) ≤ 2.7782.** Two checks:
  - float B&B `s34_bb_D4.py`: best sample 2.776968, 7.13e6 boxes;
  - **CONFIRMED by the independent libm-free `s35_interval_bb.py`** (max ub 2.778057, min invisible gap 2.69e−3).
  Artifact: `bb235/torus_L1.3_D8_p235_cert.json`.
- {7,11,13} D=8 (bundle best 1.6875–1.6880). Float B&B toward target 1.6920 was still running at commit time (best sample 1.688557, max ub 1.7177 at iteration 17).
  - If it closes: **A_cert ≤ 2.7782 + 1.692 ≈ 4.470**, giving β(700) ≈ 0.242 and β>0 down to T ≈ 2π e^{4.47} ≈ 549.
  - Status: NUMERICAL until that B&B finishes and is interval-confirmed. `bb71113/bb.log` is for S37 to harvest.

## Dead ends
- DEAD END: ball blocked-Schur at K≈740 with cond ≈ 1e62 at 1536 bits (and at 1024 bits for K=1090). Use cholres instead.
- DEAD END: `--prec` increases as a cure for rE (rE = K·Eq; raise ng instead).
- DEAD END (NUMERICAL): grid-based "true max" at NG=96 for D=8 axis frequencies up to 24 underestimates the sup by 0.028. Always B&B.
