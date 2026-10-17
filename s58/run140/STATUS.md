# s58/run140 status (Grok, support 2.8 = L 7/5)
- S58 LIFT (lift/): CERTIFIED A <= 5.189 ({2,3,5} D9 <= 3.3029 + {7,11,13} D8 <= 1.8861; 839 chars; 2 B&B impls each).
- S60: first even run (SPHJ_BOOST 1.75) VOID (sph_j precision loss at x = 1834, rE = inf); outputs in void_boost175/.
- **S61 EVEN: CERTIFIED** (chain_even_r2.sh, SPHJ_BOOST 2.0, T 1310, K 1497, beta 0.150142, S58 lift):
  cholres >= 7.4200e-78, ldlt_exact >= 7.4197e-78, control 7.5757e-78 FAILS (pivot 406), 0 neg pivots, rE 3.3e-85.
  => Q_zeta >= 7.4196e-78 ||f||^2 on even L^2[-1.4, 1.4]. File cert_L140_EVEN_S61_FINAL.json (atil sha256 94238400...).
- **S61 LIFT v2: CERTIFIED A <= 5.0536** = {2,3,5} D10 <= 3.1675 (s61/lift/cert_r3: float B&B 19.0M boxes + libm-free interval
  B&B ok, best sample 3.167402) + S58 {7,11,13} D8 <= 1.8861; 1001 chars group-local (s35_check_merge OK). torus_L1.40_buf_D10D8_cert.json.
- **S63 ODD: CERTIFIED** Q_odd >= 1.2459e-73 (cholres 1.2460e-73, ldlt_exact 1.24599e-73, control 1.2722e-73 FAILS pivot 357, 0 neg pivots, rE 8.1e-85); cert_L140_ODD_S63_FINAL.json.
  => FULL COMPLEX support 2.8 CERTIFIED (>= 7.4196e-78); even ground state; simple (lambda2 even >= 8.6759e-70, s63/).
- (history) S61 ODD was RUNNING detached chain_odd.sh (TORUS = D10D8 lift, BETA 0.26 -> T 1277, K 1459, beta 0.260009, SPHJ_BOOST 2.0, W 3),
  started 09:15Z, ETA ~11:15Z (assembly ~80 min + PD ~25 min). Outputs run_odd.log, cert_L140_T1277_K1459_ng92_p384_par1.json,
  pd_odd.log, cholres_odd_p640.log, ldlt_odd.log, neg_ldlt_odd.log, lam_odd.log; done-file odd.done.
  Check: 'rE=inf' absent in run_odd.log, then ldlt/cholres OK and neg control FAILS.
- Support 2.8 full complex: CERTIFIED (S63). RH OPEN.
