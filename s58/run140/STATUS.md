# s58/run140 status (Grok S58)
- LIFT (lift/): CERTIFIED A <= 5.189 at L = 7/5 = {2,3,5} D9 <= 3.3029 (best sample 3.29991; s34 float B&B + s35 libm-free interval B&B both ok)
  + {7,11,13} D8 <= 1.8861 (best sample 1.88349; both ok). 839 characters, group-local (s35_check_merge OK). Buffer EPS = 0.0056;
  closest characters: log(50/3) = 2.81341 (gap 0.013411, eps*T = 17.6) and gap 0.020796 for {7,11,13}. Comb now includes n = 16.
- EVEN (T = 1310, K = 1497, ng = 92, prec 384, beta = 0.150142, W = 3): RUNNING detached (lift done ~05:23Z, assembly started ~05:30Z; expect ~07:00Z) (chain_even.sh PID 54432, python PID 60944).
  Outputs: run_even.log, cert_L140_T1310_K1497_ng92_p384_par0.json, atil_L140_T1310_K1497_even.json (gitignored), pd_even.log,
  cholres_even_p640.log, ldlt_even.log, neg_ldlt_even.log, lam_even.log; done-file even.done. NOT yet certified.
- ODD: not started (needs beta ~ 0.25, i.e. T ~ 1450, K ~ 1655). RH OPEN.
