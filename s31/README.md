# Grok S31 (2026-10-16) — RH OPEN, GRH OPEN
- `s31_bb_D4.py` — rigorous sup bound for the S30 D=4 torus lift (reuses s28_torus_bb float64 B&B + rounding margin; Arb invisibility check).
  CERTIFIED (float B&B with explicit rounding margin, libm 2-ulp assumption, same standard as S28): A(L=1.0) <= 2.6329 (`torus_L1.0_D4_cert.json`, 22.4M boxes),
  A(L=1.19) <= 3.4086 (`torus_L1.19_D4_cert.json`, 6.9M boxes). Per-prime: 3.0476 / 3.9227. Floor λ_max(P_L): 1.9306 / 2.6613.
- `cert_L119_T300_torusD4.json/.log` — end-to-end test of the `--torus` path: ζ Weil positivity on support 2.38 now at T=300, K=312
  (vs T=700, K=730 in S28): even >= 2.8396e-48, odd >= 2.9199e-44 (Arb tails audit included). CERTIFIED (Arb, Lean-open, single implementation).
- `s31_lp_L125.py` — time-capped torus LP. 5-prime D=4 at L=1.25: A≈3.858 (NUMERICAL, round 4; 5D B&B infeasible in-session).
  4-prime {2,3,5,7} D=4 at L=1.25 + separate p=11 per-prime lift: `s31_merge_p11.py` → A_cert = 3.4646 + 0.74248 = 4.2071 (vs per-prime 4.5903).
- `cert_L125_T700_torus_ng70.json/.log` — support 2.5 attempt (T=700, K=730, ng=70, prec 384, β=0.5047): ODD CERTIFIED >= 2.3480e-51 (Arb); EVEN PENDING (λ_est 3.30e-55, rE 4.1e-61, ball Schur at 384 bits failed; rerun --prec 512). ng=60 run stopped (M_S 4.8e5 too large).
