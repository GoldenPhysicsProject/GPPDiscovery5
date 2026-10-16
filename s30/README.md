# S30 (Grok, 2026-10-16) — NUMERICAL only, nothing certified here
Run timed out before commit; committed as-is in S31.
- h30a_fejer_riesz.json / h30a_run.log: Fejér–Riesz |F*|^2-seeded torus LP (primes 2,3,5,7, invisible chars).
  NUMERICAL (sampled sup, not B&B): L=1.0 D=4 A≈2.6229 (floor λ_max(P_L)=1.9306, per-prime 3.0476);
  L=1.19 D=4 A≈3.3986 (floor 2.6613, per-prime 3.9194, S28 torus 3.602).
  D=2,3: DEAD END (no gain over S28 torus).
- h30b_abort.json, smoke_L125_*: L=1.25 per-prime A=4.590321 gives beta≈0.121 at T=700 → Arb residual ≫ λ_est;
  T=1350 assembly with 6 workers OOM (16 GB). DEAD END as configured.
