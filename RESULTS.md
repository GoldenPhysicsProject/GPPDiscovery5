# Results

Tables written up from `results.jsonl`, newest first. Every row says whether it is floating-point or certified.

| date | experiment | result | label | float/certified |
|---|---|---|---|---|
| 2026-10-14 | `li_certified/li_arb.py 2001 6000` | λ_1..λ_2001(F_{1/4}) all > 0 (balls exclude 0); L=10600 run for sign change near 9250 in progress | CERTIFIED | certified (Arb) |
| 2026-10-14 | `dh_phi/dh_phi.py` | DH κ = √(φ+2)−φ = 5^{1/4}√φ−φ; ε(χ)=(φ+i)/√(φ+2); φ from choice q=5 | PROVED LEMMA (algebra) + Arb consistency | certified check |
| 2026-10-14 | `ftheta_window` h=19/25 (even) | F_{1/4} −0.021969145467, ζ +0.002030094923, matches expected | NUMERICAL | float |
| 2026-10-14 | `ftheta_window` h=17/25 (odd) | F_{1/4} −0.008863529402, ζ +0.000809109424, matches expected | NUMERICAL | float |
| 2026-10-14 | `li_coefficients/li_ftheta_onset.py 0.25 600` | λ_n(F_{1/4}) is first negative near n≈9250 (250-step grid); −5.5e5 at n=10500 | NUMERICAL (ξ part asymptotic) | float |
| 2026-10-14 | `li_coefficients/li_coeffs.py 100` | λ_1..λ_100 > 0 for ξ, F_{1/4}, F_{1/20}, DH; λ_1(ξ)=0.023096 | NUMERICAL | float |
| 2026-10-14 | `dh_zero/dh_zero.py` | DH zero 0.808517182456637 + 85.699348485377592i, \|f\|≈8e-32, FE residual 2.5e-32 | NUMERICAL | float (dps≈30) |
| 2026-10-14 | `pole_parity/pole_parity.py` | Weil pole term: even +1.71887, odd −0.00528, mixed 1.69775 ≠ 2G(½)²=2.12099 (pole = 2G(½)G(−½)) | NUMERICAL | float |
