# Results

Tables written up from `results.jsonl`, newest first. Every row says whether it is floating-point or certified.

| date | experiment | result | label | float/certified |
|---|---|---|---|---|
| 2026-10-14 | `li_coefficients/li_ftheta_onset.py 0.25 600` | λ_n(F_{1/4}) is first negative near n≈9250 (250-step grid); −5.5e5 at n=10500 | NUMERICAL (ξ part asymptotic) | float |
| 2026-10-14 | `li_coefficients/li_coeffs.py 100` | λ_1..λ_100 > 0 for ξ, F_{1/4}, F_{1/20}, DH; λ_1(ξ)=0.023096 | NUMERICAL | float |
| 2026-10-14 | `dh_zero/dh_zero.py` | DH zero 0.808517182456637 + 85.699348485377592i, \|f\|≈8e-32, FE residual 2.5e-32 | NUMERICAL | float (dps≈30) |
| 2026-10-14 | `pole_parity/pole_parity.py` | Weil pole term: even +1.71887, odd −0.00528, mixed 1.69775 ≠ 2G(½)²=2.12099 (pole = 2G(½)G(−½)) | NUMERICAL | float |
| pending | `ftheta_window/ftheta_zero_sum.py <coeffs> 150` | reproduce F_{1/4} window sign (h=19/25: −0.02197, ζ +0.00203; h=17/25: −0.00886, ζ +0.00081); needs Claude Code's private coeff files | NUMERICAL | float |
