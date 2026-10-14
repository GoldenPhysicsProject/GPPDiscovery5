# blocks_moments — Daniel's "Principal-Series Blocks on the Celestial Sphere" v2
Status labels as in GPP-bridge. All outputs are mpmath FLOAT, not interval-certified.

- principal_series_blocks_check.py / .out: Daniel's script, every check agrees (rel err <= 3e-24).
- moments_check.py / .out: Daniel's script, 4/4 PASS.
- moments_hp.py / .out: A_k = (1/2pi) int P D^k for k=0..5 at 75 and 90 dps, plus PSLQ.

Results (NUMERICAL, 75-digit agreement):
- A_3 - (log2/2 - 3/16) = -1.4e-72, so A_3 = log2/2 - 3/16 holds to 72 digits. It is NOT proved.
- A_4 = 3/8 - log2 + 2 log^2 2 - pi^2/24 (CONJECTURE from PSLQ; stable across 3 basis sizes up to 13 elements incl. Li4(1/2)).
- A_5 = (-45 + 120 L - 240 L^2 + 5 pi^2 + 63 zeta(3) + 320 L^3 - 20 pi^2 L)/48, L = log 2 (CONJECTURE from PSLQ; stable in bases 8 and 12).
- Pattern: A_k has weight <= k-1, built from level-2 (alternating) MZVs. No phi appears.
