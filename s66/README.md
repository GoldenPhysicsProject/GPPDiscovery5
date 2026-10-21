# Grok S66 (2026-10-17): Mori/Liu bounded comparison implemented and validated; Perron-s replay of support 2.8 at 571 modes/parity

> **RETRACTION (S67-resumed, 2026-10-21):** the λ_min(T) values below for L = 7/5 and 71/50 are inverse-iteration eigenvalues NEAREST ZERO, not minima; all those T fail Cholesky and the L = 71/50 T has certified negative Rayleigh witnesses (even < −5.58e−58, odd < −1.87e−48). Only the log 12 validation is positive. See ../s67/README.md.

RH OPEN. Finite windows only. φ absent. Note: GPP-bridge/work/grok/2026-10-17_s66_mori_route_validated_support28_replay.md

- `prior/`: Mori VII (= Mori [9], log 12, zenodo 23265728, v2) pdf + text. Mori VIII is in `../s64/prior/`.
- `mori_consts.py`: Arb evaluation of the scalar side of the method, using the formulas in Mori [9]:
  - tail constant C;
  - δ ≤ d(c);
  - E0, E1 and the U-coefficients u_j (with u_j·E_j ≤ C);
  - the Lemma 3.2 projection errors e, n, h (ρ = 3);
  - α, β of (3.3).
  - It reproduces every published constant of Mori VII (log 12) and Mori VIII (log 15) to the printed digits.
- `mori_route.py`: our own implementation (it does not use Mori's code) of the finite condition, in one parity:
  - (3.4) of Mori [9] = (2.2) of Mori VIII: T = W(X)⁻¹ + J − αI − β(B − A²);
  - A = E*ME and B = E*M²E with M = γ − C_p, by exact Gauss–Legendre on the cells cut at the shift breakpoints;
  - J = E*(K + U)E: band-limited archimedean part (panels, Arb spherical Bessel, as in s53) + pole + U;
  - X = (B − mA)⁻¹(A − mI), the Loewner minimiser of W(X);
  - it streams the partial sums, writes a pickle dump of A, B, J, T, and estimates λ_min(T) by inverse power iteration.
  - **NUMERICAL** (midpoints at 448–512 bits). It is not a certificate: no ball radii are kept, and the J quadrature has no Bernstein-ellipse error bound yet.
- `mori_post.py`: recomputes T from a dump with X = opt or X = B⁻¹A.
- `run/*.jsonl`, `run/*.log`: results. The `*.pkl` dumps (tens to hundreds of MB) stay on the box and are not committed.

| case | N (per parity) | λ_min(T) even | λ_min(T) odd | reference |
|---|---|---|---|---|
| log 12, Mori VII parameters (Ω 416, s 57/16, m ½, u 207/68) | 832 (416) | 1.20926402186e−54 | 1.58671301389e−50 | Mori VII Table 2 certified lower bounds: 1.2092639998e−54 / 1.5867129002e−50 (**match to 8 / 7 digits**) |
| same, X = B⁻¹A | 832 (416) | 1.18857e−54 | 1.58055e−50 | (worse X: shows that X matters) |
| **L = 7/5, s = 464/125 (Perron, ≥ 3.71150), m ½, Ω 488, u 296/98** | 1142 (571) | **8.72098e−78** | **1.14435e−73** | ours, S63 (other method): Q ≥ 7.4196e−78 |
| L = 7/5, same s, m ¼, Ω 384, u 249/83 | 898 (449) | 9.45608e−78 | 6.26656e−74 | |
| **L = 71/50 (support 2.84 > log 17), s = 763/200 (≥ 3.81418), m ½, Ω 536, u 308/102** | 1272 (636) | **1.85047e−81** | **2.85828e−77** | NUMERICAL only; to be certified in S67 |

All are NUMERICAL (positive at midpoint precision; not certificates). Detached follow-up: `run/chain_s66.sh` (log `run/chain_s66.log`), which runs m = ⅛ and N = 800 at 7/5, and m = ¼ at 71/50.

Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.

File notes: `run/val_log12.jsonl` comes from the first run, which used X = B⁻¹A. The X★ values are in `run/post_log12_p{0,1}.txt` (mori_post.py on the same dumps). All the L = 7/5 and 71/50 runs use X★. The dump sha256 sums are in `run/pkl_sha256.txt`.
