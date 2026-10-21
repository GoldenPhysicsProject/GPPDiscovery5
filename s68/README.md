# Grok S68 (2026-10-21): Mori VIII log 15 replay with the fixed S67 code; the negative direction at L = 71/50 and 7/5
H68a result: BOTH parities CERTIFIED (even >= 2.395e-70, odd >= 3.9419e-66; controls fail), matching Mori VIII. RH OPEN. Finite windows only. φ absent. Note: GPP-bridge/work/grok/2026-10-21_s68_mori_viii_log15_replay_negdir.md

Code: `mori_cert.py`, `lam_from_dump.py`, `s36_cholres.py`, `s37_ldlt_exact.py`, `neg_witness.py` and `quad_bound.py` are symlinks to `../s67` (unchanged).

- `run/chain_s68a.sh` → `run/chain_s68a.log`: H68a. Mori VIII parameters: L = ½log 15, N = 1664 (832/parity), Ω 736, s 33/8, m ½, u 428/142, δ = 2^−2901, 512 bits.
  Parities run strictly one after another. Per parity: certified assembly, then the nearest-zero estimate, then s36 cholres (0.999 × est), then s37 exact LDLᵀ (0.999 × est and the control at 1.005 × est).
  - s36 stops at the first success, so the failing control comes from s37.
- `run/quad_log15.json`: CERTIFIED J quadrature bound. eps_J ≤ 4.94e−112 and eps_bb ≤ 7.49e−113 (51 panels, 96 nodes, L ≤ 13541/10000).
- `h68b_bisect_vec.py`: H68b, NUMERICAL at 640 bits.
  - Bisects λ_min with approximate Cholesky (success means T − l0 I > 0).
  - Then runs shifted inverse iteration at the success end and gives the mode profile of the eigenvector.
  - Run by `run/h68b_chain2.sh`. Outputs: `run/h68b_vec_*.json`.
- `h68b_negdir.py` / `run/h68b_negdir_p0.log`: Ritz values of the leading k×k blocks (even, L = 71/50).
  - Its shifted inverse iteration at 1.5 × Ritz₃₂₀ did NOT converge (the shift was on the wrong side). That part is superseded by `h68b_bisect_vec.py`.
- The T dumps (`run/T_log15_p*.json`, about 110 MB) stay on the box and are not committed.
