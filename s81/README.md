# s81 — Grok S81 H80c step 1: Ozaki-sliced cell products vs Arb (support 2.84, even). RH OPEN.
`sliced_cell.py`: regenerates, with the s71/mori_cert.py code (legv, shifts, cells), the Arb V, F, FW of one x-space cell job
(L=71/50, N=2376, parity 0, prec 512, chunk k=64; leading `--nm` modes), computes Cp=V·FW and Cp2=F·FW in Arb and by exact
integer slices (float64 dgemm), and checks |mid_arb − sliced| ≤ rad_arb + rad_sliced entrywise.
- fixed point: X_mid rounded to 2^−(P−e), e = ceil log2 max|X|; u = 2^−(P−e+1).
- rad_sliced = Σ_c (rV+uV)(|FW|+rFW) + |Vr|(rFW+uFW) + dropped-pair bound, in float64 × (1+1e−9).
- default: S pairs dgemm, b = 23 (k = 64). `H80C_CAT=1`: one dgemm per anti-diagonal (inner dim ≤ S·k, b = 21), faster.
- `--keep q`: drop slice pairs i+j > S−1+q, with a rigorous bound added. `H80C_CONTROL=1`: perturb V[0,0] by 2^20 units (must fail).
- limb path (int64 base-2^b carry propagation) is checked equal to the big-int recombination on the first/last row.
Logs: h80c_nm300.jsonl (9 jobs, 3 cells × 3 chunks, 300 modes), h80c_full.jsonl (Np 1188), h80c_keep2.jsonl,
h80c_timing.jsonl (Np 1188 timings), h80c_control.jsonl (control, fails as required).
Not done: full-matrix sliced assembly / dump reproduction / PSD re-check (H80c step 2).
