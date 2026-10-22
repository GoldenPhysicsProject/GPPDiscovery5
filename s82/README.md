# s82 — Grok S82 H80c step 2: speed levers for the Ozaki-sliced cell products (support 2.84 even). RH OPEN.
`levers.py` regenerates the Arb V, F, FW (and G = F·diag(sqrt w)) for one x-space cell (cell 6) and `nch` consecutive 64-node quadrature
chunks with the unchanged s71/mori_cert.py code (L=71/50, N=2376, parity 0, prec 512, leading `--nm` modes). It times Arb
`Cp = V*FW; Cp2 = F*FW` on the concatenated K = 64·nch inner dimension (ctx.threads 1 and 8), then three sliced variants per batch size k:
- `Cp`: anti-diagonal concatenated dgemm, b = largest with 2b + ceil log2(S·k) ≤ 53 (k 64: b 21 S 23; 256: b 20 S 25; 512/1024: b 19 S 26).
- `Cp2_gemm`: G·Gᵀ by the same dgemm. `Cp2_syrk`: upper triangle via dsyr2k (pairs i<j per anti-diagonal) + dsyrk (i=j).
  Gauss–Legendre weights 2hh·w_k are asserted > 0 in Arb, so sqrt splitting is valid; G carries its own Arb radius.
- amortised carry: all batch products are added into int64 base-2^b limbs (each product < 2^53, headroom 2^10), one carry at the end.
  `t_carry_if_per_chunk64` = what a per-chunk carry would have cost.
Radius (entrywise, float64 × (1+1e−9)): gemm (rA+uA)(|W|+rW)ᵀ + |A|(rW+uW)ᵀ; syrk |G|Eᵀ + E|G|ᵀ + EEᵀ with E = rG+uG.
Validation: exact dyadic |mid_arb − sliced| ≤ rad_arb + rad_sliced for every entry (both triangles for syrk, mirrored).
`--control`: V[0,0] and G[0,0] perturbed by 2^20 fixed-point units (must fail).
Logs: levers_np600.jsonl (Np 600, K 512, k 512/256/64, 1 & 8 threads, validated at 1 thread), levers_control.jsonl (Np 300, fails),
smoke_nm200.jsonl (Np 200, K 256), levers_np1188_k1024_partial.jsonl + np1188_k1024_aborted.log (Np 1188, K 1024: Arb timings and
the k=64 Cp row only; the run was stopped by our RSS guard at 3.6 GB before the other variants).
Not done: strip-wise full assembly, dump reproduction, PSD re-check.
