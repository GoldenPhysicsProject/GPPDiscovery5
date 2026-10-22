**Grok S82 (H80c step 2: speed levers, validation CHECKED, timing NUMERICAL, L ≈ 1.71 ESTIMATE). RH OPEN.**

**Comms read since S81 (6330667732):**
- Turbo 6331078651 + b86ae39 relays Daniel's new standing directive (all seats): work on what a zero is as a state of the sewn theory, and where a non-unitary mirror pair is excluded; analytic lanes continue as support and falsification. Grok acknowledges. Our Weil-positivity certificates are the finite-support shadow of exactly that reflection form (Ŵ ≥ 0), so from S83 each Grok run carries a structural slot. First structural HYPOTHESIS to test: the near-null vector of the certified support-2.84 even form, read as the boundary state closest to losing unitarity. Check whether it localises at the lowest zeros, and whether DH's failure at L ≥ 1.7105 is a mirror pair forming in that same direction (S75 showed DH fails along a ζ-stiff direction).
- Codex, Claude, Astra, Sol: nothing new. Nothing else to answer. Still waiting on a teammate replay of the support-2.84 certificates (S72/S74) and on Turbo's W(14.5) = 158.1 vs our certified W(log 2e6) = 164.5143 convention question.

**Own line: the three speed levers from S81, measured like for like against Arb on the same inputs** (even support 2.84: L = 71/50, N = 2376, prec 512, dump sha b609d6f5…; cell 6, consecutive 64-node chunks; V/F/FW from the unchanged S71 code). This time Arb was given the **same batched inner dimension** as the sliced code.
- CHECKED, entrywise ball containment with exact dyadic differences. All three levers pass:
  - (1) batching to k = 256/512 (b 20/19 so 2b + log2(S·k) ≤ 53);
  - (2) Cp2 = F diag(w) Fᵀ = G Gᵀ with G = F·diag(√w), using dsyrk + dsyr2k on the upper triangle (the Gauss–Legendre weights are > 0, asserted in Arb; radius |G|Eᵀ + E|G|ᵀ + EEᵀ);
  - (3) amortised int64 carry, once per cell.
  - Worst diff/(rad_arb + rad_s): 0.31 for Cp, 0.11 for Cp2 (Np 600, K 512; Np 200, K 256), and 0.31 for Cp at Np 1188, K 1024.
  - Control: V[0,0] and G[0,0] + 2^20 units **fails** (ratios 1863 and 3559).
- NUMERICAL, Np 600, K 512, both products:
  - Arb: **4.75 s** at 1 thread, **1.20 s** at 8 threads.
  - Best sliced: 4.24 s (**1.12×**) at 1 thread, 2.03 s (**0.59×**) at 8 threads.
  - S81 configuration (k 64, carry per chunk): 0.93× and 0.38×.
  - Np 1188, K 1024: Arb takes 29.6–31.0 s (1 thread) and 7.9 s (8 threads). Sliced Cp alone takes 18.1 s.
  - Larger k doesn't help the sliced route (b shrinks). syrk saves 15–25 % at 1 thread but loses at 8. The amortised carry removes ~10 %.
- **Correction to S81 / dead end:** the 2–3.5× in S81 came from Arb running at k = 64. One Arb mul at K = 1024 is **3.8–4.2× faster** than sixteen k = 64 jobs and scales ~4× on 8 threads. So on this CPU box the sliced float64 route gives no net speed-up. The real lever is batching Arb itself. The exactness and radius machinery stays CHECKED in case GEMM-accelerated hardware appears.
- ESTIMATE for ζ at L = 1.7105 cell work: ≈ 1e7 / 4 ≈ 2.5e6 CPU-s ≈ **3.5 days on 8 cores**, or ≈ 2.7 days if Cp2 is done as upper-triangle strips (HYPOTHESIS). This assumes the batched-Arb gain holds at Np ≈ 9.6k (untested). The out-of-core strip plan is still required (a full Np² Arb output is ≈ 12 GB).
- NOT DONE: the strip-wise full assembly, the dump reproduction and the PSD re-check (s36/s37 with pos 0.999 / control 1.005). None of these is claimed. Next: the S71 assembly with batched Arb cells (K 1024) to rebuild b609d6f5… within radius, then s36/s37.
- Fakes (F_θ, DH, backward heat): none ruled out. This is assembly engineering. Non-vacuity checks: the control fails, syrk is validated on both triangles, and the validated result is the amortised-carry one.

Note: `work/grok/2026-10-22_s82_h80c_levers.md`; scripts in GPPDiscovery5 `s82/`. No asks for other seats this run.
