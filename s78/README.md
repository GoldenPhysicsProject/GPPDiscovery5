# Grok S78 — H77a: can the high block of the Mori-route T be assembled at low precision? (RH OPEN)
Scripts: time_cell.py (one x-space cell job, Arb vs float64/BLAS), mori_ingr.py (copy of s71/mori_cert.py that perturbs
ingredient Gram blocks and measures the induced change of T by block; also the block-diagonal-X probe H78b).
Logs: run/time_cell_p{192,320,512}.json, run/sens_N500_p0.log, run/sens_N700_p1.log, run/h78b_N500_p0.log, run/h78b_N700_p1.log.
All results NUMERICAL (no new certificate; the support-2.84 certificates of S72/S74/S76 are unchanged).

1. Profile (from s71 log, odd Ω1150 N 2732, 3 workers): panels (J band) 255 s, x-space cells (Cp, Cp2) 5701 s → cells = 96 %.
2. Timing, one cell job (64 nodes, cell 5, N 2732, odd; 3 jobs concurrent): Arb prec 512: 19.8 s, 320: 14.6 s, 192: 10.95 s
   (≈5.6 s precision-independent Legendre part at +1.27N boost bits, ≈0.028 s/bit products). float64/BLAS same quadrature: 0.085 s
   (233× vs 512 bits), max |float − Arb| = 8.6e−16 (Cp), 3.0e−15 (Cp2) (NUMERICAL, no rigorous float bound yet).
3. Sensitivity (L 71/50 small instances: even N500 Ω210 NL100; odd N700 Ω300 NL200; same X, prec 320): a random symmetric
   perturbation of size eps in the HIGH block of Cp changes T by ≈1.2–6.7·eps in LL, ≈8.4–9.7·eps in LH, ≈40–45·eps in HH; Cp2 ≈0.02–0.07·eps LL,
   0.5·eps LH; Ja (J) only HH (gain 1). Linear in eps (1e−20, 1e−40). Reason: T = W(X)⁻¹ + J − αI − β(B − A²) and W⁻¹ is global.
   ⇒ the S77 T-level budget (C to 1e−4) does NOT transfer to the ingredients: high-mode Cp/Cp2 must be accurate to ≲ a/n
   (≈1e−84 even at 2.84), i.e. ≳ 290 bits. Only the J band (4 % of time) may go low precision in its HH block.
4. H78b (block-diagonal X to decouple): Loewner loss T(X_bd) − T(X*) ⪯ 0 with eigenvalues down to −5.8 and O(1) LL entries
   (X*_LH up to 0.30 / 0.023). DEAD END (numerical): an O(1) loss against a 1e−81 margin.
5. Extrapolation to L ≈ 1.71 (Ω ≈ 6.25–7.23k, N ≈ 17.9–20.7k, Np ≈ 9.6k, ~16 prime-power shifts): per cell job ≈ 450 s products
   (320 bits) + ≥ 450–900 s Legendre, × ≈ 7.5k jobs ≈ 1e7 CPU-s (≈ 40 days on 3 workers); one Np² Arb matrix at 320 bits ≈ 8–10 GB and
   ≥ 6 are live ⇒ does NOT fit 16 GB. Two-level elimination does not cut this because the whole T is needed at ~300 bits.
Fakes: none ruled out (feasibility study). A ζ PSD certificate at L ≥ 1.7105 would rule out DH; F_θ, H₋ₜ status unchanged (S63).
Next (H78c, HYPOTHESIS): cut the Arb cost without lowering accuracy — (a) Legendre/Bessel values by the stable forward/Miller
recurrence at ~320 bits with an a-priori error bound instead of the +1.27N-bit ball boost; (b) products by an error-free
Ozaki-type slicing into doubles (exact dgemm slices, rigorous). Expected 5–10× on cells; L ≈ 1.71 would still need ≫ 16 GB.
