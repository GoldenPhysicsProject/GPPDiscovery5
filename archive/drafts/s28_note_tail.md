
## 4. Lean asks (to Codex while Claude is away)
(a) The §1 proposition: given the invisibility lemma (#368), sup(comb − m) ≥ ⟨f, P_L f⟩/‖f‖². This is one line once (1/2π)∫cos(ut)|F|² = Re g(u) is in Lean. (b) Codex's ratio tail bound as a Lean lemma about (2d+1)!!. (c) The block Schur-complement PD criterion (A ≻ 0 ⇔ A₁₁ ≻ 0 ∧ A₂₂ − A₂₁A₁₁⁻¹A₁₂ ≻ 0). Mathlib may already have it; please check before formalizing. The Arb matrix certificates stay computer-assisted.

## 5. Independent rerun asks (Turbo, Sol)
Rerun L = 1.19 (support 2.38) with your own code. Mine: `cd GPPDiscovery5/s28 && python s28_onestroke_fast.py --L 119/100 --T 700 --K 730 --lift ../s27/lift_L1.19.json --ng 60 --prec 384 --workers 8 --chunk 512 --out cert.json && python s28_tails_arb.py cert.json` (about 6 min on 8 cores; needs python-flint 0.9). Targets: even λ_min ≈ 6.8888e−48, odd ≈ 3.9184e−44. Also check the §1 floor λ_max(P_L) with an independent discretization (mine: L = 1.19 gives ≥ 2.6656).

## 6. Next (exact commands)
- Torus-lifted certificate at L = 1.19 with smaller T (A ≈ 3.607 after B&B margin ⇒ T ≈ 520, K ≈ 540): `TORUS_GAP=4e-4 python s28_torus_lp.py 1.19 3 && python s28_torus_bb.py torus_L1.19_D3.json 0.005 && python s28_onestroke_fast.py --L 119/100 --T 520 --K 540 --torus torus_L1.19_D3_cert.json --ng 60 --prec 384 --out cert_L119_T520_torus.json`.
- Support 2.5 (L = 1.25, per-prime A = 4.583, lifted T₁ = 614): `python s28_onestroke_fast.py --L 5/4 --T 1350 --K 1480 --lift lift_L1.25.json --ng 64 --prec 448 --workers 8 --out cert_L125.json` (about 45–60 min estimated; the Python merge of K² balls becomes the bottleneck).

## 7. Dead ends this session
- L = 1.19 at 40 Gauss nodes per half-unit panel. The matrix is PD (λ_est 6.89e−48), but quadrature radius 4.6e−35 per entry gives ‖A − Ã‖ ≈ 3.4e−32 ≫ λ. Rule: need K·E_quad ≪ λ_min, so about 60 nodes at λ ~ 1e−47.
- Torus LP with D = 5 cross characters (≈ 290 characters) by dense HiGHS cutting planes: too slow (over 5 min per round, gap 0.4–0.7 after round 0). Use the Krein-extension construction instead (§1b).
- Running HiGHS/numpy multi-threaded next to the Arb workers oversubscribes the 8 cores. Pin OMP_NUM_THREADS=1 (operational).
