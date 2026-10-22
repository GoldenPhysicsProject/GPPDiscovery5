# Precision needed for the ζ window factorization at L = 1.7105, with a DH failing control

Goddess of Prime Number Theory, 2026-10-22. **RH is OPEN.** Daniel Toupin / Golden Physics Project, with AI assistance from Grok.

**Convention.** L is the **half-width**: f ∈ L²(−L, L), g = f * f̃ on [−2L, 2L], prime powers with log n < 2L. This is the
team's binding convention, and it is the one used by `s62/cosritz.py`, which we use for assembly. L = 1.7105 = 3421/2000.
Basis: even-sector cosine modes cos((2m+1)πu/(2L)), m = 0..K−1. This basis is nested, so leading K×K blocks are the K-mode sections.
Odd sector, where stated: sin(mπu/L).

Labels: CERTIFIED = Arb-rigorous; NUMERICAL = floating point at stated precision, checked against a higher-precision reference.
**Every statement here is about finite K-sections.** PD of a section is necessary for Weil positivity, not sufficient: there is no tail bound.
The full certificate also needs the Mori-route truncation control, which is not done here.

## 1. ζ at L = 1.7105: λ_min of the K-section (NUMERICAL, 960-bit Arb assembly, inverse iteration)
| K | 120 | 160 | 200 | 266 |
|---|---|---|---|---|
| λ_min | 7.67e-137 | 1.32e-148 | 6.44e-153 | **3.09e-153** |

From K = 200 to K = 266 the value is still drifting by a factor of about 2. −ln λ(266) = 351.2. This matches the interpolated Zhu
upper-bound curve at L = 1.7105 exactly (−ln = 351.2, `zhu_law.out`). The C ≈ 20.13 law prediction 10^(−151.8) is high by
about 1 decade. **Estimate for the true q_min(1.7105), even sector: 10^(−153.3 ± 0.7).** This is CONJECTURE-grade
extrapolation: the upper end is the K = 266 value, and the lower end allows the K drift to continue.

## 2. Precision scan (ζ), K = 200 and 266, assembly entries modelled at p bits (radius ≥ 2^{8−p}|a_ij|)
| p (bits) | 480 | 496 | 504 | 512 | **520** | 528 | 544 | 576 | 640 |
|---|---|---|---|---|---|---|---|---|---|
| float LDLᵀ (all pivots > 0) | **fails** (pivot 158) | ok | ok | ok | ok | ok | ok | ok | ok |
| CERTIFIED PD, K = 266 (‖D^{-1/2}L^{-1}EL^{-T}D^{-1/2}‖_F) | — | 4e21 ✗ | 2e14 ✗ | 1e7 ✗ | **0.61 ✓** | 1e-6 ✓ | 1e-11 ✓ | 3e-21 ✓ | 1e-40 ✓ |
| CERTIFIED PD, K = 200 | — | ✗ | ✗ | ✗ | **0.010 ✓** | ✓ | ✓ | ✓ | ✓ |
| CERTIFIED PD, K = 160 (λ = 1.3e-148) | fails | **0.075 ✓** | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| CERTIFIED PD, K = 120 (λ = 7.7e-137) | **✓** | ✓ | | | | | | | |

**Threshold: minimal working precision is 520 bits (≈ 157 decimal digits) at K = 200–266.** Below that, the factorization is not resolved.
- Float LDLᵀ breaks at ≤ 480 bits. A negative pivot appears at index 158, i.e. ζ would be wrongly reported as failing.
- At 496–512 bits the float factorization succeeds, but it cannot be certified.

The empirical rule is p_cert ≈ log2(1/λ_min) + 5…13 bits, with ‖W‖ ≈ 5. **Recommendation: 576 bits (≈ 174 digits).** That gives
≥ 20 orders of magnitude of headroom on the certificate norm. Entries must then be assembled to ≈ 2^{−568} absolute: Arb assembly
at 576–640 bits. Our assembly at 960 bits has radii ≤ 1e-284.
This **revises down** my earlier 1050–1150-bit estimate for the LDLᵀ. A (λ_min)²-type requirement does *not* occur for the congruence/residual certificate.
It would occur only for a residual bound that ignores the D-scaling.

Method (`ldl_certify.py`, CERTIFIED):
1. Run p-bit float LDLᵀ to get the exact binary L̃ and D̃.
2. In Arb, compute E = A_ball − L̃D̃L̃ᵀ and G = D̃^{−1/2}L̃^{−1}EL̃^{−T}D̃^{−1/2}.
3. If ‖G‖_F < 1, then A = L̃D̃^{1/2}(I+G)D̃^{1/2}L̃ᵀ is PD for every matrix in the ball.

A naive ball LDLᵀ (`ldl_scan.py ball`) never certifies ζ at any p ≤ 896, because of the wrapping effect. That is a dead end.

## 3. Failing control: DH at the same L (same assembly code, FORM=DH, S47 form)
- K = 266 has λ_min = **−2.31e-31** (NUMERICAL; next eigenvalues 8.4e-22, 5.7e-16). The witness is **low-frequency and nodeless**: 80% in mode 0 and 20% in mode 1.
- Ball LDLᵀ: pivots 0..239 are positive balls and pivot 240 = [−7.3303 ± 3e-5] → **CERTIFIED: the DH window form at L = 1.7105 is not PSD**. This holds from 256 bits up. At 192 bits it is unresolved.
- The ζ certifier rejects DH at every p (FLOAT_LDL_FAILED at pivot 240).
- K = 200 alone is PD for DH (λ = 6.9e-31, certified at 256 bits). **DH's failure at 1.7105 needs K ≥ 241 modes.** The negative direction is low-frequency, but it is only resolved at the 1e-31 level once high modes are included.
- So precision distinguishes the two cleanly:
  - DH's negative eigenvalue (−2.3e-31) is resolved from ~192–256 bits.
  - ζ's margin (3e-153) needs 520 bits.
  - Between 256 and 480 bits a float LDLᵀ "fails" ζ and DH alike. That is a false negative for ζ.

## 4. DH failure onset L* (half-width), converged in K
λ_min of the DH even section, K = 320 (NUMERICAL; `runs/onset/dh_margin_even_K320.jsonl`):

| L | 1.5 | 1.6052 | 1.6578 | 1.6973 | 1.7072 | 1.7103 | 1.7105 |
|---|---|---|---|---|---|---|---|
| λ_min | 4.4e-18 | 3.5e-23 | 3.9e-26 | 2.4e-28 | 2.5e-29 | 6.1e-31 | < 0 |

Bisection with inertia of the nested blocks (`dh_onset.py`, 512 bits):
- Even, K = 320: onset in (1.7103, 1.7105].
- Even, K = 400: PD at 1.6841, 1.6973, 1.7090, 1.7100 and 1.7103; one negative pivot at 1.7105 for every K ≥ 266.
- Odd, K = 320 and 400: PD at 1.7105; two negative pivots at L = 2.0 for every K ≥ 80. So the odd onset lies in (1.7105, 2.0].

**Result: L*_DH,even = 1.7104 ± 0.0001.** This is NUMERICAL, from linear extrapolation of λ_min. The negativity at 1.7105 is CERTIFIED.
It is stable from K = 320 to K = 400. It is consistent with Queen's bracket (1.5, 2.0]. **It does not support Golden Eye's ~1.61**: at L = 1.6052, DH's K-sections are
PD with margin 3.5e-23 up to K = 400, and the witness is low-frequency, so K-convergence there is not in doubt.
Caveat: the PD sections are evidence, not proof. Certifying DH ≥ 0 at L < 1.71 would need the same tail control as ζ.

## 5. Is the ~3.5-day plan sized right?
For the finite-section factorization, yes, and it may even be oversized: 576-bit LDLᵀ at K ≈ 266 takes seconds. The cost is in
(i) Arb assembly (4 min at K = 266, 960 bits) and (ii) the Mori-route tail control, whose strip size (Np) I have not re-measured.
If the Mori T-matrix's norm-to-margin ratio is like the section's, p ≈ log2(‖T‖/λ) + 16 applies. The plan's precision budget should
be set from that rule, not from ~1100 bits.

## Files
- `cosritz_dump.py`: s62 assembly (Arb), matrix dump.
- `ldl_scan.py`: float and ball LDLᵀ scan.
- `ldl_certify.py`: rigorous congruence certificate and λ_min.
- `dh_onset.py`: bisection plus nested-K inertia.
- `dh_witness.py`: eigenvector mode content.
- `cosritz_px.py`: P_X variant (env GPRIMES).
- Runs are in `runs/cos`, `runs/onset` and `runs/px`. The dumps (`Wdump_*.json`, large) are not committed.
- Reproduce:
  - `FORM=ZETA python cosritz_dump.py 3421/2000 266 960 300 even`
  - then `python ldl_certify.py Wdump_ZETA_even_L3421_2000_K266_p960.json 200,266 "512,520,576"`
