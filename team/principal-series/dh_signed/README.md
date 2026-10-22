# dh_signed: signed DH onset eigenvectors (even sector) for ⟨f₀, χ_L⟩ (New Bot; NUMERICAL; RH OPEN)
Daniel Toupin / Golden Physics Project, with AI assistance from Grok.

## Conventions
- **L is the half-width.** f lives on [−L, L].
- **Basis.** f(x) = Σ_{m=0}^{K−1} v_m cos(a_m x), with a_m = (2m+1)π/(2L). Gram = L·I, i.e. ∫_{−L}^{L} cos(a_m x)cos(a_n x)dx = L δ_mn.
- **Form.** The DH window Weil form of `s75/dhritz.py` (FORM=DH, S47 difference trick, no pole, S14 normalisation), even sector. The Arb assembly uses prec 320, and the Ritz step uses mpmath dps 60.
- **Stored vector (raw, NOT unit-normalised).** The `vec` entries are exact rationals: the mpmath Ritz vector scaled so that max|v_m| = 1, then rounded to 50 decimal digits. Hence ‖f‖²_{L²[−L,L]} = L·Σ v_m². Values are in `INDEX.json` (≈ 2.1416 at L = 1.7105). The overall sign is arbitrary, since eigsy's sign is kept as is. Normalise with f₀ = f/‖f‖.
- **RQ.** The field `RQ` is the Arb ball for vᵀWv/(L‖v‖²) on the stored rational vector. That is a rigorous Ritz upper bound on the even-sector inf, in the same normalisation. It is CERTIFIED negative at L = 3421/2000 for K = 300 and 360.

## Files
- `cosvec_DH_even_L3421_2000_K{60,120,180,240,300,360}.json` (L = 1.7105). K = 240–360 are the S75 witnesses; K = 60–180 are New Bot S2 runs of the unmodified script.
- Neighbours: `L6841_4000_K300` (1.71025), `L171_100_K360` (1.71), `L1711_1000_K240` (1.711).
- **L = 1.7104 (= 1069/625):** `cosvec_DH_even_L1069_625_K60.json` and `_K120.json` come from `dhritz_resid.py` and carry `resid_rel` and `ritz2`. K = 240 is still running; if it finishes it will land here, otherwise see PARKED.md.
- **L = 1.7105, residual reruns:** `cosvec_DH_even_L3421_2000_K{60,120}_resid.json`. These are the same vectors as the K60/K120 files to all 50 digits, and they carry the residual.
- **Residuals** (‖Wv/L − RQ·v‖/‖v‖, Arb):

  | L | K | residual | RQ | ritz2 |
  |---|---|---|---|---|
  | 1.7104 | 60 | 7.8e−50 | +2.133e−29 | 6.06e−22 |
  | 1.7104 | 120 | 1.3e−49 | +4.437e−30 | 5.30e−22 |
  | 1.7105 | 60 | 7.2e−50 | +2.066e−29 | 6.00e−22 |
  | 1.7105 | 120 | 1.4e−49 | +4.067e−30 | 5.26e−22 |

  The residual is set by rounding the stored vector to 50 digits. At these K, both L are still positive Ritz values: the certified negativity at 1.7105 needs K ≥ 300 (S75).
- **Residual.** `dhritz_resid.py` is dhritz.py plus one block. It computes ‖Wv/L − RQ·v‖/‖v‖ in Arb on the stored vector (`resid_rel`) and the second Ritz value (`ritz2`). The S75 witnesses predate this block, so their residuals are not stored. For them, the K-convergence ‖v_300 − v_360‖ = 4.2e−7 (unit vectors) and ritz2 ≈ 4.9e−22 ≫ |RQ| ≈ 1e−31 are the available accuracy data.
- `overlap.py` (float64 grid) and `ovcheck.py` (mpmath dps 40, closed form): the 1.4e−14 overlap. Outputs are in `overlap_s75.jsonl` and `overlap_moderateK.jsonl`.

## The 1.4e−14 overlap: exact definition
- φ_L(x) = e^{−iγ₀x} sinh(β₀x) on [−L, L], with γ₀ = 85.699348 and β₀ = 0.308517 (Queen's `dh_sinh.py`; DH zero 0.808517 + 85.699348i).
- φ_L = cos(γ₀x)sinh(β₀x) − i sin(γ₀x)sinh(β₀x). The real part is **odd**, so it is exactly orthogonal to the even f₀. The even part is −i·q_e, with q_e(x) = sin(γ₀x)sinh(β₀x).
- **Reported number:** |⟨f, q_e⟩| / (‖f‖·‖q_e‖) = 1.382e−14 at L = 1.7105 (K = 240–360; K = 120 gives 1.371e−14). This is the **normalised** overlap with the even part.
- Norms at L = 1.7105: ‖q_e‖² = 0.169380736, ‖Re φ_L‖² = 0.166350975, and ‖φ_L‖² = sinh(2β₀L)/(2β₀) − L = 0.335731712. At L = 1.7104 these are 0.169333732, 0.166336921 and 0.335670652.
- The overlap with the full complex φ_L, |⟨f, φ_L⟩|/(‖f‖‖φ_L‖), is 9.8e−15 = 1.382e−14·√(0.16938/0.33573).
- **Raw** (unnormalised) ⟨f₀, q_e⟩ with ‖f₀‖ = 1 is 1.382e−14 × √0.16938 = 5.69e−15.
- Mechanism: f₀ is prolate-like, with f′(±L)/‖f‖ ≈ 2e−11, so its Fourier content at γ₀ is tiny.
- Queen reports ≤ 4.5e−4 in her own measurement. If that is a bound or a different normalisation or probe (e.g. her φ_L ground state rather than the sinh pair vector), the two are compatible. Reconcile the definitions before comparing.
