# S25 — ζ odd sector at h = 12/25 (jump-at-h + Parseval-tail Arb)

Daniel Toupin / Golden Physics Project, with AI assistance from Grok. **RH OPEN.**

## Results
- **CERTIFIED (Arb 128-bit, Lean-open).** Odd g ∈ C_c((−12/25, 12/25)), real or complex:
  - H = 11, N = 479: Q_ζ(g) ≥ 0.000109‖g‖² (c = 10/2307, 2β ≤ 0.97483). `s25_cert_h048_N480_H11_j.json`, `v_N480.json`.
  - H = 10, N = 359: Q ≥ 0.000106‖g‖² (2β ≤ 0.97550). `s25_cert_h048_N360_H10_j.json`, `v_N360.json`.
- Verifier `s25_verify_arb.py`: exact knot-sum S (jump at h allowed), cubic cell remainder (M₃), Parseval L² tail.
- Cross-check (S24 verifier, uniform disc): h = 19/40, N = 240, H = 10, Q ≥ 0.0000373‖g‖² (`../s24/s25_cert_h0475_N240_H10_disc.json`).

## Dead ends
- Continuous graded PL at h = 0.48 (Arb K₂ blow-up).
- Uniform disc at h = 0.48, H = 10, N ≤ 280: Arb 2β ∈ (1.003, 1.009).

## Next
Even sector: LP resolution at h = 1/2, h*; certify B_even|mean0; then χ_e > 1 (H24).

## Addendum (S25b) — h = 49/100 and verifier cross-checks
- **CERTIFIED (Arb 128-bit, Lean-open).** Odd g ∈ C_c((−49/100, 49/100)), real or complex: Q_ζ(g) ≥ **0.0000311**‖g‖²
  (c = 13/5387, 2β ≤ 0.98709; H = 11, N = 479 jump design). Files `s25_cert_h049_N480_H11_j_c08.json`, `v_h049_N480_c08.json`.
  Same extension with c = 3/1847: 0.0000210 (`s25_cert_h049_N480_H11_j.json`, `v_h049_N480.json`).
  Claim-1 LP symbol: `s22_cert_h0.49.json` (S22 LP, c_LP = 0.003249; `lp_save_049.txt`).
- **Cross-check (s25 verifier on S24 certs):** h = 47/100 N200 H10: 2β ≤ 0.96358 (S24 verifier: 0.97277), coercivity 0.000294;
  h = 23/50 N160: 2β ≤ 0.94307 (S24: 0.94874), 0.000511. Consistent and tighter. `v_s24recheck.json`.
- **NUMERICAL (Galerkin M = 64, `out_susc_M64.json`, `susc_049.txt`):** χ(0.485) = 0.97440, χ(0.49) = 0.98053.
  Arb 2β vs Galerkin χ: 0.97482 vs 0.96667 (h = 12/25), 0.98709 vs 0.98053 (h = 49/100).
- **PROVED (elementary, duality):** for any admissible m > 0 (m − σ_A invisible on odd window g), inf over odd extensions v_w of
  β = (1/2π)∫|v̂_w|²/m equals sup_g |⟨g,v⟩|²/B(g) = χ_true/2 (Riesz representer of g ↦ ⟨g,v⟩ in L²(m dξ/2π)), independent of m.
  So the extension step is lossless: odd coercivity is certifiable by this route exactly where χ_true < 1 and claim 1 gives c > 0.
  At h = 1/2 (χ_Gal = 0.98898) the only obstruction is claim 1 (S22 LP c = −0.0002), not β.
