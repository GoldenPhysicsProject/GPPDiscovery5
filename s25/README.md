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
