# S24 — ζ odd sector at h = 23/50 and 47/100; Krein ceiling for the claim-1 LP; S23 verifier correction

Daniel Toupin / Golden Physics Project, with AI assistance from Grok. **RH OPEN.**

## Results
- **CERTIFIED (Arb 128-bit, Lean-open).** Odd g ∈ C_c((−h,h)), real or complex:
  - h = 23/50: Q_ζ(g) ≥ 0.000460‖g‖² (c = 46/5125, PL extension on [h,6], N = 160, 2β ≤ 0.94874). File `s24_cert_h046_N160.json`.
  - h = 47/100: Q_ζ(g) ≥ 0.000220‖g‖² (c = 67/8293, PL extension on [h,10], N = 200, 2β ≤ 0.97277). File `s24_cert_h047_N200_H10.json`.
  - Also verified: h=0.46 N120 (0.000387), h=0.47 N160 H10 (0.000147), h=0.47 N200 H8 (0.000205). Outputs `v3k_*.txt`, `verify_h046_N120.txt`.
- **S23 correction.** `s24_verify_arb.py` fixes two S23 verifier errors: the dyadic tail used (S/ξ)² where |S|² ≤ (J+K2/ξ)²/ξ² was meant
  (an extra 1/X² factor, i.e. tail undercounted), and the 1/ξ² coefficient was rederived (K2 = two integrations by parts with
  slope jumps). Re-run of S23's own cert (`recheck_s23_3000.txt`): still verified, 2β ≤ 0.91923, coercivity 0.000725 (S23 said 0.000737).
- **PROVED (elementary).** Krein ceiling for claim 1: any admissible symbol m (m − σ_A in span{cos uξ : u ≥ 2h} + {Ŵ : W ≡ 1 on (−2h,2h)})
  with m ≥ c obeys c ≤ λ_min(B restricted to mean-zero g ∈ L²(−h,h), BOTH parities). Proof: for g supported in (−h,h), |ĝ|² has
  inverse FT supported in [−2h,2h], so ∫|ĝ|² cos uξ = 0 and ∫|ĝ|²Ŵ = 2π(∫g)² = 0 on mean-zero g, hence B(g) = (1/2π)∫|ĝ|² m ≥ c‖g‖².
- **NUMERICAL (Galerkin upper bounds, `s24_krein_ceiling.py`).** λ_min(B_even | ∫g=0) = 0.0140 (h=1/2, M=56), 0.00099 (h*, M=56),
  decreasing with M; λ_min(B_odd) ≈ 0.0198 (h=1/2), 0.0163 (h*). The binding sector switches from odd to EVEN mean-zero near h ≈ 0.49.
  The current LP (Bernstein cap 3, one DC trapezoid) gives c = −0.0002 at h = 1/2, well below the 0.014 ceiling, so h = 1/2 is an LP
  resolution problem, not a structural barrier.
- **NUMERICAL (`s24_even_perp.py`, M=40).** Mean-zero and ⟨g,cosh(x/2)⟩=0 are almost the same constraint (coefficient cosine 0.9999993 at h*).
  λ_min(B_even | cosh-perp) = 0.253, 0.0157, 0.00128 at h = 0.4, 0.5, h*, while λ_min(Q_even) = 2.3e−4, 1.3e−6, 9.0e−8.
## Dead ends
- Long tails H = 25/2 with N = 359 knots (`s24_cert_h047.json`, `s24_cert_h048.json`): M2 and β blow up (2β ≈ 16, 10). Keep H ≤ 10.
- Raising the LP Bernstein cap to 6 with Ξ from the S22 grid: Ξ search runs past 5000 (needs ≈1.2e4); not a fix as coded.
