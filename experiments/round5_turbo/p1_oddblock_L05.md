# P1 at L = 0.5: certified finite odd block (Turbo 6889360979 item 2)

Heart (Grok / GPPDiscovery5), Round 5. Script `p1_oddblock_L05.py` (arb ball arithmetic via python-flint 0.9.0, 160-bit), results `p1_oddblock_L05_N12_results.json` (and `_N4`), independent cross-check `crosscheck_kernel.py` / `crosscheck_results.json`. No zero data.

## Setup (CHECKED by hand)
- Odd-face form B(f,g) = ½∫∫_{[0,ℓ]²}[K_S(s−u) − K_S(s+u)] f(s) g(u), ℓ = 1/4, K_S = pole − prime_{S_L} + arch exactly as in-tree.
- E_N = span{φ_j = √(2/ℓ) sin(ω_j x), ω_j = jπ/ℓ = 4πj, j ≤ N}: L²-orthonormal on [0, ℓ]; these are the first N odd modes on [−ℓ, ℓ].
- Exact reduction: B(φ_j, φ_k) = (2/ℓ){(1/π)∫₀^∞ Ψ² V S_j S_k dt − 2Ψ(i/2)² σ_j σ_k}, with
  S_j(t) = ∫₀^ℓ sin(ω_j x) sin(tx) dx = ω_j ℓ sinc((t − ω_j)ℓ)/(t + ω_j) (entire on the integration range),
  σ_j = ∫₀^ℓ sinh(x/2) sin(ω_j x) dx = −(−1)^j ω_j sinh(ℓ/2)/(1/4 + ω_j²) (the pole term is rank one and negative),
  V = Re ψ(1/4 + it/2) − log π − 2P, P(t) = Σ_{n ∈ {2,3,4,5,7,8,9,11}} Λ(n) n^{−1/2} cos(t log n). This finite set is exact: on |d| ≤ L = 0.5 only n < e^{2.5} = 12.18 have triG(log n ± d) ≠ 0.
- Integral on [0, 800]: acb.integral, with rigorous error enclosure, in 16 pieces of length 50.
- Tail on [800, ∞) (CHECKED by hand): Ψ² ≤ 16/t⁴. For t ≥ 100, |w(t)| ≤ log t + 2.5, from ψ(z) = log z − 1/(2z) − ∫₀^∞ f(u)e^{−zu}du with 0 ≤ f(u) ≤ u/12, so the integral term is ≤ 4/3 at Re z = 1/4. |P| ≤ P(0) = 4.2605. |S_j| ≤ (4/3)ω_j/t² for t ≥ 2ω_max (800 ≥ 2·48π). So |tail_jk| ≤ 1.10e-19·ω_jω_k ≤ 2.5e-15, and it is added to each entry as a ball radius.
- Certificate: LDLᵀ of M − βI carried out in ball arithmetic. Every pivot has a positive lower bound, so M − βI is positive definite.

## Result (NUMERICAL, interval-certified; not Lean)
- **N = 12: B(e, e) ≥ β‖e‖²_{L²} on E_12 with β = 1.66823e-8** (= 0.9 × float λ_min, where λ_min = 1.85359e-8). The LDLᵀ pivots of M − βI have lower bounds from 2.21e-5 down to 1.96e-9, all positive. The largest entry radius is 2.0e-14.
- The float spectrum of M_12 runs from 1.854e-8 to 2.771e-5. N = 4 gives λ_min = 1.094e-6 (certified β = 9.85e-7).
- Cross-check (independent code path): prime and pole were computed directly in the kernel domain with the in-tree triG kernel (2000² Gauss–Legendre), and arch by scipy quad with S_j. The entries (1,1), (1,2), (3,4), (6,7), (12,12) agree with the arb block to relative 2e-13 … 1e-8. The worst case is (12,12), set by the float quadrature of the kinked prime kernel.
- Reading: λ_min(E_N) decreases with N (1.09e-6 at N = 4, 1.85e-8 at N = 12), as it must for a compact form. So β is tied to N, and P2/P3 must be closed at the same N = 12. The block is nowhere near degenerate: the smallest pivot is 1.96e-9 after removing 0.9·λ_min.

## What remains
- P2 (high-mode coercivity on E_12^⊥ above the window crossover T_* = 113.19, where ω_13 = 52π = 163.4 > T_*) and P3 (cross term |B(e, r)|² ≤ (1 − η)β‖e‖²·κE(r)), in that order, per Turbo's order. The Schur step then needs κ and the cross-term bound on the same normalisation as β (L² on [0, ℓ]).
- Note for P2/P3: modes j ≥ 13 have frequency above T_*, so the window symbol is positive there (OF2). The remaining work is the remainder symbol and the coupling between E_12 and its complement.

Labels: P1 at L = 0.5, N = 12: NUMERICAL (interval-certified, β = 1.66823e-8), not Lean. Reduction formulas and tail bound: CHECKED. P2, P3: OPEN. OddFace(0.5): OPEN. H_ψ, WeilIsGram, ARP, RH: OPEN.
