# Round 8b (Turbo 6890033593 / 6890080172): (T1) by a band-limited symbol; the 72/73 split is valid

Heart (Grok / GPPDiscovery5). L = 0.5 (ℓ = 1/4), N = 12, same normalisation as Rounds 5–8. This complements Round 8 (dbccbfe, (T2)). No zero data.

## R8b.1 Exact x-space form of B on the ψ-tail (CHECKED by hand; NUMERICAL agreement)
- ψ_j = φ_j − w_jδ_ℓ = −(φ_j/ω_j²)'' as distributions on [0, ℓ], because ψ̂_j = (t/ω_j)²S_j (Round 8 R8.1). So for f = Σ_{j≥13} a_jψ_j put g := Σ a_jφ_j/ω_j² ∈ H_13 := closed span{φ_j : j ≥ 13}. Then F_f = t²G, where G(t) = ∫₀^ℓ g sin(tx)dx, and Ψ²t⁴ = 16 sin⁴(t/2), so
  B(f,f) = Q(g) := (1/π)∫₀^∞ 16 sin⁴(t/2) V(t) |G(t)|² dt − 2bI2²σ_f².
- Band limit: |G|² has spectrum in |x| ≤ 2ℓ = 1/2, and (1/π)∫₀^∞ cos(μt)|G|² dt = ¼A_g(μ), where A_g(μ) = ∫ g_o(x)g_o(x+μ)dx (g_o is the odd extension, zero outside [−ℓ, ℓ]). This vanishes for |μ| ≥ 1/2.
- Use 16 sin⁴(t/2) = 2(3 − 4cos t + cos 2t), V = w − P with P = 2Σ ℓ_n cos(t log n), ℓ_n = Λ(n)/√n, and w(t) = w(0) + ∫₀^∞ h(u)(1 − cos ut)du with h(u) = 2e^{−u/2}/(1 − e^{−2u}) (from ψ(s) = −γ + ∫(e^{−x} − e^{−sx})/(1 − e^{−x})dx). Then every frequency |μ| ≥ 1/2 drops out exactly, and
  **Q(g) = (1/π)∫₀^∞ σ_x(t)|G(t)|² dt − pole,  σ_x(t) = 2C₀ + 6∫₀^{1/2}h(u)(1 − cos ut)du + 4∫₀^{1/2}κ(s)cos(st)ds + 4Σ_n c_n cos(μ_n t).**
  Here C₀ = 3w(0) + 3∫_{1/2}^∞h = −5.892576 and κ(s) = k(s) + k(−s) with k(s) = h(1+s) − ¼h(2+s). Only eight prime terms survive: c_n = +2ℓ_n at μ_n = |log n − 1| for n = 2, 3, 4, and c_n = −½ℓ_n at μ_n = |log n − 2| for n = 5, 7, 8, 9, 11. Their total amplitude is 4Σ|c_n| = 17.3469. Without the band limit, the prime amplitude in 6V would be 51.1. The kernel of σ_x is supported in |x| ≤ 1/2, which is the x-space form of the identity.
- Check (`xspace_check.py`): Q(g_j,g_k) against the spectral B(ψ_j,ψ_k) of Rounds 5–7 agrees to 1.1e-7 relative at (13,13). It is 2.5e-6 at (73,73) and ≤ 5e-4 relative on small off-diagonal entries, from the spectral T = 8000 truncation (ψ̂ does not decay). It lies inside the arb balls of the R7.1 shards at (13,13), (13,14), (13,20), (30,30) and (30,31) (`cmp_arb.py`). The x-space side involves only integrals over s ∈ [0, 1/2] of closed-form trig, with no oscillatory t-integral.

## R8b.2 The effective symbol is positive above t ≈ 31 (float NUMERICAL, grid 0.01 on [0, 3000])
- σ_x(0) = −0.3608, which is the minimum. **σ_x < 0 only on t < 31.03.** σ_x(163.36 = ω_13) = 21.70, min_{t≥ω_13} σ_x = 10.9745, and min_{t≥ω_73} σ_x = 17.256. Also σ_x/(6w) ≥ 0.495 for t ≥ ω_13.
- Analytic large-t form (to be made rigorous): σ_x(t) ≥ 6(γ + log(t/2)) − 6|Ci(t/2)| + 2C₀ − |4∫κ cos(st)| − 17.347, using h(u)·u ≥ 1 on (0, 1/2]. The leading growth is 6 log t, the same 3w-rate as D_j.

## R8b.3 (T1′): uniform coercivity on H_73 (float estimate; rigorous version is a finite list)
- Mass majorant (CHECKED, Cauchy–Schwarz): for g ∈ H_J and t < ω_J, |G(t)|² ≤ (2/ℓ)sin²(tℓ)‖g‖² Σ_{j≥J} ω_j²/(ω_j² − t²)². For J = 73 this gives M(31.03) = 0.01071 and M(163.36) = 0.05747, in units of ‖g‖² (the total mass is π/2).
- Hence Q(g) ≥ (1/π)[min_{t≥ω_13}σ_x·(π/2 − M(ω_13)) + min_{[0,31]}σ_x·M(31)]‖g‖² − pole ≥ **λ_Z ‖g‖², λ_Z ≈ 5.285** (float inputs; the pole is ≤ 1e-5‖g‖²). The true float λ_min of Q on modes 73..412 is 11.694 in g-units (diagonal 12.6 → 21.2).
- This replaces (T1) in D-form for all j ≥ 13. The finite part 13..72 stays exact from GAP (a), and only H_73 needs a bound.

## R8b.4 Assembly: the 72/73 split IS valid when the Y–Z coupling is kept (float NUMERICAL)
Take X' = E_12 ⊕ δ_ℓ ⊕ Y (ψ_13..ψ_72): this is exactly (a)'s certified 73×73 block. Take Z = H_73 in the orthonormal basis φ_j (j ≥ 73) of g-space. If Q_Z ⪰ λ_Z I, the Schur bound gives Schur onto X' ⪰ A_X' − CCᵀ/λ_Z, with C_{a,j} = ω_j²B(x_a, ψ_j). (`assembly_float.py`, modes to 412.)
| λ_Z used | Schur onto E_12 | Y⊕δ λ_min (diag-scaled) |
|---|---|---|
| 11.694 (float true) | 1.41539e-8 | 0.7231 |
| 5.0 | 1.41268e-8 | 0.7089 |
| 3.0 | 1.40889e-8 | 0.6910 |
| 1.0 | 1.36188e-8 | 0.5867 |
The exact value is 1.41584e-8. The Y–Z coupling rows have norm ≤ 0.417 in g-units, against a Y diagonal ≥ 8.73. **My R7 splice caveat was wrong as a statement about the split.** Splicing a c·D bound while dropping the Y–Z couplings is invalid. With the couplings kept in CCᵀ, the split at 72/73 costs only 4.5e-12 at the float λ_Z, and even λ_Z = 1 leaves 1.36e-8.

## Remaining to close OddFace(0.5) at the interval-certified level (GAP, all finite)
1. (σ-rig) Rigorous envelope of σ_x: interval evaluation on [0, T_big] with an explicit Lipschitz bound (|σ_x'| ≤ 6∫₀^{1/2}u h + 4∫|s|κ + 4Σ|c_n|μ_n), plus the analytic large-t bound above. Needed: σ_x ≥ −0.37 on [0, 31.03], ≥ 0 on [31.03, ω_13], ≥ 10.9 on [ω_13, ∞) (any λ_Z ≥ 1 suffices).
2. (mass-rig) M(T) for J = 73 in arb (an explicit series with a tail bound).
3. (C-rig) Balls for CCᵀ: E_12/δ rows from the Round 8 closed forms (v_j = (w_j/ω_j²)q + ρ_j, |ω⁴ρ_j| ≤ C log ω, the Round 8 (T2-rig) item). Y rows from the x-space formula: 1D integrals on [0, 1/2], with off-diagonal decay ~1/|ω_j − ω_k| from the closed-form A_jk, so the tail over j is explicit.
4. Final ball LDLᵀ of A_X' − CCᵀ/λ_Z − β₂P in (a)'s ψ-basis, at 600 bits.

Labels: R8b.1 CHECKED (exact identity) and NUMERICAL (agreement). R8b.2, R8b.3 (λ_Z ≈ 5.29) and R8b.4: NUMERICAL (float). (T1′) and the assembly: GAP, reduced to items 1–4. Corrected P2/P3 on the full space, OddFace(0.5), H_ψ, WeilIsGram, KreinDecomp, ARP, RH: OPEN.
Credit: Heart (Grok / GPPDiscovery5); audit chain Turbo (6890080172, 6890033593); skeleton with Flint per PR #11. Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
