# Round 7 plan (Turbo 6889791941): GAP (a) then (b), OddFace at L = 0.5, N = 12

Heart (Grok / GPPDiscovery5). Written before any computation.

## Ask
(a) Fix K ≈ 60, Y = E_{12+K}. Arb-certify the Schur complement onto E_12 of the bordered Gram on Y ⊕ {δ_ℓ} is ⪰ β₁ I (L² on [0,1/4]). Report β₁ and the pivot certificate.
(b) Bound the tail j > 12+K with the explicit ω⁻³ terms of R6.1 and show the extra Schur loss is < β₁. Float prediction: ≲ 1e-10 against ~1.3e-8.

## Plan (a)
- Matrix: basis φ_1..φ_72, δ_ℓ (73×73). Entries M_jk by acb.integral on [0, T_jk] with T_jk ≥ 2 max(ω_j, ω_k), plus the Round-5 tail lemma (Ψ² ≤ 16/t⁴, |V| ≤ log t + C, |S_j| ≤ (4/3)ω_j/t²). δ column: Round-6 tail lemma, T = 4000. g, u_i reuse the Round-6 arb balls only if recomputed at the same precision; recompute to be safe.
- Certificate: ball LDLᵀ of the full 73×73 matrix minus diag(β₁ I_12, 0), eliminating Y ⊕ δ first (pivots > 0 certify D ≻ 0, i.e. the corrected P2 on Y ⊕ δ), then E_12 (pivots > 0 certify Schur − β₁ ≻ 0).
- Parallel workers (6 processes), each writes its own JSON shard; assembly script checks every entry present.
- Float run first (same reduction) to predict β₁ and check conditioning of the pivots against the entry radii.

## Plan (b) (first pass this run; full proof may remain GAP)
- Float: extra Schur loss from modes 73..12+K' for K' up to 400 against β₁, with and without δ in the bordered space.
- Structural test: R6.5 showed M22 − g wwᵀ indefinite because of the second boundary term −∂²_sG(ℓ,·)/ω³. Test whether adding the second boundary functional (∂²_s G(ℓ,·), transform −t² sin(tℓ) direction, coefficient w_j/ω_j²) to the bordered space makes the remaining high block coercive at the window-symbol rate 3w(ω)/ω⁴. If yes, the tail bound becomes: two boundary functionals (finite Gram, certified) + interior form with explicit ω⁻⁴ log ω coercivity + ω⁻³ remainders, which is summable.
- Explicit constants for R6.1: sup|∂²_sG(ℓ,·)|, ∫|∂³_sG| (arch part log-integrable, prime kink jumps a_k, pole smooth).

## Backup
- If ball LDLᵀ at 73 dims loses certification to radius growth: raise prec, tighten tails (larger T), or certify with K = 40.
- If (b)'s remainder is not coercive even with the second boundary functional: report KILLED for that route and name the next form.
