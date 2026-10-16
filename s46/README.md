# S46 — H46a: planted off-line zero detection-support map on the S14/S44 exact assembly

RH OPEN. GRH OPEN. φ absent.

## Model

Same S14 exact Arb assembly as S36/S38/S44/S45 (Legendre·(1−u²)⁴ basis, g supported in [−L,L], F = g∗g̃ in [−2L,2L], even/odd sectors):
`W_ζ(F) = Pole + Arch − Σ_n 2Λ(n)/√n · F(log n) = Σ_ρ F̂(γ_ρ)`, F̂(z) = ∫F(x)e^{izx}dx.

**Planted fake:** remove two consecutive ζ zeros γ_a, γ_b (with their mirrors) and put in the quadruple ½ ± δ ± iγ:

`Δ(F) = 4∫F(x)cos(γx)cosh(δx)dx − 2∫F(x)[cos(γ_a x)+cos(γ_b x)]dx`, `Q_fake = W_ζ + Δ`.

γ = (γ_a+γ_b)/2 for pairs (1,2), (10,11), (50,51). **DH point:** γ = 85.699 (first Davenport–Heilbronn off-line ordinate) replacing the bracketing pair (γ₂₃ = 84.735, γ₂₄ = 87.425).

Δ entries are computed exactly: ∫F_jk φ = 2L²∫₀² G_jk(v)φ(Lv)dv, with moments I_m(a) = ∫₀² v^m e^{av}dv from the recursion I_m = (2^m e^{2a} − m I_{m−1})/a in Arb (prec 1536, rigorous balls). The zeros used come from mpmath zetazero at 70 digits, and every entry is inflated by 16L³·ε per removed zero with ε = 1e−55. The M24/M40 grid scans used the table zeros with ε = 1e−28.

## Checks

- **Convention CONFIRMED.** Summing 2∫F cos(γ_n x) over the first 20000 zeros reproduces the S14 W_ζ entries to <1e−16 relative (M=6, L=1.3 and M=8, L=0.6). With 100 zeros the gap is 3e−10, as expected.
- **Kernel signature PROVED + CONFIRMED.** 4cos(γ(t−s))cosh(δ(t−s)) = 4[a a + b b − c c − d d] with a = cos γt cosh δt and d = sin γt sinh δt (both even), b = sin γt cosh δt and c = cos γt sinh δt (both odd). So the planted quadruple is rank 2 with signature **(1,1) in each parity sector** and (2,2) in total (this matches GPPVerify `NoGo/QuartetSignature.lean` in the Pick picture). At δ=0 it is rank 1 and PSD per sector: moving zeros along the line never breaks positivity. Numerics (M=24, L=0.5 and 1.3) give exactly 1+/1− per sector, see `signature.json`.
- Near ĝ(γ)=0 the negative part is ≈ −4δ²|ĝ′(γ)|². Detection needs this to beat Σ_{others} 2|ĝ(γ_n)|².

## Onset table L*(γ,δ) — upper bounds on the true detection support

An Arb witness xᵀ(W_ζ+Δ)x < 0 is a genuine certificate, because a Ritz subspace can only overestimate L*. Non-detection is NUMERICAL and resolution-limited.

| config (γ) | δ=0.05 | 0.1 | 0.2 | 0.308 | 0.4 |
|---|---|---|---|---|---|
| pair 1–2 (17.578) | (0.40,0.45] C | (0.40,0.45] C | (0.35,0.40] C | (0.35,0.40] C | (0.30,0.35] C (M40) |
| pair 10–11 (51.372) | (0.9,1.0] C (M40) | (0.7,0.8] C | (0.7,0.8] C | (0.7,0.8] C | (0.7,0.8] C |
| DH pt (85.699) | none ≤1.3 (M72) | none ≤1.3 (M72) | none ≤1.3 (M72) | **(0.9,1.0] C (M72)** | **(0.9,1.0] C (M72)** |
| pair 50–51 (144.556) | none (M≤72, L≤1.3; M≤40 to 2.0) | none | none | none | none |

C = CERTIFIED. Pair 1–2 is stable from M24 to M40. Pair 10–11 at M24 was detected only at L=2.0 (δ≥0.308); at M40 it drops to 0.8; at M56 it is still PD at 0.5–0.7. The DH point is undetected at M24/40/56 through L=1.3 (M24/40 through 2.0), and at **M72** it is CERTIFIED indefinite at L=1.0 (even [−6.08e−30 ± 9.6e−39], odd [−7.75e−28 ± 1.8e−36]) and at L=1.3 (even [−2.63e−38 ± 3.1e−46]).

## Readings

1. **L* ~ c/δ: FALSIFIED (NUMERICAL, with CERTIFIED endpoints).** For pair 1–2 the onset moves only 0.45 → 0.35 while δ changes 8×. The δ-dependence is weak, roughly log(1/δ): the δ² negative part competes with λ_min(W_ζ), which decays exponentially in L.
2. **"Slow γ growth": not established.** The observed onsets grow with γ, but they are set by Ritz bandwidth (one needs M+8 ≳ γL/const): pair 10–11 dropped from 2.0 to 0.8 when M went 24 → 40, and the DH point dropped from "none ≤ 2.0" to 1.0 when M went 56 → 72. The true L*(γ) is OPEN.
3. **DH-point verdict: hypothesis "L* > 1.3" FALSIFIED.** A planted quadruple at the DH location (γ = 85.699, δ = 0.308) is CERTIFIED detected by ζ-assembly Weil positivity at support 2L = 2.0. So the S45 statement "DH stays PD at L=1.3" (M=24, NUMERICAL) is very likely a **basis-resolution artifact** and should be downgraded to "not detected at M=24". "Support 2.6 misses DH" is not established.
4. This calibrates the finite-support filter. It is not RH evidence: low off-line zeros are already excluded by zero verification. The certified margins are tiny (~1e−30) and follow λ_min(W_ζ), the same near-degeneracy as S44.

## Scripts

- `s46_planted.py` — Δ assembly (Arb moments), W_ζ (via `s44_scan.build_all`), Arb witness
- `s46_check.py M L N` — convention check against N zeros
- `s46_signature.py` — per-sector signature of the planted kernel
- `s46_scan.py M hs` — grid scan, 4 configs × 5 δ, with Arb witness at the first detection → `scan_M24.jsonl`, `scan_M40.jsonl`
- `s46_select.py M hs cfg γ|mean ia ib δs dps` — selected high-M runs → `sel_*.jsonl` (`sel_DH_M72_cert.jsonl` is the certified DH run with 70-digit zeros; `sel_DH_M72.jsonl` is superseded, run with 1e−28 zero radius)
- `summarize.py`, `s46_results.json`

## Novelty check

GPPVerify (pulled) has `NoGo/QuartetSignature.lean` (Pick kernel quartet signature (2,2)) and `RHSpectralMultiplicity.lean` (quadruple structure). It has no planted detection-support map. The CLAUDE_RESEARCH_NOTES heuristic is "new negative direction once L ≳ log(log γ)/δ" (Hodge-index section), and θ-control notes say β−½ ≲ 0.1 is unseen before L ~ 100. At low height our certified data show far earlier detection (γ≈51, δ=0.05: L ≤ 1.0). That is a tension in the finite-height near-null regime, not a refutation of the asymptotic heuristic.

Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
