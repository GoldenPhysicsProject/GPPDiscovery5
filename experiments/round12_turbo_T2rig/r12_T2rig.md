# Round 12 (Turbo 6890743146): (T2-rig) closed. The explicit constant for |ω⁴ρ_j|

Heart (Grok / GPPDiscovery5). L = 0.5 (ℓ = 1/4), ω_j = 4jπ, same normalisation as Rounds 5–10 (V = Re ψ(1/4+it/2) − log π − 2Σ_{n≤12} Λ(n)n^{-1/2}cos(t log n)). No zero data.

## R12.1 Exact reduction (CHECKED: algebra, plus numerical agreement with Round 8 to ≤1.4e-6 abs)
Write P(t) = sin⁴(t/2) sin²(t/4) = Σ_{k=0}^{5} c_k cos(kt/2) with c = (3/16, −1/16, −1/4, 3/32, 1/16, −1/32) (exact). Define, for ω = 4mπ (m ≥ 1),
  Q(ω) := ω ∫_0^∞ V(t) P(t)/(t² − ω²) dt,   A(t) := Re ψ(1/4 + it/2).
Since F_δ S_j = √(2/ℓ)ω_j(−1)^j sin²(tℓ)/(t²−ω_j²), the Round 8 identity gives
  **ω_j⁴ρ_j(δ) = (16/π)√8 (−1)^j Q(ω_j) − (bI2²/8) sinh(1/8) sig(j)**, and by partial fractions
  **ω_j⁴ρ_j(φ_i) = (128/π)(−1)^{i+j} ω_iω_j/(ω_j²−ω_i²) · [Q(ω_j)/ω_j − Q(ω_i)/ω_i] − (bI2²/8) sig(i) sig(j)**.
Evaluating Q uses PV∫_0^∞ cos(at)/(t²−ω²)dt = −π sin(|a|ω)/(2ω), a ≠ 0, which is 0 for a = 0. At ω = 4mπ every sin(kω/2) vanishes, so the constant parts A(ω) and −log π drop out exactly. Each prime cosine gives sin((k/2 ± log n)ω) = ±sin(ω log n). The result:
  **Q(ω) = R(ω) + Π(ω),  R(ω) = ω∫_0^∞ (A(t) − A(ω)) P(t)/(t² − ω²) dt,  Π(ω) = π Σ_{n≤12} Λ(n)n^{-1/2} C_n sin(ω log n),**
with C_n = Σ_{k<2 log n} c_k: C_2 = 1/8, C_3 = C_4 = −1/8, C_5 = C_7 = −1/32, C_8 = C_9 = C_11 = 1/32.
Check (`t2rig.py`, `phicheck.py`): δ row at j = 73, 123, 173, 273, 472 gives −5.3506638, −3.1017664, −2.8002350, −12.9563027, 8.0557758, against Round 8's −5.3506638, −3.1017667, −2.8002357, −12.9563052, 8.0557870. All twelve φ rows at j = 73 agree to 9 digits. R(ω) ≈ 0.4626 at every j tested, consistent with the limit (3/16)(π²/4) = 0.46264.

## R12.2 Bounds (PROVED; the constants are Arb-evaluated in `t2rig_bound.py`)
- **Lemma A (concavity).** B(s) := A(√s) is concave on s > 0, because each series term −(n+¼)/((n+¼)² + s/4) is concave in s. So G(t) := ω(A(t) − A(ω))/(t² − ω²) (a divided difference of B in s = t²) is positive and decreasing on [0, ∞), with G(0) = (A(ω) − A(0))/ω. Float check: max second difference of B on a 2e6-point grid is −1.2e-13 ≤ 0.
- **Lemma B (derivative).** A'(t) = Σ_{n≥0} f(n+¼) with f(x) = bx/(x²+b²)², b = t/2. f is unimodal, ∫_0^∞ f = 1/t and max f = 9/(4√3 t²). Hence |A'(t) − 1/t| ≤ κ/t² with κ = (5/4)·9/(4√3) = 1.62380. Float sup of t²|A' − 1/t| is 0.654. Consequences: A(ω) ≤ log(ω/2) + κ/ω, and |e(t) − e(ω)| ≤ κ|1/t − 1/ω| for e = A − log.
- **R bound.** P = 3/16 + Σ_{k≥1} c_k cos(kt/2). Integrate the oscillatory part by parts against the monotone G: |∫G cos(kt/2)| ≤ (2/k)G(0). For the mean part: ∫_0^1 G ≤ G(0); on [1, ∞) the log part is ≤ ∫_0^∞ ω log(t/ω)/(t²−ω²) = π²/4, and the e part is ≤ (κ/ω)log(1+ω). So
  **0 ≤ R(ω) ≤ (3/16)[π²/4 + G(0) + κ log(1+ω)/ω] + 0.48125·G(0)**, which decreases in ω.
  This gives 0.507517 at ω_13 and 0.855473 at ω_1 = 4π.
- **Π bound.** |Π| ≤ π Σ Λ(n)n^{-1/2}|C_n| = 0.851515.
- **(T2-rig), δ row.** For all j ≥ 13: **|ω_j⁴ρ_j(δ)| ≤ C_δ = 19.5770**. The σ term is ≤ 7e-6. In the requested shape this is **|ω⁴ρ_j| ≤ 3.84166 log ω_j** for j ≥ 13. Float sweep max: 16.09, so C_δ is 1.22× above it.
- **(T2-rig), φ rows.** With |Q| ≤ 1.70699 for every ω_m (m ≥ 1): **|ω_j⁴ρ_j(φ_i)| ≤ 69.5490 (ω_i + ω_j)/(ω_j² − ω_i²) + (σ term ≤ 4e-6) ≤ 69.549/(4π(j − 12))**. That is 5.53 at j = 13 and 0.0907 at j = 73 (float ≤ 0.036).

## R12.3 Correction to Round 8
At L = 0.5, ρ_j is O(ω⁻⁴), not O(ω⁻⁴ log ω). The Round 8 "log ω from the arch log|d| endpoint singularity" reading is wrong at this width. The log term is A(ω)·(−π/2)Σ_a c_a sin(aω_j), which vanishes identically because every frequency of P (multiples of 1/2 = 2ℓ) times ω_j = jπ/ℓ lies in 2πℤ. For a general width (ℓ with 1/ℓ ∉ 2ℤ), the same algebra gives Q_ℓ(ω_j) = R_ℓ + Π_ℓ − (π/2)A(ω_j)Σ_a c_a sin(aω_j), and the log ω returns (CHECKED algebra; not needed here).

## For consumers
- Flint (H-Env): use C = C_δ = 19.5770 (δ row) and the φ-row bound above. Better still, use the exact form ω⁴ρ_j = 14.40506(−1)^j(R + Π) with R ∈ [0, 0.5076] and Π explicit.
- My (T3): the Round 8 envelope R_δ = 2.1 log ω is replaced by the constant 19.577 (≤ 2.1 log ω only for j ≳ 1000). R_φ = 0.045·ω_73/ω is replaced by 69.549/(ω_j − ω_12). Re-running the T3 float with these is the next step for items 3′/4.

Labels: R12.1 CHECKED. R12.2 PROVED (elementary analysis, with constants evaluated in Arb). The lemma float checks are NUMERICAL sanity checks only. (T2-rig): CLOSED at L = 0.5. Items 3′/4 interval assembly: GAP. OddFace(0.5) proof-grade: OPEN (pending (H-A), (H-ρ), and 3′/4). RH: OPEN.
Credit: Heart (Grok / GPPDiscovery5); audit chain Turbo (6890743146), Flint run 13; Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
