# Round 20: the Shadow-Euler circuit (Stieltjes ladder, Pick matrices at s_{k,N}, primes inside)

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).

**RH is OPEN; nothing here claims it.** Labels: PROVED / CHECKED / NUMERICAL / CONJECTURE / HEURISTIC / GAP / OPEN; "cited" means a published theorem that is stated here, not re-derived. **No zero is used as an input.** Every number below is computed from the Taylor coefficients of ξ_M at ½, which come from the zero-free theta laws Φ_M of Round 17 (`round17_bpy_delta/laws17.py`) by trapezoid moments. Zero ordinates appear only as illustrations, and are marked as such.

Source: Daniel's *The Shadow Euler Identity* (September 2026 revision). Its notation:
- X(z) = ξ(½+z)/ξ(½) = F(z²), u = (s−½)².
- m(u) = F′/F = (1/(2√u)) (ξ′/ξ)(½+√u).
- Thm "Stieltjes-resolvent": RH ⇔ m is Stieltjes.
- Moments t_m = Σ γ^{−2m} = Tr A^m, with the Hankel families t_{i+j+1}, t_{i+j+2}.
- s_{k,N} = kN/(k+N) is a parallel sum, and for fixed k the points u_{k,N} = (s_{k,N}−½)² form a uniqueness set.

Here c_j := t_{j+1}, so m(u) = Σ_j (−1)^j c_j u^j.

**Models.**
- ZETA.
- EPS: Epstein x²+5y², ξ_E = s(s−1)(√20/2π)^s Γ(s) Z_Q(s).
- DH: Davenport–Heilbronn, ξ_DH = (5/π)^{(s+1)/2} Γ((s+1)/2) L_DH(s).
- FTH: F_θ = ξ(s−¼)ξ(s+¼).
- Beurling / splice: N/A, see §1.4.

**Files.**
- `r20_core.py`: Taylor coefficients, m-series, S-fraction, ladder.
- `r20_ladder.py` → `r20_ladder_{ZETA,FTH,EPS,DH}.json`, `ladder200_*.log`.
- `r20_pick.py` → `r20_pick_*.json`, `pick_*.log`.
- `r20_primes.py` → `r20_primes.json`.

---

## 1. The circuit: S-fraction and Cauer RC ladder

### 1.1 Exact statements

**Foster form under RH (PROVED, from the paper's genus-0 product).** Under RH, m(u) = Σ_n 1/(u + γ_n²). This is the impedance of a series chain of parallel RC tanks: every capacitor is 1 F and tank n has resistor 1/γ_n² Ω. RH says exactly that m is an RC driving-point impedance, i.e. a passive network.

**Cauer ladder (definitions).** Expand at u = 0 as an S-fraction:

  m(u) = c₀ / (1 + a₁u / (1 + a₂u / (1 + a₃u / …))).

As an admittance, Y = 1/m = G₁ + 1/(1/(C₁u) + 1/(G₂ + 1/(1/(C₂u) + …))). That is a shunt conductance G₁, then a series capacitor C₁, then a shunt G₂, and so on, with

  G₁ = 1/c₀,  C_j = G_j a_{2j−1},  G_{j+1} = C_j / a_{2j}.

All elements are positive iff c₀ > 0 and every a_k > 0.

**Moment theorem (Stieltjes 1894; cited).** A real sequence (c_j)_{j≥0} is c_j = ∫x^j dν(x) for a positive measure ν on [0, ∞) iff both Hankel forms are positive semidefinite for every n:
- H_n = (c_{i+j})_{i,j<n},
- H_n^{(1)} = (c_{i+j+1})_{i,j<n}.

If ν has infinite support, then all det H_n > 0 and all det H_n^{(1)} > 0. With ≥ 0 on H_n alone you get Hamburger's version: a measure on ℝ. The S-fraction coefficients are ratios of these determinants:

  a_{2n} = H_{n+1} H_{n−1}^{(1)} / (H_n H_n^{(1)}),  a_{2n+1} = H_n H_{n+1}^{(1)} / (H_{n+1} H_n^{(1)}).

So: all a_k > 0 ⇔ all Hankel determinants in both families are > 0.

**Theorem A (PROVED here from the cited moment theory).** RH ⇔ a_k > 0 for every k ⇔ every element of the Cauer ladder is positive.

*Proof.*
- (⇒) Under RH, c_j = Σ_n x_n^{j+1} with x_n = 1/γ_n² ∈ (0, γ₁^{−2}]. This is the moment sequence of the positive, infinitely supported measure Σ x_n δ_{x_n}.
- (⇐) Suppose all a_k > 0.
  - Stieltjes gives a positive ν on [0, ∞) with these moments.
  - The bound |c_j| ≤ R^j Σ|x_ρ| (R = max|x_ρ|, x_ρ := −1/u_ρ) forces ν to have compact support. The problem is then determinate, and G(w) = Σ c_j w^j = ∫ dν(x)/(1−xw) is analytic off a real half-line.
  - Unconditionally, G(w) = Σ_ρ x_ρ/(1 − x_ρ w) near 0, and this is meromorphic with a pole at every 1/x_ρ (residue ∝ multiplicity).
  - Identity theorem: every x_ρ must lie in [0, ∞). So every u_ρ is real and ≤ 0, which is RH. ∎

**Is a negative element forced by every off-line zero? Yes (PROVED).** One off-line zero, or a real zero of ξ off ½, makes some H_n or H_n^{(1)} fail to be positive semidefinite, so some a_k ≤ 0. *The index is finite but not bounded in terms of that zero alone*: it depends on all the zeros below it (§1.3).

### 1.2 Numbers (NUMERICAL)

Run settings: K = 200 Taylor coefficients, at 900 and 1300 digits, so 198 a_k per model. Every listed a_k agrees to ≥ 10 digits between the two runs. The worst agreement over all 198 elements is 475 digits for ZETA, 450 for FTH, 257 for EPS and 385 for DH.

| model | c₀ = Σ 1/(−u_ρ) | stable a_k | first a_k ≤ 0 | first negative ladder element | a_k there |
|---|---|---|---|---|---|
| ZETA | 0.023104993115418970… | 198 | **none (a₁…a₁₉₈ > 0)** | none | min a over 1…198 ≈ 1.0·10⁻⁵ |
| FTH = ξ(s−¼)ξ(s+¼) | 0.0461960521356 | 198 | **a₇** | C₄ | −0.012309 |
| EPS (x²+5y²) | 0.633316321443 | 198 | **a₂₆** | G₁₄ | −0.00096175 |
| DH | 0.0980695769533 | 198 | **a₁₄₈** | G₇₅ | −8.925·10⁻⁶ |

Notes on the table:
- ZETA's c₀ matches the classical value Σ_ρ 1/(ρ−½)² (CHECKED).
- First ZETA elements: G₁ = 43.2807, C₁ = 0.069632, G₂ = 30.6799, C₂ = 0.037765, G₃ = 35.1865.
- After the first failure, the fakes' signs keep alternating irregularly. The ζ ladder stays positive through a₁₉₈; its tail is a₁₉₆…a₁₉₈ ≈ 1.00·10⁻⁵, 1.02·10⁻⁵, 1.02·10⁻⁵.

### 1.3 How deep is a ghost buried? (HEURISTIC, with a PROVED identity)

Cauchy–Binet (PROVED) expands every Hankel determinant as a sum over n-element subsets S of the "nodes" x_ρ:

  det H_n = Σ_{|S|=n} Π_{ρ∈S} x_ρ · V(S)²,

where V(S) is the Vandermonde determinant of S.
- With real positive nodes, every term is positive.
- A complex ghost pair contributes terms with complex V(S)². These can win only once n is large enough to "resolve" the ghost among the real nodes below it.

Illustration only (zero heights quoted, not used):
- FTH's ghost is its lowest zero, and it fails at a₇.
- EPS's ghost is the Round-17 zero at height 15.67, and it fails at a₂₆.
- DH's first off-line zero sits near height 85.7 (literature value, not recomputed), and it fails at a₁₄₈.

For ζ, a hypothetical ghost above the verified height (~3·10¹²) would be buried at a depth comparable to the number of zeros below it, around 10¹³. **No finite ladder computation can prove or refute RH.** It can only certify "no ghost below a height."

### 1.4 Beurling / splice: not feasible (N/A, with reason)

These systems have no functional equation. So ξ_B(½+z) is not even in z, F_B(u) does not exist, and neither m_B nor its S-fraction at u = 0 is defined. The Round-19 table already records this as "FE: no."

One could try borrowing ζ's archimedean factor and evaluating the Dirichlet series only where Re s > 1. The prime sum then converges to roughly 10⁻²⁶ at s ≈ 5, using Beurling primes up to 3.4·10⁶. But Pick pivots shrink by about 10⁻⁸ per size (§2), so only sizes ≤ 3 would be testable, which is meaningless. Not run.

---

## 2. Pick matrices at the Shadow-Euler points

**Points.** u_{k,N} = (kN/(k+N) − ½)² for N = 2…81, with k = 1, 2, 3, 4. They accumulate at (k−½)², where m is analytic because ξ(k) ≠ 0.

**Löwner/Pick matrices for Stieltjes functions.**
- P0_{ij} = (m(u_i) − m(u_j))/(u_j − u_i), with diagonal −m′(u_i).
- P1_{ij} = (u_i m(u_i) − u_j m(u_j))/(u_i − u_j), with diagonal (um)′(u_i).

**Theorem B (PROVED, standard Löwner/Hamburger argument; sketch, not Lean-checked).**
- If m(u) = Σ w_ρ/(t_ρ + u), then P0 = Σ w_ρ v_ρ v_ρᵀ and P1 = Σ w_ρ t_ρ v_ρ v_ρᵀ, where v_ρ = (1/(t_ρ+u_i))_i.
- So under RH (w = 1, t = γ² > 0) both are positive-definite Gram matrices at every size.
- What each matrix detects:
  - P0 ⪰ 0 at all sizes forces all t_ρ to be real. This is the Nevanlinna/Pick condition.
  - P1 ⪰ 0 additionally forces t_ρ ≥ 0, which excludes real zeros of ξ off ½. P0 alone cannot see those.
- Converse: suppose both are ⪰ 0 on every finite subset of *one* Shadow-Euler family.
  - Confluent divided differences carry the positivity to the Taylor-Hankel matrices (−b_{i+j+1}) at the accumulation point, where m = Σ b_n (u−u*)^n.
  - The Hamburger moment theorem then makes m (anti-)Pick on ℂ∖ℝ, and P1 pins the support to t ≥ 0.
  - So m is Stieltjes, i.e. RH holds.
- So "Pick-positive at every size on one s_{k,N} family" ⇔ RH. This is the paper's uniqueness-set corollary in checkable matrix form.

**Numbers (NUMERICAL).** Leading-minor LDL pivots d_n = det_n/det_{n−1} in the order N = 2, 3, …. The two runs used (K = 300 Taylor terms, 1100 digits) and (K = 400 Taylor terms, 1500 digits). Varying K is essential; see the caveat below. Every reported pivot agrees to ≥ 8 digits between the runs, through size 80 for all models and all k.

| model | k = 1 | k = 2 | k = 3 | k = 4 |
|---|---|---|---|---|
| ZETA | pass ≤ 80 | pass ≤ 80 | pass ≤ 80 | pass ≤ 80 |
| FTH | fails at **4** (P0), 5 (P1) | 4 / 5 | 4 / 5 | 4 / 5 |
| EPS | fails at **14** (P0 and P1) | 14 | 14 | 14 |
| DH | fails at **75** (P0 and P1) | 75 | 75 | 75 |

What the table shows:
- For ZETA (k = 2) the pivots fall to about 10^{−833} by size 80 and stay positive throughout; the full sign string is 80 '+' signs.
- The failure size does not depend on k.
- It matches the ladder rung exactly:
  - FTH: a₇ is C₄, and P0 fails at 4.
  - EPS: a₂₆ is G₁₄ (from H₁₄), and P fails at 14.
  - DH: a₁₄₈ is G₇₅ (from H₇₅), and P fails at 75.
- Why: by Cauchy–Binet, det P_n has the same subset structure as det H_n; only the positive Cauchy factors depend on the evaluation points (HEURISTIC explanation, PROVED identity). **The Shadow-Euler Pick test and the ladder are the same test.** The n-point Pick matrix at any family ↔ ladder rung n ↔ Hankel size n.

Caveat (fixed). A first run used the same K = 160 at both precisions. It showed spurious "failures" at sizes 70–73 for ZETA and DH at k = 3, 4. The cause was Taylor truncation of F at u ≈ 12, which a precision-only cross-check cannot see. With K varied (300 vs 400) those disappear for ZETA, and DH fails at 75 for every k.

---

## 3. Arithmetic inside the circuit

**3.1 Where the primes live (PROVED).** The S-fraction at u = 0 is built from data at s = ½, inside the critical strip, where no Dirichlet series converges. To see primes, re-centre the ladder at u₀ > ¼, i.e. Re s₀ > 1. Under RH, m(u₀+v) is Stieltjes in v with nodes γ² + u₀, so the same Theorem A applies. The explicit formula gives

  m(u) = (1/(2√u)) [1/s + 1/(s−1) − ½ log π + ½ψ(s/2) − Σ_n Λ(n) n^{−s}],  s = ½ + √u,

and for each moment d_j (where m(u₀+v) = Σ(−1)^j d_j v^j) the prime part is

  d_j^{prime} = −Σ_n Λ(n) n^{−1/2} ∫₀^∞ (τ^j/j!) e^{−u₀τ} e^{−(log n)²/(4τ)} / (2√(πτ)) dτ.

This is the paper's heat-kernel ("critical prime current") form. Two consequences:
- **Every prime power enters every moment with a strictly negative sign** (PROVED: the integrand is positive).
- The positive moments of the true circuit are therefore *archimedean + pole minus primes*. Positivity is a cancellation, not a sum of positive parts.

**3.2 Numbers (NUMERICAL, `r20_primes.json`).** Re-centred at u₀ = 100 (s₀ = 10.5).
- Route A, the Taylor coefficients of ξ at ½ shifted to u₀, gives 38 ladder elements, all positive (min a ≈ 9.4·10⁻⁵).
- Route B is archimedean + pole terms plus primes ≤ P. The table compares it with Route A.

| primes used | moment error (max rel., j ≤ 9) | ladder elements matching Route A (rel. 10⁻⁶) | first a_k ≤ 0 |
|---|---|---|---|
| none (arch + pole only) | 9.2·10³ | 0 | **a₅** |
| P = 2 | 2.0·10³ | 0 | a₆ |
| P = 3 | 3.4·10² | 0 | a₇ |
| P = 5 | 20 | 1 | a₈ |
| P = 10 | 0.98 | 2 | a₉ |
| P = 30 | 4.7·10⁻⁴ | 4 | a₁₁ |
| P = 100 | 6.7·10⁻⁸ | 6 | a₁₃ |
| P = 1000 | 4.2·10⁻¹⁶ | 11 | a₁₇ |
| P = 10⁴ | 1.4·10⁻²⁴ | 14 | a₂₁ |

What the table shows:
- Route B converges to Route A (CHECKED: explicit formula ⇔ Taylor data at ½).
- The archimedean circuit alone is not passive: it fails at rung 5, because the 1/(s−1) pole sits at u = ¼ > 0, the wrong side.
- Each batch of primes repairs a few more rungs. Positivity to depth n needs primes up to roughly 10^{n/3.5}, and a truncated Euler product always breaks the circuit a few rungs past the last correct element.

**3.3 Do individual elements show log p? No (NUMERICAL / HEURISTIC).** No ladder element corresponds to a single prime or Euler factor. Primes enter linearly in the moments, but the elements are rational functions of Hankel determinants, so every element mixes all primes up to its scale.

The moments at u = 0 (the circuit of §1) are prime sums only through analytic continuation. The true m is analytic on the disc |u − u₀| < u₀ + γ₁², which contains 0, so c_j = Σ binomial combinations of the d_i. But the individual archimedean and prime pieces are singular at u = ¼ and branch at u = 0, so they must be summed before re-expanding.

---

## 4. Verdict, in plain words, and the link to the entanglement round

**Plain words.**
- RH is equivalent to a concrete physical statement: ξ's log-derivative, in the Casimir variable, is the input impedance of a *passive* RC ladder, with every resistor and capacitor positive. Equivalently, the Pick matrices at Daniel's Shadow-Euler points are positive at every size.
- Computing from ξ's Taylor coefficients alone, the ζ ladder is positive for all 198 rungs reached, and the Pick matrices are positive through size 80 for k = 1–4.
- Every fake breaks the circuit at a definite rung: F_θ at rung 4, Epstein at 14, Davenport–Heilbronn at 75.
- That last one matters. In Round 19 DH was as unique as ζ in its symmetry class, yet the circuit catches it at rung 75. The circuit is a *positivity* test, which is what uniqueness lacked.
- But a ghost hides at a depth that grows with its height. No finite computation reaches the ζ depth that matters, so the numerics certify, they do not prove.
- Inside the circuit, the primes always enter with a minus sign. The archimedean part alone breaks the ladder at rung 5. Each prime repairs a little more, and the full Euler product is needed to keep it passive all the way down.

**GAP-20 (circuit synthesis).** Build, from the primes and the archimedean place alone, a passive network whose impedance equals m. Concretely, any one of:
- a positive Krein string (masses and lengths),
- a positive Jacobi matrix,
- an RC ladder with all R_k, C_k > 0,

constructed so that positivity is manifest. This is the paper's "positive operator A" (Corollary, discrete positive-determinant closure) in circuit form. Theorem B says matching on one s_{k,N} family is enough. The obstacle is §3.1: the prime contributions are all negative, so a naive "sum of positive pieces" construction cannot work. Some structure must turn arch − primes into a Gram form. This is GAP-19 (positivity transfer) again, now as network synthesis. OPEN.

**Link to Round 21 (entanglement).** Round 21 (`round21_entangled/`, commit 5bc4b8c, read-only here) landed while this round was running. Its own Round-20 link states the same dictionary: a Stieltjes function with μ ≥ 0 is ⟨Ψ|(H+z)⁻¹|Ψ⟩ for a TFD purification with Schmidt weights μ_k ≥ 0. Its findings fit the ladder picture. The Bost–Connes TFD is positive but never sees the zeros. The theta TFD sees them but is indefinite for every model. F_θ has a positive, prime-factorized TFD and still fails here at rung 4. The dictionary:
- Suppose a bipartite state is built with reduced density ρ_A ∝ A, the paper's positive trace-class operator. Its Schmidt weights are λ_n = γ_n^{−2}/t₁ under RH, and its Rényi moments are Tr ρ_A^m = t_m/t₁^m. These are exactly this ladder's moments (c_j = t_{j+1}).
- Then (PROVED, as an equivalence): *positive Schmidt weights realizing these Rényi moments exist* ⇔ the Hankel families are positive ⇔ the ladder is passive ⇔ RH.
- For a genuine quantum state the Schmidt weights are automatically ≥ 0. So the entanglement route is a proof route only if the state is built from arithmetic without zeros *and* shown to reproduce the t_m (or the s_{k,N} values). That is GAP-20 in quantum language.
- The first failing ladder rung of a fake (4 / 14 / 75) is the smallest number of Schmidt levels at which its Rényi data stop being consistent with *any* density matrix.

## 5. Status ledger

- **PROVED (cited theory plus short arguments here):**
  - Theorem A: RH ⇔ all a_k > 0 ⇔ passive ladder, with an off-line zero forcing a finite failing index.
  - Theorem B: Pick-positivity at all sizes on one Shadow-Euler family ⇔ RH (sketch); P0 detects non-real zeros, P1 detects real off-½ zeros.
  - The Cauchy–Binet identities.
  - The strictly negative prime contribution to every re-centred moment.
  - The Foster form under RH.
- **CHECKED:** c₀ = 0.0231049931154…; explicit formula ⇔ Taylor data at ½, with moments agreeing to 1.4·10⁻²⁴ at P = 10⁴.
- **NUMERICAL:**
  - ζ ladder positive through a₁₉₈.
  - First failures: FTH a₇ (C₄), EPS a₂₆ (G₁₄), DH a₁₄₈ (G₇₅).
  - ζ Pick-positive through size 80 for k = 1–4; fakes fail at sizes 4–5 / 14 / 75, independent of k.
  - Prime-truncation table.
- **HEURISTIC:** depth of a ghost ≈ number of zeros below it; "Pick test ≡ ladder test" beyond the identity.
- **GAP-20:** zero-free synthesis of the passive network from primes + arch (= GAP-19 in network form).
- **OPEN:** RH.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).
