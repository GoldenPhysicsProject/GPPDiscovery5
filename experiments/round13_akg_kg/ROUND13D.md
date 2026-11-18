# Round 13d: non-symmetric arithmetic input (Λ ≥ 0, Euler product, 3-4-1) as a source of normalizability

Heart of the Cosmos (Grok / GPPDiscovery5), 2026-11-18 (Toronto). This follows `ROUND13B.md` (`fb95f79`). No PR post.
Labels: PROVED / CHECKED / NUMERICAL / CONJECTURE / GAP / OPEN. AKG, H_ψ, ARP, RH: **OPEN**.
Scripts:
- `nonsym.py` → `nonsym_results.json`, `nonsym.log`
- `beurling_toy.py` → `beurling_toy_results.json`

The off-line quartets ρ = 1/2 ± β ± 28i are hypothetical toys. No zeta zero is a construction input.

**Bottom line.** None of the three candidates beats κ = 1/2 on unconditional input. Nothing here implies a zero-free strip. The literature check below explains why. It agrees with the numbers.

## (a) Λ ≥ 0 through Landau's lemma: RIGHT, but it does not beat κ = 1/2
- **Setup.** The KG norm in Laplace form is N(σ) = ∫ e^{−σu} Φ_g(u)² du. The integrand is nonnegative, so Landau's lemma (Laplace-integral form; Montgomery–Vaughan MNT I, Thm 15.1 generalisation of Thm 1.7) says the abscissa of convergence is a real singularity of N.
- **What it gives.**
  - For the true primes the singularity is at σ = 0 (CHECKED, from R13.1 under RH): N(σ) ≈ W/σ.
  - For a field with an off-line quartet it is at σ = 2β.
  - Landau's lemma is satisfied in both cases. It *locates* the divergence at a real point, which is the Ω-direction (Landau's oscillation theorem, MV Thm 15.2). It gives no upper bound on N, and normalizability needs an upper bound. RIGHT as stated, but no forcing.
- **Numbers (NUMERICAL; L = 0.5, Laplace norm / W; ratio of the U = 13 to U = 9 truncation of the Laplace integral):**

  | σ | true primes | toy (β = 0.1) | toy / true at U = 13 |
  |---|---|---|---|
  | 0.05 | 1.476 | 1.493 | 1.052 |
  | 0.1 | 1.359 | 1.372 | 1.048 |
  | 0.2 | 1.201 | 1.208 | 1.041 |
  | 0.3 | 1.110 | 1.114 | 1.035 |

  Both are consistent with their predicted abscissae, 0 for the true primes and 0.2 for the toy. At U ≤ 13 the toy's excess is a few percent, so Landau's real singularity is far too weak a signal to build a bound on.
- The "Dirichlet series in prime pairs" version (the C-Pair expansion of N) has coefficients Λ(m)Λ(n)(mn)^{−1/2}g g. The diagonal and off-diagonal parts are nonnegative only before the main-term subtraction. After it they are signed, so Landau's lemma no longer applies to the part that matters.
- **Verdict: right, does not beat κ = 1/2.**

## (b) Euler-product constraints on an injected quartet: the naive test is WRONG (an artifact); a consistent model shows positivity cannot force anything
1. **Injecting the quartet on top of the true primes** means dψ_fake = dψ − Σ_q x^{ρ_q−1}dx.
   - **Λ ≥ 0 between integers is violated.** The fake prime density there runs from −1.20 to +1.27 (β = 0.1) and from −1.50 to +1.55 (β = 0.3), and is negative on 49.7% of [10, 10⁴].
   - **Λ(n) ≤ log n is not violated**, because the atoms are unchanged.
   - **The 3-4-1 combination 3D(σ) + 4Re D(σ+it) + Re D(σ+2it), with D = −ζ′/ζ − Σ_q 1/(s − ρ_q), goes negative for every β tested:**

     | Quartet | min E | at (σ, t) |
     |---|---|---|
     | true ζ (no quartet) | +1.85 | — |
     | β = 0.1 | −2.43 | (1.6, 28) |
     | β = 0.3 | −2.95 | (1.6, 28) |
     | β = 0.47 | −5.21 | (1.3, 28) |
     | β = 0.49 | −9.85 | (1.1, 28) |

     NUMERICAL.
   - **But the inference "true primes obey these, the fake breaks them, so no off-line zero" is WRONG.** Adding zeros to ζ while keeping the rest of −ζ′/ζ fixed is not a possible state of the world. It is ζ·P for a polynomial P, which has no Euler product, so its log-derivative has signed coefficients by construction. A genuine off-line zero would come with a globally different prime measure. The violation is an artifact of the injection.
2. **A consistent counter-model (CHECKED analytically + NUMERICAL).** Take continuous Beurling primes dΠ(x) = 1_{x≥20}(1 − Σ_q x^{ρ_q−1})dx/log x with ρ_q = 1/2 ± 0.1 ± 28i.
   - Then −ζ_C′/ζ_C(s) = 20^{1−s}/(s − 1) − Σ_q 20^{ρ−s}/(s − ρ). The residue at ρ is −1 (checked: −1.0000000017), so ζ_C has simple zeros exactly at the quartet.
   - The prime measure is nonnegative: density ≥ 0.127 on [20, 10⁸].
   - The 3-4-1 combination is ≥ 0, with min 0.057 over σ ∈ [1.001, 2], t ∈ [0.5, 60]. It must be, because 3 + 4cos θ + cos 2θ = 2(1 + cos θ)² ≥ 0.
   - So **Λ ≥ 0 plus Euler-product (Beurling) structure is fully compatible with a zero at Re ρ = 0.6.**
3. **Literature, which settles it (CHECKED):**
   - Diamond–Montgomery–Vorhauer, *Beurling primes with large oscillation*, Math. Ann. 334 (2006). They build discrete Beurling primes with integer counting N(x) = κx + O(x^θ) whose zeta has infinitely many zeros on σ = 1 − a/log t. So the de la Vallée Poussin region is the best that positivity and multiplicative structure alone can give. In their words, any wider zero-free region "will need to exploit … the additive structure of the integers". Vinogradov–Korobov does exactly that, and still gives no strip.
   - Davenport–Heilbronn (1936) and Potter–Titchmarsh: Epstein zeta functions with class number > 1 have nonnegative integer-supported coefficients and a Riemann-type functional equation, yet have zeros with σ > 1. So FE plus coefficient positivity, without an Euler product, also fails.
   - The 3-4-1 inequality's actual reach is the classical region σ > 1 − c/log t. In the width language that is a defect O(e^{L/2 − c√L}), with **κ still 1/2** and only a subexponential gain.
- **Verdict:**
  - Λ ≥ 0 and Λ ≤ log n: right facts, but they cannot force normalizability.
  - The 3-4-1 inequality: right, and it beats e^{L/2} only by the classical subexponential factor (de la Vallée Poussin), not in the exponent.
  - The injection test: WRONG as evidence. It is an artifact.

## (c) The signed sum (prime − pole cancellation): exponent 0 numerically, but only because the true primes cancel; any proof is RH-equivalent; unconditionally κ = 1/2
- **Exact cancellation (CHECKED).** pole − prime = 2∫A e^{−t/2}dt − 2∫₀^L A(t)e^{−t/2} dR(e^t), with R = ψ − x. The e^{t/2} parts of the pole and the prime main term cancel identically, so W = arch + 2∫A e^{−t/2} − 2∫A e^{−t/2} dR(e^t).
- **Signed separate bound λ_min(arch) − λ_max(prime − pole), Galerkin with M = max(24, 10L′) sines (NUMERICAL):**

  | L′ | λ_max(prime), absolute | λ_max(prime − pole), signed | λ_min(arch) | signed bound |
  |---|---|---|---|---|
  | 2 | 1.80 | 1.49 | −2.59 | −4.08 |
  | 4 | 7.78 | 1.56 | −3.56 | −5.13 |
  | 6 | 22.0 | 1.55 | −4.10 | −5.65 |
  | 8 | 57.9 | 1.60 | −4.44 | −6.04 |

  The signed prime operator stays flat at about 1.5–1.6 while the absolute one grows like e^{L′/2}. The bound tends to C_∞ − 1.6 ≈ −7. **Numerically the exponent drops from 1/2 to 0.**
- **Why this is not a result (CHECKED):**
  - The flat 1.6 is the true primes' cancellation. Since W ≥ 0 numerically, prime − pole = arch − W ≤ arch, and λ_max(arch) ≈ 1.6 on this basis. Proving λ_max(prime − pole) ≤ C uniformly in L′ would give a bounded defect, hence κ = 0, hence RH by R13.3. **Proving the signed bound is RH-equivalent.**
  - The Galerkin basis only resolves frequencies |r| ≲ Mπ/L′ ≈ 31. At very large r, Kronecker/Bohr alignment of p^{−ir} makes |Σ_{p<e^{L′}} log p · p^{−1/2−ir}| approach the full absolute sum, about 2e^{L′/2}. There, only the arch growth (log r) can compete. The flat 1.6 is a statement about low and moderate frequencies.
  - The best unconditional signed inputs (de la Vallée Poussin, Vinogradov–Korobov error terms for ψ(x) − x, explicit versions by Platt–Trudgian 2021) improve the defect only to e^{L′/2 − c L′^{3/5}(log L′)^{−1/5}}. Even that holds only where the transfer from the envelope to the operator symbol works: |r| ≲ e^{O(√L′)} by partial summation. **κ = 1/2, with a subexponential gain.**
- **Verdict: right identity, beats the absolute bound numerically, does not beat κ = 1/2 on unconditional input.**

## Did anything beat κ = 1/2? Was the error found first?
- **No unconditional candidate beats κ = 1/2.** Signed (c) beats it only when it is fed the true primes, and then it is equivalent to the target. A claim of κ < 1/2 would be a zero-free strip Re ρ ≤ 1 − δ.
- That is not known. The best known region is Vinogradov–Korobov, σ > 1 − c(log t)^{−2/3}(log log t)^{−1/3}, which shrinks to 1. DMV show positivity and multiplicativity alone cannot even reach Vinogradov–Korobov.
- So any apparent strip from these inputs must contain an error. The one I found and removed is the injection artifact in (b): a quartet added on top of the true ζ "violates" 3-4-1 even at β = 0.1. Taken at face value that would "prove" a strip. It is wrong, because the injected object is not an Euler product. The consistent Beurling model shows this.

## O9 after Round 13d
The forcing must come from what DMV identify as missing in Beurling theory: the **additive structure of the integers** (the n-ary structure behind exponential-sum bounds). The ordering positivity of Λ, the multiplicative structure, and the functional equation all fail as sources, each with an explicit counter-model:
- Beurling (DMV, and the toy above): positivity and Euler-product structure, with off-line zeros.
- Epstein (Davenport–Heilbronn): FE and positive coefficients, with off-line zeros.
- Off-line quartets: FE symmetry and S-unitarity hold for them (13c).

**Next step:**
- Define the quantity that separates ζ from the DMV and Epstein counter-models inside the KG framework. The test: compute the Cesàro/KG field Φ_g and the signed bound (c) for (i) a DMV-type *discrete* Beurling system and (ii) an Epstein zeta function (class number 2, e.g. Q = x² + 5y²), using their own "primes" or coefficients.
- Then find which term of the two-space decomposition differs from ζ's. That term is where integer additivity must enter.
- This is a falsifier harness: any proposed normalizability source must fail on (i) and (ii). It is not a proof route.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5). Audit: Turbo.
