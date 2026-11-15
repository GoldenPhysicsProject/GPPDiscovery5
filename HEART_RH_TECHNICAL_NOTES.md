# HEART_RH_TECHNICAL_NOTES — the RH-strength step: state, exact open statements, plans

Daniel Toupin / Golden Physics Project … This seat: Heart of the Cosmos / Grok (GPPDiscovery5). Audit: Turbo.

Written 2026-11-14 (America/Toronto). Daniel directed (t316u, his own words) that this seat keep working past Turbo's usage hold. This file goes first, so nothing is lost if the budget runs out.

**Status line.** H_ψ, ARP and RH are **OPEN**. KreinDecomp is **OPEN** and **ON HOLD**. Nothing below claims any of them. Labels: PROVED (Lean, green CI, std axioms only), BANKED (PROVED plus Turbo PASS), NUMERICAL, SKETCH, OPEN.

---

## 0. The objects (exact, from `KreinHPsi.lean` / `WeilGramSemilocal.lean`)

- ψ(t) = max(0, 1 − |t|): the triangle on [−1, 1]. g = ψ⋆ψ = `triG`, supported on [−2, 2].
- Ψ(t) = (sin(t/2)/(t/2))². `bigPsiI2` = Ψ(i/2) = 16 sinh²(1/4).
- `poleTerm d` = 2Ψ(i/2)² cosh(d/2).
- `primeTerm d` = Σ_n Λ(n) n^{−1/2} [g(log n + d) + g(−log n + d)].
- `archTerm d` = (1/π) ∫_0^∞ Ψ(t)² cos(td) (Re ψ_Γ(1/4 + it/2) − log π) dt.
- `Kpsi` = poleTerm − primeTerm + archTerm. This is Kψ(d) = B(τ_0ψ, τ_dψ) on the zero-free side.
- `kform K x c` = Σ_i Σ_j c_i c_j K(x_i − x_j).
- ker ℓ_c: Σ c_i cosh(x_i/2) = 0. Window L: |x_i − x_j| ≤ L.
- `primeTermS S` is the prime term with n kept iff minFac n ∈ S. `KS S` = pole − primeTermS S + arch. `primesBelowExp L` = {p prime, p < e^{2+L}}.

## 1. Banked chain (all Turbo-banked, GPPVerify, none merged; heart/workbench stays at 4dce49e)

| Branch @ SHA | File | Build | Content (labels) |
|---|---|---|---|
| heart/hpsi @ 0e720c1 | ExpLaplaceMeroProof.lean | 43275270503 | ExpLaplaceMero PROVED (BANKED) |
| heart/hpsi @ 9cc5ed7 | HPsiSplit.lean | 43282044617 | 18 thms. `hPsi_iff_split`: HPsi ↔ RegKerPos ∧ RegUnitBound. `hPsiRegChain`: RegKerPos → RegUnitBound → KreinDecomp → RH (BANKED, Turbo 6833474105/6833489276) |
| heart/hpsi-plan @ 20050db | HPsiPlan.lean | 43277261882 | `not_primeSupportVanish` (ψ on [−2, 2] already sees log 2), `hPsi_iff_plan` (window doubling), poleRankOne, kform_poleTerm_of_ker, abelKernelPosDef (BANKED) |
| heart/weilgram @ 99f3b3a | WeilGramSemilocal.lean | 43285069393 | `kS_eq_Kpsi` (S covering L gives the exact Kψ on window L), `hPsi_iff_semilocalWeilIsGram` (BANKED) |
| heart/weilgram-attack @ f9f7676 | WeilGramContraction.lean | 43312478043 | 11 thms: kform_add, kform_sub, kform_KS_split, `kform_KS_of_ker`, `semilocalPos_iff_contraction`, semilocalContraction_iff_archDom, archDom_prime_le_arch, `semilocalPos_congr_cover`, `semilocalPos_anti_window`, hPsi_iff_semilocalContractionAll, `hPsi_iff_semilocalArchDomAll` (BANKED, Turbo 6837243506) |
| heart/archpd @ 5fb6c8d | WeilGramArchPD.lean | 43317360576 | 8 thms. `archNormBound_iff`: opnorm contraction ↔ ArchDom ∧ LowerDom, so it is strictly stronger than H_ψ and is not a route. `hPsi_iff_archDomAll_of_posDef` (BANKED, Turbo 6838006693) |

Earlier PROVED: CoeffVanish, ExpLaplaceMero, IdentityStep.

**NUMERICAL (Turbo, Arb).** Pole + archimedean alone is negative on ker ℓ_c. Adding the primes makes Kψ positive. So the primes carry the positivity, and every "no-prime base" route is dead.

**Literature.** Connes–Consani 2021 (archimedean/semilocal positivity) covers only widths < log 2, where no prime is visible. Connes' S-local trace formula is an identity, not a positivity statement.

## 2. Exact OPEN statements

All of the following are **OPEN**. The first five are each proved equivalent to H_ψ, so they are RH-strength.

```lean
def HPsi : Prop := ∀ (n : ℕ) (x c : Fin n → ℝ), ∑ i, c i * Real.cosh (x i / 2) = 0 → 0 ≤ kform Kpsi x c
def RegKerPos : Prop := ∀ n x c, coshSum x c = 0 → sinhSum x c = 0 → 0 ≤ kform Kreg x c
def RegUnitBound : Prop := ∀ n x c, coshSum x c = 0 → sinhSum x c = 1 → 2 * bigPsiI2 ^ 2 ≤ kform Kreg x c
def SemilocalPos (S) (L) : Prop := ∀ n x c, (∀ i j, |x i - x j| ≤ L) → ∑ i, c i * cosh (x i / 2) = 0 → 0 ≤ kform (KS S) x c
-- SemilocalWeilIsGram := ∀ L > 0, SemilocalPos (primesBelowExp L) L      (↔ HPsi, PROVED)
def SemilocalArchDom (S) (L) : Prop := ∀ n x c, (∀ i j, |x i - x j| ≤ L) → ∑ i, c i * cosh (x i / 2) = 0 →
  kform (primeTermS S) x c + 2 * bigPsiI2 ^ 2 * (∑ i, c i * sinh (x i / 2)) ^ 2 ≤ kform archTerm x c
def SemilocalArchDomAll : Prop := ∀ L > 0, SemilocalArchDom (primesBelowExp L) L   -- MY TARGET (↔ HPsi, PROVED)
def SemilocalContractionAll : Prop := ∀ L > 0, SemilocalContraction (primesBelowExp L) L  -- (↔ HPsi, PROVED)
def ArchPosDef (L) : Prop := ∀ n x c, window L → ker ℓ_c → c ≠ 0 → 0 < kform archTerm x c   -- OPEN, not RH-strength by itself
def SemilocalArchNormBound S L / SemilocalArchLowerDom S L                                  -- OPEN, stronger than H_ψ (not a route)
def KreinDecomp : Prop := AtMostOneNegSq Kpsi → ∃ Q, Q ≠ 0 ∧ Q.natDegree ≤ 1 ∧ ∃ h r, ExpPart Q h ∧ Continuous r ∧ PolyBounded r ∧ ∀ d, Kpsi d = h d + r d
   -- OPEN, ON HOLD (source check: Krein 1959 Dokl. 125:31–34; Krein–Langer 1977 Math. Nachr. 77:187–236; Stewart 1972)
```

In words, my target says: on every window, the S_L-prime energy plus the single pole correction is at most the archimedean energy, on the cosh hyperplane.

## 3. Relation to the board

**GPPVerify PR #597 (merged, Daniel's account / Claude Code).** `weilForm_autocorr_nonneg_l036_unconditional`: for every C^∞ g with compact support in [−9/25, 9/25], weilForm(g ⋆ g̃) ≥ 0. Only the prime 2 is visible (0.72 > log 2, and 0.72 < log 3). PROVED, unconditional.

What it gives my target:
- **Scale obstruction.** Every vector in my space F_ψ = span{τ_x ψ} has support length ≥ 2, since each translate has support length 2 and the window only adds to it. #597 needs support length ≤ 18/25 = 0.72 < 2. So #597 implies **no** instance of SemilocalPos for my fixed ψ, at any L ≥ 0. Formalised as `not_hasSupportIn_tau` in Phase 2 below.
- **Regularity obstruction.** IsTestFunction is C^∞, and ψ is only Lipschitz. A mollification or continuity step is needed in any case.
- **What would transfer.** Rebuild the semilocal machinery on ψ_a (triangle of half-width a), smoothed. Then #597 gives the window statement for L + 2a ≤ 18/25, with prime 2 visible. This is the base case for plan B1 below. It needs an identification weilForm(autocorr(Σc_i τ_{x_i}ψ_a)) = kform Kψ_a x c, which is OPEN (SKETCH: pole ∫g·2cosh(u/2) = 2Ψ_a(i/2)²cosh(d/2) termwise; the prime and archimedean sides match by Fourier/Parseval).
- **Method.** #597 is a certificate method: a finite Legendre block LDLᵀ plus a minorant m ≤ Ψ₂ − β₂ plus a tail. That is exactly a template for certifying SemilocalArchDom on one fixed window. It is not a template for all L.

**Flint's obligation C2 (CONVERSATION.md, eighth run, Turbo-audited).** "Continuity of the residue functional at the singular labels for the WeilIsGram approximating families, with the extracted-class norm formula." This sits on the **spectral (zero) side**. The completed pairing 𝒫(α) = −(ξ′/ξ)(1/2 + α), and the state at ρ is the residue at α = ρ − 1/2, transported by the Cesàro factor ε²/(ε² − (2σ_ρ − 1)²) (residue first, ε → 0 second).

How it relates to my target. C2 and SemilocalArchDomAll are the two faces of one Weil explicit formula:
- C2: the zero-side Gram/GNS representation exists, with norms m_ρ on the line and 0 off it.
- Mine: positivity of the zero-free (prime + pole + arch) Toeplitz form.

H_ψ ⇔ SemilocalArchDomAll is PROVED (zero-free). Flint's Gram-form WeilIsGram ⇒ HPsi is her construction, and C2 is the obligation inside it. They are distinct named targets with distinct owners (Turbo, t310u). Neither feeds the other as a proof input; in particular my side uses no zeros.

## 4. Plans

### Primary plan P — two-point slice, then window induction, with numerics guiding every step

1. **(Lean, cheap)** Exact reduction of the n = 2 slice. On ker ℓ_c with two points, c = t·(b, −a) where a = cosh(x₁/2) and b = cosh(x₂/2). For even K:
   kform K = t²[(a² + b²)K(0) − 2abK(x₁ − x₂)].
   So the 2-point H_ψ is the scalar inequality 2ab·Kψ(d) ≤ (a² + b²)Kψ(0). Prove it for general even K, and prove Kψ even.
2. **(NUMERICAL)** Evaluate Kψ on [0, L] with mpmath (the prime sum is finite: log n ≤ 2 + |d|). Check the 2-point inequality, and λ_min of the cosh-hyperplane compression of [Kψ(x_i − x_j)] on random and adversarial point sets. Report the margin as a function of L. Falsifier: a negative λ_min, which would mean a sign or normalisation error, since RH predicts ≥ 0.
3. **(Lean)** Base-window certificate in the #597 style. Fix L₀ (for example 0.72) and certify SemilocalArchDom for **all** n at window L₀. Two sub-routes:
   - (i) a continuous form: the kernel restricted to |d| ≤ L₀, written as an operator on L²-type measures with a Legendre block;
   - (ii) a Bochner route: show K_{L₀} := Kψ·1_{|d| ≤ L₀} extends to a positive-definite kernel modulo the cosh constraint.
   This is the honest base case. It is not RH.
4. **(OPEN, the RH step)** Window doubling: from SemilocalPos(2L) ⇐ SemilocalPos(L) plus a prime-block estimate. `hPsi_iff_plan` gives the doubling skeleton. The missing inequality is a uniform bound on the cross block between the two half-windows in terms of the diagonal blocks. Falsifier: any L where the cross-block ratio exceeds 1 numerically while H_ψ still holds, which would show the doubling inequality is false even though RH may be true. Check numerically before writing Lean.

### Backup B1 — rescaled smooth ψ_a, using #597 as the base case
Restate the semilocal chain for ψ_a (triangle of half-width a, mollified). Prove (a) H_{ψ_a} ⇔ the same RH chain (port HPsiSplit/HPsiPlan with a as a parameter), and (b) the weilForm/kform identification. Then #597 ⇒ SemilocalPos_a(L) for L ≤ 18/25 − 2a: a genuine prime-visible base case. Falsifier: the identification fails (normalisation 2π, or a log π shift in the arch term). Numeric check: compare weilForm(autocorr f) from #597's `primeSide`/`archimedeanTerm` against kform Kψ_a for one sample, to 1e-12.

### Backup B2 — Toeplitz/Carathéodory–Fejér in cosh-weighted coordinates
On ker ℓ_c, substitute c_i = u_i / cosh(x_i/2). The constraint becomes Σu_i = 0, and the kernel becomes K̃(x, y) = Kψ(x − y)/(cosh(x/2)cosh(y/2)). Positivity on Σu = 0 is conditional positive definiteness of K̃. This is Schoenberg's characterisation: K̃ is CPD iff exp(tK̃) is PD for all t > 0. Next Lean step: state the Schoenberg equivalence for finite matrices (PROVED-able). Numerics: test PD of exp(tK̃) on grids. Falsifier: CPD failing for small t on some grid.

### Backup B3 — Krein κ ≤ 1 / Pontryagin (Flint's #593/#594 layer), kept separate from KreinDecomp
#594 says the Toeplitz forms of an even AtMostOneNegSq kernel split off at most one negative square. H_ψ says that on ker ℓ_c the single negative square of Kψ (carried by the pole, rank one) is exactly removed. Route: prove AtMostOneNegSq Kψ (κ ≤ 1) and identify the negative direction with the cosh vector. Then H_ψ follows from "the negative eigenvector is not orthogonal to cosh". This is RH-strength at the κ ≤ 1 step. Do not touch KreinDecomp (on hold). Numeric: count the negative eigenvalues of [Kψ(x_i − x_j)] on grids. The expected count is exactly 1; 2 or more would falsify the plan.

### Backup B4 — semilocal scattering, Daniel's S-matrix
The compressed S_L-local Connes cutoff operator on ker ℓ_c: contraction ⇔ ArchDom, with the opnorm version strictly stronger (PROVED, archpd). The remaining unitarity route is to show the S_L scattering matrix is an inner function on the window, where the p-adic factors (1 − p^{−s})^{−1} for p ∈ S_L give the commutator with the archimedean cutoff. Literature check first: Connes–Consani 2021 extended to the p-adic places. Falsifier: commutator norm > 1 at L just above log 2 (NUMERICAL).

## 5. Rules for every step
No sorry, admit or custom axioms. No bare `import Mathlib`. Native stub audit before push. Green GitHub Build at the exact head. No force-push, no merge, no PR. Don't touch KreinDecomp, Flint's items 1–4/C2, Golden Eye's #567, or other seats' repos. No zeros data as proof inputs. Honest labels throughout.
