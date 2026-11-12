# H_ψ ⇒ RH, Lean-typed (Turbo PR #11 6823553961)

Source: `lean/s-eq-delta/KreinHPsiChain.lean`, byte-identical to GPPVerify `heart/hpsi` 9b67b3b
(`GppVerify/RiemannHypothesis/KreinHPsiChain.lean`, sha256 d56e5281…c22c1). Compiled locally
(`lake env lean`, Mathlib 0df444a3, Lean 4.33.1): standard three axioms only, no warnings. Python
gates (import graph, sorries, stub naming, vacuity, blueprint refs, landing claims 1281/120/0)
pass locally. The receipt is the branch Build run.

| Link | Statement | Label |
|---|---|---|
| 1 | `atMostOneNegSq_of_hPsi : HPsi → AtMostOneNegSq Kpsi` (hPsiEquiv + bridge + fit_translates) | PROVED |
| 2 | `PoleTransfer` (AtMostOneNegSq Kpsi → ∃ Q ≠ 0, deg ≤ 1, every γ ∈ Zgamma with Im γ < 0, Ψγ ≠ 0 is a root of Q or Q̄) | named premise; paper = hEF + Krein P≤1 (secondary) + (L) + (G1) + identity thm |
| 3a | `Zgamma_neg` (functional equation) | PROVED |
| 3b | `Zgamma_negconj` (riemannZeta_conj) | PROVED |
| 3c | `Zgamma_re_ne_zero` from `NoRealZeroStrip` (ζ(σ) ≠ 0, 0<σ<1) | PROVED mod named premise (classical, not formalised) |
| 4 | `Zgamma_real` via `tri_RH_of_pole_transfer` | PROVED mod premises |
| 5 | `riemannHypothesis_of_hPsi`, `hPsiChain : HPsi → PoleTransfer → NoRealZeroStrip → RiemannHypothesis` (Mathlib's RH, via riemannHypothesis_of_ne_zero) | PROVED (conditional) |
| — | `poleTransfer_of_RH : RiemannHypothesis → PoleTransfer` | PROVED (honesty check: the premise carries the analytic load) |

H_ψ, ARP, RH: OPEN. Krein–Langer primary still owed (non-blocking).
Remaining: split `PoleTransfer` into `hEF` and a Krein-representation premise as exact Props
(needs a zero-sum with multiplicities); discharge `NoRealZeroStrip` (η(σ) > 0) if wanted.

## PoleTransfer split (Turbo PR #11 6824324169)

Source: `lean/s-eq-delta/KreinPoleSplit.lean`, byte-identical to GPPVerify `heart/hpsi` 3c8fd2a
(`GppVerify/RiemannHypothesis/KreinPoleSplit.lean`, 144 lines, sha256 aa06ad46…a98b551). Compiled
locally: standard three axioms, no warnings. The receipt is the branch Build run.

| Item | Statement | Label |
|---|---|---|
| mult | `mult γ = zeroMult (rho γ) = analyticOrderNatAt riemannZeta (1/2 + iγ)` | definition |
| one_le_mult | every γ ∈ Zgamma has `1 ≤ mult γ` (no zero has weight 0) | PROVED (from `one_le_zeroMult`) |
| zero sum | `ZeroSummable w := Summable (γ : Zgamma ↦ m_γ w γ)`, `zeroSum w := ∑' γ : Zgamma, m_γ w γ` | definition |
| HEF | ∀ d, ZeroSummable and `Kψ(d) = zeroSum (γ ↦ Ψ(γ)Ψ(−γ)e^{iγd})` | named premise (CHECKED secondary, NUMERICAL 1.02e−8) |
| KreinRep | ∀ a, Kψ has expansion `Σ m_γ a(γ) e^{iγd}` → `AtMostOneNegSq Kpsi` → ∃ Q ≠ 0, deg ≤ 1, every γ ∈ Zgamma, Im γ < 0, a γ ≠ 0 is a root of Q or Q̄ | named premise (Krein P≤1 via Stewart 1972 + (L) + (G1) + identity thm) |
| assembly | `poleTransfer_of_split : HEF → KreinRep → PoleTransfer` (uses `triPsi_neg`) | PROVED |
| chain | `hPsiSplitChain : HPsi → HEF → KreinRep → NoRealZeroStrip → RiemannHypothesis` | PROVED (conditional) |
| honesty | `kreinRep_of_RH : RiemannHypothesis → KreinRep` | PROVED |

H_ψ, ARP, RH: OPEN. Next (Turbo item 2): discharge `NoRealZeroStrip`.
