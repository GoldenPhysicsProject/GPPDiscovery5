# Round 13 — AKG via the Pasterski–Shao Klein–Gordon construction: plan, backups, kill criteria

Heart of the Cosmos (Grok / GPPDiscovery5), 2026-11-18 (Toronto). Written BEFORE any Round 13 numerics, per Daniel's instruction.
Labels: PROVED / CHECKED / NUMERICAL / CONJECTURE / GAP / OPEN. AKG, H_ψ, ARP, RH: **OPEN**. Zeros (Platt tables in GPP Drive "Zeros", mpmath) are **test data only**, never a construction input.

## 0. Target (unchanged from Round 11)
**AKG (OPEN, RH-equivalent).** A Hilbert space H and a linear map g ↦ Φ_g ∈ H built only from primes, the pole and the archimedean place, with
‖Φ_g‖²_H = W(g) := pole − prime + arch  for every real test g supported in [−L/2, L/2], for every L at once.
Normalisation used throughout this round (same as the D2/Verify kernel form, checked against the zero side in `akg_kg.py`):
  W(g) = C_∞ A(0) + 2∫₀^∞ w(t)[A(0) − A(t)] dt + 2∫ A(t) cosh(t/2) dt − 2 Σ_n Λ(n) n^{−1/2} A(log n),
  A(t) = ∫ g(x) g(x+t) dx, w(t) = e^{−t/2}/(1 − e^{−2t}), C_∞ = ψ(1/4) − log π;  zero side Σ_ρ |ǧ(γ_ρ)|², ǧ(r) = ∫ g(y) e^{iry} dy.

Constraints any candidate must meet (Round 11/11c): L1 (dilation-respecting forms put off-line pairs in a (1,1) plane); O2 (a translation-invariant Q = A*A is forced to be the zero measure); O7/O8 (B = arch₋ − prime + pole = Z − arch₊; the norm must be null to exponential order on the zero-free band |t| < 14.13 and do all its work above it).

## 1. The PS mechanism, stripped to what transfers
Pasterski–Shao (1705.01027 §3.2, §4): the KG product (3.6) is positive on positive-energy solutions *because the solutions already form a Hilbert space* (plane waves, eq. 3.7). Conformal primaries φ_Δ = Mellin transforms in ω of plane waves. Their KG products are ∝ δ(i(Δ₁ + Δ₂* − d)): **δ-normalizable iff Re Δ = d/2**, otherwise divergent at one end of the ω-integral. In the 1+1 chiral case the KG product in frequency space is ∫ ψ̃₁ ψ̃₂* dω (positive frequency), and the Mellin transform diagonalises it: H_KG ≅ L²(Re Δ = 1/2). So the transferable content is: **(positive space given in advance) + (Mellin modes) ⇒ only the principal line is δ-normalizable.**

Arithmetic dictionary (u = log x is the "bulk time", ρ plays Δ):
- bulk field: Φ(u) = e^{−u/2}[ψ(e^u) − e^u + log 2π + ½ log(1 − e^{−2u})] = −Σ_ρ e^{(ρ−1/2)u}/ρ  (von Mangoldt). Built from primes, pole (e^u) and arch (the trivial-zero / constant terms) only.
- smeared field (detector): Φ_g(u) := Σ_ρ ǧ(γ_ρ) e^{iγ_ρ u}, which the explicit formula applied to g(· − u) computes from primes + pole + arch with no zeros. For u > L/2:
  Φ_g(u) = e^{−u/2} ǧ(i/2) + e^{u/2} ǧ(−i/2) − Σ_n Λ(n) n^{−1/2} g(log n − u) − ∫ w(x) g(x − u) dx.
- KG/δ-normalizability ↔ Besicovitch B² (Cesàro mean square): a mode e^{(ρ−1/2)u} has Cesàro norm 1 if Re ρ = 1/2 and ∞ if Re ρ > 1/2 (Daniel, t339u: "with Cesàro it's 1 on the line and diverges to infinity anywhere else").

## 2. Candidates, in order (each gets: statement, test at L = 0.5, 1, 2, verdict)
- **C-B2 (lead): Cesàro/Besicovitch KG norm.** H = B²(ℝ₊) (Besicovitch almost-periodic space, mean ⟨F,G⟩ = lim (1/U)∫₀^U F Ḡ). Map g ↦ Φ_g (prime/pole/arch built). Claim to test: ‖Φ_g‖²_{B²} = W(g).
  Expected: under RH, ‖Φ_g‖² = Σ_γ m_γ² |ǧ(γ)|² (almost-periodic Parseval), = W(g) iff zeros simple where ǧ ≠ 0; with an off-line zero, ‖Φ_g‖² = ∞. Normalizability ⇔ RH (Landau). Positivity free; O8 automatic (no spectrum below 14.13).
  Test: Hann-windowed mean of Φ_g² on u ∈ [u₀, U] (U ≈ 12–13, primes to e^{U+L/2}) vs W(g) (prime/arch side) vs zero side.
  Kill criteria: mismatch beyond the window-leakage estimate; or the identity only being reachable through a statement that is itself RH-equivalent (then it is a reformulation, not a construction — record exactly where).
- **C-Abel (backup 1): weighted-KG norm.** N_ε(g) = 2ε ∫₀^∞ e^{−2εu} Φ_g(u)² du, ε ↓ 0 (the PS ω-integral cut off by ω^{−2ε}). Same limit as C-B2; tests whether a one-sided Abel regulator changes anything (it should not).
- **C-Window (backup 2): finite-U KG norm without the limit.** Q_U(g) = ∫ σ_U(u) Φ_g(u)² du for a fixed smooth window. Positive, prime-built, all L. Measure Q_U − W at L = 0.5, 1, 2; it is the leakage Σ_{γ≠γ'} ǧ(γ) ǧ(γ')* σ̂_U(γ − γ') + edge terms. Dies if Q_U − W is not controllable from the prime side (expected: it is indefinite, so Q_U ≥ 0 says nothing about W at finite U).
- **C-Mult (backup 3): fix the multiplicity mismatch** (m² vs m) if C-B2 survives otherwise; candidates: cross-mean of Φ_g with the unsmeared field, or B² of the logarithmic field log ζ rather than ζ′/ζ.
- **C-Pair (backup 4, if all above reduce to normalizability):** write the B² identity on the prime side — a two-point sum over prime pairs (m, n) with |log(m/n)| ≤ L — and check which averaged Hardy–Littlewood statement it is.

## 3. What would count as progress, plainly
- A construction where Re ρ = 1/2 is *forced* by normalizability needs an a-priori reason that Φ_g ∈ H (PS have it: solutions of KG form a Hilbert space before any Mellin transform). If the a-priori bound on Φ_g is itself RH-equivalent, the candidate is a correct reformulation and an exact diagnostic, not a proof. That is the expected verdict for C-B2; the job is to pin the precise residual statement.

## 4. Deliverables
`akg_kg.py` (all numerics, JSON out), `README.md` (results, labels, obstructions, next step), `handoff_checks/` (Heart-arm checks of the D2 Mellin twist and Woodbury identities for the PR #14 hand-off to Claude).

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5). Audit: Turbo.
