# Round 11: all-L identity. Lead candidate: the Klein–Gordon/Mellin (PSS) inner product on zero modes

Heart (Grok / GPPDiscovery5), PR #11. Source of the lead argument: Daniel Toupin, "The Reason Riemann's Hypothesis Is True: Re(Δ) = d/2" (GPP note, 2026-11-11).
Template: Pasterski–Shao arXiv:1705.01027 §3.2, §4 (eqs. 4.7–4.9); Pasterski–Shao–Strominger arXiv:1701.00049 (App. A); secondary inputs arXiv:2310.04932 (Pasterski, *A Chapter on Celestial Holography*), arXiv:2202.11127 (*Goldilocks modes*), arXiv:2307.16801 (Hu–Pasterski, *Detector operators*).
Labels: PROVED / CHECKED / NUMERICAL / CONJECTURE / GAP / OPEN. H_ψ, WeilIsGram, ARP, RH: **OPEN**. No zero data is used as a proof input; s1/s2 use zeros only as comparison targets or as an explicitly hypothetical toy.

## 0. Objects (in-tree, unchanged)
- Ψ̂(t) = (sin(t/2)/(t/2))⁴, Ψ = cubic B-spline `triG` on [−2,2]; Ψ̂(i/2) = Ψ(i/2)² = 256 sinh⁴(1/4).
- K = pole − primeS + arch, pole(d) = 2Ψ̂(i/2)cosh(d/2), primeS(d) = Σ_{n ≤ e^{2+L}} Λ(n) n^{−1/2}(Ψ(log n + d) + Ψ(−log n + d)), arch(d) = (1/π)∫₀^∞ Ψ̂(t)[Re ψ(1/4 + it/2) − log π] cos(td) dt.
- kform(K)(x,c) = Σ_{ij} c_i c_j K(x_i − x_j); ker ℓ_c: C := Σ c_i cosh(x_i/2) = 0; S := Σ c_i sinh(x_i/2).
- SemilocalArchDomAll(L) ⇔ kform(K) ≥ 0 on ker ℓ_c, |x_i − x_j| ≤ L (since kform(pole) = 2Ψ̂(i/2)(C² − S²), see P1).
- Mode amplitude at s = 1/2 + iz: F_c(z) = Σ c_i e^{i z x_i}, ψ̂(z) = (sin(z/2)/(z/2))² (so Ψ̂ = ψ̂²), G_c(s) := ψ̂(z)F_c(z). Since z_{τρ} = conj z_ρ and c is real: G_c(ρ)·conj G_c(τρ) = Ψ̂(z_ρ)F_c(z_ρ)F_c(−z_ρ) (CHECKED, one line).
- τ(s) = 1 − s̄ (shadow Δ ↦ d − Δ̄ at d = 1).

## 1. The KG/Mellin pairing at d = 1 (template, from 1705.01027 eq. 4.7–4.8)
In log variable u = log ω, a weight-Δ mode is e_Δ(u) = e^{(Δ − 1/2)u}. The KG/L² pairing is
  ⟨e_a, e_b⟩ = ∫ e^{(a − 1/2)u} conj(e^{(b − 1/2)u}) du = ∫ e^{(a + b̄ − 1)u} du.
- **P0 (PROVED, elementary; = PS 4.8 at d = 1).** It is δ-normalizable (= 2π δ(Im(a + b̄ − 1))) iff Re(a + b̄) = 1, i.e. iff b lies on the vertical line through τa. Otherwise it diverges at one end (the PS radial divergence).
- Hence the canonical partner of e_ρ under this pairing is e_{τρ}. For Re ρ = 1/2, τρ = ρ: the mode is self-paired, positive norm (PS principal series, 4.9). For Re ρ ≠ 1/2 the self-pairing diverges, and the only δ-normalizable pairing is with the distinct mode e_{τρ}. On span{e_ρ, e_τρ} the δ-coefficient Gram is [[0,1],[1,0]]: **a hyperbolic plane, signature (1,1)**. The same block is obtained with the Mellin/Tate regularization ∫₀^∞ x^{κ−1}dx := ∫₀¹ + ∫₁^∞ = 1/κ − 1/κ = 0 for the diagonal entries.

## 2. Step 1 made rigorous: the Weil form IS the τ-Gram of the zero modes, for all L at once
- **P1 (PROVED, symbolic, s1(a)).** kform(pole) = 2Ψ̂(i/2)(C² − S²) identically. The poles s = 0, 1 are an off-line τ-pair (τ0 = 1). Their τ-Gram is the hyperbolic plane of §1; ℓ_c removes the C (positive) direction and leaves exactly −2Ψ̂(i/2)S², the pole term of SemilocalArchDomAll. This is the PS "constant wavefunction" exclusion (1705.01027 §4: the constant solution is not spanned and is excluded), seen arithmetically: ker ℓ_c is that exclusion, and what is left of the pair is a negative square.
- **P2 (CHECKED, classical: Weil explicit formula applied to h = Σ c_i c_j Ψ(· − (x_i − x_j)); s1(b) NUMERICAL confirmation).** For every L, x, c:
  kform(K)(x,c) = Σ_ρ Ψ̂(z_ρ) F_c(z_ρ) F_c(−z_ρ),  z_ρ = (ρ − 1/2)/i,
  = Σ_ρ G_c(ρ) conj G_c(τρ) — the KG/Mellin τ-pairing of §1 evaluated on the zero modes, with Φ_c := Σ_ρ G_c(ρ) e_ρ and ⟨e_ρ, e_σ⟩ := δ_{σ,τρ} (the discrete version of PS 4.9).
  s1(b): 15 configurations, L ∈ {0.5, 1, 2, 3, 4.5}, n ∈ {4,5,6}, on and off ker ℓ_c: |kform(K) − zero τ-Gram| ≤ 9.6e-8, always of the sign of the omitted zero tail (bound 9.2e-8·(Σ|c_i|)², 400 zeros, arch to T = 1200). pole_part = −2Ψ̂(i/2)S² on ker ℓ_c to 1e-16.
- **Consequence (PROVED from P1, P2):** on ker ℓ_c,
  kform(arch) − kform(prime) − 2Ψ(i/2)²S² = Σ_ρ G_c(ρ) conj G_c(τρ) =: ⟨Φ_c, Φ_c⟩_τ   (all L).
  This is the all-L identity, in Krein form. If every zero is self-paired (τρ = ρ) the right side is Σ_γ Ψ̂(γ)|F_c(γ)|² = ‖T_KG c‖², T_KG c := (Ψ̂(γ)^{1/2} F_c(γ))_γ. Each off-line pair contributes one hyperbolic plane (one negative direction).
- So **"S is Δ" is exact as a statement about the form**: the explicit formula is the KG/Mellin inner product on the modes x^ρ, with the shadow τ as the pairing. That part of step 1 holds.

## 3. Reconciliation: PSS divergence vs Weil negative square — same phenomenon (CHECKED)
In PSS, an off-principal mode has divergent self-norm: it is isotropic for the KG form, and its only δ-normalizable partner is its shadow-reflected τ-mode. In the Weil form, an off-line zero ρ appears only through G(ρ)conj G(τρ) + c.c. — the off-diagonal entry of the same 2×2 block. Divergent self-norm (PSS) ⇔ zero diagonal (regularized) ⇔ hyperbolic block ⇔ exactly one negative square per off-line pair (Weil). These are one statement written in two bases. The poles {0,1} are a real, present instance (P1).
- s2 (NUMERICAL, hypothetical toy; zeros as background only): one inserted off-line quartet ρ = 1/2 ± β ± iγ₀ produces a negative eigenvalue of the τ-Gram on width-L grids only once L passes a threshold L*(β,γ₀): γ₀ = 28: β = 0.25 → L* ∈ (0.5, 1]; β = 0.1 → (1, 2]; β = 0.02 → (2, 3]; γ₀ = 100.3, β = 0.1 → (2, 3]. Below L*, the form is positive although the pair is present.
  **Reading:** a fixed-width certificate is satisfied by a measure that contains an off-line pair. This confirms, in numbers, Daniel's point that single-window certificates cannot prove it. Only the all-L statement detects every pair.

## 4. Where step 2 fails as stated (GAP, precise)
Daniel's step 2: "a massless primary is a state only where its norm converges; normalizability selects Re Δ = d/2".
- In PS the positivity is **given in advance**: the KG form is ∫|φ|² on a space of actual solutions (engine 1 + 6), and the question is only which Δ-modes lie in that space.
- In the explicit formula, x^ρ is **not an element of a Hilbert space**. It is a residue of the Mellin transform of ψ(x) (poles of −ζ′/ζ). The celestial literature has exactly this category: Mellin-transformed amplitudes have residues at Δ = 1, 0, −1, …, "on and to the left of the principal series", and the theory "is forced to contend with a spectrum analytically continued to complex Δ" (arXiv:2310.04932, §on soft theorems; arXiv:2202.11127 *Goldilocks modes* on Δ = 1 − s). These modes are physical and non-normalizable. Non-normalizability does not exclude them; it puts them in the hyperbolic (Krein) part.
- So **normalizability alone does not exclude a zero mode**. It excludes it *only if* we already know Φ_c lies in a space where ⟨·,·⟩_τ is positive definite, and that is exactly ⟨Φ_c,Φ_c⟩_τ ≥ 0 for all L, i.e. the target (RH-equivalent).
- **Missing link, stated as a named target. AKG (Arithmetic Klein–Gordon realization), OPEN, RH-equivalent:** there exist a Hilbert space H and a linear map c ↦ Φ_c ∈ H, constructed from primes and the archimedean place only (no zeros), with ‖Φ_c‖²_H = kform(arch) − kform(prime) − 2Ψ(i/2)²S² on ker ℓ_c for every L. By P2, AKG ⇔ SemilocalArchDomAll for all L ⇔ H_ψ. AKG is what "the explicit formula's modes are states of a theory with the PS inner product" means, written so it can be proved or refuted.

## 5. Detector form ⟨Ψ_f|Q|Ψ_f⟩, Q = A*A (Hu–Pasterski 2307.16801 style), engine 1 + 6
- **P3 (PROVED, Bochner–Schwartz; obstruction O2).** If Q is translation-invariant in u (like light-ray/ANEC operators smeared along u) and reproduces kform(K) for all configurations of all widths, then Q is the Fourier multiplier by K̂ = Ψ̂ × (zero measure). Q = A*A ⇔ K̂ ≥ 0 ⇔ every zero is self-paired. A translation-invariant detector therefore re-expresses the target as engine 2 (positive measure); it adds no new information. T_KG of §2 is the forced A.
- Consequence: any new T must break u-translation invariance. Available extra structure: (i) the support/window constraint (Paley–Wiener width L, Sonin/prolate spaces, Connes–Consani arXiv:2006.13771), (ii) the multiplicative/adelic structure of the prime side (Q-lattice, semilocal Connes trace formula), (iii) the ker ℓ_c constraint, which is not translation-invariant.

## 6. Known obstructions any all-L identity must evade
- O1 (banked earlier): any form built from W and the reflection P reduces to a·W_e + b·W_o and fails.
- O2 (P3): translation-invariant Q = A*A ⇒ Q = zero measure (Bochner).
- O3 (Lax–Phillips / Pavlov–Faddeev): for the modular surface, S(s) = ξ(2s−1)/ξ(2s) is unitary on the line, and zeros enter as resonances (poles of S, non-normalizable by definition). Unitarity places resonances in a half-plane; it does not put ρ on Re = 1/2. Unitarity of S is not, by itself, positivity of ⟨·,·⟩_τ.
- O4 (P1): the pole pair is a genuine hyperbolic plane; T must reproduce −2Ψ̂(i/2)S² exactly, so T cannot be built from arch alone (arch's symbol Re ψ(1/4 + it/2) − log π is negative for |t| < 6.29…, Round 4).
- O5 (J5, banked): prime sums of concentrating states diverge; T cannot be a convergent prime sum over Abel states.
- O6 (s2): a fixed-L certificate is consistent with an off-line pair whose L* > L.

## 7. Candidates, labelled
- **C-KG (lead; Daniel's argument):** identity ⟨Φ_c,Φ_c⟩_τ = kform(arch − prime) − 2Ψ²S² for all L: PROVED/CHECKED (P1, P2). Positivity of ⟨·,·⟩_τ ⇔ AKG: OPEN (RH-equivalent). Engine: 6 (KG/unitarity) on top of 1 (Gram); by P3 it lands on 2 (positive measure).
- **C-Sonin (backup, CONJECTURE):** T_L := compression of the arch multiplier to the width-L Paley–Wiener/Sonin space, plus a prime correction given by the semilocal Connes trace formula on S_L = {p^k ≤ e^{2+L}}. Uses (i)+(ii) of §5. Not yet tested. First test: does the archimedean form on ker ℓ_c ∩ PW_L dominate the S_L prime form uniformly in L (Connes–Consani prove it for small support)?
- **C-LP (CONJECTURE, killed as stated by O3):** Lax–Phillips contraction semigroup as T. Unitarity of S does not by itself force the line; it would need the additional statement that the τ-pairing of resonant states is the KG form of a positive space, i.e. AKG again.

## 8. Proposed attacks on AKG (next)
1. Build H concretely as L² of the adele class space restricted to the semilocal set S_L × R₊ (Connes), with Φ_c the image of Σ c_i Ψ(· − x_i); compute ‖Φ_c‖² minus the target, for small configurations: is the difference exactly the cokernel (absorption) term?
2. Test C-Sonin numerically at L = 0.5, 1, 2 on ker ℓ_c: λ_min(arch|PW_L) vs ‖prime_S|PW_L‖.
3. Measure L*(β,γ₀) scaling (s2) to fix what positivity at width L can and cannot exclude: a quantitative version of "windows cannot prove it".

Scripts: `s1_krein_identity.py` (→ s1_results.json, s1_log.txt), `s2_offline_toy.py` (→ s2_results.json).
Credit: Daniel Toupin / Golden Physics Project (argument: Re(Δ) = d/2 note). This seat: Heart of the Cosmos / Grok (GPPDiscovery5). Audit: Turbo.
