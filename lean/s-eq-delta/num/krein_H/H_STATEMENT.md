# H for ζ, triangle ψ — exact statement (Turbo PR #11 6821018808, item 5)

Heart of the Cosmos (Grok, GPPDiscovery5), 2026-11-12. H, ARP, RH: **OPEN**. No proof claim here.

## Data (every argument of `GppKreinBridge.HAllWindows B S lc ls`)

* **ψ** (triangle): `ψ(t) = max(0, 1 − |t|)`, supp ψ = [−1, 1] (so `a = 1` in `fit_translates`),
  `Ψ(z) = ∫ ψ(t) e^{izt} dt = (sin(z/2)/(z/2))²`. `τ_x ψ := ψ(· − x)`.
* **F** := `F_ψ = span_ℝ { τ_x ψ : x ∈ ℝ }` (real finite combinations `v = Σ cᵢ τ_{xᵢ} ψ`). This is the
  minimal space `bridge` consumes. (Any larger real space of compactly supported Lipschitz
  functions gives a stronger H that implies this one.)
* **S L** := `F_ψ ∩ window L`, `window L = {f : f(t) = 0 for |t| > L}` (KreinBridge.window).
* **lc L, ls L** (independent of L):
  `ℓ_c(v) = ∫ v(t) cosh(t/2) dt`, `ℓ_s(v) = ∫ v(t) sinh(t/2) dt`.
  On `F_ψ`: `ℓ_c(v) = Ψ(i/2) Σ cᵢ cosh(xᵢ/2)`, `ℓ_s(v) = Ψ(i/2) Σ cᵢ sinh(xᵢ/2)`,
  `Ψ(i/2) = Ψ(−i/2) = 16 sinh²(1/4) = 1.0210077216…`.
* **B** (Weil form, defined zero-free; real, symmetric, bilinear), with `(u⋆w)(y) = ∫ u(x+y) w(x) dx`,
  `U(z) = ∫ u e^{izx}`:

      B(u,w) = 2 ℓ_c(u)ℓ_c(w) − 2 ℓ_s(u)ℓ_s(w)
               − Σ_{n≥2} Λ(n) n^{−1/2} [ (u⋆w)(log n) + (u⋆w)(−log n) ]
               + (1/2π) ∫_ℝ U(t) W(−t) ( Re ψ_Γ(1/4 + it/2) − log π ) dt,

  ψ_Γ = Γ'/Γ. So `B = A − P + 2ℓ_c² − 2ℓ_s²` in the ToeplitzSewing pieces
  (A = archimedean integral, P = prime sum). On `F_ψ`, `B(τ_xψ, τ_yψ) = K(x − y)` with the prime sum finite
  (`u⋆w` supported in `[x−y−2, x−y+2]`) and the archimedean integrand `O(t^{−4} log t)`.

## H (the statement attacked)

    H_ψ :≡ HAllWindows B S lc ls
        ⟺ ∀ L, ∀ v ∈ F_ψ ∩ window L, ℓ_c(v) = 0 → 0 ≤ B(v, v)
        ⟺ ∀ n, ∀ x ∈ ℝⁿ, ∀ c ∈ ℝⁿ, Σ cᵢ cosh(xᵢ/2) = 0 → Σᵢⱼ cᵢ cⱼ K(xᵢ − xⱼ) ≥ 0,

  K(d) := B(τ_0ψ, τ_dψ) = 2Ψ(i/2)² cosh(d/2) − Σ Λ(n)n^{−1/2}[g(log n + d) + g(−log n + d)]
          + (1/π) ∫_0^∞ Ψ(t)² cos(td)(Re ψ_Γ(1/4+it/2) − log π) dt,
  g = ψ∗ψ: g(y) = 2/3 − y² + |y|³/2 (|y| ≤ 1), (2 − |y|)³/6 (1 ≤ |y| ≤ 2), 0 else.

  First ⟺: the `2ℓ_s²` terms cancel identically (claim_iff_hyperplane, PROVED earlier). Second ⟺: every v ∈ F_ψ
  lies in some window (fit_translates, PROVED), so the window quantifier is vacuous on F_ψ. H_ψ is a statement
  about primes and Γ only — no zeros.

## Named premises of the chain H_ψ ⇒ RH (labels unchanged)

* **hEF** (Weil's explicit formula on translates): `K(d) = Σ_ρ m_ρ Ψ(γ_ρ)Ψ(−γ_ρ) e^{iγ_ρ d}`, ρ = 1/2 + iγ.
  CHECKED (secondary: Connes arXiv:1509.05576 §2.1.1 citing Weil). **New: NUMERICAL** — prime/Γ side vs zero side
  (first 400 zeros, T = 679.74) agree at d ∈ {0, .25, .5, 1, 1.5, 2, 3, 4, 5, 6} to ≤ 1.02e−8 (float64; zero-side
  tail bound ≈ 5e−7). `hef_tri.py`, `hef_tri.json`. This also fixes the normalisation of B above.
* **hfit**: PROVED (fit_translates, a = 1).
* bridge ⇒ `AtMostOneNegSq K`: PROVED mod hEF (KreinBridge.bridge with F = F_ψ).
* Krein: CHECKED (secondary, Stewart 1972); primary Krein–Langer, Math. Nachr. 77 (1977) 187–236: **owed**.
* (L), (G1), hPT: CHECKED on paper. Assembly ⇒ every zero real: PROVED mod hPT. Chain H_ψ ⇒ RH: SKETCH.

## Cheapest falsifier (item 5.2)

* A kill is a window L and v ∈ F_ψ ∩ window L with Σ cᵢ cosh(xᵢ/2) = 0 and Σ cᵢcⱼK(xᵢ−xⱼ) < 0. By hEF, if all
  zeros are on the line then `B(v,v) = Σ_ρ m_ρ |V(γ_ρ)|² ≥ 0`, so **a kill of H_ψ is an ARP failure and forces an
  off-line zero** (paper, from hEF; PROVED earlier as arp_failure_of_not_claim in the abstract form).
* Search run (zero-free side only, so not circular): grids xᵢ = x₀ + (i − (n−1)/2)/4, n = 2…25, shifts
  x₀ ∈ {0, 1/4, …, 3}, windows L ≤ 7; min of B on ker ℓ_c ∩ span via projected eigenvalues (`scan2.py`, `scan2.json`).
  **No kill.** Minimum over all runs = 1.55e−7 > 0 (n = 25, L = 7); for n ≤ 13 (L ≤ 5.5) the minimum is ≥ 7.9e−6,
  far above the estimated K error (~1e−8 per entry). NUMERICAL, float64, not certified. The full unconstrained
  Toeplitz matrices are also positive definite on every grid (smallest 1.5e−7), so these families show no
  negative square at all; the ker ℓ_c constraint is not doing the work.
* Sensitivity (SKETCH estimate): RH is verified for |t| ≤ 3·10¹² (Platt–Trudgian 2021), so any off-line ρ has
  |Ψ(γ)Ψ(−γ)| ≤ (2cosh(1/4)/|t|)⁴ ≈ 2·10⁻⁴⁹, times e^{δ·span}, δ < 1/2. A float64 scan cannot see it unless the
  span is ≳ 180 (L ≳ 90). The cheapest falsifier of H_ψ is therefore blind at accessible windows; it checks the
  statement and its normalisation, not H.
