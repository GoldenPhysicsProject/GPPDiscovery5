# Round 22: attacking GAP-20, synthesizing the passive ζ network from primes + the archimedean place

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).

**RH is OPEN, and nothing here claims it.**
- Labels: PROVED / CHECKED / NUMERICAL / CONJECTURE / HEURISTIC / GAP / OPEN. "Cited" means a published theorem stated but not re-derived.
- **No zero is a proof input.** Zeros appear only in two places: (i) code validation, where the prime-side and zero-side are compared, and (ii) illustrations. Both are marked.
- Frame: u = (s−½)², with m(u) = F′/F = (1/(2√u)) ξ′/ξ(½+√u). RH ⇔ m is Stieltjes (Shadow Euler paper; Round 20).
- Round 21 (`round21_entangled/ROUND21.md`) was read. Its TFD/Schmidt link says a positive ladder is a TFD purification with positive Schmidt weights at real energies γ². It also notes that the Bost–Connes TFD is positive only for β > 1. Section 4 shows that this is the same frontier as everything else here.

## Short answer

1. **Positivity by construction, piece by piece, is impossible (PROVED).**
   - Any finite family of prime pieces contributes a *mean-zero almost-periodic* density to m's spectral measure, so it can never be ≥ 0 on its own.
   - The archimedean piece alone is *not* positive: θ′(t) < 0 on (0, 6.2898), and the arch-only network fails.
   - Γ-completed finite Euler products put poles of F on u > 0.
   - So there is no regrouping into passive per-prime sections plus a passive arch section. Positivity can only live in the *whole sum*.
2. **The construction that does work (the "Fejér network").**
   - Smooth the explicit formula with the Fejér kernel at scale 1/log X. This gives, for each X, an explicit Stieltjes function m_X built from **primes ≤ X + Γ + pole only**.
   - m_X is passive ⇔ a one-variable density ρ_X(γ) ≥ 0, and m_X → m.
   - RH ⇔ m_X is passive for every X (PROVED ⇒; sketch ⇐).
   - **NUMERICAL headline:** ζ's network is passive at *every* truncation tested, X = 10, 30, 100, …, 10⁶ (γ ≤ 200). The fakes break: F_θ at X = 30, Epstein at X = 100 (near γ ≈ 14), Davenport–Heilbronn at X = 100 (near γ ≈ 84).
   - Compare plain truncation, which can never be passive (pole at u = ¼) and whose ladder breaks at a₂₁ even with P = 10⁴.
   - **But passivity is *checked*, not built in.** It is still arch + pole against primes: Weil positivity restricted to Fejér test functions.
3. **The unconditional positive network exists, but only up to the Euler-product line.**
   - By de Branges, E_h(z) = Ξ(z + ih) is Hermite–Biehler (a lossless passive scattering network) for every h ≥ ½, unconditionally, *because of the Euler product*.
   - RH ⇔ the same for every h > 0. The positivity from primes stops exactly at Re s = 1. That is the same line as R20's convergence region u > ¼ and R21's Bost–Connes Hagedorn point β = 1.
   - F_θ (Euler product) satisfies the h ≥ ½ condition. EPS and DH (no Euler product) fail it: they have zeros in Re s > 1 (cited).
4. **GAP-22 (sharpened GAP-20).** Prove the single explicit inequality
   arch_X(γ) + pole_X(γ) ≥ (2/π) Σ_{n≤X} Λ(n) n^{−½} (1 − log n/log X) cos(γ log n) for all γ and all X.
   Equivalently, push de Branges positivity from h = ½ down to h = 0. The first step, any h < ½ uniformly, is already the open zero-free-strip problem.

Files:
- `r22_models.py`: explicit-formula data and Λ_M tables.
- `r22_fejer.py` → `r22_fejer_{ZETA,FTH,EPS,DH}.json`, `r22_fejer_ZETA_ext.json`, `fejer_*.log`.
- `r22_fejer_validation.json`.
- `r22_jacobi.py` → `r22_jacobi.json`.
- `r22_debranges.py` → `r22_debranges_*.json`, `db_*.log`.

---

## 1. Network composition: what cannot work (PROVED)

**1a. Single prime-power pieces are not passive.** In the region u > ¼ the explicit formula gives

  m(u) = (1/(2√u)) [ 1/s + 1/(s−1) − ½ log π + ½ψ(s/2) − Σ_n Λ(n) n^{−s} ],  s = ½ + √u.

Each prime-power term is c·e^{−L√u}/(2√u) with L = log n. That is completely monotone in u (a Laplace transform of a positive heat kernel), but it is not Stieltjes. Its boundary density on u = −γ² is ∝ cos(γL)/γ, which changes sign.

The full Euler-factor piece (log p) p^{−s}/(1 − p^{−s}) has density ∝ Re[q e^{−iφ}/(1 − q e^{−iφ})], with q = p^{−½} and φ = γ log p. This equals +q/(1−q) at φ = 0 and −q/(1+q) at φ = π. So neither the piece nor its negative is passive.

**1b. Any finite family of primes is mean-zero.** Let D(γ) = Σ_{n≤N} c_n cos(γ log n). It is a Bohr almost-periodic function with mean value 0. A nonnegative almost-periodic function with mean 0 is identically 0. Therefore **no finite set of prime pieces, however regrouped, is passive on its own.** This holds for Euler-factor sections, Möbius-inverted forms (μ has mixed signs anyway), shadow pairings p^{−s} with p^{−(1−s)}, and parallel sums.

On parallel sums specifically: log-derivatives are additive over Euler factors, while a parallel sum is harmonic-additive, so the two operations do not match. Positivity of m is a property of the infinite sum only.

**1c. The archimedean piece is not passive either (PROVED, with values CHECKED).** The archimedean density of states is θ′(t)/π, where θ′(t) = ½ Re ψ(¼ + it/2) − ½ log π.
- θ′(0) = −2.686.
- θ′ < 0 on (0, 6.28984) (the root was computed), and θ_min = −3.53097.

So the "Γ-factor string" has negative mass at low γ, and the primes must *supply* positive mass there. There is no "positive arch string + sign-controlled prime perturbation" decomposition. Other measurements agree:
- Round 20's arch-only ladder fails at a₅.
- The arch + pole Fejér network with no primes has min ρ = −0.38, −1.40, −4.69, −15.0 at X = 10, 10², 10³, 10⁴ (`r22_jacobi.json`).

**1d. Γ-completed truncated Euler products are never passive.** Consider G(½+w)P_X(½+w) + G(½−w)P_X(½−w), or Gonek-type χ-completions. In u these have Γ/χ poles at u = (2k+½)² > 0. Poles of F give m negative-residue poles on (0, ∞), so m is not Stieltjes for any X. Gonek's ζ_X (Trans. AMS 364, 2012; cited) satisfies RH except for finitely many zeros, but it is *non-analytic*, so it has no Stieltjes or network meaning.

**1e. Plain truncation** (Round 20) keeps the 1/(s−1) pole at u = ¼ > 0, so it is never passive for finite P. Its ladder breaks at a₂₁ with P = 10⁴.

## 2. The Fejér network: a construction that is passive at every tested truncation

**Definition.** Fix L = log X, and let K(r) = (L/2π) sinc²(rL/2) (K ≥ 0, ∫K = 1, Fourier transform supported in [−L, L]). For γ₀ ≥ 0 put h(r) = K(r − γ₀) + K(r + γ₀). The explicit formula, derived in the `r22_models.py` docstring by the contour + FE argument, gives

  ρ_X(γ₀) := Σ_ρ h(γ_ρ) = 2 Σ_{poles z∈(½,c)} h(−i(z−½)) + (1/π) ∫ h(r)[log Q + Re G′/G(½+ir)] dr − (2/π) Σ_{n≤X} Λ_M(n) n^{−½} (1 − log n/L) cos(γ₀ log n).

This is a **finite** expression in the primes ≤ X. The archimedean integral is evaluated exactly through ψ(z) = −γ_E + ∫(e^{−t} − e^{−zt})/(1−e^{−t}) dt, which has compact support in t. The network is

  m_X(u) := ∫₀^∞ ρ_X(γ) dγ / (u + γ²).

**Theorem C.**
- (i) m_X is Stieltjes ⇔ ρ_X ≥ 0 (PROVED; ρ_X grows only like log γ, so the integral converges).
- (ii) ρ_X → the zero measure as X → ∞, so m_X → m (PROVED, distributionally).
- (iii) RH ⇒ ρ_X = Σ_γ (Fejér bumps) ≥ 0 for every X (PROVED).
- (iv) Converse, sketch, not Lean-checked:
  - Take an off-line zero ½ + δ + iγ* with maximal δ in a window.
  - It contributes 2 Re K(γ* − γ₀ − iδ). This oscillates in γ₀ with period 2π/L and amplitude ≍ X^δ/(L(γ₀−γ*)²).
  - The real zeros contribute O(L log γ₀).
  - So ρ_X < 0 somewhere once X is large. Hence RH ⇔ every m_X is passive.
  - This is Weil's criterion restricted to Fejér (box ∗ box) test functions.

**Code validation (CHECKED; zeros used only here).** For ζ the prime side matches the zero side summed over the first 300 zeros. The difference is a constant offset of 9.0·10⁻⁴, 4.6·10⁻⁴, 3.0·10⁻⁴ (X = 10, 100, 1000), which is the positive tail from zeros above 300 (it decays like 1/L). See `r22_fejer_validation.json`.

**Numbers (NUMERICAL).** min ρ_X over a γ grid (step 0.01):

| X | ζ (γ ≤ 100; 200 in ext run) | F_θ (γ ≤ 60) | EPS (γ ≤ 60) | DH (γ ≤ 100) |
|---|---|---|---|---|
| 10 | +0.0106 | +0.0181 | +0.302 | +0.0444 |
| 30 | +0.0066 | **−0.025** (first neg. γ 12.2) | +0.066 | +0.0270 |
| 100 | +0.0047 | −0.221 | **−0.374** (first neg. γ 14.25) | **−0.134** (first neg. γ 84.3) |
| 300 | +0.0039 | −0.350 | −1.07 | −0.364 |
| 10³ | +0.0034 | −0.731 | −2.35 | −0.875 |
| 3·10³ | +0.0024 | −1.10 | −4.47 | −1.52 |
| 10⁴ | +0.0028 | −1.65 | −8.42 | −2.72 |
| 10⁵ (ext) | +0.0022 | | | |
| 10⁶ (ext) | +0.0017 | | | |

What the table shows:
- **ζ is passive at every truncation, from 10 to 10⁶ primes-cutoff.**
- The fakes fail at a definite truncation, and once failed, never recover. Their negativity grows like X^δ, consistent with (iv): δ = ¼ (F_θ), 0.433 (EPS), 0.31 (DH).
- The places where fakes go negative sit next to their off-line zeros (illustration only: EPS ghost 15.67, DH ghost ≈ 85.7).
- ζ's margin shrinks slowly, about 1/L (0.0106 → 0.0017). ζ is passive *at the edge*, the same criticality as Λ_dBN = 0 under RH.

**Honest classification.** This is the best synthesis found. It is explicit, finite in the primes, uses only arch + pole + primes, is passive at every tested truncation for ζ, and fails for every fake. But by §1b it is **not** "positive by construction": at each X positivity results from arch + pole beating the prime cosine sum. What the construction achieves is a *reduction*. Passivity of the whole network becomes one explicit real inequality per X (GAP-22), instead of a ladder depth.

## 3. Spectral / Jacobi / Krein route (Round 20 data, no zeros)

- **Jacobi matrix.** From Round 20's 198 S-fraction coefficients (Taylor data at ½): α₀ = a₁, α_n = a_{2n} + a_{2n+1}, β_n = a_{2n−1}a_{2n}. This is a 98×98 Jacobi matrix.
  - α: 1.6089e−3, 3.5006e−3, 1.7983e−3, 1.2172e−3, …
  - β: 3.6515e−6, 1.3212e−6, 4.8622e−7, …
  - log β_n has slope ≈ −2.59 in log n (n ≥ 10).
  - **No arithmetic pattern is visible in the entries** (NUMERICAL/HEURISTIC). Like the ladder elements of R20, the entries are rational in the Hankel determinants and mix all primes.
- **The circuit's resonances are the zeros (illustration).** The Gauss nodes of this Jacobi matrix, 1/√(eigenvalue), are 14.134725141734695, 21.022039638771556, 25.010857580145690, 30.424876125859512, …. They agree with ζ's first 25 zero ordinates to 10⁻¹⁶ relative error, computed from ξ's Taylor coefficients alone. So the Krein string synthesized from Taylor data *is* the zero spectrum. This is a statement about spectrum, not about positivity: a positive Jacobi matrix arises only from moments known to be a positive sequence, which is RH again.
- **Archimedean string.** §1c: the Γ-factor density θ′/π is negative on (0, 6.29). The Γ-only object is not a positive string, so "primes as a rank-one or Schur-complement positive update of a positive arch string" has no positive base to start from. The primes must create positive low-energy mass, which is exactly the cancellation in §2.

## 4. de Branges / Lagarias check

**Theorem D (PROVED, elementary).** Let E_h(z) := Ξ(z + ih) = ξ(½ − h + iz).
- E_h is Hermite–Biehler (|E_h(z̄)| < |E_h(z)| for Im z > 0) whenever h > sup_ρ |β − ½|. Proof: Hadamard product Ξ(z) = Ξ(0) Π(1 − z²/ρ_Ξ²); each paired factor has its zeros at Im = ±Im ρ_Ξ − h < 0, so it is HB; the product converges absolutely.
- If some zero has |β − ½| > h, then E_h has a zero in the upper half-plane, so it is not HB.
- Hence RH ⇔ E_h is HB for every h > 0. By FE and real coefficients, HB is equivalent to |ξ(½ + a + ix)| < |ξ(½ + b + ix)| for |a| < b: monotonicity of |ξ| away from the line (cited: Matiyasevich–Saidak–Zvengrowski).
- **For ζ, h ≥ ½ holds unconditionally, because the Euler product makes ζ(s) ≠ 0 on Re s ≥ 1** (Hadamard–de la Vallée Poussin). The canonical system / lossless network of H(E_h) then exists with a positive Hamiltonian (de Branges' theorem; existence, not explicit).

**Numbers (NUMERICAL, `r22_debranges_*.json`).** Grid: σ step 0.05, x step 0.1. HB violated?

| h | 0.6 | 0.5 | 0.45 | 0.4 | 0.35 | 0.3 | 0.25 | 0.2 | ≤ 0.15 |
|---|---|---|---|---|---|---|---|---|---|
| ζ (x ≤ 100) | no | no | no | no | no | no | no | no | no (down to 0.05) |
| F_θ (x ≤ 60) | no | no | no | no | no | no | no | **yes** (x = 48) | yes |
| EPS (x ≤ 60) | no | no | no | **yes** (x = 30) | yes | yes | yes | yes | yes |
| DH (x ≤ 100) | no | no | no | no | no | **yes** (x = 85.7) | yes | yes | yes |

Observed frontiers in these windows:
- F_θ: h₀ ∈ (0.20, 0.25], and the true value is ¼ (zeros at ρ ± ¼).
- EPS: h₀ ∈ (0.40, 0.45], consistent with the Round-17 ghost δ = 0.433.
- DH: h₀ ∈ (0.30, 0.35], consistent with δ ≈ 0.3085 near 85.7.

EPS and DH also have zeros in Re s > 1 higher up: Davenport–Heilbronn 1936 for Epstein zeta with class number > 1; Titchmarsh §10.25 for the DH function (cited). So for them even h = ½ fails globally. **With no Euler product there is no unconditional de Branges positivity at all.** F_θ has an Euler product, but only for Re s > 5/4. Its Euler-product guarantee is h ≥ ¾, and its true frontier is ¼.

**Relation to the other rounds (PROVED as equivalences).** All of the following are the same statement:
- E_h is HB for all h > 0;
- m is Stieltjes;
- m_X is passive for all X;
- the R20 ladder is passive;
- RH.

The *unconditional* positivity that the primes provide lives exactly at h ≥ ½, i.e. Re s ≥ 1:
- the region where R20's prime sums converge (u > ¼);
- the region where R21's Bost–Connes TFD is a positive state (β > 1, the Hagedorn point at β = 1).

**Primes give positivity for free down to the Hagedorn line Re s = 1, and not one step further.**

## 5. Verdict, in plain words

- We tried to build ζ's passive circuit out of primes so that positivity is automatic. **It cannot be done piece by piece.** Every finite batch of primes, however grouped, contributes an oscillation that averages to zero, and the Γ-factor part on its own has negative mass at low frequency. Positivity can only be a property of the whole sum.
- What *does* work is the **Fejér network**. Smoothing the explicit formula at resolution 1/log X gives, for every X, an explicit circuit made of the primes up to X plus the Γ-factor and the pole. For ζ that circuit is passive at every truncation we tried, up to X = 10⁶, while F_θ, Epstein and Davenport–Heilbronn break at X = 30, 100 and 100. That is a clean separation, and a sharp reduction of GAP-20 to one explicit inequality per X.
- But it is still a cancellation (Weil positivity with Fejér test functions), not positivity by construction.
- The only *unconditional* positive network, de Branges', comes from the Euler product and reaches exactly the line Re s = 1. That is the same frontier as Round 20's convergence region and Round 21's Hagedorn point. To get to the critical line you must cross half the strip with no structural help. Even crossing *any* fixed fraction of it is the open zero-free-strip problem.

**GAP-22.** Prove, for all X ≥ 2 and all γ ≥ 0,

  ρ_X(γ) = [Γ-term + pole-term]_X(γ) − (2/π) Σ_{n≤X} Λ(n) n^{−½}(1 − log n/log X) cos(γ log n) ≥ 0.

Equivalently (Theorem D), push de Branges positivity from h = ½ to h = 0. A first, still-open milestone: any uniform h < ½, i.e. a zero-free strip. OPEN.

## 6. Status ledger

- **PROVED:**
  - 1a–1e: no passive per-prime or arch pieces; mean-zero almost-periodic obstruction; Γ-completed Euler products have poles at u > 0; plain truncation has the u = ¼ pole.
  - Theorem C (i)–(iii), and (iv) as a sketch.
  - Theorem D and the h ≥ ½ unconditional statement for ζ.
  - The equivalences in §4.
- **CHECKED:**
  - Fejér prime side = zero side, up to the explained tail.
  - θ′ root 6.28984, θ_min −3.53097.
  - EPS genus identity, to 3·10⁻⁸ (lattice truncation).
- **NUMERICAL:**
  - ζ Fejér network passive for X = 10 … 10⁶ (γ ≤ 100/200); fakes fail at X = 30 / 100 / 100.
  - de Branges frontiers: ζ none down to h = 0.05 (x ≤ 100); F_θ (0.20, 0.25]; EPS (0.40, 0.45]; DH (0.30, 0.35].
  - Jacobi data; Gauss nodes = zeros (illustration).
- **HEURISTIC:** no arithmetic pattern in the Jacobi/ladder entries; "ζ passive at the edge" (margin ~ 1/L).
- **GAP-22:** as above. **OPEN:** RH.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).
