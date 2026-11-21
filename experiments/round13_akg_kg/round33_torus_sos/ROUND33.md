# Round 33 — the torus sum-of-squares (Heart of the Cosmos / D5)

Daniel's go: "Go" (2026-11-21). Labels: PROVED / CHECKED / NUMERICAL / CONJECTURE / HEURISTIC / GAP / OPEN.
No zeta zeros are proof inputs (300 zeros appear only as test data in `r33_validate_zeros.py`). Weil positivity is an
output, never the target. **RH is OPEN.**

## 0. Verdict

- **The SOS exists and is exact (PROVED, CHECKED).** On the prime torus T = ∏_p S¹, push the Fejér kernel forward along the
  Kronecker line t ↦ (p^{−it})_p to get a positive measure ν_L. Then:
  - the R32 Möbius floor is the inner product ⟨f_c, 1⟩_{ν_L};
  - a one-parameter square gives **(m(L) − cL)² ≤ L·Q_c(L)**, where Q_c(L) = ∫₀^L |M_{½}(e^u) − c|² du and c = 1/b(½);
  - **the √X of Prop 32.3 is gone from the bound.** It has been moved into a second-moment quantity Q_c.
- **The cost is a second moment.** For ζ, Q_c grows linearly with slope 0.0291, flat to X = 10⁹. Each control fails by its own mechanism:
  - F_θ, EPS and DH grow exponentially;
  - ζ² grows polynomially (double zeros);
  - Beurling is PROVED to fail eventually, but it is not separated in range.
  - **New:** DH, which R32's pointwise floor could not separate, is separated here (Q_c/L goes 0.07 → 0.71).
- **The counted-once part is the Haar diagonal**, and it is polylog for ζ (PROVED). The work is done entirely by the
  **off-diagonal remainder R_c**, the Kronecker-line (archimedean) correlation. Beurling has the same Haar diagonal and
  must fail in R_c.
- **Equivalences.**
  - The ε-form GAP-33 is Littlewood's criterion in L²/Fejér form, and is **equivalent to RH**.
  - The linear form is the σ = ½ weak Mertens hypothesis, which gives RH + simple zeros.
  - Attacking R_c directly needs two-point (log-)Chowla correlations Σ μ(n)μ(n+h)/n with square-root
    cancellation *uniformly in all shifts h < X*. The Montgomery–Vaughan / Hilbert bound is sharp for general
    coefficients and gives back exactly √X.
- **ψ side (PROVED identity):** D(L) = κ₀ − ⟨f_Λ, 1⟩_{ν_L} + O(e^{−L/2}), with κ₀ = −8 + ψ′(¼)/4 − ζ′(½) = 0.2219784278.
  - The same square gives only the two-sided |D − κ₀| ≤ √(L·Q_Λ). A Cauchy–Schwarz square cannot see the sign of D; that sign is the Weil functional at t = 0.
  - Q_Λ's slope (0.04623 ± 0.0014) equals R₀ = (log ξ)″(½) = 0.0462100, which is computed from integers and Γ only. Under RH that is the simple-zero signature Σm_ρ² = Σm_ρ (NUMERICAL/HEURISTIC).

## 1. The torus and the Kronecker-line Fejér measure

- **The torus.** T = Hom(ℚ_{>0}^×, S¹) = ∏_p S¹, with characters χ_q(z) = ∏ z_p^{v_p(q)} for q ∈ ℚ_{>0}, and Haar measure dz.
- **(C1) Haar orthogonality (PROVED):** ∫ χ_q d(Haar) = [q = 1]. Distinct n ∈ ℕ have distinct exponent vectors, so
  ⟨χ_m, χ_n⟩_Haar = [m = n]. This is **unique factorization**.
- **The Fejér measure.** The Kronecker line k(t) = (p^{−it})_p is dense and Haar-null. Define
  ν_L := k_*( (1/2π)·4 sin²(Lt/2)/t² dt ), a positive measure. Its Fourier coefficients are
  **ν̂_L(q) = (L − |log q|)_+**, and ν_L(T) = L.
- **(PROVED, Fejér approximate identity):** ν_L/L → δ_1 weakly, where δ_1 is the point mass at the identity of T (t = 0, all z_p = 1).
- **(PROVED, Kronecker–Weyl / Bohr):** time averages of |A(½+it)|² converge to the Haar norm.
- **The torus contains both extremes.** Haar is the counted-once average; δ_1 is the identity point, where the floor lives.
  ν_L/L interpolates between them.

## 2. The SOS identity (Thm 33.1, PROVED; CHECKED)

**Ingredients.** For a finite f = Σ_{n≤X} a_n χ_n, write F(s) = Σ a_n n^{−s}. Write W_f(u) = Σ_{u−L < log n ≤ u} a_n
(a window of length L).

**The identities:**
```
‖f‖²_{ν_L} = Σ_{m,n} a_m ā_n (L − |log(m/n)|)_+
           = (1/2π) ∫ |F(it)|² · 4 sin²(Lt/2)/t² dt
           = ∫_ℝ |W_f(u)|² du ,
⟨f, 1⟩_{ν_L} = Σ_n a_n (L − log n)_+ = ∫_0^L W_f(u) du .
```

**The model vector.** Take f_c = Σ_{n≤X} μ_M(n) n^{−1/2} χ_n − c χ_1, with c = 1/b_M(½) (Arb constants from R32).
- ⟨f_c, 1⟩ = m(L) − cL = m̃(L) + β_M, where β_M is R32's intercept.
- For every real λ: 0 ≤ ‖1 − λ f_c‖²_{ν_L} = L − 2λ(m − cL) + λ² G_c.
- Hence, also using the support of the window:
```
(m(L) − cL)²  ≤  L · Q_c(L)  ≤  L · G_c(L),     Q_c(L) = ∫_0^L |M_{½}(e^u) − c|² du,
G_c = 2 Q_c + L·M_c(X)² − 2 M_c(X)(m − cL),     M_c = M_{½} − c,  M_{½}(x) = Σ_{n≤x} μ_M(n) n^{−1/2}.
```

**Diagonal and remainder.** Split G_c = H_c + R_c.
- The Haar diagonal is H_c = L·‖f_c‖²_Haar = L·((1−c)² + Σ_{2≤n≤X} |μ_M(n)|²/n).
- The remainder is R_c = Σ_{m≠n} a_m ā_n (L − |log m/n|)_+ = 2 Σ_{1≤h<X} Σ_n a_n a_{n+h} (L − log(1+h/n)), with a_1 = 1 − c.

**Exact step integrals.** These come from one segmented pass to 10⁹ (`r33_pass.py`), with no grid quadrature:
```
∫_0^{log X} M(e^u)^2 du = Σ_{n<X} M(n)^2 log(1+1/n) + M(⌊X⌋)^2 (L − log⌊X⌋)
```

**CHECKED** (`r33_identities.json`):
- the quadratic form equals the step integral to 4·10⁻¹⁵ (L = 3, 5, 7);
- the Kronecker-line t-integral equals the quadratic form to 4·10⁻⁷ (L = 3, T = 4000, with a diagonal tail estimate);
- the time average at T = 4000 gives the Haar diagonal: 4.6981 vs 4.6987;
- the SOS inequality holds on every grid point for all 7 models;
- m̃ reproduces R32 to ≤ 4·10⁻¹⁰.

## 3. What the square proves (Thm 33.2, PROVED via Landau + Plancherel)

**Laplace transform.** Let H(s) = (1/b(s+½) − c)/s be the Laplace transform of M_c(e^u).

**(a)** If Q_c(L) ≤ B e^{εL} for all L, then b ≠ 0 on Re s > ½ + ε/2.
- *Why:* M_c(e^u) ∈ L²(e^{−2δu} du) for δ > ε/2, so its Laplace transform converges absolutely on Re s > δ (Cauchy–Schwarz), and Landau gives the zero-free region.

**(b)** If Q_c(L) ≤ B·L for L ≥ L₀, then RH holds for b and every zero on the line is simple.
- Plancherel gives 2ε∫Q_c e^{−2εu}du = (1/2π)∫|H(ε+it)|²dt ≤ B/(2ε).
- A zero of order k at ½+iγ contributes ≍ ε^{1−2k}, which is too large for k ≥ 2.
- Sketch: the same comparison bounds Σ_ρ 1/(|ρ−½|²|b′(ρ)|²) ≲ B.

**(c) CLASSICAL, conditional.**
- RH ⇒ M(x) ≪ x^{½+ε} (Littlewood) ⇒ Q_c ≪ e^{εL}.
- RH + J₋₁(T) ≪ T (Gonek–Hejhal type) ⇒ Q_c ~ κ_μ L, with κ_μ = Σ_ρ 1/(|ρ−½|²|ζ′(ρ)|²) (the Ng 2004 mechanism).

**Controls unbounded (PROVED).**
- A zero at ½+δ+iγ with δ > 0 forces Q_c ≠ O(e^{(2δ−ε)L}) by (a). This covers:
  - F_θ (δ = θ);
  - EPS and DH (zeros in σ > 1: Davenport–Heilbronn; Titchmarsh §10.25);
  - the idealized DMV Beurling system (δ = 0.1).
- ζ² has double zeros, so it fails (b). It is consistent with (a)'s form: under RH Q_c ~ L³.

## 4. Numbers (NUMERICAL; each model with its own μ; L ≥ 3)

| model | X | c = 1/b(½) | Q_c(L_max) | Q_c/L at L_max | growth of Q_c (last third of range) | H_c (Haar) at L_max | R_c/H_c at L_max | R_c sign | verdict |
|---|---|---|---|---|---|---|---|---|---|
| ζ | 1e9 | −0.6848 | 2.994 | 0.1445 | **linear**, slope 0.0291 (4-window 0.02907 ± 0.0009); κ_μ from 300 zeros (test data) 0.02872 + tail | 320.7 | −0.982 | < 0 everywhere (max −8.15) | bounded slope: passes |
| ζ² | 1e9 | 0.4689 | 112.3 | 5.42 | polynomial (d log Q/dL = 0.155 ≈ 3/L) | 2.86e4 | −0.992 | < 0 | passes ε-form, fails linear form (double zeros) |
| F_0.1 | 1e9 | 0.4513 | 271.1 | 13.1 | exponential, rate 0.227 (2θ = 0.2) | 6.87e4 | −0.992 | < 0 | fails |
| F_0.25 | 1e9 | 0.3573 | 1.24e4 | 599 | exponential, rate 0.4996 (2θ = 0.5) | 3.55e6 | −0.993 | < 0 | fails |
| EPS | 1e7 | −0.8703 | 4367 | 271 | exponential, rate 0.81 | 1.29e6 | −0.982 | **> 0 at L ≈ 3–3.6** (R/H up to +0.72) | fails |
| DH | 1e7 | 1.2116 | 11.49 | 0.713 | exponential, rate 0.45 (visible from L ≈ 10) | 2.33e4 | −0.994 | < 0 | **fails (new: R32's pointwise floor did not separate DH)** |
| Beurling | 1e8 | −25.95 | 4831 | 262 | saturating transient (ζ_B(½) = −0.0385 makes c huge) | 1.34e4 | −0.290 | ≤ 0 (max 5·10⁻¹³) | not separated in range; PROVED to fail eventually |
| ζ, ψ side | 1e9 | c_Λ = −1.22574 | 0.899 | 0.0434 | **linear**, slope 0.04623 ± 0.0014 vs **R₀ = 0.0462100** | 4104 | −0.9996 | < 0 | two-sided control of D only |

- **SOS efficiency for ζ.** η = (m−cL)²/(L·Q_c) falls from 0.22 (L = 6) to 0.056 (L = 20.7). Cauchy–Schwarz loses a factor
  ≈ L against the truth, because m − cL ≈ β = 1.84 stays bounded.
- **Bound comparison** (fig. 1d): √(L·Q_c) ≈ 7.9 at X = 10⁹, against the termwise transfer cost 2√X ≈ 6.3·10⁴.

## 5. Where counted-once enters, and where it stops

- **(C1) Unique factorization ⇔ Haar orthogonality of the χ_n.** This is what makes the ν_L-norm's q = 1 Fourier
  mode equal to the diagonal H_c = L·Σ|a_n|².
- **(C2) |μ| ≤ 1** ⇒ H_c ≤ L(L + 1 + (1−c)²) for ζ: polylog (PROVED).
  - For ζ², |μ∗μ| ≤ d(n) ⇒ H ≪ L⁵: polylog, and correctly so, since RH holds for ζ² whenever it holds for ζ.
  - For F_θ, EPS and DH the diagonal is already exponential. **These controls die at the counted-once level.**
- **Beurling satisfies C1 and C2** on its own torus (a free abelian group on generalized primes), so H_B is polylog. Yet
  Q_B must grow exponentially (DMV zeros + Thm 33.2a). Hence **R_B > 0 eventually (PROVED, idealized DMV system)**.
  - So for ζ, GAP-33 is *not* a consequence of C1 + C2. It lives entirely in R_c: the off-diagonal Fourier modes
    q = m/n ≠ 1 of ν_L.
  - Those modes are the Kronecker-line (archimedean) coupling: how the integers sit on the line log n, the pairs (n, n+h).
  - This is R32's (P2) floor-function structure, seen from the torus.
- **CONJECTURE-33⁻ (for ζ): R_c(L) ≤ 0 for all L**, i.e. ‖f_c‖²_{ν_L} ≤ ν_L(T)·‖f_c‖²_Haar ("near the identity, the
  Kronecker line is no hotter than the Haar average").
  - It is sufficient: C2 + Q_c ≤ G_c ⇒ Q_c polylog ⇒ RH by 33.2a.
  - CHECKED for ζ on L ∈ [3, 20.72], with R_c/H_c → −0.98.
  - It fails for EPS at small L, and is PROVED to fail eventually for Beurling. So it is an ℕ-specific statement, not a torus generality.

## 6. GAP-33 (exact) and the named reductions

**The GAP in three forms:**
```
GAP-33 (ε-form):  ∀ε>0 ∃B:  Q_c(L) = ∫_0^L | Σ_{n≤e^u} μ(n) n^{-1/2} − 1/ζ(½) |² du  ≤  B e^{εL}  for all L ≥ 0.
GAP-33₁ (linear): ∃B, L₀:   Q_c(L) ≤ B·L  for L ≥ L₀.
Remainder form (ζ, given C1–C2): GAP-33 ⇔ R_c(L) ≤ B e^{εL};  sufficient: CONJECTURE-33⁻ R_c ≤ 0.
```

**How the forms relate:**
- GAP-33 ⇔ RH, by Thm 33.2(a) and the classical converse 33.2(c).
- GAP-33₁ ⇒ RH + simple zeros.
- GAP-33₁ ⇒ |m̃| = O(L) only, so it does not give GAP-32 (m̃ = O(1)). GAP-32 does not bound the L² mass either. They are independent strengths, and both ⇒ RH + simple zeros.
- ψ side: GAP-33ψ₁ (Q_Λ = O(L)) is a Cramér-type statement (Cramér 1922: under RH ∫(ψ(x)−x)²x^{−2}dx ≍ log X). It gives |D − κ₀| = O(L) only.

**Named reductions:**
- **Littlewood's criterion** (M(x) ≪ x^{½+ε} ⇔ RH) in Fejér/L² form: GAP-33.
- **The σ = ½ weak Mertens hypothesis** (Titchmarsh §14.29; Ng 2004 for the constant): GAP-33₁.
- **Plancherel form:** ∫ |1/ζ(½+ε+it) − 1/ζ(½)|² dt/(ε²+t²) ≪ ε^{−1}, a negative second moment of ζ near the line.
- **Nyman–Beurling / Báez-Duarte (sister, not equivalent).**
  - BBLS: d_N² = inf_A ∫ |1 − ζA|² dt/(2π(¼+t²)).
  - Its Cauchy measure pushes forward to T with ν̂(q) = min(q, 1/q)^{½}: the same Kronecker-line construction, but with the
    Cauchy weight, *fixed in width*, applied to 1 − ζA.
  - Ours uses a Fejér weight that *shrinks to the identity* (ν_L/L → δ_1), applied to M_X with no ζ factor.
  - Báez-Duarte needs a global distance → 0; ours needs polylog growth at the identity.
- **Li coefficients:** no direct link found.

**What a direct attack needs.** This is the exact failing step, not just "equivalent".

1. **The general off-diagonal bound gives back √X.** Montgomery–Vaughan / Hilbert's inequality bounds
   |R_c| ≲ Σ_n n|a_n|² = Σ_{n≤X} μ² ≍ X, so |m − cL| ≲ √(LX). This is the √X again, and it is sharp for
   general coefficients (take a_n ≥ 0, which is exactly the ψ-side / Prop 32.3 situation).
2. **So μ itself must be used.** One needs Σ_{1≤h<X} |Σ_n μ(n)μ(n+h) w_{L,h}(n)/n| ≪ X^{ε}: two-point logarithmic Chowla
   with square-root-type cancellation **uniformly in every shift h < X**.
   - Known: Tao (2016), Σ_{n≤x} μ(n)μ(n+h)/n = o(log x) for each fixed h.
   - Known: Matomäki–Radziwiłł–Tao, averaged over short ranges of h, with o(1) savings.
   - Even full Chowla (o(1) per shift) does not reach X^{ε}.
3. **The mollifier route stops in the same place.** Levinson/Conrey mean values make ∫|ζM − 1|² small on average
   (M of length up to T^{4/7}), but turning that into |1/ζ| bounds fails on the exceptional set, where
   (size × sup|1/ζ|²) is uncontrolled.
   - A fixed-level mollifier M_y (y < X) makes ⟨·, 1⟩ an explicit, zero-free truncation (R32 lesson), so y = X is forced.

## 7. The ψ side and the √X, on the torus (Prop 33.3, PROVED; CHECKED)

**The ψ-side vector.** f_Λ = Σ_{n≤X} (Λ(n) − 1) n^{−1/2} χ_n − c_Λ χ_1, with c_Λ = −ζ′/ζ(½) − ζ(½) = −1.2257372008.

**Identity:** D(L) = κ₀ − ⟨f_Λ, 1⟩_{ν_L} + O(e^{−L/2}), with κ₀ = −8 + ψ′(¼)/4 − ζ′(½) = 0.2219784278359.
- *Derivation:* the R31 closed form for P, plus A(Y) = 4√Y + ζ(½) log Y + ζ′(½) + …, plus c1 = −ζ′/ζ(½) (FE).
- CHECKED: D + ⟨f_Λ, 1⟩ has mean 0.2219784283 and spread 7·10⁻⁷ over L ∈ [8, 20.7].

**Consequences:**
- **The SOS bound is two-sided:** |D − κ₀| ≤ √(L·Q_Λ). D ≥ 0 is the one-sided *linear* statement ⟨f_Λ, 1⟩_{ν_L} ≤ κ₀ + O(e^{−L/2}). That is the Weil functional at t = 0; it is not targeted.
- **Where the √X sits.** On T, Dirichlet convolution is pointwise multiplication, and f_Λ ≈ −Z′_X·M_X − Z_X. The
  multiplier Z_X = Σ n^{−1/2}χ_n has positive coefficients, so its sup is at the identity: Z_X(1) = 2√X + ….
  - Since ν_L/L → δ_1, transferring Möbius ↔ ψ through the product costs ‖Z_X‖_{L^∞(ν_L/L)} → 2√X.
  - **Prop 32.3's √X is the value of the integers' own generating function at the identity of T.**
  - The SOS avoids it by never multiplying: each side carries its own square and its own second moment.
- **ψ slope (NUMERICAL):** 0.04623 ± 0.0014 (4-windows; long baseline 0.04595–0.04607), against R₀ = Σ_ρ m_ρ/|ρ−½|² = 0.04620999 (integer + Γ side, Arb/mpmath).
  - Under RH the mean square of the ψ window is Σ_distinct m_ρ²/|ρ−½|², so equality ⇔ all zeros simple.
  - This is HEURISTIC evidence, the same signature as Burnol's m_ρ². It aligns with R32's R₀ and Round 13g's 0.046.

## 8. Lean candidates (proposals; claim on the bridge first; NOT built)

- **(L8) Finite Fejér-torus SOS.** For finitely supported a : ℕ → ℝ and L ≥ 0:
  - Σ_{m,n} a_m a_n (L − |log m − log n|)_+ = ∫ |Σ_n a_n 1_{[log n, log n + L)}(u)|² du ≥ 0;
  - corollary: (Σ_n a_n (L − log n)_+)² ≤ L·∫_0^L |Σ_{log n ≤ u} a_n|² du.
  - Finite sums plus interval integrals; Mathlib-friendly.
- **(L9) Landau step.** If ∫_0^∞ (Σ_{n≤e^u} μ(n)n^{−1/2} − 1/ζ(½))² e^{−2δu} du < ∞, then riemannZeta s ≠ 0 for re s > ½ + δ.
  It needs the Laplace ↔ Dirichlet-series identity and a Landau/analytic-continuation argument (heavier).
- Spare **(L10)**: χ_m = χ_n on ∏_p S¹ ⇒ m = n (Nat.factorization injective): Haar orthogonality as unique factorization.

## 9. Files

- `r33_pass.py` → `r33_pass_1e09.json`: ζ, ζ² (θ = 0), F_0.1, F_0.25 and the ψ-side vector; exact step integrals; 344 s.
- `r33_controls.py` → `r33_controls_1e07.json`: EPS and DH, each with its own Dirichlet-inverse μ.
- `r33_beur_dump.c` + `r33_beurling.py` → `r33_beurling.json`: 18,645,213 squarefree Beurling integers ≤ 10⁸, primes from R32.
  The 298 MB dump `beur_sqf.bin` and the binary are not committed.
- `r33_identities.py` → `r33_identities.json`: T1–T4.
- `r33_validate_zeros.py` → `r33_validate_zeros.json`: test data only.
- `r33_analysis.py` → `r33_results.json`, `r33_series.npz`.
- `r33_figs.py` → **`fig1_torus_sos.png`** (main figure):
  - (a) Q_c/L for all models;
  - (b) R_c/H_c;
  - (c) windowed slopes against R₀ and κ_μ;
  - (d) the SOS bound against 2√X.

## 10. Labels and credit

- **PROVED:** Thm 33.1 (identities, the SOS inequality, the H + R split); Thm 33.2(a, b); controls unbounded; Beurling R_B > 0 eventually (idealized DMV); Prop 33.3; (C1), (C2).
- **CHECKED:** T1–T4; κ₀ to 5·10⁻¹⁰; m̃ vs R32.
- **NUMERICAL:** the table and the slopes.
- **CLASSICAL-conditional:** 33.2(c), Cramér.
- **CONJECTURE:** GAP-33₁ and CONJECTURE-33⁻ for ζ.
- **HEURISTIC:** the R₀ = slope simplicity signature.
- **GAP:** GAP-33.
- **OPEN:** RH, simplicity, uniform Chowla.

Credit: Bohr, Kronecker–Weyl, Fejér, Landau, Littlewood, Cramér, Titchmarsh, Montgomery–Vaughan, Ng, Gonek–Hejhal,
Nyman, Beurling, Báez-Duarte, Balazard–Landreau–Saias, Burnol, Tao, Matomäki–Radziwiłł, Levinson, Conrey,
Davenport–Heilbronn, Diamond–Montgomery–Vorhauer. Project: Alpöge (window architecture), Codex, Turbo (R31/R32 audits;
R32 wording corrections carried: the Selberg agreement is per-y 2–21%, not a uniform "~10%"; "ζ² passes GAP-31" holds only under RH),
Claude Code.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).
