# T3a_plunge_mellin — Queen's (T3a) (low-Mellin mass of the trimmed plunge states): the stated rate fails; the true scale is (log Λ)/Λ²

Spectrum & Measure seat, GPP RH team. Daniel Toupin / Golden Physics Project, with AI assistance from Grok. **RH is OPEN.**
This file is RH-silent. Numerics test one step only and are labelled NUMERICAL.

Source, used verbatim: team/queen 47e63a4, ARP_queen_Hkappa.md §H.9 (table row (T3a)), quoted here:
> (T3a) ω^{arch}_Λ(S) := sup over unit φ ∈ span{ψ_n : n ∈ I_Λ \ I^{(K)}} of ∫_{|s|≤S}|Mφ|² / ∫|Mφ|² ≤ Λ^{−B−2−δ} for S = log Λ
> (any polynomial rate faster than Λ^{−B−2} suffices). HEURISTIC (plunge states live in phase space near |x| ≈ Λ, |ξ| ≈ Λ,
> i.e. at Mellin height s ≍ 2πΛ², H.1b).

Conventions (Queen §H.0/H.9): ψ_n are even prolates on [−Λ, Λ] with band Λ; c = 2πΛ² (Slepian), M = 4Λ²;
Mφ(s) = ∫_0^Λ φ(x)x^{−½+is}dx; ∫|Mφ|² = 2π∫_0^∞|φ|². The modulus |Mφ| is dilation-invariant up to a constant factor:
φ(x) ↦ φ(Λx) multiplies Mφ by Λ^{−½−is}. So ω^{arch} may be computed on [0, 1] with c = 2πΛ².

## 1. The phase-space premise is wrong for x (HEURISTIC, standard WKB) — and why it matters

- The prolate operator W = −∂(1 − x²)∂ + c²x² has classical symbol (1 − x²)p² + c²x².
- Transition-band states have χ_n ≈ c², whose level curve (1 − x²)p² + c²x² = c² is exactly **p = ±c for all |x| < 1**, together with the edges x = ±1.
- So plunge states are spread along the horizontal edges |ξ| ≈ Λ of the phase-space rectangle, **at all x, including x ≈ 0**. They are not concentrated near |x| ≈ Λ.
- The Mellin height is the dilation symbol s = 2πxξ. On the horizontal edge it ranges over all of [0, 2πΛ²]. The portion with |s| ≤ S is |x| ≤ S/(2πΛ), a fraction ≈ S/(2πΛ²) of the orbit.
- Hence the expected size is **ω^{arch}_Λ(S) ≍ S/Λ²**, i.e. (log Λ)/Λ² at S = log Λ. It is polynomially small, but at the rate Λ^{−2}, not Λ^{−B−2−δ}.

**Model computation (PROVED, exact).**
- For φ(x) = cos(cx) on [0, ∞): M(s) = Γ(½ + is)cos(π(½ + is)/2)c^{−½−is}.
- So |M(s)|² = (π/cosh πs)·(cosh πs/2)/c = **π/(2c)**, flat in s.
- A unit-amplitude wave at frequency c therefore has low-Mellin mass ≈ 2S·π/(2c) = πS/c, against total mass π·(length). On [0, 1] this is a fraction ≈ S/c ≍ S/Λ².
- Truncation to [0, 1] adds an endpoint term of relative size O(1/(c|½ + is|)) = O(Λ^{−2}). It does not change the order (HEURISTIC).
- Plunge states have amplitude ψ_n(0)² ≍ ∫_0^1ψ_n² (WKB: amplitude ∝ ((1 − x²)k(x))^{−½} with k ≈ c, roughly flat). They behave like this model near x = 0.

## 2. Numerical test of the step (NUMERICAL)

- **Method.** Even prolates from the Legendre–Galerkin matrix of W (exact tridiagonal-in-steps-of-2 entries). λ_n = (c/2π)|μ_n|² from the integral equation. Low-Mellin fraction = ∫_{|s|≤log Λ}|Mψ_n|²/(2π∫_0^1ψ_n²), via x = y² and 4000-point Gauss–Legendre. Plunge = {Λ^{−K}-type: 10^{−12} < 1 − λ_n, λ_n ≥ ½}.
- **Script.** /workspace/sm/t3a/t3a_test.py on the box. Its logic is fully described above. Commit policy here is one file, so it is not committed. It can be committed on request.

| Λ | c | plunge states tested | low-Mellin fraction ω (range) | ω·Λ²/log Λ | ψ_n(0)²/∫_0^1ψ_n² |
|---|---|---|---|---|---|
| 2 | 25.1 | n = 6…14 | 2.0e−2 … 7.3e−2 | 0.12 … 0.42 | 0.70 … 1.57 |
| 3 | 56.5 | n = 22…34 | 1.1e−2 … 2.9e−2 | 0.09 … 0.24 | 0.60 … 1.14 |
| 4 | 100.5 | n = 48…62 | 5.9e−3 … 1.7e−2 | 0.07 … 0.20 | 0.54 … 0.95 |
| 5 | 157.1 | n = 84…98 | 5.1e−3 … 8.7e−3 | 0.08 … 0.14 | 0.51 … 0.82 |

- Each single plunge state has ω ≈ (0.1–0.2)·(log Λ)/Λ². The sup over their span is at least the largest of these.
- Λ^{−B−2−δ} with B ≥ 0 would require ω·Λ²/log Λ → 0 like Λ^{−B−δ}/log Λ. The data show no such decay over Λ = 2…5.
- Every ψ_n(0)² is O(1): the plunge states are not edge-localized in x.
- Kept states (1 − λ ≤ 10^{−12}) include ψ_0, whose low-Mellin fraction is O(1) (0.76 … 0.96), as expected for a slowly varying profile.

## 3. Conclusions (labels)

| # | Statement | Label |
|---|---|---|
| T1 | (T3a) as stated (ω^{arch}_Λ(log Λ) ≤ Λ^{−B−2−δ}, any B ≥ 0) | **expected FALSE**: CONJECTURE (negation), supported by §1 (HEURISTIC) and §2 (NUMERICAL) |
| T2 | Corrected size: ω^{arch}_Λ(S) ≍ S/Λ² for transition-band states (upper and lower) | CONJECTURE (HEURISTIC + NUMERICAL) |
| T3 | Model |M cos(c·)|² = π/(2c) on [0, ∞) | PROVED |
| T4 | Consequence in Queen's reduction (§H.9 "Reduction of ω_Λ"): ω ≤ CΛ^B log Λ·ω^{arch}(log Λ) ≍ Λ^{B−2}log²Λ, so **(T3) ⇐ (T2) iff B < 2** (with r_K = O(K log²Λ), one needs Λ^{B−2}log²Λ = o(1/log²Λ)) | PROVED given T2 (algebra) |
| T5 | The Λ^B loss is the worst-case m_Λ. Plunge states carry their |Mφ|² at heights up to 2πΛ², where |ζ(½ + is)|² has mean ≍ log Λ. A route weighting |ζ|² against the actual profile, not by m_Λ, would replace Λ^B by (log Λ)^{O(1)} | HEURISTIC (suggestion for (T3b), Queen) |

**Bottom line for the lead and Queen.**
- (T3a) at the requested rate is not provable because it is (by all evidence) false. The plunge states have low-Mellin fraction ≍ (log Λ)/Λ², since they oscillate at frequency ≈ Λ uniformly in x, including near x = 0.
- (T3) still closes along Queen's reduction if κ_Λ ≤ CΛ^B with **B < 2**, or if (T3b) avoids the m_Λ loss (T5).
- First-principles upper bound ω^{arch} ≤ CS/Λ² for the trimmed span: **GAP**. Natural route: a WKB/Plancherel–Rotach bound sup_{|x|≤δ}|ψ_n(x)|² ≤ C∫ψ_n² uniformly on the transition band (KRD/Bonami–Karoui-type, not yet located verbatim), plus the model of §1. Not claimed here.
