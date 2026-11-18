# Round 13g — Nyman–Beurling / Báez-Duarte distance d_N

Code: `nyman_beurling.py`. Data: `nyman_beurling.json`.

## 1. ζ: an integer-only finite criterion, computed
- d_N² = inf ‖χ_(0,1] − Σ_{k≤N} c_k {1/(kx)}‖²_{L²(0,∞)}. **RH ⇔ d_N → 0** (Báez-Duarte 2003, *Rend. Mat. Acc. Lincei* 14; NB with natural dilations).
- Gram matrix, integer-only (Vasyunin 1995; Báez-Duarte–Balazard–Landreau–Saias 2005). For coprime h, k:
  - ⟨e_h, e_k⟩ = (log 2π − γ)/2 · (1/h + 1/k) + (k−h)/(2hk) · log(h/k) − π/(2hk) · (V(h,k) + V(k,h));
  - V(h,k) = Σ_{m<k} {mh/k} cot(πm/k);
  - ⟨e_{dh}, e_{dk}⟩ = ⟨e_h, e_k⟩/d;
  - ⟨χ, e_k⟩ = (log k + 1 − γ)/k.
  - The formula was checked against direct quadrature at (1,1), (1,2), (2,3), (3,5) to 1e-6, the quadrature truncation level. CHECKED.
- One Cholesky factorization at 60 digits gives every d_n², n ≤ 300 (prefix property). Results (NUMERICAL):

| N | 1 | 2 | 5 | 10 | 20 | 50 | 100 | 200 | 300 |
|---|---|---|---|---|---|---|---|---|---|
| d_N² | 0.8582 | 0.2900 | 0.03610 | 0.02377 | 0.01650 | 0.01187 | 0.01019 | 0.00892 | 0.00799 |
| d_N² log N | — | 0.201 | 0.058 | 0.0547 | 0.0494 | 0.0464 | 0.0469 | 0.0473 | 0.0456 |

  - Báez-Duarte / Burnol conjecture: d_N² log N → 2 + γ − log 4π = 0.04619 (Burnol's lower bound gives ≥ this under simple zeros). The data sit at 0.0456–0.047 from N = 50 on. Consistent, not a proof.
- **Falsifier reading.** Each d_N is finite and integer-only. The criterion, though, is the limit N → ∞, so no fixed N decides anything. That is the same structure as the finite-width W_L.

## 2. The counter-models
- **Necessity half is generic** (CHECKED sketch, the NB form on L²(0,1) / H²(Re s > 1/2) via Mellin).
  - Let F be any Dirichlet series whose generators have Mellin transforms vanishing at a zero ρ of F with Re ρ > 1/2.
  - Then point evaluation at ρ gives d² ≥ (2β − 1)/|ρ|². That is the reciprocal of the H² reproducing-kernel norm, |1/ρ|²·(2β−1).
  - Epstein, ρ = 0.933 + 15.668i: **d² ≥ 3.5·10⁻³**.
  - Spliced Beurling, ρ = 0.6 + 28i: **d² ≥ 2.5·10⁻⁴**.
  - With ζ's observed d_N² ≈ 0.046/log N, ζ only drops below 3.5e-3 at **N ≈ e^{13} ≈ 5·10⁵**, and below 2.5e-4 at **N ≈ e^{181}**. So no computable N separates ζ from the Beurling splice. Separating it from Epstein would need N ~ 10⁶, which is out of reach of this dense 60-digit Gram code.
- **Epstein: the analogue exists in form.** The generators are the lattice-count remainders P_Q(1/(kx)), where P_Q(y) = #{v ∈ ℤ²∖0 : Q(v) ≤ y}/2 − πy/(2√5), with Mellin transform −Z_Q(s)k^{−s}/s. They are integer-valued counts.
  - Their Gram matrix would need ℤ²-lattice analogues of Vasyunin's cotangent sums. I know of no closed form, and I did not compute it (OPEN).
  - The sufficiency half (RH ⇒ d_N → 0) fails structurally. Báez-Duarte's proof mollifies with c_k from μ(k), i.e. 1/ζ as a Dirichlet series with bounded multiplicative coefficients, which needs the Euler product. 1/Z_Q has unbounded, non-multiplicative coefficients, and Z_Q has off-line zeros, so d_N^{(Q)} does not tend to 0.
- **Spliced Beurling: no integer-only analogue.** The generators would count Beurling integers N_B(1/(kx)). Above 55 those are not integers, and there is no cotangent-sum (additive-ℤ) structure. Vasyunin's {mh/k} cot(πm/k) is exactly the additive structure of ℤ / kℤ.
- **Verdict.** RIGHT as a criterion, and the numerics are consistent with the conjecture. It is not a finite certificate: like every criterion so far it is an N → ∞ limit, and the separating scale is astronomically far away for the Beurling splice. The integer-only feature (cotangent sums) is precisely the additive ℤ-structure the 13e sentence asks for. Its sufficiency proof also uses μ, i.e. the Euler product. Labels: Gram formula CHECKED; d_N table NUMERICAL; the fake-model bounds CHECKED sketch; Epstein Gram OPEN; RH OPEN.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5). Audit: Turbo.
