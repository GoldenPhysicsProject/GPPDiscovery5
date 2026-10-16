# S47 — H47a: ACTUAL Davenport–Heilbronn Weil form on the S14 exact Arb assembly; H47b: carrier-modulated planted-quartet onsets

RH OPEN. GRH OPEN. φ only in DH κ = √(φ+2) − φ, which comes from the INPUT CHOICE q = 5 (genuine for conductor-5 DH; not forced).

## H47a — DH on the unmodulated S14 basis (Legendre·(1−u²)⁴, g ⊂ [−L,L], support 2L)

**Assembly (new, exact):** W_DH = (Pole+Arch)_ζ − Pole_ζ + log5·F(0) + ∫₀^{2L} F(x) sech(x/2) dx − Σ_{2≤n<e^{2L}} 2b(n)/√n F(log n),
b = (a·log)∗a⁻¹, a(n) = Re χ(n) + κ Im χ(n), χ mod 5 with χ(2)=i.
Key identity (PROVED, elementary): ψ(3/4+ir/2) − ψ(1/4+ir/2) = ∫₀^∞ e^{−irx}/cosh(x/2) dx, so Re = π sech(πr) = FT of ½sech(x/2).
The arch difference is therefore a smooth sech kernel. All moments use acb.integral (rigorous); entry radii are ≤ 4e−175.
- **Validation CONFIRMED:** entrywise agreement with the independent S45 float digamma-quadrature form. The difference is 5e−14 at L = 0.7 and 1.0, and 6e−11 at L = 1.3 (float r-cutoff). b(2) = 0.1969085883 = κ log 2, which matches Claude's dh_b_coefficients.

**λ_min (Ritz, plain coefficient norm, mp dps 100–120) — NUMERICAL** (entries are Arb-rigorous, the eigen solve is not):

| M | L | DH even | DH odd | ζ even | ζ odd |
|---|---|---|---|---|---|
| 24 | 1.0 / 1.3 | 6.7e−11 / 1.6e−11 | 5.4e−11 / 1.7e−11 | 3.5e−20 / 2.1e−24 | 2.1e−18 / 2.6e−22 |
| 56 | 1.0, 1.1, 1.2, 1.3 | 8.0, 7.8, 6.5, 5.1 e−16 | 7.2, 6.9, 6.3, 5.3 e−16 | 2.3e−28 … 5.6e−38 | 1.2e−25 … 1.1e−35 |
| 72 | 1.0, 1.1, 1.2, 1.3 | 2.07, 2.04, 2.05, 1.68 e−17 | 1.86, 1.89, 1.74, 1.60 e−17 | 5.5e−30 … 5.1e−44 | 4.1e−27 … 5.9e−41 |
| 88 | 1.0 / 1.3 | 1.05e−18 / 9.5e−19 | 8.9e−19 / 8.6e−19 | 3.9e−30 / 2.4e−49 | 3.2e−27 / 9.1e−47 |

**Answer:** generic support ≤ 2.6 does **NOT** detect actual DH at M ≤ 88 (NUMERICAL). There is no negative eigenvalue in either sector, so there is no Arb witness to certify.
- DH λ_min is about 13–30 orders of magnitude **larger** than λ_min(W_ζ) on the same space.
- DH λ_min is nearly **L-independent**, but it falls with M (×~40 and ×~20 per +16), so it is **not converged**. Ritz λ_min only decreases in M (nested basis), so the true infimum could still be ≤ 0. This is OPEN, not a gap.
- **Contrast with the S46 planted DH point** (CERTIFIED detected at L=1.0, M=72). The planted model sits on the ζ background (λ_min ~1e−30), so a ~1e−30 off-line negative wins. The real DH background is the conductor-5 zero set plus its own Euler-less prime side, and its positive floor ~1e−17 swamps it. **The planted quartet on ζ is NOT a faithful proxy for DH** (NUMERICAL). The S46 reading "S45 DH-PD was an M=24 artifact" is therefore only half right: M matters (λ_min drops 1e−11 → 1e−18), but at M ≤ 88 DH is still PD.

## H47b — carrier-modulated planted-quartet onsets (NUMERICAL, upper bounds on true L*)

Q_fake via zero-sum (5000 ζ zeros, removed pair a,b, plus quartet ½±δ±iγ). The convention is CONFIRMED vs S14 Arb W_ζ, with relative error 4.7e−15 at M=6, L=0.6.
Basis: 8 unmodulated p_j plus 12 p_j·{cos,sin}(γx), L²-orthonormalized on the numerically independent span. Detection threshold: λ < −1e−9.

| γ | δ=0.1 | δ=0.308 | ½log(γ/2π) | best CERTIFIED (S46, unmodulated) |
|---|---|---|---|---|
| 17.58 | (0.4,0.5] | (0.3,0.4] | 0.51 | ≤0.40 / ≤0.40 |
| 51.37 | (1.0,1.2] | (1.0,1.2] | 1.05 | ≤0.8 |
| 85.699 (DH pt) | (1.2,1.5] | (1.2,1.5] | 1.31 | ≤1.0 (δ=0.308); none δ=0.1 |
| 144.56 | **(1.5,2.0]** | **(1.5,2.0]** | 1.57 | none (M≤72, L≤1.3) |

- **New:** γ = 144.6 is detected for the first time, at L = 2.0 (NUMERICAL).
- With this carrier basis the onset is essentially δ-independent and tracks L ≈ ½ log(γ/2π) = π/(mean zero spacing): the window must resolve one local spacing.
- **Grok conjecture "L* depends on γ only via local density": SUPPORTED for the carrier-basis threshold (NUMERICAL), but NOT as the true L*.** At γ = 85.7 the CERTIFIED unmodulated bound ≤ 1.0 beats the carrier 1.2–1.5 and the ½log(γ/2π) = 1.31 law. So the law is a basis property, not intrinsic. The true L*(γ,δ) is OPEN. Both bases are resolution-limited, so neither gives the true onset.

## Scripts / data
- `s47_dh.py M Ls [dps] [prec]` — Arb W_DH + W_ζ, sector λ_min, Arb witness if negative → `dh_M{24,56,72,88}.jsonl`
- `s47_validate.py` — entrywise check vs S45 float DH form
- `s47_carrier.py K K0 Ls` (env NZ) — carrier-modulated planted scan → `carrier_K12_K08.jsonl`
- `s14_exact.py`, `s44_scan.py`, `s46_planted.py` (copied unchanged); `s47_results.json`

## Novelty check
- GPPVerify: DH appears only in NoGo/Ledger and RH criteria/no-go files (DH as a counterexample for Fisher/shifted kernels). There is no DH Weil-form Ritz on polynomial windows.
- Claude `2026-10-13_dh_modulated_window`: carrier-adapted DH detection at h = 2.5 (float) and certificates at h = 3 and 3.5. That is a different basis.
- The sech-kernel difference form W_DH − W_ζ is not in the ledgers (grep "sech"/"DH" in CLAUDE_*/CODEX_*).

Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
