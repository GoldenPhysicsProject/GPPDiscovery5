# S45 — H45a/H45b: F_θ & DH falsifiers of support-window positivity (+ p=3 spot)

RH OPEN. GRH OPEN. φ appears only in DH κ (from choice q=5).

## Goal (Claude task 3 remainder)

Test support-≤2.6 ζ positivity against shifted F_θ and Davenport–Heilbronn on the **same** S14 exact Arb Legendre·(1−u²)⁴ assembly used for ζ support-2.6 certs (S36/S38) and S44 p=2 windows.

## H45a — F_θ / DH

### F_θ prime-weight deformation (Arb-CERTIFIED)

```
W_θ = Pole_ζ + Arch_ζ − Σ_n 2Λ(n)/√n · (n^θ + n^{−θ}) · F(log n)
```

**GAP:** arch and pole stay at ζ. Full F_θ needs shifted digamma and poles at ±½±θ. CERTIFIED claims are for this single-Euler-style cosh-weight deformation only (still useful for Claude task 3 beyond p=2).

| θ | L* onset (even) | First CERTIFIED witness | At L=1.3 (M=40) |
|---|---|---|---|
| 1/4 | (0.36, 0.38] | L=0.38: xᵀWx ∈ [−1.70e−9 ± 1e−29] | even ≈ −0.684 (CERTIFIED) |
| 1/20 | ≤0.40 (scan); CERT from 0.53 | L=0.53: even ≈ −0.00213 | even ≈ −0.655 (CERTIFIED) |

ζ stays PD on the same Ritz spaces. Binding sector: **even**.

### Full F_θ (NUMERICAL, float Ritz)

Matches S4: θ=1/4 fails at h*∈(0.44,0.46]; θ=1/20 at (0.50,0.52]. Arch shift delays failure vs the prime-only deformation.

### DH (NUMERICAL)

Still PD at L=1.3 (even λ_min ≈ 3.9e−10, M=24). Tiny negatives at L≥1.6 are float-noise (~1e−14 floor). Claude CERTIFIED failure at h=3. κ=√(φ+2)−φ with **φ from choice q=5** (INPUT-CHOICE; genuine for conductor-5).

### Claude task (3) reading

Support-≤2.6 ζ positivity **excludes F_θ** (prime-def CERTIFIED by L≤0.38; full NUMERICAL by ~0.46) but **does not exclude DH** (PD through L=1.3). Partial RH-strength filter only.

## H45b — p=3 coupling spot (CERTIFIED outer edges, M=40)

c_true = 2 log 3 / √3 ≈ 1.268568201195128. Same gen-eig as S44 with T₃.

| L | even outer [−δ₋,+δ₊] | odd |
|---|---|---|
| 0.55 | near entry; T₃ tiny → gaps huge/fragile | — |
| 0.70 | [−2.35e−5, +5.10e−8] | [−6.03e−7, +3.14e−4] |
| 1.0 | [−2.61e−18, +1.31e−21] | [−1.50e−20, +1.27e−16] |
| 1.3 | [−2.04e−26, +1.27e−28] | [−1.07e−28, +2.00e−24] |

Even = upper cap, odd = lower floor (same as p=2). **Rigidity ≠ RH evidence** (S44 DEAD END).

## Scripts

- `s14_exact.py` — copied from s44 (unchanged)
- `s45_ftheta_scan.py M hs θ [dps] [prec]` — Ritz λ_min for ζ vs F_θ-def
- `s45_ftheta_cert.py M hs θ [dps] [prec]` — Arb rational-vector xᵀWx < 0 witnesses
- `s45_dh_float.py M hs` — NUMERICAL DH + full F_θ float Ritz
- `s45_p3_spot.py M hs [dps] [prec]` — p=3 coupling CERTIFIED outer gaps

## Outputs

`scan_*.jsonl`, `cert_*.jsonl`, `dh_ftheta_float.json`, `p3_cert_M40.jsonl`, `s45_results.json`, `pr11_comment_s45.md`

## Novelty

- GPPVerify main `1e3c8f6`: ShiftedFisher/Hagedorn/ShiftedProductControl exist; **no** Arb S14-assembly F_θ-def / DH support-window witnesses. Do not duplicate CfsClosedChain #372.
- S4: FLOAT F_θ/DH; S43: NUMERICAL Pick F_θ; Claude: CERTIFIED F_θ windows elsewhere (different pipeline).
- **NEW here:** Arb-CERTIFIED negativity on the **same** S14/S44 assembly.

Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
