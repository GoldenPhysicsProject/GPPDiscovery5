# S19 — Independent Arb audit of Codex’s χ₄ first-prime-window ground-state certificate

Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
**RH OPEN. Dirichlet GRH OPEN.** Finite Mellin window only: `0 < h ≤ h* = ½ log 3`.

Audits Codex’s notes and scripts (PR #11 comments 6187965201, 6188137220):
- `GPP-bridge/work/codex/2026-10-15_chi4_positive_ground_state_first_prime_window.md`
- Discovery2 `chi4_groundstate_interval_certificate.py` (mpmath.iv) and
  `chi4_groundstate_pure_rational_certificate.py` (Fraction) — **reference only; not imported**.

GPPVerify: `#362 ConductorPrimeCoercivity` MERGED. `#363 GroundStatePicone` OPEN
(Blueprint green; finite edge algebra only — does **not** yet certify infinite Picone +
interval bound). Arb audit remains novel.

| file | what | label |
|---|---|---|
| `s19_chi4_arb.py` → `out_chi4_arb.txt` | Independent 15-mode / 100-cell Arb certificate (`ctx.prec=512`): kh*<π/2, s-sign pattern, monotone cell lowers, δ=c₁−log4, gap vs 0.025/0.02 | CERTIFIED (Arb) finite check; Picone + j≥15 remainder ANALYTIC / Lean-open |
| `s19_chi4_ritz_crosscheck.py` → `out_chi4_ritz.txt` | Ritz λ_min of true Q_χ4 at h* (M=32,48), reuse S17/S18 | NUMERICAL (upper on λ_min; not a lower bound) |

Reproduce:
```
/workspace/venv/bin/python s19_chi4_arb.py | tee out_chi4_arb.txt
/workspace/venv/bin/python s19_chi4_ritz_crosscheck.py 32 48 | tee out_chi4_ritz.txt
```

## Claim audited (labels)

For primitive odd χ₄ mod 4: f=4, ε=1, χ₄(2)=0. On `0<h≤h*` only n=2 could enter,
but χ₄(2)=0, so
`Q_χ4(g) = E_odd(g) + (log4 − c₁)||g||² = E_odd(g) − δ||g||²`
with `δ = c₁ − log4 = log(2π)+γ−π/2 ≈ 0.84429640451598…`.

Ground-state comparison with `φ=cos(kx)`, `k=19/10`, 15 modes, 100 cells:
`min C_i ≈ 0.87958346…` (cell 64) `> δ`, hence
`Q_χ4(g) > 0.025 ||g||²` (mpmath) / `> 0.02 ||g||²` (rational).

**Non-claims:** RH, Dirichlet GRH, Lean-certified Picone, Lean-certified infinite-mode
remainder. φ does not appear.

Full write-up: `GPP-bridge/work/grok/2026-10-15_s19_chi4_groundstate_arb_audit.md`.
