# S18 — Audit of Codex's high-conductor / first-active-prime Dirichlet Weil coercivity

Daniel Toupin / Golden Physics Project, with AI assistance from Grok. **RH OPEN. Dirichlet GRH OPEN.**
Audits `GPP-bridge/work/codex/2026-10-15_high_conductor_uniform_window_GRH_coercivity.md` (PR #11 asks 6187022122, 6187131901).
Builds on S17 (`../s17`: `Tarb`, `consts`, `build`, `forms`, `E_comb`). Python: `/workspace/venv/bin/python` (python-flint 0.9, prec 256).

| file | what | label |
|---|---|---|
| s18_common.py, s18_ritz_common.py | shared helpers (h*, (log2)/√2, primitive characters through flint `dirichlet_char`, Ritz plus Arb Rayleigh quotient) | — |
| s18_gaps_arb.py → out_gaps_arb.txt | B(h) ≤ 8he^h on a grid, plus looseness; T0/T1/c0/c1/(log2)/√2; Codex gaps; E_comb(h*); per-character analytic lower bounds for the residual set; improved universal thresholds | CERTIFIED (Arb) |
| s18_disjoint_ritz.py → out_disjoint_ritz.txt | sup \|F_g(log2)\|/‖g‖² by Ritz (M=16/32/48) with Arb Rayleigh quotients, h ∈ {0.36,…,1.0} | NUMERICAL + CERTIFIED witnesses |
| s18_low_cond.py → out_low_cond.txt, low_cond_32_48.json | Ritz λ_min of the true completed real-g form at h* for every residual primitive class, plus ζ and F_{1/4} controls | NUMERICAL (upper), CERTIFIED upper witnesses |

Reproduce: `python s18_gaps_arb.py; python s18_disjoint_ritz.py 16 32 48; python s18_low_cond.py 32 48`. The full write-up is
`GPP-bridge/work/grok/2026-10-15_s18_first_prime_conductor_audit.md`.
