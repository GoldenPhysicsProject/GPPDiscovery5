# S19b (χ₇): Arb audit of Astra's cubic mod-7 signed ground-state block certificate

Daniel Toupin / Golden Physics Project, with AI assistance from Grok. **RH OPEN. GRH OPEN.** Fixed window h ≤ H = ½log3, **real g only**.
This is separate from the χ₄ S19 files (`README.md`, `s19_chi4_*`), which belong to the sibling audit.

Audits `GPP-bridge/work/2026-10-15_astra_chi7_signed_groundstate_block.md` (bridge f4e383f) and
`GPPDiscovery3/experiments/haar_parity_shift/verify_chi7_signed_block.py` (Discovery3 5788008, Fraction intervals).

| file | what | label |
|---|---|---|
| s19_chi7_arb.py → out_chi7_arb.txt | [1] characters through flint `dirichlet_char` plus the explicit generator-3 table; [2] T(H), α, c0, d, b, r, β, v; [3] 2×2 certificate at q=1/200 plus the optimal λ_min(M); reproduces Astra's table; [4] T by quadrature, c0 by digamma, E0(e) Codex form = Astra split; [5] S17-basis test that Q(g) ≥ LB(g) on 403 real vectors, Ritz minimisers included; [6] complex-g caveat | CERTIFIED (Arb, prec 256), except that [5] is a finite-vector sanity test |

Reproduce: `/workspace/venv/bin/python s19_chi7_arb.py` (< 1 s). Uses `../s17` (`Tarb`, `consts`, `build`), `../s14` (basis) and `../s18` (Ritz helper).
Write-up: `GPP-bridge/work/grok/2026-10-15_s19_chi7_signed_groundstate_arb_audit.md`.
