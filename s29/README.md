# Grok S29 — Codex audits + H29a Krein/eigenfunction lift

RH OPEN. GRH OPEN. Seat: Grok / Loverbot under Daniel Toupin.

## Audits of Codex asks (after S28)
| Script | Target | Verdict |
|---|---|---|
| `audit_codex_rational_tails.py` + `audit_rational_tails_arb_inputs.json` | Exact-rational L=1 tails (6204544406) | CONFIRMED (conditional; Arb single-impl) |
| `audit_one_positive_mode.py` → `.json` | Continuum K_L one-positive-mode (6204606897) | CONFIRMED (minor exposition GAP) |
| `audit_pnt_sobolev.py` → `.json` | PNT Sobolev \|P_L−C_L\| bound (6204562225) | CONFIRMED |

## H29a
- `s29_krein_lift.py` / `h29a_krein_lift.json` — **DEAD END** (finite-t Chebyshev overfits).
- `s29_krein_torus.py` / `h29a_krein_torus.json` — eigen-seeded torus LP; matches S28 A≈3.60 at L=1.19 D=3; does not beat floor gap.

## Run
```bash
source /workspace/venv/bin/activate
python audit_codex_rational_tails.py
python audit_one_positive_mode.py
python audit_pnt_sobolev.py
OMP_NUM_THREADS=1 python s29_krein_torus.py   # slow at L=1.19
```

Bridge note: `GPP-bridge/work/grok/2026-10-15_s29_codex_audits_krein_lift.md`.
