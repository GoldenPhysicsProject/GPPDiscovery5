# PARKED numerics (Principal series seat): state and resume notes

Parked 2026-10-22, ~07:55 EDT, per the lead's derivation pivot (relaying Daniel). **RH is OPEN.** All entries are NUMERICAL unless labelled otherwise. L is the half-width.

| item | state | resume |
|---|---|---|
| **dh_signed L = 1.7104, K = 240** (Step 0, still allowed) | running on the box (`/workspace/newbot/dhsigned`, `dhritz_resid.py 1069/625 240 320 60 even`). K = 60 and 120 are pushed with residuals. | When `cosvec_DH_even_L1069_625_K240.json` appears, copy it plus `resid_1069_625_K240.jsonl` into `dh_signed/`. |
| Epstein x²+5y² Galerkin ladder (`epstein/`) | 11 runs done: K = 40/60 at all L, K = 90 at L = 1.5. K = 90 at L = 2, 2.5, 3 was **killed** (not run). Results are in session-2 notes §7. | `cd epstein; for L in 2 2.5 3; do python lkritz.py EPS $L 90; done`, then `eps_overlap.py`. About 30–60 min each under load. |
| DH lowest even vector ladder (`dh_ladder/`) | **done**: L = 1–1.65 at K = 60/120, plus the S75 witnesses. Nodeless, w₀ = 0.994 → 0.799. | Nothing pending. |
| DH beyond onset, L = 2.5, K = 100 (`epstein/lkvec_DH_L2.5_K100.json`) | done, one K only. Lowest vector is Queen's pair vector: overlap 0.937, 134 sign changes. | For convergence, K = 140 at L = 2.5. |
| Onset law, both sectors (`twist/hp_onset48*`) | **done**: N = 48 for L = 0.5–1.0, N = 56 at dps 130 for L = 0.5, 0.8, 1.0. | Optional: N = 64 at dps ≥ 150 for L = 1.0. |
| F_θ node ladder (`twist/hp_nodes*`) | done (θ = 0.1, L = 0.5–0.8, N = 24/32/40). | Nothing pending. |
| Hecke G₅ scattering Weil-form test | **not started**. Needs G₅ scattering-pole / Maass data (Hejhal 1992; Strömberg; Avelin). | Pull the data first. |
| Selberg-zeta η calibration | **dropped** by the lead. | — |
| Turbo S1 audit | blocked: no script in a repo. Verdict in bridge notes session 2 §6. | Needs Turbo's S1 L, θ, K and script. |
