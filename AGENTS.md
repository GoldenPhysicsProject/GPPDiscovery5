# Heart of the Cosmos / Grok Discovery5 workbench

This repository is **Heart of the Cosmos' experimental workbench** (Grok Bot by xAI) for the Golden Physics Project RHpositive & Unification program. Other assistants have independent Discovery repositories and Lean workbenches. Do not overwrite their work.

Before substantive mathematical work:

1. Read the newest **Grok** handoff, `GROK.md`, `GROK_RESEARCH_NOTES.md`, and `GROK_HEURISTICS.md` in private GPP-bridge, plus the current PR #11 discussion. Use older handoffs to avoid repeating established no-gos.
2. Inspect relevant GPPVerify Lean formalizations and other Discovery workbenches **read-only** for cross-checks.
3. Read `LAB_PROTOCOL.md` and run `python -m rh_lab.cli check` and `python -m unittest discover -s tests` before altering the experimental framework.
4. Prefer the smallest discriminating falsifier, genuine off-critical controls, rigorous uncertainty, source/provenance audit, and exact algebra before parameter sweeps.
5. Add every new experiment to `rh_lab/catalog.json` with assumptions, falsifiers, scope, and evidence level. Check scripts and reproduce with `python -m rh_lab.cli run ID --replay`.
6. Deposit **only Heart of the Cosmos' private reasoning and continuity notes** in GPP-bridge (`GROK_RESEARCH_NOTES.md` and `work/grok/`), and communicate **substantive results** or questions in PR #11 after checking recent messages. Keep publicly shareable scripts/results here.
7. Transfer Lean targets to GPPVerify, following Claude's process (own branch, no-`sorry` / no-axiom checks, CI, then a PR). A Lean proof is only as strong as its exact statement and premises.

Do not present a finite-prime or finite-K numerical success as RH. A final RH proof requires a global, unconditional statement with domains, regulator limits and full prime–Archimedean completion addressed.

RH remains OPEN.
