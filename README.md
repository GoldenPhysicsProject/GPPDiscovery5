# GPPDiscovery5

Public mathematical workspace for the Golden Physics Project, operated by the guest seat (Grok).

This repository contains exploratory mathematical code, formal statements, and reproducible computational work. Detailed working notes, unpublished strategy, session records, and research history are maintained in the private GPP bridge (`work/guest/`, `GUEST_RESEARCH_NOTES.md`) and mirrored to the project database.

Results remain research-stage unless explicitly promoted to the canonical verification repository (`GoldenPhysicsProject/GPPVerify`). Promotion goes through Claude Code, which audits each hand-off and formalizes its exact parts in Lean 4.

Public credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.

## Layout (mirrors GPPDiscovery2)

| Path | Contents |
|---|---|
| `discovery/<thread>/` | Longer threads: write-ups, maps and syntheses, each labelled THEOREM / PROVED LEMMA / CONDITIONAL / NUMERICAL / OPEN / CONJECTURE. |
| `research/` | Dated research notes and checkpoints (`<YYYY-MM-DD>_<topic>.md`). |
| `experiments/` | Self-contained runnable scripts, each with saved output. |
| `scripts/` | Reusable checkers and helpers. |
| `RESULTS.md` + `results.jsonl` | Scan results written up as tables, plus the raw rows. The workflow appends to `results.jsonl`. |
| `.github/workflows/run-experiment.yml` | A `workflow_dispatch` runner. It runs any script under `experiments/` or `scripts/` on GitHub Actions and commits its JSON output to `results.jsonl`. |
| `requirements.txt` | Python dependencies. |

## Standing rules

- RH is OPEN. Nothing here proves it.
- A numerical check is evidence, not a proof. Every output says whether it is floating-point or certified (interval arithmetic).
- No credentials anywhere in this repo.
- Before calling anything a positivity mechanism, test it against the fakes F_θ = ξ(s−θ)ξ(s+θ) and Davenport–Heilbronn.
