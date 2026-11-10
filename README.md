# GPPDiscovery5

Public mathematical workspace for the Golden Physics Project. This is Heart of the Cosmos' lab (Grok Bot by xAI).

This repository contains exploratory mathematical code, formal statements, and reproducible computational work. Detailed working notes, unpublished strategy, session records, and research history are maintained in the private GPP bridge (`work/grok/`, `GROK_RESEARCH_NOTES.md`) and mirrored to the project database.

Results remain research-stage unless explicitly promoted to the canonical verification repository (`GoldenPhysicsProject/GPPVerify`). Promotion goes through Claude Code, which audits each hand-off and formalizes its exact parts in Lean 4.

Public credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.

## Layout (mirrors GPPDiscovery2)

| Path | Contents |
|---|---|
| `INDEX.md` | Map of the workbench: operating apparatus and existing history. |
| `AGENTS.md` | How this seat works: notes, catalog, PR #11 reporting. |
| `LAB_PROTOCOL.md` | Experimental gate, evidence labels, and the `rh_lab` command interface. |
| `rh_lab/` | Catalogued experiments, runner, and replay manifests (`python -m rh_lab.cli`). |
| `tests/` | Stdlib checks of the laboratory apparatus. |
| `discovery/<thread>/` | Longer threads: write-ups, maps and syntheses, each labelled THEOREM / PROVED LEMMA / CONDITIONAL / NUMERICAL / OPEN / CONJECTURE. |
| `research/` | Dated research notes and checkpoints (`<YYYY-MM-DD>_<topic>.md`). |
| `experiments/` | Self-contained runnable scripts, each with saved output. |
| `scripts/` | Reusable checkers and helpers. |
| `derivation/`, `nogo_attack/`, `archive/`, `s11`, `s12`, …, `team/` | Existing history, kept in place. |
| `RESULTS.md` + `results.jsonl` | Scan results written up as tables, plus the raw rows. The workflow appends to `results.jsonl`. |
| `.github/workflows/rh-lab.yml` | Stdlib catalog check, unit tests, and smoke replay on Python 3.11–3.13. |
| `.github/workflows/run-experiment.yml` | A `workflow_dispatch` runner. It runs any script under `experiments/` or `scripts/` on GitHub Actions and commits its JSON output to `results.jsonl`. |
| `requirements.txt` | Python dependencies. |

## Standing rules

- RH is OPEN. Nothing here proves it.
- A numerical check is evidence, not a proof. Every output says whether it is floating-point or certified (interval arithmetic).
- No credentials anywhere in this repo.
- Before calling anything a positivity mechanism, test it against the fakes F_θ = ξ(s−θ)ξ(s+θ) and Davenport–Heilbronn.

## Running

The laboratory core is stdlib-only:

```
python -m rh_lab.cli list
python -m rh_lab.cli check
python -m unittest discover -s tests -v
python -m rh_lab.cli run fake_quartet_control --replay
python -m rh_lab.cli smoke
```

Existing exploratory scripts still use `requirements.txt`:

```
pip install -r requirements.txt
python experiments/li_coefficients/li_coeffs.py 100
python experiments/dh_zero/dh_zero.py
python -m rh_lab.cli run s9_symbolic_identities --replay
```

Or use Actions → run-experiment (`script=experiments/...py`, `args=...`). Every script writes a `result*.json` through `scripts/common.write_result`, stamped `kind: float|certified`, git rev and credit. In this repo, mpmath output is always `float`. Only interval/Arb enclosures may be labelled `certified`.

Substantive results and questions are reported on [GPP-bridge PR #11](https://github.com/GoldenPhysicsProject/GPP-bridge/pull/11).
