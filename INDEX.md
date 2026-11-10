# Index of GPPDiscovery5

A map of this workbench. Nothing listed here is a claimed proof of RH. For conventions, start with `GUIDE.md` and `LAB_PROTOCOL.md`.

## Operating apparatus (parity with GPPDiscovery2)

| Path | Role |
|---|---|
| `README.md` | Public introduction. |
| `GUIDE.md` | Layout and relation to the rest of the project. |
| `AGENTS.md` | How this seat works: read notes, catalog experiments, report to PR #11. |
| `LAB_PROTOCOL.md` | Experimental gate, evidence labels, command interface. |
| `rh_lab/` | Catalog, runner, replay manifests (`python -m rh_lab.cli`). |
| `tests/` | Stdlib tests of the laboratory apparatus. |
| `RESULTS.md` + `results.jsonl` | Human-readable ledger and raw rows. |
| `.github/workflows/rh-lab.yml` | Stdlib catalog check, unit tests, and smoke replay on Python 3.11–3.13. |
| `.github/workflows/run-experiment.yml` | On-demand runner for any `experiments/` or `scripts/` script. |

## Research layout (existing history, kept in place)

| Path | Contents |
|---|---|
| `discovery/` | Longer labelled threads. |
| `experiments/` | Self-contained runnable scripts with saved output. |
| `research/` | Dated research notes. |
| `derivation/` | Derivation threads, including `ARP.md`. |
| `nogo_attack/` | No-go attacks. |
| `archive/` | Retired material. |
| `scripts/` | Shared helpers (`scripts/common.py`). |
| `s11`, `s12`, … | Per-session workspaces and their result files. Do not rewrite. |
| `team/` | Sub-team threads (celestial-cft, prime-number-theory, principal-series, queen, spectrum-measure). |
| `GROK_RESEARCH_NOTES.md` | Public copy of selected research notes; the rolling private record lives in GPP-bridge. |

## Reporting

Substantive results and questions go to [GPP-bridge PR #11](https://github.com/GoldenPhysicsProject/GPP-bridge/pull/11). Private continuity notes go to GPP-bridge `GROK_RESEARCH_NOTES.md` and `work/grok/`. Lean targets go to [GPPVerify](https://github.com/GoldenPhysicsProject/GPPVerify).

Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
