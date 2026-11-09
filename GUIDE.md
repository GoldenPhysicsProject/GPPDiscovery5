# Guide to GPPDiscovery5

GPPDiscovery5 is the lab of Heart of the Cosmos (Grok Bot by xAI): a public mathematical workspace of the Golden Physics Project holding runnable numerics and finished math from that seat, including the sessions filed under the `s*` folders. Its layout mirrors GPPDiscovery2.

This guide is the map: where things are, how results are recorded, and how
the repository relates to the rest of the Golden Physics Project. It
describes the layout as of 2026-11-09. The workbench moves; the map is
revised when the layout changes.

## Start here

- `README.md`: the repository's own introduction.
- `RESULTS.md`: the human-readable results ledger.
- `results.jsonl`: the raw result rows behind the ledger.
- `GROK_RESEARCH_NOTES.md`: the seat's research notes.
## Layout

- `discovery/<thread>/`: longer threads: write-ups, maps, and syntheses, each labelled THEOREM, PROVED LEMMA, CONDITIONAL, NUMERICAL, OPEN, or CONJECTURE.
- `research/`: dated research notes and checkpoints.
- `experiments/`: self-contained runnable scripts, each with saved output.
- `derivation/`, `nogo_attack/`, `archive/`: derivation threads, no-go attacks, and retired material.
- `s11`, `s12`, …: per-session workspaces with their own result files.
- Results are produced by the GitHub Actions workflow in `.github/workflows/`, which runs scripts and commits their JSON output to `results.jsonl`.

## Conventions

Results carry epistemic labels: THEOREM, PROVED LEMMA, CONDITIONAL,
NUMERICAL, OPEN, KILLED. Numerical work is evidence, not proof; what
hardens is promoted to
[GPPVerify](https://github.com/GoldenPhysicsProject/GPPVerify), the
kernel-checked Lean 4 library, through the build seat's audit. Candidate
positivity mechanisms are tested against the controls F_θ = ξ(s−θ)ξ(s+θ)
and the Davenport–Heilbronn function before they are taken seriously.
Working notes, heuristics, and internal strategy live in the private
GPP-bridge repository; this repository is public by design and is the
publication of the workbench's math, code, and results.

## The family

- [GPPVerify](https://github.com/GoldenPhysicsProject/GPPVerify): the formalization library.
- [GPPDiscovery](https://github.com/GoldenPhysicsProject/GPPDiscovery): Claude Code's workbench (Anthropic).
- [GPPDiscovery2](https://github.com/GoldenPhysicsProject/GPPDiscovery2): ChatGPT Sol's workbench (OpenAI, chat mode).
- [GPPDiscovery3](https://github.com/GoldenPhysicsProject/GPPDiscovery3): ChatGPT Astra's workbench (OpenAI, work mode).
- [GPPDiscovery4](https://github.com/GoldenPhysicsProject/GPPDiscovery4): Turbo's workbench (Muse by Meta).
- [GPPDiscovery5](https://github.com/GoldenPhysicsProject/GPPDiscovery5): Heart of the Cosmos' lab (Grok Bot by xAI).

License: Apache-2.0. Citation: see `CITATION.cff`. Author: Daniel Toupin,
Golden Physics Project.
