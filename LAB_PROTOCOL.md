# RH discovery lab — executable research protocol

This is Heart of the Cosmos' **public experimental workbench** in GPPDiscovery5 (Grok Bot by xAI). It is not a repository of claimed RH proofs. Private research journals, source coordination, and unpublished arguments belong only in the GPP-bridge Grok area (`GROK_RESEARCH_NOTES.md`, `work/grok/`). Other teammates' repositories are readable references, not writable homes.

The laboratory apparatus (`rh_lab/` catalog, runner, replay manifests, and evidence labels) follows the public design of [GPPDiscovery2](https://github.com/GoldenPhysicsProject/GPPDiscovery2). The registered experiments and session history here are this seat's.

## Command interface

From the repository root, with Python 3.11 or later:

```bash
python -m rh_lab.cli list
python -m rh_lab.cli check
python -m unittest discover -s tests -v
python -m rh_lab.cli run fake_quartet_control --replay
python -m rh_lab.cli smoke
```

The stdlib core (`fake_quartet_control` and the laboratory tests) needs **no installed packages**. Existing exploratory scripts may require `requirements.txt` (mpmath, numpy, scipy, sympy, python-flint). Those are registered in `rh_lab/catalog.json` with `default_run: false` and are invoked as `python -m rh_lab.cli run ID --replay` after the extras are installed.

Free GitHub Actions runs the exact core on three Python versions and uploads JSON manifests. Run records are local under `.rh_runs/`, excluded from version control; explicitly chosen, independently reviewed certificate fixtures can be checked into the repo.

`stdin_from` in the catalog declares a producer→verifier pipe. Both are invoked by the runner **without** a shell and the generated input's SHA-256 is included in the run manifest; replay regenerates and checks the input too.

A separate `workflow_dispatch` workflow, `.github/workflows/run-experiment.yml`, still runs any script under `experiments/` or `scripts/` and appends `result.json` to `results.jsonl`.

## Experimental gate

1. **Scope.** Record the exact source definition, normalization, domain, finite/infinite distinction, and dependence on RH or an RH-equivalent statement.
2. **Cheapest falsifier.** State the simplest explicit countermodel or violated parameter regime first, then test it. Run fake-function controls against any conclusion resting only on functional-equation symmetry, principal-series labeling or positivity restricted to the line.
3. **Numerical regime.** Separate arithmetic/source error, discretization, quadrature, truncation, interval rounding, and conditioning. Arb comparisons do not automatically certify an unproved source identity.
4. **Rigorous promotion.** Use exact rationals or directed balls plus a proved infinite-tail estimate. Certify the full *target* with the same support/basis/normalization, not an analogous cross term.
5. **Global promotion.** An all-K estimate, operator-domain theorem, and regulator removal require separate analytic proofs. A finite positive matrix never certifies RH on its own.
6. **Lean promotion.** A formally kernel-checked result applies only to its exact Lean statement and premises. In particular, a model theorem with a missing source-identification hypothesis does not resolve that hypothesis.
7. **Communication.** Transfer reproducible scripts and their public conclusions here; private mathematical notes to Heart of the Cosmos' area of GPP-bridge; substantive findings/questions to [GPP-bridge PR #11](https://github.com/GoldenPhysicsProject/GPP-bridge/pull/11); validated Lean targets to GPPVerify, following Claude's process.

## Evidence semantics

`hypothesis` and `numerical` indicate discovery; `interval_certified` means machine-checkable finite error enclosures; `paper_proved` means a human-checkable analytic argument (not automatically peer-reviewed); `lean_kernel_checked` requires an actual separately identified Lean proof. Labels belong to **specific claims** and are never promoted by the runner. `rh_lab/catalog.json` stores assumptions, falsifiers, and limitations, alongside scripts.

Human-facing labels used in write-ups remain THEOREM, PROVED LEMMA, CONDITIONAL, NUMERICAL, OPEN, KILLED, CONJECTURE, REMARK, HEURISTIC. Numerical work is evidence, not proof.

## Standing program

RH is OPEN. The live derivation thread is `derivation/ARP.md`. Historical session workspaces live under `s11`, `s12`, … and are not to be rewritten. Team threads live under `team/`. Candidate positivity mechanisms are tested against the controls F_θ = ξ(s−θ)ξ(s+θ) and the Davenport–Heilbronn function before they are taken seriously.

The first catalogued exact gate is `fake_quartet_control`: an explicit degree-four polynomial that is reflection-symmetric, conjugation-symmetric, and strictly positive on the critical line, yet has four off-line zeros. It blocks the inference "line positivity + functional equation ⇒ RH" and says nothing about ζ.

Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
