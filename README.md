# GPPDiscovery5

Discovery lab of the **guest seat (Grok)** for the Golden Physics Project. Set up by Daniel Toupin on 2026-10-14.

Public credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.

This repo is **public**. It holds runnable numerics and finished, checkable math. RH strategy, unpublished ideas and proof attempts belong in the private coordination repo, `GoldenPhysicsProject/GPP-bridge`, under `work/guest/` and `GUEST_RESEARCH_NOTES.md`.

## Layout

| Path | Contents |
|---|---|
| `research/<YYYY-MM-DD>_<topic>.md` | Finished write-ups. Each carries a status label (THEOREM / PROVED LEMMA / CONDITIONAL / NUMERICAL / OPEN / CONJECTURE) and a "what this does NOT claim" section. |
| `experiments/<topic>/` | Runnable scripts, each with a README (how to run, dependencies) and its saved output. Every output says whether it is floating-point or certified (interval arithmetic). |
| `notes/` | Scratch. Not citable. |
| `scripts/` | Shared helpers. |
| `requirements.txt` | Python dependencies for the experiments. |

## How work moves

1. Explore here, and keep the private strategy in GPP-bridge.
2. Write an exact result up as a hand-off: statement, proof sketch or computation, status label, reproducible script, and what it does NOT claim.
3. Claude Code audits the hand-off and formalizes the exact parts in Lean 4 in `GoldenPhysicsProject/GPPVerify` (no `sorry`, no axioms, CI green).

## Rules

- **RH is OPEN.** Nothing here proves it.
- **A numerical check is evidence, not a proof.**
- **No credentials** anywhere in this repo.
- **Before calling anything a positivity mechanism,** test it against the two fake L-functions (the shifted product F_θ and Davenport–Heilbronn).
