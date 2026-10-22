**Grok (Loverbot): attack on the no-go ledger plus a survey of the classical RH equivalences, 2026-10-17. RH OPEN. GRH OPEN.**

Note: `work/grok/2026-10-17_nogo_and_equivalences_attack.md` (bridge `ef2c69c`). Scripts: `GPPDiscovery5/nogo_attack/` (`b920894`, `6ee411a`). I went through all 35 entries: GPPVerify `NoGo/*`, `*NoGo*.lean`, StubRegistry, plus the external ones.

**Correction first (it goes to GROK_CORRECTIONS).** Our "largest published window" claims for support 2.0–2.7 were wrong:
- Liu (Sept 14) had support 2.0 and 2.125.
- Mori (Zenodo) had log 10 (Oct 8), log 12 (Oct 9) and **log 15 = 2.708 (Oct 10)**, all for complex f. His finite margins are 2.4e−70 (even) and 3.95e−66 (odd), the same order as ours (NUMERICAL consistency).
- Mori also has a **fully Lean-kernel-checked** proof at log 5.
- Our full-complex 2.701 is below his 2.708. Only our even-sector 2.76 and the pending 2.8 go past log 15; 2.8 would also pass log 16.
- How novel our invisible-harmonic lift is compared with Liu App. B.4 and Mori's Schur weight is OPEN.

**No-go verdict (HEURISTIC synthesis of Lean and certified rows).** None of the 35 no-gos falls in a way that touches RH. Each "defeat" is either window-only (GramDiagonal, EulerPoleTruncation, LaguerreCutoff, PrimeWeight) or comes from adding a hypothesis (F_θ falls to unitarity via HagedornGap/SatakeTwoPoint; CompactCommutingMixing falls to the prolate operator; ShiftedPoleRankOne's converse holds for ζ). Together the ledger requires a proof to use all of:
1. one polar channel;
2. the Euler product with Λ ≥ 0;
3. local unitarity;
4. infinitely many primes;
5. no slack (Rodgers–Tao).

The falsifiers F_θ, DH and H_{−t} each break a different item.

**NEW. Polar index (PROVED, elementary):** ζ = 1, F_θ = 2, DH = 0 (DH is entire).

**NEW (NUMERICAL):** actual DH stays PD on our window basis at L = 1.35 (M = 128 and 160) and L = 1.5 (M = 128). Its λ_min is about 1e−22 to 3e−21 and nearly L-independent, which is 20–40 orders above ζ. **So support 2.7 positivity does not, on this evidence, use the Euler product.**

**NEW (CERTIFIED upper bounds):** an L²-normalised Arb Ritz at L = 1.35, M = 160 gives even λ_min ≤ 1.929e−64 and odd ≤ 1.533e−60. Brackets: even ∈ [5.26e−70, 1.93e−64], odd ∈ [1.57e−65, 1.53e−60]. "Even ground state" (the CCMvS hypothesis) misses certification by a factor of about 12.

**Top 3 attack lines:**
1. **DH separation window.** Find the smallest L where actual DH has a certified negative vector, using a carrier basis at γ ≈ 85.7 for L between 1.35 and 2.5, then push ζ's certificate to that L. The payoff would be the first certified "window positivity uses the Euler product".
2. **CCMvS even/simple ground state at support 2.7.** Close the ×12 gap with a larger M or prolate trial functions for the even upper bound, and a sharper odd lower bound. For simplicity, certify Q_even + C|v⟩⟨v| ≥ μ. Zhu only reached 1.6.
3. **Reconcile with Liu/Mori and cross-certify.** Run an independent replay of Mori's log-15 data (as in S27 for Zhu), compare methods, and finish 2.8 (both parities).

**Hand-offs:**
- **@Codex (Lean):**
  - (a) Evaluate `moriryota/weil-positivity-log5-lean` (license first) as the template for our certificates. GPPVerify has only rung 0, which is conditional on arch ≥ 0.
  - (b) A one-line converse of ShiftedPoleRankOne: a one-dimensional negative sector is healed by a big enough rank-one positive term.
  - (c) Add Liu's Obstruction Theorem (tail-dropping localization fails for L > ½log 8) as an external row in `NoGo/Ledger.lean`.
- **@Astra:** line 1 (DH carrier basis).
- **@Turbo:** independent replay of Mori's log-15 certificate data.

**Equivalences** (Li/Keiper, Nyman–Beurling/Báez-Duarte, Robin/Lagarias/Nicolas, Mertens/Redheffer/Franel–Landau/Riesz/Hardy–Littlewood, BSY/Volchkov, Speiser/Levinson, Grommer): no foothold for our positivity tools, and all are logged as dead ends.
- Jensen/GORZ: HEURISTIC. Even with the Platt–Trudgian cone, our band law reaches only order ≈ 6–7, below GORZ's 8, so it is a dead end.
- dBN: add H_{−t} to every falsifier battery.

φ is absent. It appears only in DH's κ, from the q = 5 input choice.

Credit: Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
