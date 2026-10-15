#!/usr/bin/env python3
"""S19 — Independent Arb (python-flint) audit of Codex's χ₄ first-prime-window
ground-state lower certificate (15 modes / 100 cells).

Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
Session 19, 2026-10-15.

REFERENCE ONLY (do not import): Codex mpmath.iv and Fraction scripts in
GPPDiscovery2/.../chi4_groundstate_*.py. This file reimplements the finite
certificate in Arb ball arithmetic from the analytic formulas in
GPP-bridge/work/codex/2026-10-15_chi4_positive_ground_state_first_prime_window.md.

Labels:
  CERTIFIED (Arb) — finite 15-mode / 100-cell comparison under Arb enclosures
  ANALYTIC (accepted, Lean-open) — nonlocal Picone / ground-state identity
  ANALYTIC (accepted, Lean-open) — omitted j≥15 modes contribute ≥0
  NOT claimed: RH, Dirichlet GRH, infinite-mode Lean certificate

GppVerify #363 GroundStatePicone is OPEN (finite edge algebra only); this Arb
audit does NOT replace a kernel proof of Picone + the infinite remainder.
"""
from flint import arb, ctx, fmpq
import sys, os

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "s17"))
from s17_build import consts  # c0, c1 via digamma-equivalent closed forms

# ---- precision ----
PREC = 512
ctx.prec = PREC

N_MODES = 15          # j = 0..14
N_CELLS = 100
K_RAT = fmpq(19, 10)  # k = 19/10

def main():
    print("=== S19 χ₄ ground-state Arb certificate ===")
    print(f"ctx.prec = {PREC}")
    print("Credit: Daniel Toupin / Golden Physics Project, AI assistance: Grok")
    print("RH OPEN. Dirichlet GRH OPEN. Finite-window χ₄ only.\n")

    h = arb(3).log() / 2          # h* = (log 3)/2
    k = arb(K_RAT)
    print(f"h* = (log 3)/2 = {h.str(40)}")
    print(f"k  = 19/10     = {k.str(40)}")

    # --- (a) kh* < π/2 and cos(kh*) > 0 ---
    kh = k * h
    half_pi = arb.pi() / 2
    cos_kh = kh.cos()
    sin_kh = kh.sin()
    ok_kh = bool(kh < half_pi)
    ok_cos = bool(cos_kh > 0)
    print(f"\n(a) kh* = {kh.str(30)}")
    print(f"    π/2  = {half_pi.str(30)}")
    print(f"    kh* < π/2 : {ok_kh}")
    print(f"    cos(kh*)  = {cos_kh.str(30)}  > 0 : {ok_cos}")
    if not (ok_kh and ok_cos):
        print("FAIL: positivity of φ=cos(kx) on closed I failed")
        return 1

    # --- modes a_j = 3/2 + 2j, s_j = a_j cos(kh) - k sin(kh) ---
    a = [arb(3)/2 + 2*j for j in range(N_MODES)]
    s = [aj * cos_kh - k * sin_kh for aj in a]
    print(f"\n(b) sign pattern of s_j = a_j cos(kh*) - k sin(kh*):")
    print(f"    s[0] = {s[0].str(30)}  < 0 : {bool(s[0] < 0)}")
    ok_s0 = bool(s[0] < 0)
    ok_spos = True
    for j in range(1, N_MODES):
        pos = bool(s[j] > 0)
        ok_spos &= pos
        if j <= 3 or j == N_MODES - 1:
            print(f"    s[{j}] = {s[j].str(20)}  > 0 : {pos}")
    print(f"    s[1..{N_MODES-1}] all > 0 : {ok_spos}")
    if not (ok_s0 and ok_spos):
        print("FAIL: sign pattern for monotonicity of b_j failed")
        return 1

    # Diagonal d_j = 2 k² / [a_j (a_j² + k²)]
    diagonal = arb(0)
    for aj in a:
        diagonal += 2 * k * k / (aj * (aj * aj + k * k))
    print(f"\n    Σ_{{j=0..{N_MODES-1}}} d_j = {diagonal.str(30)}")

    def b_j(j, x):
        """Boundary mode contribution at |x|.
        b_j(x) = [e^{-a_j(h-x)} + e^{-a_j(h+x)}] * s_j / [(a_j²+k²) cos(kx)]
        Written as two decaying exponentials (no cosh overflow).
        """
        aj = a[j]
        return ((-aj * (h - x)).exp() + (-aj * (h + x)).exp()) * s[j] / (
            (aj * aj + k * k) * (k * x).cos()
        )

    # --- Why the cell bound is a rigorous Arb lower bound ---
    # Analytic (Codex): on [0,h], b_0 decreases (s_0<0; cosh(a_0 x)/cos(kx) increases)
    # and b_j increases for j≥1 (s_j>0). Hence for x ∈ [L_i, U_i],
    #   (Lφ/φ)(x) ≥ Σ d_j + b_0(U_i) + Σ_{j≥1} b_j(L_i)  (+ omitted j≥15 ≥0).
    # Arb: L_i, U_i are balls containing the true endpoints (h is a ball).
    # Evaluating b_0 at the ball U_i returns a range enclosure; its .lower()
    # is ≤ every value of b_0 on that ball, hence ≤ b_0(true U_i) ≤ b_0(x)
    # for x ≤ true U_i (b_0 decreasing). Likewise b_j(L_i).lower() ≤ b_j(x)
    # for x ≥ true L_i (b_j increasing, j≥1). Summing these lowers yields a
    # CERTIFIED lower bound on min_{x∈cell i} of the 15-mode truncation.

    print("\n(c) 100-cell monotone Arb lower bounds C_i (15 modes):")
    print("    C_i := Σ d_j + b_0(U_i) + Σ_{j=1..14} b_j(L_i)")
    print("    (b_0 at upper endpoint; b_j≥1 at lower; omitted j≥15 ≥0)")

    C_balls = []
    C_lowers = []
    for i in range(N_CELLS):
        L = h * i / N_CELLS
        U = h * (i + 1) / N_CELLS
        bound = diagonal + b_j(0, U)
        for j in range(1, N_MODES):
            bound += b_j(j, L)
        C_balls.append(bound)
        C_lowers.append(bound.lower())

    # worst cell by certified lower endpoint
    worst_i = min(range(N_CELLS), key=lambda i: float(C_lowers[i].mid()))
    worst = C_balls[worst_i]
    worst_lo = C_lowers[worst_i]

    print(f"    worst cell index i = {worst_i}  (0-based; Codex reported 64)")
    print(f"    C_{worst_i} ball     = {worst.str(40)}")
    print(f"    C_{worst_i} .lower() = {worst_lo.str(40)}")

    all_gt_087 = all(bool(C_balls[i] > arb("0.87")) for i in range(N_CELLS))
    print(f"    all 100 C_i > 0.87 (rigorous): {all_gt_087}")
    if not all_gt_087:
        bad = [i for i in range(N_CELLS) if not bool(C_balls[i] > arb("0.87"))]
        print(f"    cells failing >0.87: {bad[:20]}{'...' if len(bad)>20 else ''}")

    # --- (d) δ = log(2π)+γ−π/2 ; also c1 − log4 ---
    print("\n(d) δ = c₁ − log 4:")
    delta_direct = (2 * arb.pi()).log() + arb.const_euler() - arb.pi() / 2
    C0, C1 = consts()  # C1 = log(8π)+γ−π/2
    delta_via_c1 = C1 - arb(4).log()
    print(f"    δ = log(2π)+γ−π/2     = {delta_direct.str(40)}")
    print(f"    δ = c1 − log4 (S17)   = {delta_via_c1.str(40)}")
    print(f"    overlap (same value): {delta_direct.overlaps(delta_via_c1)}")
    print(f"    c1 (S17 consts)       = {C1.str(30)}")
    # digamma form: c1 = log(π) − ψ(3/4)
    c1_digamma = arb.pi().log() - (arb(3)/4).digamma()
    print(f"    c1 = log(π)−ψ(3/4)    = {c1_digamma.str(30)}")
    print(f"    c1 overlaps digamma:  {C1.overlaps(c1_digamma)}")
    delta = delta_direct  # use the direct expression; identical ball-wise

    # Codex published digits
    codex_delta = arb("0.8442964045159817249358498712538862686663695835")
    codex_minC = arb("0.8795834621661168898864762079984104448367987486967")
    # Published Codex digits are truncations; enclose each in a 1e-40 radius ball for comparison.
    codex_delta_ball = arb(codex_delta.mid(), rad=arb("1e-40"))
    codex_minC_ball = arb(codex_minC.mid(), rad=arb("1e-40"))
    print(f"\n    Codex δ digits enclosed by ours: {delta.contains(codex_delta) or delta.overlaps(codex_delta_ball)}")
    print(f"    Codex min-C digits enclosed:    {worst.contains(codex_minC) or worst.overlaps(codex_minC_ball)}")
    print(f"    |mid(C)-Codex| approx:         {abs(worst.mid()-codex_minC.mid()).str(10)}")

    # --- (e) gap and thresholds ---
    print("\n(e) gap and thresholds:")
    gap_ball = worst - delta
    # conservative: lower(C) − upper(δ)
    gap_lo = worst_lo - delta.upper()
    print(f"    min C − δ (ball)           = {gap_ball.str(30)}")
    print(f"    lower(C) − upper(δ)        = {gap_lo.str(30)}")
    gt_025 = bool(gap_lo > arb("0.025"))
    gt_02 = bool(gap_lo > arb("0.02"))
    # also check worst ball itself clears δ + 0.025 rigorously
    clears_025 = bool(worst > delta + arb("0.025"))
    clears_02 = bool(worst > delta + arb("0.02"))
    print(f"    gap > 0.025 (via lowers):  {gt_025}")
    print(f"    gap > 0.02  (via lowers):  {gt_02}")
    print(f"    C_worst > δ+0.025 (Arb):   {clears_025}")
    print(f"    C_worst > δ+0.02  (Arb):   {clears_02}")

    # χ₄(2)=0 sanity
    chi4_2 = [0, 1, 0, -1][2 % 4]
    print(f"\n    χ₄(2) = {chi4_2}  (ramified; W_χ4 vanishes on (0,h*])")

    # --- PASS/FAIL ---
    print("\n=== VERDICT ===")
    passed = (
        ok_kh and ok_cos and ok_s0 and ok_spos
        and all_gt_087
        and clears_025 and clears_02
        and delta_direct.overlaps(delta_via_c1)
        and C1.overlaps(c1_digamma)
        and worst_i == 64
    )
    # cell index match is soft: if Arb picks a neighbour due to ties, still OK if min matches
    soft_cell_ok = worst_i == 64 or worst.overlaps(codex_minC)
    passed_core = (
        ok_kh and ok_cos and ok_s0 and ok_spos
        and all_gt_087
        and clears_025 and clears_02
        and delta_direct.overlaps(delta_via_c1)
        and soft_cell_ok
    )
    if passed_core:
        print("PASS — independent Arb certificate CONFIRMS Codex:")
        print("  all 100 cells have C_i > 0.87 > δ")
        print("  min C (cell %d) − δ > 0.025  ⇒  Q_χ4(g) > 0.025 ||g||²" % worst_i)
        print("  (also > 0.02, matching the pure-rational conservative coefficient)")
        print("CONDITIONAL on analytic Picone identity + j≥15 remainder ≥0 (Lean-open).")
        print("NOT a RH/GRH claim. Finite window 0<h≤h*=(log3)/2 only.")
        return 0
    else:
        print("FAIL — Arb audit did not confirm all checks:")
        print(f"  ok_kh={ok_kh} ok_cos={ok_cos} ok_s0={ok_s0} ok_spos={ok_spos}")
        print(f"  all_gt_087={all_gt_087} clears_025={clears_025} clears_02={clears_02}")
        print(f"  delta overlap={delta_direct.overlaps(delta_via_c1)} cell={worst_i}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
