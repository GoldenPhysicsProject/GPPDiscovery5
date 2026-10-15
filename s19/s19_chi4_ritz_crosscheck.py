#!/usr/bin/env python3
"""S19 optional Ritz sanity check for Q_χ4 at h*=(log3)/2.

NUMERICAL: Ritz λ_min is an UPPER bound on the true Rayleigh infimum
(finite subspace). It is NOT a lower bound on coercivity. Use only as a
sanity check that the form looks positive and that λ_min ≫ 0.02.

Daniel Toupin / Golden Physics Project, AI assistance: Grok. RH/GRH OPEN.
"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "s18"))
from s18_ritz_common import *
from s18_common import forms

def main(Ms):
    print("=== S19 χ₄ Ritz cross-check at h* (NUMERICAL; upper on λ_min) ===")
    print("Credit: Daniel Toupin / Golden Physics Project, AI assistance: Grok")
    print("RH OPEN. GRH OPEN. Not a coercivity lower bound.\n")
    C0, C1 = consts()
    delta = C1 - arb(4).log()
    print(f"h* = {HSTAR.str(20)}")
    print(f"δ = c1−log4 = {delta.str(20)}")
    print(f"target coefficients 0.025 / 0.02 (from ground-state certificate)\n")
    for M in Ms:
        B = build(M, HSTAR)
        F = forms(B, HSTAR)
        # forms() key for chi4
        Q = F["chi4 (mod 4, odd)"]
        N = B["N"]
        print(f"--- M={M} ---")
        for sn, I in SECT(M):
            lam, cert, fx = ritz_min(Q, N, I, extra=B["T"].get(2))
            print(f"  sector {sn:4s}: Ritz λ_min = {mp.nstr(lam, 12)}  "
                  f"Arb Rayleigh upper = {cert.str(12)}  "
                  f"F(log2)/N = {fx.str(8) if fx is not None else 'n/a'}")
            # sanity: Ritz should sit well above 0.02 if certificate is true
            above = bool(cert > arb("0.02"))
            print(f"           Arb Rayleigh > 0.02: {above}  (sanity only; not a lower proof)")
    print("\nLabel: NUMERICAL Ritz / CERTIFIED Rayleigh upper witnesses.")
    print("Ground-state Arb lower certificate is in s19_chi4_arb.py.")

if __name__ == "__main__":
    Ms = list(map(int, sys.argv[1:])) or [32, 48]
    main(Ms)
