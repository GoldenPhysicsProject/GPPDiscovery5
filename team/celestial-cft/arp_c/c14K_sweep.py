"""C.14.6 sweep driver (NUMERICAL). See c14K_gamma_test.py. RH is OPEN.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
import numpy as np, sys
from c14K_gamma_test import gamma, gstar, a
ells = [0.01, 0.03, 0.06, 0.1, 0.119]
Ts = [10.0, 14.0, 20.0]
for sector in ('even', 'odd'):
    for T in Ts:
        print(f"== sector {sector}, T = {T}, gamma* = {gstar(T):.4f}"); sys.stdout.flush()
        for lt in ells:
            row = []
            for ns0 in (2, 4, 8):
                M0 = round(a*2/lt) + 2
                ns, M = ns0, M0*ns0//2
                r = gamma(M, ns, T, sector)
                row.append(r)
            print(f"  ell={row[0]['ell']:.4f}  " + "  ".join(f"[ns={2*2**i:d},M={round(a*2/lt)*2**i+2*2**i:d}] g={r['gamma']:.6f} (condU {r['condU']:.1e})" for i, r in enumerate(row)))
            sys.stdout.flush()
