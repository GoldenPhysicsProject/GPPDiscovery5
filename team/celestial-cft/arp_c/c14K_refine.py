"""C.14.6 refinement study at T = 14 (and T = 10), n_s = 2..32 (NUMERICAL). RH is OPEN.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
import sys
from c14K_gamma_test import gamma, gstar, a
for T in (14.0, 10.0):
    for sector in ('even', 'odd'):
        print(f"== T={T} sector={sector} gamma*={gstar(T):.4f}"); sys.stdout.flush()
        for lt in (0.01, 0.02, 0.03, 0.06, 0.1, 0.119):
            M0 = round(a*2/lt) + 2
            gs = []
            for ns in (2, 4, 8, 16, 32):
                M = M0*ns//2
                if M > 2400: break
                gs.append(gamma(M, ns, T, sector)['gamma'])
            ell = 2*a/(M0-2)
            d = [gs[i+1]-gs[i] for i in range(len(gs)-1)]
            # geometric extrapolation from last three
            ext = gs[-1] + d[-1]**2/(d[-2]-d[-1]) if len(d) >= 2 and d[-2] > d[-1] > 0 else float('nan')
            print(f"  ell={ell:.4f} gammas(ns=2,4,..)=" + ", ".join(f"{g:.4f}" for g in gs) + f"  extrap={ext:.4f}"); sys.stdout.flush()
