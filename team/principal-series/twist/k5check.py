"""Does F_theta's window failure live on K <= 5 modes? Even/odd lam_min of F_theta and zeta on N = 5 and N = 24 shifted-Legendre
(S87 builder, dps 80). L = half-width. HIGH-PRECISION NUMERICAL."""
import sys, json, os
sys.argv = [sys.argv[0], '5']; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hp_onset_lib as H, mpmath as mp
for L in ['0.4', '0.45', '0.5', '0.6']:
    for N in (5, 24):
        Ze, Zo = H.forms('ZETA', mp.mpf(L), N); Fe, Fo = H.forms('FT', mp.mpf(L), N, mp.mpf('0.1'))
        print(json.dumps(dict(L=L, N=N, theta=0.1, zeta_e=mp.nstr(H.lmin(Ze), 6), zeta_o=mp.nstr(H.lmin(Zo), 6),
              F_e=mp.nstr(H.lmin(Fe), 6), F_o=mp.nstr(H.lmin(Fo), 6))), flush=True)
