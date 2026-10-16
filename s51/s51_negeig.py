# Grok S51: magnitude of the negative eigenvalue of H_N (unscaled) at/after detection depth; #neg eigenvalues (NUMERICAL)
import mpmath as mp
mp.mp.dps = 300
gam = [mp.mpf(l.split()[0]) for l in open('zeros_arb_3000.txt')]
def H(k, delta, N, plant=True):
    xs = [1/g**2 for g in gam]
    w = 1/(-(mp.mpc(delta, gam[k-1]))**2)
    pw = []
    for p in range(2*N+1):
        v = mp.fsum(x**p for x in xs)
        if plant: v += -xs[k-1]**p + 2*mp.re(w**p)
        pw.append(v)
    return mp.matrix([[pw[i+j+2] for j in range(N)] for i in range(N)])
for k, d, Ns in [(4, mp.mpf('0.1'), [11, 12, 13, 16]), (4, mp.mpf('0.4'), [9, 10, 11]), (10, mp.mpf('0.1'), [22, 23, 24])]:
    for N in Ns:
        E = sorted(mp.eigsy(H(k, d, N), eigvals_only=True))
        E0 = sorted(mp.eigsy(H(k, d, N, False), eigvals_only=True))
        print('k=%d delta=%s N=%d  lambda_min planted=%s  (#neg=%d)  unplanted lambda_min=%s' % (k, d, N, mp.nstr(E[0], 6), sum(1 for e in E if e < 0), mp.nstr(E0[0], 6)), flush=True)
