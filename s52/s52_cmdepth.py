# Grok S52 H52b(2): complete-monotonicity detection depth for a planted off-line quartet (NUMERICAL).
# K_pl(tau) = sum_{j != j0} e^{-a_j tau} + 2 Re e^{w tau},  a_j = gamma_j^2, w=(delta+i gamma_j0)^2  (quartet REPLACES zero j0,
# same convention as S51 Hankel depth).  Order-k CM test: D_k(tau) = (-1)^k K^(k) = sum a^k e^{-a tau} + 2 Re (-w)^k e^{w tau}.
# k* = least k with min_tau D_k < 0.  Log-domain numpy; truncation to 3000 zeros only drops positive terms.
import sys, numpy as np
g = np.array([float(l.split()[0]) for l in open('../s51/zeros_arb_3000.txt')])
def Dmin(k, j0, delta, ngrid=None):
    g0 = g[j0-1]; w = complex(delta**2 - g0**2, 2*delta*g0); mw = -w
    lam0 = g0**2
    # tau window: Post-Widder centre lambda=k/tau over lam0*[e^-s, e^s], s = 12/sqrt(k)+ phase room
    s = min(3.0, 12/np.sqrt(max(k, 1)) + np.pi*g0/(2*delta*max(k, 1)))
    taus = np.linspace(k/(lam0*np.exp(s)), k/(lam0*np.exp(-s)), 1) if k == 0 else None
    tlo, thi = k/(lam0*np.exp(s)), k/(lam0*np.exp(-s))
    dt = 0.05/(2*delta*g0)
    n = int(min(4e5, max(2000, (thi-tlo)/dt)))
    T = np.linspace(tlo, thi, n)
    a = np.delete(g**2, j0-1)
    lw = k*np.log(abs(mw)) + w.real*T            # log|planted|
    ph = k*np.angle(mw) + w.imag*T
    best = (np.inf, None)
    for c in range(0, n, 2000):
        Tc = T[c:c+2000]
        L = k*np.log(a)[None, :] - a[None, :]*Tc[:, None]
        m = np.maximum(L.max(axis=1), lw[c:c+2000])
        pos = np.exp(L - m[:, None]).sum(axis=1)
        neg = 2*np.exp(lw[c:c+2000] - m)*np.cos(ph[c:c+2000])
        r = (pos + neg)/(pos + 2*np.exp(lw[c:c+2000]-m))
        i = r.argmin()
        if r[i] < best[0]: best = (r[i], Tc[i])
    return best
def spacing(x): return 2*np.pi/np.log(x/(2*np.pi))
def kstar(j0, delta, kmax=10**7):
    # geometric scan then bisection on the first sign change (reports first k on grid found negative)
    ks = [0, 1, 2, 4] + [int(round(4*1.25**i)) for i in range(1, 200) if 4*1.25**i <= kmax]
    prev = None
    for k in ks:
        r, t = Dmin(k, j0, delta)
        if r < 0:
            lo, hi = (prev if prev is not None else -1), k
            while hi - lo > max(1, hi//200):
                mid = (lo+hi)//2
                if Dmin(mid, j0, delta)[0] < 0: hi = mid
                else: lo = mid
            r, t = Dmin(hi, j0, delta)
            return hi, r, t
        prev = k
    return None, None, None
if __name__ == '__main__':
    targets = [int(x) for x in sys.argv[1].split(',')]
    DL = sys.argv[2].split(',')
    print('j0 gamma spacing | delta: k*  (min ratio, tau*, k*/gamma^2*delta^2, lam=k/tau / gamma^2)')
    for j0 in targets:
        g0 = g[j0-1]; sp = spacing(g0); row = []
        for lab in DL:
            d = float(lab[3:]) if lab.startswith('abs') else float(lab[:-2])*sp
            k, r, t = kstar(j0, d)
            if k is None: row.append('%s:none'%lab); continue
            row.append('%s(d=%.3g):k*=%d r=%.2e tau=%.4g c=%.3g lamrel=%.4f' % (lab, d, k, r, t, k*d*d/g0**2, (k/t)/g0**2))
        print(j0, '%.4f'%g0, '%.4f'%sp, '|', '  '.join(row), flush=True)
