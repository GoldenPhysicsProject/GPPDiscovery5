"""New Bot S2: even-sector Galerkin of the Queen/Celestial Levy-Khintchine window form (L = HALF-width, f on [-L,L]).
  Q(f) = Om0||f||^2 + int_0^inf(||f||^2 - G(u)) nu(u) du - 2 sum_{n<e^{2L}} c(n) n^{-1/2} G(log n) + pole,  G = f*f~ on [0,2L].
Basis cos(a_m x), a_m=(2m+1)pi/(2L) (same as s75/dhritz.py; Gram = L I; closed-form F_mn from dhritz docstring).
Forms: EPS (x^2+5y^2; Celestial's Om0, nu=1/sinh(u/2), c(n) via her rQ/logderiv, pole as zeta), DH (Queen's data, no pole), ZETA.
mpmath dps 50; NUMERICAL. Usage: lkritz.py FORM L K"""
import sys, json, math, mpmath as mp, numpy as np
mp.mp.dps = 50
FORM, L, K = sys.argv[1], mp.mpf(sys.argv[2]), int(sys.argv[3])
NN = int(math.exp(2*float(L))) + 2
if FORM == 'EPS':
    src = open('epstein_sinh_lib_celestial_copy.py').read().split('lam = np.zeros')[0].replace('N = 170000', 'N = %d' % NN)
    g = {}; exec(src, g); r = g['rQ']((1, 0, 5)); coef = g['logderiv'](r / r[1])
    Om0 = mp.log(mp.mpf(20)/(4*mp.pi**2)) + 2*mp.digamma(mp.mpf(1)/2); nu = lambda u: 1/mp.sinh(u/2); PF = 1
elif FORM == 'DH':
    src = open('/workspace/newbot/queen/dh_sinh_queen_copy.py').read().split('lam = np.zeros')[0].replace('N = 400000', 'N = %d' % NN)
    g = {}; exec(src, g); coef = g['c']
    Om0 = mp.digamma(mp.mpf(3)/4) - mp.log(mp.pi/5); nu = lambda u: mp.e**(-u/2)/mp.sinh(u); PF = 0
else:
    coef = np.zeros(NN+1)
    for n in range(2, NN+1):
        for p in range(2, n+1):
            if n % p == 0:
                m = n
                while m % p == 0: m //= p
                if m == 1: coef[n] = math.log(p)
                break
    Om0 = mp.digamma(mp.mpf(1)/4) - mp.log(mp.pi); nu = lambda u: mp.e**(u/2)/mp.sinh(u); PF = 1
def mp_logderiv(aa):   # same recursion as Queen/Celestial logderiv, in mpmath (float64 c(n) caps eigenvalues at ~1e-16)
    c = [mp.mpf(0)]*(NN+1)
    for d in range(2, NN+1):
        c[d] = aa[d]*mp.log(d) - mp.fsum(c[e]*aa[d//e] for e in range(2, d) if d % e == 0)
    return c
if FORM == 'EPS': coef = mp_logderiv([mp.mpf(int(t))/int(r[1]) for t in r])
if FORM == 'DH':
    s5 = mp.sqrt(5); xi = (mp.sqrt(10 - 2*s5) - 2)/(s5 - 1); per = {0: 0, 1: mp.mpf(1), 2: xi, 3: -xi, 4: mp.mpf(-1)}
    coef = mp_logderiv([mp.mpf(0)] + [per[n % 5] for n in range(1, NN+1)])
if FORM == 'ZETA': coef = [mp.mpf(t) if t == 0 else mp.log(round(math.exp(t))) for t in coef]
ns = [n for n in range(2, NN+1) if coef[n] != 0 and mp.log(n) < 2*L]
w = {n: 2*coef[n]/mp.sqrt(n) for n in ns}; lg = {n: mp.log(n) for n in ns}
a = [(2*m+1)*mp.pi/(2*L) for m in range(K)]
pts = lambda f: mp.linspace(0, 2*L, int(4*f*L/mp.pi) + 8)
tail = mp.quad(nu, [2*L, mp.inf])
Is = []; Dg = []
for m in range(K):
    am = a[m]; P = pts(am)
    Is.append(mp.quad(lambda u: mp.sin(am*u)*nu(u), P))
    Dg.append(mp.quad(lambda u: (L - (2*L-u)*mp.cos(am*u)/2 - mp.sin(am*u)/(2*am))*nu(u), P))
def Fmn(m, n, u):
    if m == n: return (2*L-u)*mp.cos(a[m]*u)/2 + mp.sin(a[m]*u)/(2*a[m])
    return (-1)**(m+n)*(a[m]*mp.sin(a[n]*u) - a[n]*mp.sin(a[m]*u))/(a[m]**2 - a[n]**2)
cp = [PF*2*a[m]*(-1)**m*mp.cosh(L/2)/(a[m]**2 + mp.mpf(1)/4) for m in range(K)]
W = mp.matrix(K, K)
for m in range(K):
    for n in range(m, K):
        if m == n: arch = Om0*L + L*tail + Dg[m]
        else: arch = -(-1)**(m+n)*(a[m]*Is[n] - a[n]*Is[m])/(a[m]**2 - a[n]**2)
        W[m, n] = W[n, m] = arch - mp.fsum(w[k]*Fmn(m, n, lg[k]) for k in ns) + 2*cp[m]*cp[n]
ev, U = mp.eigsy(W/L); o = sorted(range(K), key=lambda t: ev[t])
v = [U[i, o[0]] for i in range(K)]
json.dump(dict(form=FORM, L=str(L), K=K, ev=[mp.nstr(ev[o[i]], 12) for i in range(min(6, K))], vec=[mp.nstr(t, 30) for t in v]),
          open('lkvec_%s_L%s_K%d.json' % (FORM, sys.argv[2], K), 'w'))
print(json.dumps(dict(form=FORM, L=sys.argv[2], K=K, ev=[mp.nstr(ev[o[i]], 10) for i in range(min(4, K))])), flush=True)
