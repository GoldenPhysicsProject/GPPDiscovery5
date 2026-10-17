"""Grok S64 H64c stage 3: window test for the backward-heat fake H_t (t<0), zero side (NUMERICAL).
Q_t(f) = sum over zeros gamma of H_t (gamma = z/2, all four of +-gamma, +-conj gamma) Ff(gamma) conj(Ff(conj gamma)),
the same normalisation as Mori (1.1)/our S14 for zeta (t = 0). Zhu trig basis, even: cos(a_m x), a_m=(2m+1)pi/(2L);
odd: sin(b_m x), b_m = m pi/L. Ff closed forms: even (-1)^m cos(zL) 2a/(a^2-z^2); odd i(-1)^(m+1) sin(zL) 2b/(b^2-z^2).
Zeros: real zeros of H_t with x<=300 (scan), non-real zeros (Newton), and for gamma>150 the zeta zeros as PROXY
(H_t zeros above 150 not computed; their positions drift by O(|t|) - caveat). lam_min = min eig of G/L (Gram = L I).
Tiers: gamma<=GMP in mpmath (dps), rest in float64 (all real => PSD).
usage: ht_window.py t L K parity [GMP]"""
import sys, json, math
import mpmath as mp, numpy as np
t, Ls, K, par = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]; GMP = float(sys.argv[5]) if len(sys.argv) > 5 else 1500
mp.mp.dps = 80; L = mp.mpf(Fr) if False else mp.mpf(eval(Ls))
zz = [float(l.split()[1]) for l in open('/workspace/GPPVerify/discovery/cutkosky_weil/zeta_zeros_100k.txt')]
if t == '0':
    real = [g for g in zz if g <= 150]; cplx = []
else:
    real = [x/2 for x in json.load(open(f'scan_t{t}.json'))['real_zeros_x']]
    cplx = [complex(a, b) for a, b in json.load(open(f'cplx_t{t}.json'))['gamma_form']]
prox = [g for g in zz if g > 150]
if par == 'even': fr = [(2*m+1)*mp.pi/(2*L) for m in range(K)]
else: fr = [m*mp.pi/L for m in range(1, K+1)]
def vec(z):
    z = mp.mpc(z)
    if par == 'even': return [(-1)**m*mp.cos(z*L)*2*a/(a*a-z*z) for m, a in enumerate(fr)]
    return [1j*(-1)**(m+2)*mp.sin(z*L)*2*b/(b*b-z*z) for m, b in enumerate(fr)]
G = mp.matrix(K, K)
def add(z, mult):   # adds mult * u u'^H-like term: Ff(z) conj Ff(conj z) as real symmetric form on real coeffs
    u = vec(z); w = vec(mp.conj(mp.mpc(z)))
    for i in range(K):
        for j in range(K):
            G[i, j] += mult*mp.re(u[i]*mp.conj(w[j]))
for g in real + [g for g in prox if g <= GMP]: add(g, 2)          # +-gamma
for g in cplx:                                                      # +-gamma, +-conj gamma
    add(g, 2); add(g.conjugate(), 2)
# float tier
Gf = np.zeros((K, K)); frf = np.array([float(a) for a in fr]); Lf = float(L)
pf = np.array([g for g in prox if g > GMP])
for chunk in np.array_split(pf, 50):
    z = chunk[:, None]
    if par == 'even': V = (-1.0)**np.arange(K)[None, :]*np.cos(z*Lf)*2*frf[None, :]/(frf[None, :]**2-z**2)
    else: V = (-1.0)**(np.arange(K)[None, :]+2)*np.sin(z*Lf)*2*frf[None, :]/(frf[None, :]**2-z**2)
    Gf += 2*V.T @ V
tail_est = float(np.trace(Gf))
for i in range(K):
    for j in range(K): G[i, j] += mp.mpf(Gf[i, j])
Gs = (G + G.T)/2/L
E, Q = mp.eigsy(Gs)
lam = min(E)
# pair (non-real zeros) contribution along the minimising vector
k = [i for i in range(K) if E[i] == lam][0]; v = [Q[i, k] for i in range(K)]
def form(zlist):
    s = 0
    for g in zlist:
        for zc in (g, g.conjugate()):
            u = vec(zc); w = vec(mp.conj(mp.mpc(zc)))
            s += 2*mp.re(mp.fsum(u[i]*v[i] for i in range(K))*mp.conj(mp.fsum(w[i]*v[i] for i in range(K))))
    return s/L
pc = form(cplx) if cplx else 0
r = {"t": t, "L": Ls, "K": K, "parity": par, "lam_min": mp.nstr(lam, 8), "pair_part_on_minvec": mp.nstr(pc, 8),
     "n_real": len(real), "n_cplx_pairs": len(cplx), "float_tier_trace": tail_est}
print(json.dumps(r)); open('window.jsonl', 'a').write(json.dumps(r)+'\n')
