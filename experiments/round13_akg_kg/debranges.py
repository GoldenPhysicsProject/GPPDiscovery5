"""Round 13f: de Branges / Lagarias test of the 'projection curvature' R = dd^c log K(z,z).
E_h(z) = xi_F(1/2 + h - i z), s = 1/2 - i z, y = Im z = sigma - 1/2, Theta = E#/E.
 ZETA, EPS (functional equation): Theta = xi_F(s-h)/xi_F(s+h);  HB defect D_h = log|xi(s+h)| - log|xi(s-h)|  (K(z,z) > 0 <=> D_h > 0)
 BEUR (spliced, no FE): E# = xi_B(1+h-s) (real coefficients), xi_B := s(s-1) Gamma_R(s) zeta_B(s) with zeta's gamma borrowed;
   zeta_B only continued to Re s > 0 here, so the kernel is evaluated only where 1+h-sigma >= 0.25.
Curvature (exact identity): dd^c log K(z,z) = 1/(4y^2) - |Theta'|^2/(1-|Theta|^2)^2 on {|Theta| < 1};
 Schwarz-Pick ratio SP = 2y|Theta'|/(1-|Theta|^2); curvature >= 0 <=> SP <= 1.  Zeros are NOT used here."""
import mpmath as mp, numpy as np, json, sys, math
mp.mp.dps = 15
which = sys.argv[1]
SIG = [0.51, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2]; HS = [0.05, 0.1, 0.2, 0.3, 0.45, 0.5]
TS = np.round(np.arange(0.5, 60.0001, 0.1), 4)
SIGC = [0.52, 0.55, 0.6, 0.7, 0.8, 1.0]; HC = [0.1, 0.3]; TC = np.round(np.arange(0.5, 60.0001, 0.25), 4)

if which in ('ZETA', 'EPS'):
    if which == 'ZETA':
        def logxi(s): return mp.log(s*(s-1)/2) - s/2*mp.log(mp.pi) + mp.loggamma(s/2) + mp.log(mp.zeta(s))
        def dlogxi(s): return 1/s + 1/(s-1) - mp.log(mp.pi)/2 + mp.digamma(s/2)/2 + mp.zeta(s, derivative=1)/mp.zeta(s)
    else:
        src = open('epstein_zeros.py').read().split('# coefficient check')[0]; ns = {}; exec(src, ns); Z = ns['Z']
        A = mp.sqrt(20)/(2*mp.pi)
        def logxi(s): return mp.log(s*(s-1)) + s*mp.log(A) + mp.loggamma(s) + mp.log(Z(s))
        def dlogxi(s): return 1/s + 1/(s-1) + mp.log(A) + mp.digamma(s) + mp.diff(Z, s)/Z(s)
    cache = {}
    def lx(sig, t):
        k = (round(sig, 4), float(t))
        if k not in cache: cache[k] = float(mp.re(logxi(mp.mpc(sig, t))))
        return cache[k]
    D = {}
    for h in HS:
        for sig in SIG:
            D[(h, sig)] = np.array([lx(sig+h, t) - lx(sig-h, t) for t in TS])
    def theta_parts(sig, t, h):
        sp, sm = mp.mpc(sig+h, t), mp.mpc(sig-h, t)
        logT = logxi(sm) - logxi(sp); T = abs(mp.exp(logT)); dT = T*abs(dlogxi(sm) - dlogxi(sp))
        return float(T), float(dT)
else:
    import harness as h_
    Zm = h_.zeta_model(); B = h_.beurling_model(x0=55.0)
    m = Zm['n'] < 55
    n = np.concatenate([Zm['n'][m], B['n']]); c = np.concatenate([Zm['c'][m], B['c']])
    wts = c/np.log(n); logn = np.log(n); ell = math.log(n.max()); q = B['quartet']
    from scipy.special import exp1, loggamma, digamma
    def E1c(w): return complex(mp.e1(w))
    def logzB(s):  # sum_{atoms<=X} n^-s/m + int_X^inf x^-s dPi_C (analytic tail)
        v = np.sum(wts*np.exp(-s*logn))
        return v + E1c((s-1)*ell) - sum(E1c((s-r)*ell) for r in q)
    def dlogzB(s):
        v = -np.sum(c*np.exp(-s*logn))
        return v - np.exp(-(s-1)*ell)/(s-1) + sum(np.exp(-(s-r)*ell)/(s-r) for r in q)
    def logxi(s): return complex(np.log(s*(s-1)/2) - s/2*math.log(math.pi) + loggamma(s/2)) + logzB(s)
    def dlogxi(s): return 1/s + 1/(s-1) - math.log(math.pi)/2 + digamma(s/2)/2 + dlogzB(s)
    SIG = [0.51, 0.55, 0.58, 0.6, 0.62, 0.65, 0.7, 0.8]; HS = [0.05, 0.08, 0.1, 0.2, 0.3]
    SIGC = [0.52, 0.55, 0.58, 0.6, 0.65, 0.7]; HC = [0.05, 0.3]
    D = {}
    for h in HS:
        for sig in SIG:
            if 1+h-sig < 0.25: continue
            D[(h, sig)] = np.array([logxi(complex(sig+h, t)).real - logxi(complex(1+h-sig, t)).real for t in TS])
    def theta_parts(sig, t, h):
        sp, sm = complex(sig+h, t), complex(1+h-sig, -t)   # E#(z) = xi_B(1+h-s); |.| same at conj argument
        logT = logxi(sm).real - logxi(sp).real; T = math.exp(logT)
        # d/ds log Theta = -(xi'/xi)(1+h-s) - (xi'/xi)(s+h)
        dT = T*abs(-np.conj(dlogxi(sm)) - dlogxi(sp))
        return T, dT
out = dict(model=which, D={}, curv={})
for (h, sig), arr in D.items():
    i = int(np.argmin(arr)); neg = TS[arr <= 0]
    out['D'][f'h{h}_s{sig}'] = dict(min=float(arr[i]), t_at_min=float(TS[i]), n_nonpositive=int(len(neg)),
                                     t_range_nonpos=[float(neg.min()), float(neg.max())] if len(neg) else None)
for h in HC:
    for sig in SIGC:
        if which == 'BEUR' and 1+h-sig < 0.25: continue
        y = sig-0.5; worst = (-1.0, None); Kneg = []
        for t in TC:
            T, dT = theta_parts(sig, t, h)
            if T >= 1: Kneg.append(float(t)); continue
            sp = 2*y*dT/(1-T*T)
            if sp > worst[0]: worst = (sp, float(t))
        out['curv'][f'h{h}_s{sig}'] = dict(max_SP_ratio=worst[0], t_at_max=worst[1], n_K_nonpositive=len(Kneg),
                                           t_K_nonpositive=[min(Kneg), max(Kneg)] if Kneg else None)
    print(h, json.dumps({k: v for k, v in out['curv'].items() if k.startswith(f'h{h}_')}), flush=True)
json.dump(out, open(f'debranges_{which}.json', 'w'), indent=1)
print(json.dumps(out['D']))
