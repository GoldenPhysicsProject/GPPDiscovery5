"""Round 17, part C: Laguerre inequality (k=1) on the real line and de Bruijn-Newman heat flow, per model.
  L_M(z) = Xi'(z)^2 - Xi(z) Xi''(z)   (>= 0 for all real z is necessary for Xi in Laguerre-Polya, i.e. RH_M)
  H_t(z) = int e^{t u^2} Phi_M(u) cos(z u) du,   Lambda_M = inf{t : H_t has only real zeros}.
Zero labels (from ../epstein_zeros_lowband.json, DH literature values, zeta's first zeros) are used ONLY as
Newton starting points to follow off-line pairs to their collision time t* (=> Lambda_M >= t*). NUMERICAL.
Output r17_heat.json"""
import sys, os, json, time
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import mpmath as mp
from laws17 import *
DPS = 60
H = mp.mpf('0.01')

def phigrid(name):
    mp.mp.dps = DPS
    if name in ('ZETA', 'EPS', 'DH'): return grid(Model(name), H)
    us, ph = grid(Model('ZETA'), H)
    uZ = [-u for u in us[:0:-1]] + us; pZ = ph[:0:-1] + ph; c0 = len(uZ)//2
    par = mp.mpf('0.2') if name == 'FTH' else mp.mpf(3)
    nv = int(3.4/float(H)); vs = []; vals = []
    for j in range(nv + 1):
        v = mp.mpf(j)*H; s = []
        for i in range(len(uZ)):
            kk = j - (i - c0) + c0
            if 0 <= kk < len(uZ):
                u = uZ[i]; w = pZ[i]*pZ[kk]
                w *= (mp.exp(par*(v - 2*u)) if name == 'FTH' else mp.cos(par*(v - 2*u)))
                s.append(w)
        vs.append(v); vals.append(H*mp.fsum(s))
    return vs, vals

ZMAX = dict(ZETA=60, EPS=60, DH=100, FTH=40, FTAU=40)
LABELS = dict(EPS=[(0.9329696974854141, 15.668249531278473), (0.9376669067003994, 29.983395235155516),
                   (0.6969270453122914, 36.37406368638513), (0.8231873164780866, 44.000113180236895),
                   (0.6349076019151867, 46.7584086570051), (0.5974373015976651, 52.743266408321766),
                   (0.8111720812060284, 58.507737327731896)],
              DH=[(0.808517, 85.699348), (0.650830, 114.163343)],
              FTH=[(0.7, 14.134725), (0.7, 21.022040), (0.7, 25.010858)])

def task(name):
    mp.mp.dps = DPS
    t0 = time.time()
    us, ph = phigrid(name)
    def X(z, t=0, d=0): return Xi_heat(us, ph, H, z, t, d)
    rec = {}
    # Laguerre scan
    step = mp.mpf('0.05'); z = mp.mpf(0); neg = []; minrel = mp.inf; zmin = None
    while z <= ZMAX[name]:
        a, b, c = X(z), X(z, 0, 1), X(z, 0, 2)
        L = b*b - a*c; rel = L/(b*b + abs(a*c))
        if rel < minrel: minrel, zmin = rel, z
        if L < 0: neg.append(float(z))
        z += step
    # group negative z into intervals
    iv = []
    for x in neg:
        if iv and x - iv[-1][1] < 0.051: iv[-1][1] = x
        else: iv.append([x, x])
    rec['laguerre'] = dict(z_range=[0, ZMAX[name]], step=0.05, negative_intervals=iv, min_rel=mp.nstr(minrel, 6),
                           at=float(zmin), passes=not iv)
    # heat-flow tracking of off-line pairs
    tr = []
    for beta, gam in LABELS.get(name, []):
        zz = mp.mpc(gam, -(beta - 0.5)); t = mp.mpf(0); hist = []; dt = mp.mpf('0.002'); ok = True
        while True:
            for _ in range(40):
                f = X(zz, t); fp = X(zz, t, 1); dz = f/fp; zz -= dz
                if abs(dz) < mp.mpf(10)**-25: break
            if abs(dz) > 1e-8: ok = False; break
            hist.append((float(t), float(zz.real), float(abs(zz.imag))))
            if t == 0: start = (float(zz.real), float(zz.imag))
            if abs(zz.imag) < 0.02 or t > 1.0: break
            # adapt step: Im^2 roughly linear in t
            if len(hist) >= 2:
                (t1, _, y1), (t2, _, y2) = hist[-2], hist[-1]
                slope = (y2*y2 - y1*y1)/(t2 - t1) if t2 > t1 else -1
                if slope < 0: dt = min(mp.mpf(0.05), max(mp.mpf('1e-4'), mp.mpf(0.5*y2*y2/(-slope))))
            t += dt
        if len(hist) >= 2:
            (t1, _, y1), (t2, _, y2) = hist[-2], hist[-1]
            tstar = t2 + y2*y2*(t2 - t1)/(y1*y1 - y2*y2) if y1 > y2 else None
        else: tstar = None
        tr.append(dict(label=[beta, gam], H0_zero=start if hist else None, converged=ok, steps=len(hist),
                       last=hist[-1] if hist else None, t_star=tstar))
    rec['heat_tracks'] = tr
    rec['Lambda_lower_bound'] = max([x['t_star'] for x in tr if x['t_star'] is not None], default=None)
    rec['seconds'] = round(time.time() - t0, 1)
    print(name, json.dumps(rec)[:600], flush=True)
    return name, rec

if __name__ == '__main__':
    with Pool(5) as p: res = dict(p.map(task, ['ZETA', 'EPS', 'DH', 'FTH', 'FTAU']))
    json.dump(res, open(os.path.join(HERE, 'r17_heat.json'), 'w'), indent=1)
