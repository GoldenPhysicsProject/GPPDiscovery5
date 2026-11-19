"""Round 17, part A: the BPY law of each model. Normalisation checks, positivity, log-concavity (PF_2),
PF_3 / PF_4 minors, BPY / Kuiper identification for zeta. Integer data only; zeros never used.
Output: r17_laws.json"""
import sys, os, json, random, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from laws17 import *
mp.mp.dps = 40
H = mp.mpf('0.01')
out = {}

def full(us, ph):
    """even grid -> full line arrays"""
    return [-u for u in us[:0:-1]] + us, ph[:0:-1] + ph

G = {}
for name in ['ZETA', 'EPS', 'DH']:
    G[name] = grid(Model(name), H)

# --- FTH / FTAU on the same grid by trapezoid convolution of Phi_Z (spectrally accurate) ---
uZ, pZ = full(*G['ZETA'])
idx0 = len(uZ)//2
def conv(v_index, kind, par):
    s = []
    for i in range(len(uZ)):
        k = v_index - (i - idx0)          # v - u index offset (both on same lattice)
        kk = k + idx0
        if 0 <= kk < len(uZ):
            u = uZ[i]; v = mp.mpf(v_index)*H
            w = pZ[i]*pZ[kk]
            if kind == 'FTH': w *= mp.exp(-par*u)*mp.exp(par*(v - u))
            else: w *= mp.cos(par*(v - 2*u))
            s.append(w)
    return H*mp.fsum(s)
TH, TAU = mp.mpf('0.2'), mp.mpf(3)
nv = int(3.0/float(H))
for kind, par in [('FTH', TH), ('FTAU', TAU)]:
    us = [mp.mpf(j)*H for j in range(nv + 1)]
    G[kind] = (us, [conv(j, kind, par) for j in range(nv + 1)])

def xi_from(us, ph, s):
    w = s - mp.mpf(1)/2
    return trap_even(ph[0], [ph[j]*mp.cosh(w*us[j]) for j in range(1, len(us))], H)

ref = {'ZETA': lambda s: s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s),
       'DH': lambda s: (5/mp.pi)**((s+1)/2)*mp.gamma((s+1)/2)*mp.dirichlet(s, [0, 1, kappa(), -kappa(), -1])}
ref['FTH'] = lambda s: ref['ZETA'](s - TH)*ref['ZETA'](s + TH)
ref['FTAU'] = lambda s: ref['ZETA'](s - 1j*TAU)*ref['ZETA'](s + 1j*TAU)

def eps_direct3():
    import numpy as np
    N = 1500; m = np.arange(-N, N + 1); Q = (m[:, None]**2 + 5*m[None, :]**2).astype(float); Q[N, N] = np.inf
    Z = float(np.sum(Q**-3.0)); s = 3
    return s*(s-1)*(mp.sqrt(20)/(2*mp.pi))**s*mp.gamma(s)*Z
ref_eps3 = eps_direct3()

for name, (us, ph) in G.items():
    rec = {}
    chk = []
    for s in [mp.mpf(3), mp.mpc(0.5, 2), mp.mpc(0.7, 5)]:
        a = xi_from(us, ph, s)
        b = ref[name](s) if name in ref else (ref_eps3 if s == 3 else None)
        chk.append(dict(s=str(s), from_law=mp.nstr(a, 14), reference=(mp.nstr(b, 14) if b is not None else None),
                        relerr=(float(abs(a - b)/abs(b)) if b is not None else None)))
    rec['xi_check'] = chk
    rec['xi0'] = mp.nstr(xi_from(us, ph, mp.mpf(0)), 15)
    neg = [float(us[j]) for j in range(len(us)) if ph[j] < 0]
    rec['phi0'] = mp.nstr(ph[0], 12)
    rec['min_phi'] = mp.nstr(min(ph), 8)
    rec['positive_on_grid'] = not neg
    rec['negative_u_interval'] = [-max(neg), max(neg)] if neg else None
    # log-concavity on the positive part (second differences of log phi)
    lcv = []
    for j in range(1, len(us) - 1):
        if ph[j-1] > 0 and ph[j] > 0 and ph[j+1] > 0:
            d2 = mp.log(ph[j+1]) - 2*mp.log(ph[j]) + mp.log(ph[j-1])
            if d2 > 0: lcv.append(float(us[j]))
    rec['log_concave_where_positive'] = not lcv
    rec['log_convex_points'] = lcv[:10]
    # PF_k minors det[Phi(x_i - y_j)], x,y increasing, sampled on the lattice (k=2,3,4)
    uf, pf = full(us, ph); c0 = len(uf)//2
    def P(k):
        kk = c0 + k
        return pf[kk] if 0 <= kk < len(pf) else mp.mpf(0)
    random.seed(17)
    pfres = {}
    for k in [2, 3, 4]:
        worst = None; nneg = 0; ntry = 4000
        for _ in range(ntry):
            span = random.choice([20, 40, 80, 150])
            xs = sorted(random.sample(range(-span, span), k)); ys = sorted(random.sample(range(-span, span), k))
            Mx = mp.matrix([[P(x - y) for y in ys] for x in xs])
            d = mp.det(Mx)
            scale = mp.fprod([P(xs[i] - ys[i]) for i in range(k)]) or mp.mpf(1)
            r = d/abs(scale) if scale != 0 else d
            if d < -mp.mpf(10)**(-30)*abs(scale):
                nneg += 1
                if worst is None or r < worst[0]: worst = (r, [float(x*H) for x in xs], [float(y*H) for y in ys])
        pfres[f'PF{k}'] = dict(tries=ntry, negative_minors=nneg,
                               worst=None if worst is None else dict(rel=mp.nstr(worst[0], 5), x=worst[1], y=worst[2]))
    rec['PF_minors'] = pfres
    out[name] = rec
    print(name, json.dumps({k: rec[k] for k in ['phi0', 'min_phi', 'positive_on_grid', 'negative_u_interval', 'log_concave_where_positive']}),
          {k: v['negative_minors'] for k, v in pfres.items()}, [c['relerr'] for c in chk], flush=True)

# --- BPY identification: P(Y <= y) from the law vs Kuiper (Brownian-bridge range) CDF at r = y sqrt(pi/2) ---
uf, pf = full(*G['ZETA'])
dens = [2*pf[i]*mp.exp(-uf[i]/2) for i in range(len(uf))]
bpy = []
for y in [mp.mpf('0.6'), mp.mpf('0.8'), mp.mpf(1), mp.mpf('1.3'), mp.mpf(2)]:
    ly = mp.log(y)
    # trapezoid of density on lattice up to ly, plus partial cell via quad of model phi
    cdf = mp.quad(lambda u: 2*Model('ZETA').phi(u)*mp.exp(-u/2), [-3, -1, 0, ly]) if ly > -3 else None
    r = y*mp.sqrt(mp.pi/2)
    kuiper = 1 - 2*mp.nsum(lambda n: (4*n*n*r*r - 1)*mp.exp(-2*n*n*r*r), [1, mp.inf])
    bpy.append(dict(y=float(y), law_cdf=mp.nstr(cdf, 15), kuiper_cdf=mp.nstr(kuiper, 15), diff=float(abs(cdf - kuiper))))
out['BPY_kuiper_check'] = bpy
print('BPY', bpy)
mp.mp.dps = 40
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'r17_laws.json'), 'w'), indent=1)
