"""Clue probe A: Hermite-Kakeya-Obreschkoff (HKO) interlacing for the two fakes that are pencils.
EPS: Xi_Q1 = A + B,  A = completed zeta_K (K=Q(sqrt-5)) = zeta*L(chi_-20), B = completed L(chi_-4)L(chi_5).
DH : Xi_DH ~ A + r B, A = Z_chi, B = Z_chibar (chi mod 5, chi(2)=i), r computed (expected +-1).
HKO: A + lam*B is real-rooted for ALL real lam iff zeros of A and B are real and interlace.
We list the heights where the on-line zeros of A and B fail to interlace and compare with the fakes'
known off-line zeros (test data only). Zeros never enter any construction."""
import mpmath as mp, json, numpy as np
mp.mp.dps = 30
def chi_from(q, table):  # table: dict n mod q -> value
    return [table.get(n % q, 0) for n in range(q)]
chi_m4 = chi_from(4, {1: 1, 3: -1})
chi_5 = chi_from(5, {1: 1, 4: 1, 2: -1, 3: -1})
leg = {}
def chi_m20_f(n):
    return chi_m4[n % 4] * chi_from(5, {1:1,4:1,2:-1,3:-1})[n % 5] * 1  # chi_-20 = chi_-4 * chi_5 ? (check below)
# chi_-20 is the Kronecker symbol (-20/.) ; = chi_-4 * chi_5 (both primitive, conductors 4,5) -> conductor 20, odd
chi_m20 = [chi_m4[n % 4] * chi_5[n % 5] for n in range(20)]
chi_c = chi_from(5, {1: 1, 2: 1j, 4: -1, 3: -1j})      # chi(2)=i, odd
chi_cb = [complex(v).conjugate() for v in chi_c]
def L(s, chi): return mp.dirichlet(s, chi)
def eps_A(t):
    s = mp.mpf(0.5) + 1j*t
    g = (mp.sqrt(20)/(2*mp.pi))**s * mp.gamma(s)
    ph = g/abs(g)
    return ph*mp.zeta(s)*L(s, chi_m20), ph*L(s, chi_m4)*L(s, chi_5)
def dh_AB(t):
    s = mp.mpf(0.5) + 1j*t
    g = (5/mp.pi)**(s/2) * mp.gamma((s+1)/2)
    ph = g/abs(g)
    return ph*L(s, chi_c), ph*L(s, chi_cb)
# root-number rotation for DH: Lam(1/2+it,chi) = sqrt(eps)*real
a0, b0 = dh_AB(mp.mpf(3.3))
rotA = mp.sqrt(a0/mp.conj(a0)); rotB = mp.sqrt(b0/mp.conj(b0))
def zeros_by_sign(f, T0, T1, h):
    ts = np.arange(T0, T1, h); vals = [float(f(t)) for t in ts]
    z = []
    for i in range(len(ts)-1):
        if vals[i] == 0 or vals[i]*vals[i+1] < 0:
            z.append(float(mp.findroot(lambda x: f(x), (ts[i], ts[i+1]), solver='bisect' if False else 'anderson')))
    return z, ts, vals
def interlace_defects(za, zb):
    ev = sorted([(x, 'A') for x in za] + [(x, 'B') for x in zb])
    bad = []
    for i in range(len(ev)-1):
        if ev[i][1] == ev[i+1][1]:
            bad.append((round(ev[i][0], 4), round(ev[i+1][0], 4), ev[i][1]))
    return bad
out = {}
T1 = 60.0
# EPS
fa = lambda t: mp.re(eps_A(t)[0]); fb = lambda t: mp.re(eps_A(t)[1])
imA = max(abs(float(mp.im(eps_A(t)[0]))) for t in (5.0, 20.0, 40.0)); imB = max(abs(float(mp.im(eps_A(t)[1]))) for t in (5.0, 20.0, 40.0))
za, _, _ = zeros_by_sign(fa, 0.05, T1, 0.02); zb, _, _ = zeros_by_sign(fb, 0.05, T1, 0.02)
# check pencil value at known off-line EPS zero (test data) : A+B real form vanishes there? (analytic continuation check)
out['EPS'] = dict(imag_residual_A=imA, imag_residual_B=imB, zeros_A=za, zeros_B=zb, defects=interlace_defects(za, zb))
# DH
fa2 = lambda t: mp.re(dh_AB(t)[0]/rotA); fb2 = lambda t: mp.re(dh_AB(t)[1]/rotB)
imA2 = max(abs(float(mp.im(dh_AB(t)[0]/rotA))) for t in (5.0, 20.0, 40.0)); imB2 = max(abs(float(mp.im(dh_AB(t)[1]/rotB))) for t in (5.0, 20.0, 40.0))
# DH combination per 14F: DH = alpha L(chi) + beta L(chibar), alpha,beta=(1 -+ i kappa)/2
kap = (mp.sqrt(10-2*mp.sqrt(5))-2)/(mp.sqrt(5)-1)
al, be = (1-1j*kap)/2, (1+1j*kap)/2
r = (be*rotB)/(al*rotA)
T2 = 130.0
za2, _, _ = zeros_by_sign(fa2, 0.05, T2, 0.02); zb2, _, _ = zeros_by_sign(fb2, 0.05, T2, 0.02)
out['DH'] = dict(imag_residual_A=imA2, imag_residual_B=imB2, pencil_ratio_r=[float(mp.re(r)), float(mp.im(r))],
                 zeros_A=za2, zeros_B=zb2, defects=interlace_defects(za2, zb2))
# DH check: alpha L(chi)+beta L(chibar) at the 14F off-line zero
rho = mp.mpc('0.808517182457', '85.6993484854')
out['DH']['abs_DH_at_offline_zero'] = float(abs(al*L(rho, chi_c)+be*L(rho, chi_cb)))
rhoE = mp.mpc('0.932969697485', '15.6682495313')
out['EPS']['abs_EPS_at_offline_zero'] = float(abs(mp.zeta(rhoE)*L(rhoE, chi_m20) + L(rhoE, chi_m4)*L(rhoE, chi_5)))
json.dump(out, open('hko_interlace.json', 'w'), indent=1)
for k in out:
    print(k, 'nA', len(out[k]['zeros_A']), 'nB', len(out[k]['zeros_B']), 'defects', out[k]['defects'], {x: out[k][x] for x in out[k] if 'residual' in x or 'abs' in x or 'ratio' in x})
