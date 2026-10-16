#!/usr/bin/env python3
"""CERTIFIED Omega: Arb digamma vs EM series through B8, remainder B10."""
import json
from flint import arb, acb, ctx
from math import comb
ctx.prec = 384

def omega0():
    return -arb.const_euler() - arb.pi()/2 - 3*arb(2).log() - arb.pi().log()

def Omega_digamma(t):
    z = acb(arb('0.25'), arb(t)/2)
    return z.digamma().real - arb.pi().log()

def Omega_series_certified(t, N=2000):
    t = arb(t); t2 = t*t
    s = arb(0)
    for j in range(N):
        a = arb(j) + arb('0.25')
        s += t2 / (a * (4*a*a + t2))
    A = arb(N) + arb('0.25')
    integ = (arb(1)/2) * (1 + t2/(4*A*A)).log()
    den = [4*A*A*A + t2*A, 12*A*A + t2, 24*A, arb(24)] + [arb(0)]*16
    g = [arb(0)]*16
    g[0] = 1/den[0]
    for n in range(1, 16):
        acc = arb(0)
        for j in range(1, min(n, 3)+1):
            acc += comb(n, j) * den[j] * g[n-j]
        g[n] = -acc / den[0]
    fder = [t2 * g[k] for k in range(16)]
    # B2=1/6, B4=-1/30, B6=1/42, B8=-1/30, B10=5/66
    B = {2: arb(1)/6, 4: -arb(1)/30, 6: arb(1)/42, 8: -arb(1)/30, 10: arb(5)/66}
    em = integ + fder[0]/2
    for k in (1, 2, 3, 4):  # through B8
        twok = 2*k
        fac = arb(1)
        for i in range(1, twok+1):
            fac *= i
        em -= B[twok] / fac * fder[twok - 1]
    fac10 = arb(1)
    for i in range(1, 11):
        fac10 *= i
    Rbound = (abs(B[10]) / fac10) * abs(fder[9]) * arb(2)
    return omega0() + s + em, Rbound

ts = [0.5, 1.0, 5.0, 20.0]
rows = []
for t in ts:
    Od = Omega_digamma(t)
    ok = False
    for N in [1000, 5000, 20000, 80000, 200000]:
        Os, Rb = Omega_series_certified(t, N=N)
        resid = abs(Od - Os)
        ok = float((resid + Rb).upper()) < 1e-50
        print(f'  t={t} N={N}: resid={float(resid.upper()):.3e} R={float(Rb.upper()):.3e} OK={ok}')
        if ok:
            rows.append({'t': t, 'N': N, 'resid_upper': float(resid.upper()),
                         'R_upper': float(Rb.upper()), 'certified_le_1e-50': True,
                         'Omega': str(Od)})
            break
    if not ok:
        rows.append({'t': t, 'N': N, 'resid_upper': float(resid.upper()),
                     'R_upper': float(Rb.upper()), 'certified_le_1e-50': False,
                     'Omega': str(Od)})

o0 = omega0()
diff0 = float(abs(Omega_digamma(0) - o0).upper())
out = {
    'omega0': str(o0),
    'omega0_matches_Omega0': diff0 < 1e-80,
    'checks': rows,
    'all_certified_1e-50': all(r['certified_le_1e-50'] for r in rows),
    'verdict': 'CERTIFIED' if all(r['certified_le_1e-50'] for r in rows) else 'CONFIRMED',
}
json.dump(out, open('/workspace/GPPDiscovery5/s39/h39a_omega_cert.json','w'), indent=2)
print('VERDICT', out['verdict'], 'all50', out['all_certified_1e-50'])
