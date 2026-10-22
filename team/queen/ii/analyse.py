"""Queen S1d (ii) analysis: tables, k-convergence, power-law extrapolation, eigenvector cosines vs step-(c) Legendre vectors.
NUMERICAL. Cosine basis/assembly: Goddess of Prime Number Theory / team/prime-number-theory, assembly from Grok S62."""
import json, glob, os
import mpmath as mp
mp.mp.dps = 60
from fractions import Fraction
Ls = ['1','11_10','6_5','13_10','7_5','3_2','8_5','17_10']
def lval(s): return float(Fraction(s.replace('_','/')))
res = {}
for sec in ('even','odd'):
    for L in Ls:
        f = f'r_{sec}_{L}.json' if L != '1' else f'r_{sec}_L1_K266_p576.json'
        if os.path.exists(f): res[(sec, L)] = json.load(open(f))
leg = {float(e['L']): e for e in json.load(open('/workspace/queen/GPPDiscovery5/team/queen/step_c/lambda2_Re_eigenvectors.json'))['vectors']}
legvals_o = {1.0: '8.1749e-21'}
def extrap(rows):
    r = rows[-3:]
    if len(r) < 3: return None, None
    (k1,l1),(k2,l2),(k3,l3) = [(mp.mpf(x['k']), mp.mpf(x['lam'])) for x in r]
    d1, d2 = l1-l2, l2-l3
    if d1 <= 0 or d2 <= 0: return None, None
    f = lambda p: (k1**-p - k2**-p)/(k2**-p - k3**-p) - d1/d2
    try: p = mp.findroot(f, 2)
    except Exception: return None, None
    A = d2/(k2**-p - k3**-p); return l3 - A*k3**-p, p
# function evaluation for cosines on [0,L]
def cosfun(vec, L, sec):
    v = [mp.mpf(x) for x in vec]
    if sec == 'even': return lambda u: sum(v[m]*mp.cos((2*m+1)*mp.pi*u/(2*L)) for m in range(len(v)))
    return lambda u: sum(v[m]*mp.sin((m+1)*mp.pi*u/L) for m in range(len(v)))
def legfun(coef, L):
    c = [mp.mpf(x) for x in coef]
    return lambda s: sum(c[j]*mp.sqrt((2*j+1)/L)*mp.legendre(j, 2*s/L-1) for j in range(len(c)))
def fcos(f, g, L, n=400):
    xs, ws = zip(*[(L*(x+1)/2, w*L/2) for x, w in zip(*gl(n))])
    fv = [f(x) for x in xs]; gv = [g(x) for x in xs]
    ip = sum(w*a*b for w,a,b in zip(ws,fv,gv)); nf = sum(w*a*a for w,a in zip(ws,fv)); ng = sum(w*b*b for w,b in zip(ws,gv))
    return abs(ip)/mp.sqrt(nf*ng)
_gl = {}
def gl(n):
    if n not in _gl:
        import numpy as np
        x, w = np.polynomial.legendre.leggauss(n); _gl[n] = ([mp.mpf(a) for a in x], [mp.mpf(b) for b in w])
    return _gl[n]
mp.mp.dps = 30
out_md = ["| L | sector | k=200 | 266 | 320 | 400 | k-extrap (p) | last-step rel. change | residual(k_max) | Weyl(k_max) |", "|---|---|---|---|---|---|---|---|---|---|"]
vecs = []; summ = []
for L in Ls:
    for sec in ('even','odd'):
        rows = res.get((sec, L))
        if not rows: continue
        Lf = lval(L); d = {x['k']: x for x in rows}
        cells = [mp.nstr(mp.mpf(d[k]['lam']), 5) if k in d else '—' for k in (200,266,320,400)]
        ex, p = extrap(rows)
        last = rows[-1]; prev = rows[-2]
        rel = float((mp.mpf(prev['lam'])-mp.mpf(last['lam']))/mp.mpf(last['lam']))
        out_md.append(f"| {Lf:g} | {sec} | " + " | ".join(cells) + f" | {mp.nstr(ex,4) if ex else '—'} ({mp.nstr(p,3) if p else '—'}) | {rel:.2%} | {last['residual']:.1e} | {last['weyl']:.1e} |")
        # cosine nested
        f_last = cosfun(last['vec'], Lf, sec); f_prev = cosfun(prev['vec'], Lf, sec)
        c_nest = fcos(f_last, f_prev, Lf)
        c_leg = None
        if sec == 'even' and Lf in leg and Lf <= 1.4:
            c_leg = fcos(f_last, legfun((lambda c: eval(c) if isinstance(c,str) else c)(leg[Lf]['coefficients']), Lf), Lf)
        summ.append(dict(L=Lf, sector=sec, k=last['k'], lam=last['lam'], lam_prev_k=prev['k'], lam_prev=prev['lam'], extrap=(mp.nstr(ex,6) if ex else None), extrap_p=(mp.nstr(p,4) if p else None),
                         rel_change_last_step=rel, residual=last['residual'], weyl=last['weyl'], one_minus_abs_cos_nested=float(1-c_nest),
                         one_minus_abs_cos_vs_legendre_N64=(float(1-c_leg) if c_leg is not None else None),
                         legendre_N64_lam=(leg[Lf]['eigenvalue_lam2_Re'] if sec=='even' and Lf in leg else legvals_o.get(Lf))))
        vecs.append(dict(L=str(Lf), K=str(last['k']), sector=sec, prec_bits=str(last['prec']), eigenvalue=last['lam'],
                         eigenvalue_name=('lambda2_Re' if sec=='even' else 'lambda1_Ro'), converged_flag_5pct=str(abs(rel) < 0.05),
                         abs_cos_with_K_prev=mp.nstr(c_nest, 12), K_prev=str(prev['k']), coefficients=str(last['vec'])))
        print(out_md[-1], summ[-1]['one_minus_abs_cos_nested'], summ[-1]['one_minus_abs_cos_vs_legendre_N64'], flush=True)
D = '/workspace/queen/GPPDiscovery5/team/queen/ii/'
open(D+'table_ii.md','w').write("\n".join(out_md)+"\n")
json.dump(summ, open(D+'summary_ii.json','w'), indent=1)
json.dump(dict(description="Queen S1d (ii): eigenvector of R_e (lambda_2) / R_o (lambda_1) for zeta, half-width L, in Goddess of Prime Number Theory's cosine basis (team/prime-number-theory precision/cosritz_dump.py, assembly from Grok S62): even c_m(u)=cos((2m+1)pi u/(2L)), odd s_m(u)=sin(m pi u/L), m=1..K, u in [-L,L]; Gram = L*I, so divide |coef|^2 by... unit coefficient vector has L^2[-L,L] norm sqrt(L). Coefficients: unit Euclidean norm, sign arbitrary, 40 digits. R = Q - pole (R_e = W_e - 2cp cp^T, R_o = W_o + 2 sp sp^T). Eigenvalue = Rayleigh quotient / L in Arb from midpoint inverse iteration (HIGH-PRECISION NUMERICAL, not certified: index of eigenvalue not certified). Ritz values are upper bounds in K. RH OPEN.", vectors=vecs), open(D+'eigvec_ii.json','w'))
