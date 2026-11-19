"""Task 2: Stieltjes-Pick (Loewner) matrices of m_M at the Shadow-Euler points u_{k,N} = (kN/(k+N) - 1/2)^2, N = 2..Nmax.
P0[i,j] = (m(u_i) - m(u_j))/(u_j - u_i)   (diag -m'(u_i))       -> for Stieltjes m = int dmu(t)/(t+u): Gram of 1/(t+u)
P1[i,j] = (u_i m(u_i) - u_j m(u_j))/(u_i - u_j) (diag (u m)'(u_i)) -> Gram of sqrt(t)/(t+u)
Both PSD for all sizes <=> m Stieltjes (on a set with an interior accumulation point) <=> RH for the model.
We report LDL pivots d_n = det_n/det_{n-1} (leading minors in the order N = 2, 3, ...), the first n with d_n <= 0,
and stability between two precisions. m is evaluated from the Taylor coefficients of F (no zeros)."""
import sys, json, time
import mpmath as mp
from r20_core import *
name = sys.argv[1]; K = int(sys.argv[2]); D1, D2 = int(sys.argv[3]), int(sys.argv[4]); NMAX = int(sys.argv[5])
K2 = int(sys.argv[6]) if len(sys.argv) > 6 else K   # second run uses more Taylor terms: checks truncation too
KOF = {D1: K, D2: K2}
KS = [1, 2, 3, 4]
def ldl_pivots(P):
    n = P.rows; A = P.copy(); piv = []
    for k in range(n):
        d = A[k, k]; piv.append(d)
        if d == 0: break
        for i in range(k + 1, n):
            l = A[i, k]/d
            for j in range(k + 1, i + 1):
                A[i, j] -= l*A[j, k]
                A[j, i] = A[i, j]
    return piv
res = {}
for D in (D1, D2):
    f = F_coeffs(name, KOF[D], D); mp.mp.dps = D
    fp = [(k + 1)*f[k + 1] for k in range(len(f) - 1)]; fpp = [(k + 1)*fp[k + 1] for k in range(len(fp) - 1)]
    ev = lambda c, u: mp.fsum(c[k]*u**k for k in range(len(c)))
    res[D] = {}
    for k in KS:
        us = [(mp.mpf(k*N)/(k + N) - mp.mpf(1)/2)**2 for N in range(2, NMAX + 2)]
        Fv = [ev(f, u) for u in us]; F1 = [ev(fp, u) for u in us]; F2 = [ev(fpp, u) for u in us]
        m = [a/b for a, b in zip(F1, Fv)]; mprime = [(c*b - a*a)/(b*b) for a, b, c in zip(F1, Fv, F2)]
        n = len(us); P0 = mp.matrix(n, n); P1 = mp.matrix(n, n)
        for i in range(n):
            for j in range(n):
                if i == j:
                    P0[i, i] = -mprime[i]; P1[i, i] = m[i] + us[i]*mprime[i]
                else:
                    P0[i, j] = (m[i] - m[j])/(us[j] - us[i]); P1[i, j] = (us[i]*m[i] - us[j]*m[j])/(us[i] - us[j])
        res[D][k] = (ldl_pivots(P0), ldl_pivots(P1), us)
    print(name, D, 'done', flush=True)
out = {}
mp.mp.dps = D2
for k in KS:
    row = {}
    for idx, lab in ((0, 'P0'), (1, 'P1')):
        a, b = res[D1][k][idx], res[D2][k][idx]
        dig = [float(-mp.log10(abs(x - y)/abs(y))) if x != y else float(D1) for x, y in zip(a, b)]
        ok = 0
        while ok < len(dig) and dig[ok] >= 8: ok += 1
        neg = next((n + 1 for n, x in enumerate(b[:ok]) if x <= 0), None)
        row[lab] = dict(n_stable=ok, first_failing_size=neg, log10_abs_pivots=[round(float(mp.log10(abs(x))), 2) for x in b[:ok]],
                        pivot_signs=''.join('+' if x > 0 else '-' for x in b[:ok]))
    row['u_limit'] = float((k - mp.mpf(1)/2)**2)
    out[f'k={k}'] = row
    print(name, 'k', k, {l: (row[l]['n_stable'], row[l]['first_failing_size']) for l in ('P0', 'P1')}, flush=True)
out['settings'] = dict(K=K, K2=K2, D1=D1, D2=D2, NMAX=NMAX)
json.dump(out, open(f'r20_pick_{name}.json', 'w'), indent=1)
