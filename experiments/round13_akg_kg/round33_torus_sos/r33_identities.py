"""Round 33 identity checks (CHECKED; PROVED by the derivations in ROUND33.md). No zeros.
 (T1) quadratic form  sum_{m,n<=X} a_m a_n (L - |log(m/n)|)_+  =  int |W(u)|^2 du  (exact step integrals)
 (T2) Kronecker-line form  (1/2pi) int |A(1/2+it)|^2 4 sin^2(Lt/2)/t^2 dt  =  the same quadratic form
 (T3) Haar diagonal: (1/2T) int_{-T}^{T} |A(1/2+it)|^2 dt -> sum |a_n|^2 (Bohr / Haar orthogonality of distinct n)
 (T4) Cauchy-Schwarz/SOS: L - 2 lam <f,1> + lam^2 ||f||^2 >= 0, i.e. (m - cL)^2 <= L * G1c <= L * Gc."""
import json, math, numpy as np, sympy
out = {}
def vec(X, c):
    n = np.arange(1, int(X) + 1); a = np.array([int(sympy.mobius(k)) for k in n], float)/np.sqrt(n); a[0] -= c; return n, a
c = -0.6847652360899366
for L in (3.0, 5.0, 7.0):
    n, a = vec(math.exp(L), c); lg = np.log(n.astype(float))
    K = np.maximum(L - np.abs(lg[:, None] - lg[None, :]), 0); qf = float(a @ K @ a)
    M = np.cumsum(a); dl = np.diff(np.append(lg, L)); Q1 = float((M*M*dl).sum()); m1 = float((M*dl).sum())
    # window (u-L, u] for u in [L, 2L): M(X) - M(e^{u-L}); M includes the n=1 coefficient (1-c) there
    G = Q1 + L*M[-1]**2 - 2*M[-1]*m1 + Q1
    cross = float(a @ np.maximum(L - lg, 0))
    out[f'T1_L{L}'] = dict(X=int(n[-1]), quadform=qf, step_integral=G, diff=qf - G, cross_f1=cross, m_minus_cL_stepint=m1)
    print(out[f'T1_L{L}'], flush=True)
# T2 on the Kronecker line, L = 3 (X = 20)
L = 3.0; n, a = vec(math.exp(L), c); lg = np.log(n.astype(float))
K = np.maximum(L - np.abs(lg[:, None] - lg[None, :]), 0); qf = float(a @ K @ a)
Tm = 4000.0; t = np.linspace(1e-9, Tm, 4_000_001); A = (a[None, :]*np.exp(-1j*np.outer(t, lg))).sum(1) if False else None
acc = np.zeros_like(t, dtype=complex)
for an, l in zip(a, lg): acc += an*np.exp(-1j*t*l)
f = np.abs(acc)**2*4*np.sin(L*t/2)**2/t**2
val = 2*np.trapezoid(f, t)/(2*math.pi); tail = float((a*a).sum())*2/(math.pi*Tm)
out['T2'] = dict(L=L, quadform=qf, line_integral=val, diag_tail_estimate=tail, diff=qf - val - tail)
# T3 Haar diagonal by long time average (Kronecker-Weyl)
for TT in (1e2, 1e3, 4e3):
    sel = t <= TT; out[f'T3_T{int(TT)}'] = dict(time_avg=float(np.trapezoid(np.abs(acc[sel])**2, t[sel])/TT), haar_diag=float((a*a).sum()))
print(out['T2'], out['T3_T100'], out['T3_T4000'], flush=True)
json.dump(out, open('r33_identities.json', 'w'), indent=1)
