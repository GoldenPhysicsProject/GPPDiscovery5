"""C.14.K rigorous refutation at ell = 0.1155..., T = 14, even sector (Arb ball arithmetic, python-flint).
Profile pair (xi_b, xi_s): the (M, n_s) extremal pair (argv; runs: M=28,n_s=4 -> ell = a/6 = 0.11552; M=28,n_s=2 -> ell = a/13 = 0.05332) of c14K_gamma_test.gamma, rounded to 6 decimals,
symmetrised (even), xi_b linear on the strips, xi_s supported in the open strips. Node j at x_j = -L + j h, h = 2L/M.
P(u, v) = sum_{j,k} u_j v_k P(hat_j, hat_k), hat Toeplitz entries:
  P(hat_0, hat_k) = sum_{p,q in {-1,0,1}} w_p w_q K((k + p - q) h),  w = (1, -2, 1)/h,
  K(d) = sum_{n>=0} e^{-c_n|d|}(c_n^{-2} - c_n^{-4}/4) - |d| zeta(3, 1/4)/32,  c_n = 2n + 1/2
  (closed form of K = Phi1/2 + Psi of c14K_gamma_test.py; tail of the n-sum bounded by e^{-c_N d}/(c_N^2(1 - e^{-2d}))).
P_{<T}(hat_0, hat_k) = (1/pi) int_0^T (1/4 + t^2)(Omega(t) - Omega(0)) h^2 sinc^4(th/2) cos(k h t) dt  (acb.integral),
Omega(t) = (psi(1/4 + it/2) + psi(1/4 - it/2))/2 - log pi.  P_{>T} = P - P_{<T}.
Certifies |2 P_{>T}(xi_b, xi_s)| / [P_{>T}(xi_b) + P_{>T}(xi_s)] > gamma*(14) = Omega(14)/(Omega(14) + |Omega(0)|).
Daniel Toupin / Golden Physics Project, with AI assistance from Grok. RH is OPEN."""
import numpy as np
from flint import arb, acb, ctx
from c14K_gamma_test import gamma as gfloat, a as af
ctx.prec = 128
import sys
M, ns, T = int(sys.argv[1]), int(sys.argv[2]), 14
r = gfloat(M, ns, float(T), 'even', retvec=True)
xb = np.round(r['xi_b']/np.max(np.abs(r['xi_b'])), 6); xs = np.round(r['xi_s']/np.max(np.abs(r['xi_b'])), 6)
n = M - 1
xb = np.round((xb + xb[::-1])/2, 6); xs = np.round((xs + xs[::-1])/2, 6)
# enforce structure: xi_s zero off strips (nodes 1..ns-1 and M-ns+1..M-1), xi_b linear on strips
strip = [j for j in range(1, n+1) if j < ns or j > M - ns]
for j in range(1, n+1):
    if j not in strip: xs[j-1] = 0.0
alpha = xb[ns-1]
for j in strip:
    xb[j-1] = alpha*j/ns if j < ns else alpha*(M-j)/ns
ln2 = arb(2).log(); L = (ln2 + arb(ns)*ln2/arb(M-ns))/2        # L = a M /(2(M - ns)) exactly
h = 2*L/M; ell = 2*L - ln2
pi = arb.pi(); z3 = acb(3).zeta(acb(0.25)).real
def K(d):
    d = abs(d)
    if d == 0:
        return acb(2).zeta(acb(0.25)).real/4 - acb(4).zeta(acb(0.25)).real/64
    s = arb(0); Nn = 4000
    for k in range(Nn):
        c = arb(2*k) + arb(0.5); s += (-c*d).exp()*(1/c**2 - 1/(4*c**4))
    cN = arb(2*Nn) + arb(0.5); tail = (-cN*d).exp()/(cN**2*(1 - (-2*d).exp()))
    s += arb(0, 1)*tail   # enclose the n >= Nn tail
    return s - d*z3/32
Kv = [K(arb(k)*h) for k in range(M+3)]
w = [arb(1)/h, arb(-2)/h, arb(1)/h]
colP = []
for k in range(M+1):
    s = arb(0)
    for p in range(3):
        for q in range(3):
            s += w[p]*w[q]*Kv[abs(k + p - q)]
    colP.append(s)
Om0 = acb(0.25).digamma().real - pi.log()
def lowcol(k):
    hk = h; dk = arb(k)*h
    def g(t, analytic):
        Om = (acb(0.25) + acb(0, 0.5)*t).digamma()/2 + (acb(0.25) - acb(0, 0.5)*t).digamma()/2 - acb(pi.log())
        sn = (t*acb(hk)/2).sinc()
        return (acb(0.25) + t*t)*(Om - acb(Om0))*acb(hk)**2*sn**4*(t*acb(dk)).cos()/acb(pi)
    return acb.integral(g, 0, T).real
colL = [lowcol(k) for k in range(M+1)]
G = [[None]*n for _ in range(n)]
for i in range(n):
    for j in range(n):
        G[i][j] = colP[abs(i-j)] - colL[abs(i-j)]
def form(u, v):
    s = arb(0)
    for i in range(n):
        if u[i] == 0: continue
        for j in range(n):
            if v[j] == 0: continue
            s += arb(repr(float(u[i])))*arb(repr(float(v[j])))*G[i][j]
    return s
Pbb, Pss, Pbs = form(xb, xb), form(xs, xs), form(xb, xs)
OmT = (acb(0.25) + acb(0, 7)).digamma().real - pi.log()
gstar = OmT/(OmT + abs(Om0))
ratio = abs(2*Pbs)/(Pbb + Pss)
print("ell =", ell, " L =", L, " h =", h)
print("xi_b nodes j=1..M-1:", list(map(float, xb)))
print("xi_s nodes j=1..M-1:", list(map(float, xs)))
print("P_{>T}(xi_b,xi_b) =", Pbb); print("P_{>T}(xi_s,xi_s) =", Pss); print("P_{>T}(xi_b,xi_s) =", Pbs)
print("ratio |2P(b,s)|/(P(b)+P(s)) =", ratio); print("gamma*(14) =", gstar)
print("CERTIFIED ratio > gamma*:", bool(ratio > gstar), "; P_{>T}(xi_b) > 0:", bool(Pbb > 0), "; P_{>T}(xi_s) > 0:", bool(Pss > 0))
