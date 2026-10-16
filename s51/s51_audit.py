# Grok S51 H51a: independent audit of Codex (B) finite xi Hankel coercivity H3,H4,H5.
# Labels: CERTIFIED = Arb ball arithmetic / exact Fractions; NUMERICAL = high-precision float.
import sys, math
from fractions import Fraction as Q
import mpmath as mp
import flint
from flint import arb, acb, arb_mat, ctx
ctx.prec = 1400
K = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
G = [arb(l.split()[0], '1e-120') for l in open('zeros_arb_%d.txt' % K)]
H = 3 * 10**12
print('== (iii-a) first-five brackets vs Arb-certified zeros (acb.zeta_zeros, on line) ==')
br = [("14.13","14.14"),("21.02","21.03"),("25.01","25.02"),("30.42","30.43"),("32.93","32.94")]
for j,(a,b) in enumerate(br):
    ok = (G[j] > arb(a)) and (G[j] < arb(b))
    # independent Hardy Z sign change at bracket endpoints
    def Z(t):
        t = arb(t); th = acb(0.25, t/2).lgamma().imag - t/2*arb.pi().log()
        return (acb(0, th).exp() * acb(0.5, t).zeta()).real
    za, zb = Z(a), Z(b)
    print(j+1, G[j].mid().str(25), 'in', (a,b), 'CERTIFIED' if ok else 'FAIL', '| Z(a)=%s Z(b)=%s signchange=%s' % (za.str(5), zb.str(5), bool((za*zb) < 0)),
          '| lo-gap %.4g hi-gap %.4g' % (float(G[j]-arb(a)), float(arb(b)-G[j])))
print('== (iii-b) Codex Lagrange/Frobenius bound recomputed independently (exact Fractions) ==')
glo = [Q(a) for a,b in br]; ghi = [Q(b) for a,b in br]
zlo = [1/x**2 for x in ghi]; zhi = [1/x**2 for x in glo]
def frob_upper(n):
    tot = Q(0)
    for j in range(n):
        den = zlo[j]
        for k in range(n):
            if k == j: continue
            den *= (zlo[k] - zhi[j]) if k < j else (zlo[j] - zhi[k])
        # coefficients of prod_{k!=j}(t - z_k): |c_m| = e_m(z_k) (all z>0), majorize by z_hi
        e = [Q(1)]
        for k in range(n):
            if k == j: continue
            e = [ (e[i] if i < len(e) else 0) + (zhi[k]*e[i-1] if i >= 1 else 0) for i in range(len(e)+1)]
        tot += sum(c*c for c in e) / den**2
    return tot
def frob_true(n):  # ||V^-1||_F^2 at the true zeros, Arb
    z = [1/(g*g) for g in G[:n]]
    V = arb_mat(n, n, [z[j]**(i+1) for i in range(n) for j in range(n)])
    Vi = V.inv(); return sum(Vi[i,j]**2 for i in range(n) for j in range(n))
def lmin_mp(M):
    E = mp.eigsy(M, eigvals_only=True); return min(E)
mp.mp.dps = 120
for n in (3,4,5):
    fu = frob_upper(n); ft = frob_true(n)
    zt = [1/mp.mpf(G[j].mid().str(120, radius=False))**2 for j in range(n)]
    Gm = mp.matrix(n, n)
    for i in range(n):
        for k in range(n): Gm[i,k] = sum(z**(i+k+2) for z in zt)
    print('N=%d  Codex-style upper ||V^-1||_F^2 = %.6e  (true %.6e, ratio %.3f) ; true lambda_min(G_N, first N zeros) = %s ; 1/||V^-1||_F^2 = %.4e'
          % (n, float(fu), float(ft), float(fu)/float(ft), mp.nstr(lmin_mp(Gm), 8), 1/float(ft)))
print('Codex asserts: inv4<5e23, inv5<2e31 ->', frob_upper(4) < 5*10**23, frob_upper(5) < 2*10**31)
print('== (ii) tail sum_{gamma>H} gamma^-4 ==')
lH = math.log(H)
print('log H = %.6f' % lH)
print('2/H^2 = %.4e ; (4logH/3+4/9)/H^3 = %.4e ; 41/H^3 = %.4e' % (2/H**2, (4*lH/3+4/9)/H**3, 41/H**3))
# Rosser-type count check N(t) <= t log t for t>=H (worst case is t=H; ratio decreasing)
ross = lambda t: t/(2*math.pi)*math.log(t/(2*math.pi*math.e)) + 7/8 + 0.137*math.log(t) + 0.443*math.log(math.log(t)) + 1.588
print('Rosser upper N(H) = %.6e ; H log H = %.6e ; H^2 = %.3e' % (ross(H), H*lH, H**2.0))
# near-true value of the tail with Riemann-von Mangoldt main term (NUMERICAL, any zero placement)
tv = mp.quad(lambda t: mp.log(t/(2*mp.pi))/(2*mp.pi)*t**-4, [H, mp.inf])
print('main-term tail int_H^inf (1/2pi)log(t/2pi) t^-4 dt = %s (= ~%.2f/H^3)' % (mp.nstr(tv, 6), float(tv*H**3)))
print('needed: H5 margin 1/(2e31)-4e-32 = 1e-32 -> need 5*S <= 1e-32, S<=2e-33; N(t)<=t^2 gives 2.2e-25 (FAILS for H5), t^1.5 gives', '%.2e' % (1.6/H**2.5))
print('== (iv) actual h_n: Phi-moment (zero-free) vs zero sums ==')
mp.mp.dps = 80
def Phi(t):
    s = 0
    for n in range(1, 30):
        e = mp.pi*n*n*mp.exp(4*t)
        if e > 400: break
        s += (2*mp.pi**2*n**4*mp.exp(9*t) - 3*mp.pi*n*n*mp.exp(5*t))*mp.exp(-e)
    return s
NM = 22
mom = [mp.quad(lambda t: Phi(t)*t**(2*k), [0, 0.25, 0.5, 1, 1.5, 2.5]) for k in range(NM)]
a = [mom[k]*4**k/mp.factorial(2*k)/mom[0] for k in range(NM)]   # xi(1/2+w) ∝ ∫Phi cosh(2wt)
q = []
for n in range(NM-1):
    q.append((n+1)*a[n+1] - sum(a[k]*q[n-k] for k in range(1, n+1)))
hPhi = [(-1)**n*q[n] for n in range(NM-1)]
zs = [mp.mpf(g.mid().str(80, radius=False)) for g in G]
T = zs[-1]
def tailint(p):  # smooth RvM density tail beyond last zero (NUMERICAL, assumes on-line density)
    return mp.quad(lambda t: mp.log(t/(2*mp.pi))/(2*mp.pi)*t**(-2*p), [T, mp.inf])
hZ = [sum(z**(-2*(n+1)) for z in zs) + tailint(n+1) for n in range(NM-1)]
for n in range(8):
    print('h_%d  Phi: %s   zeros+tail: %s   rel diff %.2e' % (n, mp.nstr(hPhi[n], 20), mp.nstr(hZ[n], 20), float(abs(hPhi[n]/hZ[n]-1))))
def Hmat(h, N):
    M = mp.matrix(N, N)
    for i in range(N):
        for j in range(N): M[i,j] = h[i+j+1]
    return M
print('det H3 (Phi) = %s ; Codex 3.11996028999615e-31' % mp.nstr(mp.det(Hmat(hPhi,3)), 12))
print('h1*h5-h3^2 (Phi) = %s ; Codex 1.50575915582118e-19' % mp.nstr(hPhi[1]*hPhi[5]-hPhi[3]**2, 12))
mp.mp.dps = 80
for N in range(3, 11):
    l1 = lmin_mp(Hmat(hPhi, N)); l2 = lmin_mp(Hmat(hZ, N))
    print('N=%2d lambda_min(H_N): Phi-moment %s | zero-sum %s' % (N, mp.nstr(l1, 10), mp.nstr(l2, 10)))
print('== CERTIFIED extension: G_K (K Arb-certified on-line zeros) > c I by Sylvester on G_K - cI ==')
# sharper tail: |R_ij| <= S_{i+j+2} <= S4*H^{-2(i+j)}, so ||R||_op <= ||R||_F <= S4*(1+1e-40);
# S4 = sum_{gamma>H} gamma^-4 <= 4 int_H^inf N(t) t^-5 dt with N(t) <= (t/2pi) log t (t>=H; Rosser),
# = (1/2pi)(4 log H/3 + 4/9)/H^3.
S4 = (4*arb(H).log()/3 + arb(4)/9)/(2*arb.pi())/arb(H)**3
print('S4 bound via N(t)<=(t/2pi)log t:', S4.str(6), '; via N(t)<=t log t: 41/H^3 =', (arb(41)/arb(H)**3).str(6))
x = [1/(g*g) for g in G]
pw = [sum((xx**p for xx in x), arb(0)) for p in range(2, 2*12+1)]
for N in range(3, 11):
    M = [[pw[i+j] for j in range(N)] for i in range(N)]
    # numerical lambda_min of G_K
    Mm = mp.matrix(N, N)
    for i in range(N):
        for j in range(N): Mm[i,j] = mp.mpf(M[i][j].mid().str(100, radius=False))
    lam = lmin_mp(Mm)
    c = arb(mp.nstr(lam*mp.mpf('0.9'), 30))
    ok = all(arb_mat(k, k, [M[i][j] - (c if i == j else 0) for i in range(k) for j in range(k)]).det() > 0 for k in range(1, N+1))
    tail = S4*(1+arb('1e-40'))
    print('N=%2d lambda_min(G_K)~%s ; certified G_K > %s I: %s ; minus tail %.3e -> H_N > %s I : %s'
          % (N, mp.nstr(lam, 8), c.mid().str(6), ok, float(tail), (c - tail).mid().str(6), bool(ok and (c - tail) > 0)))
