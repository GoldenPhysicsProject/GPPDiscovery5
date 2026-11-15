"""C1 face-split discriminator (Heart, Turbo 6888439939 item 2). NUMERICAL.
archTerm by period-wise Gauss-Legendre (96 nodes per 2*pi, t <= 3200*pi; tail < 1e-11), cross-checked at d=0 by mp.quadosc (agree to 4e-12).
Kernels exactly as in-tree KreinHPsi.lean: archTerm, primeTermS(S_L), bigPsiI2; no zero data.
Configuration: Flint's clustered0.1L, L=0.5, n=6: x_k = (k-2.5)*0.01, k=0..5 (P-symmetric).
D = A - P_S - 2 Psi(i/2)^2 s s^T, s_i = sinh(x_i/2). Reversal J commutes with D, so ker l_c splits
into even (c = Jc, l_c(c)=0, s.c=0) and odd (c = -Jc, l_c(c)=0 automatically; pole term lives here)."""
import mpmath as mp, json, time
mp.mp.dps = 40
def triG(y):
    ay = abs(y)
    if ay <= 1: return mp.mpf(2)/3 - y*y + ay**3/2
    if ay <= 2: return (2-ay)**3/6
    return mp.mpf(0)
def primepow(n):
    for p in range(2, n+1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return p if m == 1 else None
L = mp.mpf('0.5')
Smax = mp.e**(2+L)
def primeS(d):
    d = mp.mpf(d); s = mp.mpf(0)
    for n in range(2, int(mp.e**(2+abs(d)))+3):
        p = primepow(n)
        if p is None or p >= Smax: continue
        ln = mp.log(n)
        s += mp.log(p)/mp.sqrt(n)*(triG(ln+d)+triG(-ln+d))
    return s
def arch(d):
    d = mp.mpf(d)
    f = lambda t: (mp.mpf(1) if t == 0 else (mp.sin(t/2)/(t/2))**2)**2*mp.cos(t*d)*(mp.re(mp.digamma(mp.mpc(mp.mpf('0.25'), t/2))) - mp.log(mp.pi))
    pts = [0] + [2*mp.pi*k for k in range(1, 41)] + [mp.inf]
    return mp.quad(f, pts)/mp.pi
bI2 = 16*mp.sinh(mp.mpf('0.25'))**2
n = 6
x = [(mp.mpf(k) - mp.mpf('2.5'))*mp.mpf('0.01') for k in range(n)]
def lag(i, j): return int(round(abs(float(x[i]-x[j]))*100))
lags = sorted(set(lag(i, j) for i in range(n) for j in range(n)))
t0 = time.time()
def arch_flint(d):
    d = mp.mpf(d)
    f = lambda t: (mp.mpf(1) if t == 0 else (mp.sin(t/2)/(t/2))**2)**2*mp.cos(t*d)*(mp.re(mp.digamma(mp.mpc(mp.mpf('0.25'), t/2))) - mp.log(mp.pi))
    return (1/mp.pi)*mp.quad(f, [0, mp.inf])
A_ = {k: arch_flint(mp.mpf(k)/100) for k in lags}
P_ = {k: primeS(mp.mpf(k)/100) for k in lags}
s = [mp.sinh(v/2) for v in x]; ell = [mp.cosh(v/2) for v in x]
A = mp.matrix(n, n); P = mp.matrix(n, n); D = mp.matrix(n, n)
for i in range(n):
    for j in range(n):
        A[i,j] = A_[lag(i,j)]; P[i,j] = P_[lag(i,j)]
        D[i,j] = A[i,j] - P[i,j] - 2*bI2**2*s[i]*s[j]
def orth(cols):
    Z = []
    for v in cols:
        w = mp.matrix(v)
        for q in Z: w -= q*(q.T*w)[0]
        Z.append(w/mp.sqrt((w.T*w)[0]))
    return Z
def proj(M, Z):
    m = len(Z); R = mp.matrix(m, m)
    for i in range(m):
        for j in range(m): R[i,j] = (Z[i].T*M*Z[j])[0]
    return sorted(mp.eigsy(R, eigvals_only=True))
def e(k):
    v = mp.matrix(n, 1); v[k] = 1; return v
lc = lambda v: sum(ell[i]*v[i] for i in range(n))
Zfull = orth([e(k) - e(0)*(ell[k]/ell[0]) for k in range(1, n)])
ev_full = proj(D, Zfull)
U = [e(k)+e(n-1-k) for k in range(3)]
Zeven = orth([U[k] - U[0]*(lc(U[k])/lc(U[0])) for k in (1, 2)])
Mnp = A - P
ev_even_AP = proj(Mnp, Zeven); ev_even_D = proj(D, Zeven)
Zodd = orth([e(k)-e(n-1-k) for k in range(3)])
ev_odd = proj(D, Zodd); ev_odd_np = proj(Mnp, Zodd)
chk_lc_odd = max(abs(lc(z)) for z in Zodd); chk_s_even = max(abs(sum(s[i]*z[i] for i in range(n))) for z in Zeven)
ns = lambda v: mp.nstr(v, 15)
res = {"label": "NUMERICAL", "config": [ns(v) for v in x], "L": "0.5", "S_L": [p for p in range(2, 13) if primepow(p)==p and p < Smax],
  "full_kerlc_eigs_D": [ns(v) for v in ev_full],
  "even_face_eigs_arch_minus_prime": [ns(v) for v in ev_even_AP], "even_face_eigs_D": [ns(v) for v in ev_even_D],
  "odd_face_eigs_D": [ns(v) for v in ev_odd], "odd_face_eigs_without_pole": [ns(v) for v in ev_odd_np],
  "check_lc_on_odd": ns(chk_lc_odd), "check_s_on_even": ns(chk_s_even),
  "archTerm": {str(k/100): ns(A_[k]) for k in lags}, "primeTermS": {str(k/100): ns(P_[k]) for k in lags},
  "bigPsiI2": ns(bI2), "seconds": round(time.time()-t0,1)}
print(json.dumps(res, indent=1)); json.dump(res, open('c1_flintquad_results.json','w'), indent=1)
