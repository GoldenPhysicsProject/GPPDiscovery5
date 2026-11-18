"""Round 13g: Baez-Duarte distance d_N^2 = inf || chi_(0,1] - sum_{k<=N} c_k {1/(kx)} ||^2_{L^2(0,inf)}  (RH <=> d_N -> 0).
Integer-only Gram (Vasyunin): for coprime h,k
  <e_h,e_k> = (log 2pi - gamma)/2 (1/h+1/k) + (k-h)/(2hk) log(h/k) - pi/(2hk) (V(h,k)+V(k,h)),  V(h,k)=sum_{m<k} {mh/k} cot(pi m/k)
  <e_{dh},e_{dk}> = <e_h,e_k>/d ;  <chi, e_k> = (log k + 1 - gamma)/k.
One Cholesky gives every d_n^2, n <= N (prefix property). Formula checked against direct quadrature for small (h,k)."""
import mpmath as mp, json, math, sys
mp.mp.dps = 60
N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
def V(h, k):
    return mp.fsum(mp.frac(mp.mpf(m*h)/k)*mp.cot(mp.pi*m/k) for m in range(1, k))
def gram_coprime(h, k):
    return (mp.log(2*mp.pi) - mp.euler)/2*(mp.mpf(1)/h + mp.mpf(1)/k) + mp.mpf(k-h)/(2*h*k)*mp.log(mp.mpf(h)/k) - mp.pi/(2*h*k)*(V(h, k) + V(k, h))
def gram(h, k):
    d = math.gcd(h, k); return gram_coprime(h//d, k//d)/d
# direct quadrature check: <e_h,e_k> = int_0^inf {y/h}{y/k} dy/y^2
def direct(h, k, periods=4000):
    P = h*k//math.gcd(h, k); pts = sorted(set([0] + [h*i for i in range(1, periods*P//h+1)] + [k*i for i in range(1, periods*P//k+1)]))
    tot = mp.mpf(0)
    for a, b in zip(pts[:-1], pts[1:]):
        A, B = mp.mpf(a//h), mp.mpf(a//k)  # floors on (a,b)
        f = lambda y: (y/h - A)*(y/k - B)/y**2
        tot += (mp.mpf(1)/(h*k))*(b-a) if a == 0 else mp.quad(f, [a, b])
    Y = pts[-1]; mean = mp.quad(lambda y: mp.frac(y/h)*mp.frac(y/k), mp.linspace(0, P, 2*P+1))/P
    return tot + mean/Y
out = dict(check={})
mp.mp.dps = 20
for h, k in [(1, 1), (1, 2), (2, 3), (3, 5)]:
    out['check'][f'{h},{k}'] = dict(vasyunin=float(gram(h, k)), quadrature=float(direct(h, k, 300)))
print(out['check'], flush=True)
mp.mp.dps = 60
G = mp.matrix(N, N)
for i in range(1, N+1):
    for j in range(i, N+1):
        G[i-1, j-1] = G[j-1, i-1] = gram(i, j)
b = mp.matrix([(mp.log(k) + 1 - mp.euler)/k for k in range(1, N+1)])
L = mp.cholesky(G)
z = mp.zeros(N, 1); s = mp.mpf(0); d2 = {}
for n in range(N):
    z[n] = (b[n] - mp.fsum(L[n, j]*z[j] for j in range(n)))/L[n, n]
    s += z[n]**2; d2[n+1] = 1 - s
C0 = 2 + mp.euler - mp.log(4*mp.pi)
out['C0_2+gamma-log4pi'] = float(C0)
out['d2'] = {n: float(d2[n]) for n in d2}
out['d2_logN'] = {n: float(d2[n]*mp.log(n)) for n in d2 if n > 1}
out['min_cholesky_pivot'] = float(min(L[i, i] for i in range(N)))
json.dump(out, open('nyman_beurling.json', 'w'), indent=1)
for n in (1, 2, 5, 10, 20, 50, 100, 150, 200, 250, 300):
    if n <= N: print(n, float(d2[n]), float(d2[n]*mp.log(n)) if n > 1 else None, flush=True)
