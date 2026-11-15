"""Round 5 (Turbo 6889360979 item 2): P1 at L = 0.5, certified finite odd block. NUMERICAL (interval arithmetic, arb via python-flint).
Form on the odd face: B(f,g) = 1/2 ∫∫_{[0,l]^2} [K_S(s-u) - K_S(s+u)] f(s) g(u) ds du, l = L/2 = 1/4,
K_S = pole - prime_{S_L} + arch exactly as in-tree (pole = 2 Psi(i/2)^2 cosh(d/2), Psi(i/2) = 16 sinh^2(1/4)).
Basis E_N: phi_j = sqrt(2/l) sin(w_j x), w_j = j pi / l, j = 1..N (L^2-orthonormal on [0,l]; = first N odd modes on [-l,l]).
Exact reduction (CHECKED by hand):
  B(phi_j,phi_k) = (2/l) { (1/pi) ∫_0^∞ Psi(t)^2 V(t) S_j(t) S_k(t) dt  -  2 Psi(i/2)^2 sigma_j sigma_k }
  S_j(t) = ∫_0^l sin(w_j x) sin(t x) dx = w_j l sinc((t-w_j) l)/(t+w_j),   sigma_j = ∫_0^l sinh(x/2) sin(w_j x) dx = -(-1)^j w_j sinh(l/2)/(1/4+w_j^2),
  V(t) = Re psi(1/4+it/2) - log pi - 2 P(t),  P(t) = sum_{n in {2,3,4,5,7,8,9,11}} Lambda(n) n^{-1/2} cos(t log n)
  (only these prime powers of S_L have triG(log n ± d) != 0 for |d| <= L = 0.5, since e^{2.5} = 12.18).
Rigorous parts: acb.integral on [0,T] (rigorous error bounds), explicit tail bound on [T,∞), LDL^T of M - beta*I in ball arithmetic.
Tail bound (CHECKED by hand): Psi^2 <= 16/t^4; |w(t)| <= log t + 2.5 for t >= 100 (psi(z) = log z - 1/(2z) - ∫ f(u) e^{-zu} du, 0 <= f(u) <= u/12, Re z = 1/4 gives |.| <= 4/3);
  |P| <= P(0); |S_j| <= w_j/(t^2-w_j^2) <= (4/3) w_j / t^2 for t >= 2 w_max. So |tail_jk| <= (1/pi) 16 (16/9) w_j w_k T^-7 [(log T + C)/7 + 1/49], C = 2.5 + 2 P(0).
No zero data."""
from flint import acb, arb, ctx
import json, time, sys
ctx.prec = 160
N = int(sys.argv[1]) if len(sys.argv) > 1 else 12
T = 800
l = arb(1)/4; pi = arb.pi(); lp = pi.log()
NS = [(2,2),(3,3),(4,2),(5,5),(7,7),(8,2),(9,3),(11,11)]
W = [j*pi/l for j in range(1, N+1)]
bI2 = 16*(l).sinh()**2
def S(j, t):
    w = W[j]; return w*l*((t-w)*l).sinc()/(t+w)
def V(t):
    i = acb(0,1)
    wv = ((acb(l) + i*t/2).digamma() + (acb(l) - i*t/2).digamma())/2 - lp
    return wv - 2*sum(arb(p).log()/arb(n).sqrt()*(t*arb(n).log()).cos() for n,p in NS)
def integrand(j, k):
    return lambda t, a: (t/2).sinc()**4 * V(t) * S(j,t) * S(k,t) / pi
P0 = sum(arb(p).log()/arb(n).sqrt() for n,p in NS)
C = 2.5 + 2*P0
TT = arb(T)
assert T >= 2*float(W[-1].mid()) and T >= 100
tailc = 16*arb(16)/9/pi * TT**-7 * ((TT.log() + C)/7 + arb(1)/49)
sigma = [-(-1)**(j+1)*W[j]*(l/2).sinh()/(arb(1)/4 + W[j]**2) for j in range(N)]
t0 = time.time()
cuts = list(range(0, T+1, 50))
M = [[None]*N for _ in range(N)]
for j in range(N):
    for k in range(j, N):
        f = integrand(j, k)
        I = sum((acb.integral(f, cuts[m], cuts[m+1]) for m in range(len(cuts)-1)), acb(0))
        assert abs(float(I.imag.mid())) < 1e-30
        val = I.real + arb(0, float((tailc*W[j]*W[k]).upper())) - 2*bI2**2*sigma[j]*sigma[k]
        M[j][k] = M[k][j] = (2/l)*val
    print("row", j, round(time.time()-t0,1), flush=True)
import mpmath as mp
mp.mp.dps = 40
Mm = mp.matrix([[mp.mpf(M[j][k].mid().str(40, radius=False)) for k in range(N)] for j in range(N)])
ev = sorted(mp.eigsy(Mm, eigvals_only=True))
beta_try = mp.mpf(ev[0])*mp.mpf('0.9')
beta = arb(mp.nstr(beta_try, 6))
# LDL^T of M - beta I in ball arithmetic
A = [[M[j][k] - (beta if j == k else 0) for k in range(N)] for j in range(N)]
Lm = [[arb(0)]*N for _ in range(N)]; D = [arb(0)]*N
for j in range(N):
    D[j] = A[j][j] - sum((Lm[j][m]**2*D[m] for m in range(j)), arb(0))
    assert D[j] > 0, ("pivot not certified", j, D[j])
    for i in range(j+1, N):
        Lm[i][j] = (A[i][j] - sum((Lm[i][m]*Lm[j][m]*D[m] for m in range(j)), arb(0)))/D[j]
out = dict(label="NUMERICAL (interval-certified; not Lean)", L="0.5", N=N, T=T, prec_bits=ctx.prec,
           prime_powers=[n for n,_ in NS], P0=P0.str(15), tail_coeff_per_wjwk=tailc.str(5),
           beta_certified=beta.str(10), lambda_min_float=mp.nstr(ev[0], 12), eigs_float=[mp.nstr(e, 8) for e in ev],
           LDL_pivots_lower=[D[j].lower().str(6) for j in range(N)],
           max_entry_radius=max(float(M[j][k].rad()) for j in range(N) for k in range(N)),
           diag=[M[j][j].str(12) for j in range(N)], M=[[M[j][k].mid().str(20, radius=False) for k in range(N)] for j in range(N)],
           pole_rank1_sigma=[s.str(12) for s in sigma], seconds=round(time.time()-t0,1))
print(json.dumps({k: v for k, v in out.items() if k != "M"}, indent=1))
json.dump(out, open(f"p1_oddblock_L05_N{N}_results.json", "w"), indent=1)
