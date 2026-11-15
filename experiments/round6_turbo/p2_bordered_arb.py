"""Round 6: certified leading-order Schur condition (NUMERICAL, interval-certified via arb; not Lean).
Bordered Gram of the odd-face form on {phi_1..phi_12, delta_l}, l = 1/4:  g = B(delta_l, delta_l) = G(l,l),  u_i = B(phi_i, delta_l),
with delta_l's transform F(t) = sin(t l), sigma = sinh(l/2). Claim certified: g > 0 and M11 - u u^T/g - beta0 I is positive definite.
M11 = the Round 5 certified block (midpoints from p1_oddblock_L05_N12_results.json, each widened by its filed max radius + 1e-20 for
the 20-digit rounding). Integral on [0,T] by acb.integral; tail on [T,inf) (CHECKED by hand, same lemmas as Round 5):
 Psi^2 <= 16/t^4, |V| <= log t + C (C = 2.5 + 2 P(0)), |sin| <= 1, |S_i| <= (4/3) w_i / t^2 (t >= 2 w_12):
 |tail_g| <= (16/pi) T^-3 [(log T + C)/3 + 1/9],  |tail_u_i| <= (16/pi)(4/3) w_i sqrt(2/l) T^-5 [(log T + C)/5 + 1/25]."""
from flint import acb, arb, ctx
import json, time
ctx.prec = 160
T = 4000; N = 12
l = arb(1)/4; pi = arb.pi(); lp = pi.log()
NS = [(2,2),(3,3),(4,2),(5,5),(7,7),(8,2),(9,3),(11,11)]
W = [j*pi/l for j in range(1, N+1)]; c2 = (2/l).sqrt()
bI2 = 16*l.sinh()**2
def S(j, t): w = W[j]; return c2*w*l*((t-w)*l).sinc()/(t+w)
def V(t):
    i = acb(0,1)
    return ((acb(l)+i*t/2).digamma()+(acb(l)-i*t/2).digamma())/2 - lp - 2*sum(arb(p).log()/arb(n).sqrt()*(t*arb(n).log()).cos() for n,p in NS)
P0 = sum(arb(p).log()/arb(n).sqrt() for n,p in NS); C = 2.5 + 2*P0; TT = arb(T)
sd = (l/2).sinh()
sig = [-c2*(-1)**(j+1)*W[j]*sd/(arb(1)/4+W[j]**2) for j in range(N)]
cuts = list(range(0, T+1, 50)); t0 = time.time()
def integ(f): return sum((acb.integral(f, cuts[m], cuts[m+1]) for m in range(len(cuts)-1)), acb(0)).real
tg = 16/pi*TT**-3*((TT.log()+C)/3+arb(1)/9)
g = integ(lambda t, a: (t/2).sinc()**4*V(t)*(t*l).sin()**2/pi) + arb(0, float(tg.upper())) - 2*bI2**2*sd**2
print("g", g, round(time.time()-t0,1), flush=True)
u = []
for i in range(N):
    tu = 16/pi*arb(4)/3*W[i]*c2*TT**-5*((TT.log()+C)/5+arb(1)/25)
    u.append(integ(lambda t, a, i=i: (t/2).sinc()**4*V(t)*S(i,t)*(t*l).sin()/pi) + arb(0, float(tu.upper())) - 2*bI2**2*sig[i]*sd)
    print("u", i+1, u[-1], round(time.time()-t0,1), flush=True)
R5 = json.load(open('../round5_turbo/p1_oddblock_L05_N12_results.json'))
rad = R5['max_entry_radius'] + 1e-20
M = [[arb(R5['M'][j][k]) + arb(0, rad) for k in range(N)] for j in range(N)]
assert g > 0
A0 = [[M[j][k] - u[j]*u[k]/g for k in range(N)] for j in range(N)]
import mpmath as mp; mp.mp.dps = 40
ev = sorted(mp.eigsy(mp.matrix([[mp.mpf(A0[j][k].mid().str(40, radius=False)) for k in range(N)] for j in range(N)]), eigvals_only=True))
beta0 = arb(mp.nstr(mp.mpf(ev[0])*mp.mpf('0.9'), 6))
A = [[A0[j][k] - (beta0 if j == k else 0) for k in range(N)] for j in range(N)]
Lm = [[arb(0)]*N for _ in range(N)]; D = [arb(0)]*N
for j in range(N):
    D[j] = A[j][j] - sum((Lm[j][m]**2*D[m] for m in range(j)), arb(0))
    assert D[j] > 0, ("pivot not certified", j, D[j])
    for i in range(j+1, N):
        Lm[i][j] = (A[i][j] - sum((Lm[i][m]*Lm[j][m]*D[m] for m in range(j)), arb(0)))/D[j]
out = dict(label="NUMERICAL (interval-certified; not Lean)", T=T, g=g.str(15), g_tail_radius=tg.str(5), u=[x.str(12) for x in u],
           lam_min_float_schur_leading=mp.nstr(ev[0], 12), beta0_certified=beta0.str(10),
           LDL_pivots_lower=[D[j].lower().str(6) for j in range(N)], seconds=round(time.time()-t0,1))
print(json.dumps(out, indent=1)); json.dump(out, open('p2_bordered_arb_results.json','w'), indent=1)
