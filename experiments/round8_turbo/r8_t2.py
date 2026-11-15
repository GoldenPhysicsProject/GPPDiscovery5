"""Round 8 (T2): closed-form q and accurate rho_j (NUMERICAL test of CHECKED identities). Float64, chunked Gauss-Legendre to T, + mean tail."""
import numpy as np, json, sys
from scipy import special
l = 0.25; NS = [(2,2),(3,3),(4,2),(5,5),(7,7),(8,2),(9,3),(11,11)]
bI2 = 16*np.sinh(l)**2; N = 12
om = lambda j: j*np.pi/l
w_ = lambda j: np.sqrt(2/l)*(-1)**(j+1)/om(j)
sig = lambda j: -np.sqrt(2/l)*(-1)**j*om(j)*np.sinh(l/2)/(0.25+om(j)**2)
def S(j, t):
    w = om(j); return np.sqrt(2/l)*w*l*np.sinc((t-w)*l/np.pi)/(t+w)
JJ = [73, 123, 173, 273, 472]
sx = np.array([sig(i) for i in range(1, N+1)] + [np.sinh(l/2)])
x0, w0 = np.polynomial.legendre.leggauss(40)
# mean of sin^4(t/2) sin^2(t/4) over its period 4 pi
tt = np.linspace(0, 4*np.pi, 400001)[:-1]; m = float(np.mean(np.sin(tt/2)**4*np.sin(tt/4)**2))
def run(T, h=2.0, chunk=50000):
    qI = np.zeros(N+1); rI = np.zeros((N+1, len(JJ))); vI = np.zeros((N+1, len(JJ)))
    npan = int(T/h); a0 = 0
    while a0 < npan:
        a1 = min(npan, a0+chunk); ed = np.arange(a0, a1)*h
        t = (ed[:,None] + h/2*(x0[None,:]+1)).ravel(); wt = np.tile(h/2*w0, a1-a0)
        V = special.digamma(0.25+0.5j*t).real - np.log(np.pi) - 2*sum(np.log(p)/np.sqrt(n)*np.cos(t*np.log(n)) for n,p in NS)
        s4 = 16*np.sin(t/2)**4                                    # Psi^2 t^4
        Psi2 = np.where(t > 1e-3, s4/np.maximum(t,1e-300)**4, 1.0)
        Fx = np.vstack([S(i, t) for i in range(1, N+1)] + [np.sin(t*l)])
        base = wt*V/np.pi
        qI += (Fx*(base*Psi2*t*t*np.sin(t*l))).sum(1)
        for k, j in enumerate(JJ):
            Sj = S(j, t)
            rI[:, k] += (Fx*(base*s4*Sj)).sum(1)
            vI[:, k] += (Fx*(base*Psi2*(t/om(j))**2*Sj)).sum(1)
        a0 = a1
    # analytic mean tails (delta row only; phi_i rows decay faster): V ~ log(t/2pi), sin(tl)S_j ~ sqrt(2/l)(-1)^j om sin^2/t^2
    tail = (np.log(T/(2*np.pi))+1)/T
    qT = qI.copy(); qT[N] += (16/np.pi)*m*tail
    rT = rI.copy()
    for k, j in enumerate(JJ): rT[N, k] += (16/np.pi)*np.sqrt(2/l)*(-1)**j*om(j)*m*tail
    q = qT + 0.5*bI2**2*sx*np.sinh(l/2)
    r4 = rT - (bI2**2/8)*np.outer(sx, [sig(j) for j in JJ])          # omega^4 rho_j
    v = vI - 2*bI2**2*np.outer(sx, [sig(j)-w_(j)*np.sinh(l/2) for j in JJ])
    return q, r4, v
out = dict(mean_sin4_sin2=m, runs={})
for T in [8000.0, 1e5, 4e5, 1e6]:
    q, r4, v = run(T)
    ident = [float(np.max(np.abs(v[:,k] - (w_(j)/om(j)**2*q + r4[:,k]/om(j)**4)))/np.max(np.abs(v[:,k]))) for k,j in enumerate(JJ)]
    out['runs'][str(T)] = dict(q=q.tolist(), omega4_rho_maxabs=[[j, float(np.max(np.abs(r4[:,k]))), float(r4[N,k])] for k,j in enumerate(JJ)],
                               omega3_v_maxabs=[[j, float(np.max(np.abs(v[:,k]))*om(j)**3)] for k,j in enumerate(JJ)], identity_relerr=ident,
                               omega4_rho_rows=r4.tolist())
    print(T, 'q_delta', q[N], 'w4rho', [round(float(np.max(np.abs(r4[:,k]))),5) for k in range(len(JJ))], 'ident', max(ident), flush=True)
json.dump(out, open('r8_t2_results.json','w'), indent=1)
