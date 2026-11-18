"""Round 14D face 3 (b): the global Plancherel norm of the Connes-Consani map,
   N_F(g) := ||E_F g||^2 = (1/2pi) int |F(1/2+it)|^2 |g^(t)|^2 dt   (exact for g with g^ vanishing at F's poles; here the
   form itself, from F on the line = integer data via Hurwitz sums, no zeros), versus the Weil form W_F(g) (prime side).
Question: W = c N - R with R of definite sign?  Pencil W v = mu N v:  mu_min = sup{c : W >= c N}, mu_max = inf{c : W <= c N}.
Also the 'identity' check W - N (sign of R = N - W)."""
import sys, os, math, json, numpy as np, scipy.linalg as sl
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'round14b_euler')); sys.path.insert(0, os.path.join(HERE, '..'))
import weil_galerkin as wg, models14b as md
def phihat(L, M, t):
    k = (np.arange(1, M+1)*np.pi/L)[:, None]; j = np.arange(1, M+1)[:, None]; t = t[None, :]
    den = k**2 - t**2; den = np.where(np.abs(den) < 1e-9, 1e-9, den)
    return np.exp(-1j*t*L/2)*math.sqrt(2/L)*k*(1 - (-1.0)**j*np.exp(1j*t*L))/den
def Nmat(L, M, line):
    t, F2 = line; P = phihat(L, M, t); wt = F2.copy(); wt[0] *= 0.5; wt[-1] *= 0.5
    dt = t[1] - t[0]; return (P*wt) @ P.conj().T * dt / np.pi
out = {}
models = {'ZETA': md.zeta(), 'EPS': md.eps(), 'DH': md.dh(), 'FTH': md.fth(0.2)}
lines = {k: np.load(os.path.join(HERE, f'line_{k}.npy')) for k in models}
# Plancherel sanity (F = 1): N should be the identity up to the t > 200 tail
for L in [1.0, 4.0]:
    M = int(math.ceil(45*L/math.pi))+8; t = lines['ZETA'][0]; I = Nmat(L, M, (t, np.ones_like(t))).real
    out[f'plancherel_check_L{L}'] = dict(max_offdiag=float(np.abs(I-np.diag(np.diag(I))).max()), diag_min=float(np.diag(I).min()), diag_max=float(np.diag(I).max()))
    print(out[f'plancherel_check_L{L}'], flush=True)
Ls = [math.log(2)-0.05, math.log(2)+0.05, 1.0, 2.0, 3.0, 4.0]
for name, mo in models.items():
    rows = []
    for L in Ls:
        M = max(16, int(math.ceil(45*L/math.pi))+8)
        a, p, q = wg.operators(mo, L, M); W = a+p-q; N = Nmat(L, M, lines[name]).real; N = 0.5*(N+N.T)
        mu = sl.eigh(W, N, eigvals_only=True); R = N - W; evR = np.linalg.eigvalsh(R); evN = np.linalg.eigvalsh(N)
        rows.append(dict(L=L, M=M, mu_min=float(mu[0]), mu_max=float(mu[-1]), lam_min_N=float(evN[0]), lam_max_N=float(evN[-1]),
                         R_min=float(evR[0]), R_max=float(evR[-1]), lam_min_W=float(np.linalg.eigvalsh(W)[0])))
        print(name, rows[-1], flush=True)
    out[name] = rows
json.dump(out, open(os.path.join(HERE, 'pencil.json'), 'w'), indent=1)
