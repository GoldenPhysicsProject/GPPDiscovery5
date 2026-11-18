"""Round 14B part 1 - Daniel's No-Ghosts theorem (idea: Daniel Toupin).
Reflected pairing on the integer state space: <f,f>_Theta := W(f * Theta f), Theta f(x) = x^{-1} conj f(1/x)
(time reflection x -> 1/x, unitary on L^2(dx), combined with conjugation). Computed from integer data only
(prime side c(n), own arch, own poles; weil_galerkin.py). Ghost = negative-norm direction of the pairing.
(A) Minkowski/Weil-reflected pairing: lam_min, ghost count, ghost frequency, Rayleigh parts, per model and width L.
(B) Euclidean OS pairing  <Theta f, f>_OS = int int f(u) f(v) k(u+v), f supported at positive time u in (0,T):
    the literal Osterwalder-Schrader form of the same kernel. Zeros never used (they appear only as labels)."""
import sys, os, math, json, time, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))
import weil_galerkin as wg, models14b as md
OUT = os.path.join(HERE, os.environ.get('NG_OUT', 'noghosts.json'))
res = json.load(open(OUT)) if os.path.exists(OUT) else {}
def save(): json.dump(res, open(OUT, 'w'), indent=1)

def run_A(model, Ls, fmax):
    rows = res.setdefault('A', {}).setdefault(f"{model['name']}_f{fmax}", [])
    done = {r['L'] for r in rows}
    for L in Ls:
        if L in done: continue
        M = int(math.ceil(fmax*L/math.pi)) + 8; t0 = time.time()
        ar, po, pr = wg.operators(model, L, M); W = ar + po - pr
        ev, V = np.linalg.eigh(W); sc = np.abs(W).max(); neg = np.nonzero(ev < -1e-7*sc)[0]
        gh = [dict(lam=float(ev[k]), freq=wg.peak_freq(L, V[:, k], fmax=fmax+10),
                   arch=float(V[:, k]@ar@V[:, k]), pole=float(V[:, k]@po@V[:, k]), prime=float(V[:, k]@pr@V[:, k])) for k in neg[:12]]
        rows.append(dict(L=L, M=M, lam_min=float(ev[0]), scale=float(sc), n_ghosts=int(len(neg)), ghosts=gh, sec=round(time.time()-t0, 1)))
        print(model['name'], fmax, L, M, ev[0], len(neg), [round(g['freq'], 1) for g in gh], flush=True); save()

def run_B(model, Ts, fmax=45.0, NQ=900):
    rows = res.setdefault('B_OS', {}).setdefault(model['name'], [])
    done = {r['T'] for r in rows}
    XG, WG = np.polynomial.legendre.leggauss(NQ)
    for T in Ts:
        if T in done: continue
        M = int(math.ceil(fmax*T/math.pi)) + 8
        u = 0.5*T*(XG+1); wu = 0.5*T*WG
        phi = lambda x: np.where((x >= 0) & (x <= T), math.sqrt(2/T)*np.sin(np.arange(1, M+1)[:, None]*np.pi*x/T), 0.0)
        P = phi(u)*wu; S = u[:, None] + u[None, :]
        K = -model['w'](S) if model.get('C') is not None else np.zeros_like(S)
        for a, wt in model.get('poles', []): K = K + wt*2*np.cosh(a*S)
        Mos = P @ K @ P.T
        x2, w2 = np.polynomial.legendre.leggauss(200)
        m = np.log(model['n']) < 2*T
        for n, cc in zip(model['n'][m], model['c'][m]):
            l = math.log(n); lo, hi = max(0.0, l-T), min(T, l)
            if hi <= lo: continue
            uu = 0.5*(hi-lo)*(x2+1)+lo; ww = 0.5*(hi-lo)*w2
            H = (phi(uu)*ww) @ phi(l-uu).T
            Mos -= cc/math.sqrt(n)*0.5*(H+H.T)
        ev = np.linalg.eigvalsh(Mos); sc = np.abs(Mos).max()
        rows.append(dict(T=T, M=M, lam_min=float(ev[0]), lam_max=float(ev[-1]), n_neg=int(np.sum(ev < -1e-7*sc)),
                         n_pos=int(np.sum(ev > 1e-7*sc))))
        print('OS', model['name'], T, M, ev[0], rows[-1]['n_neg'], rows[-1]['n_pos'], flush=True); save()

if __name__ == '__main__':
    part = sys.argv[1] if len(sys.argv) > 1 else 'all'
    Ls = [1.0, 2.0, 3.0, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0]
    if part in ('A', 'all'):
        Z = md.zeta(); run_A(Z, Ls, 45.0)
        E = md.eps(); run_A(E, Ls, 45.0)
        B = md.beur2(); run_A(B, Ls, 45.0)
        F = md.fth(0.2); run_A(F, Ls, 45.0)
        D = md.dh(); run_A(D, Ls, 125.0)
        run_A(Z, [2.0, 4.0, 6.0, 8.0], 125.0)
    if part in ('B', 'all'):
        Z = md.zeta(); run_B(Z, [0.3, 0.5, 1.0, 2.0, 3.0, 4.0])
        run_B(md.eps(), [1.0, 2.0, 4.0]); run_B(md.dh(), [1.0, 2.0, 4.0])
