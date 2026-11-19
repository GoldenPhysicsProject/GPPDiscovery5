"""Round 16 Part 3a': the Round 15 one-shift Hasse block in matrix form (all g in the l = 0.6 window, M = 17 sines)
extended to the Round 16 cutoffs: s <= 14.0 for every model (fakes' data end at e^{14.65}).
rho(s) = || W^{-1/2} Bsym_s W^{-1/2} ||, Hasse block <=> rho <= 1 (Round 15 s2.2).
Arch + pole from round15 hasse_core.Bs (model with no prime data); prime part vectorised here (closed-form X_ij)."""
import sys, os, math, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import prime_window as PW
hc = PW.hc
L = 0.6; M = hc.Mof(L, 45)
def Xvec(t):
    """sum-ready X_ij(t) for t >= 0 array: shape (M, M, n)."""
    a = np.arange(1, M+1)*np.pi/L; ai = a[:, None, None]; aj = a[None, :, None]; t = t[None, None, :]; Y = L - t
    d = ai - aj; s = ai + aj; same = np.abs(d) < 1e-14; dd = np.where(same, 1.0, d)
    T1 = np.where(same, Y*np.cos(aj*t), (np.sin(d*Y - aj*t) - np.sin(-aj*t))/dd)
    T2 = (np.sin(s*Y + aj*t) - np.sin(aj*t))/s
    return (T1 - T2)/L
def prime_mat(u, w):
    P = np.zeros((M, M))
    for sgn in (1, -1):
        m = (u >= 0) if sgn == 1 else (u < 0); uu = np.abs(u[m]); ww = w[m]
        acc = np.zeros((M, M))
        for k in range(0, len(uu), 20000): acc += np.einsum('ijn,n->ij', Xvec(uu[k:k+20000]), ww[k:k+20000])
        P += acc if sgn == 1 else acc.T
    return P
def main():
    D = PW.data(); out = {}
    for nm in ('ZETA', 'FTAU', 'FTH', 'BEUR2', 'EPS', 'DH'):
        base = hc.MODELS[nm](); base = dict(base, n=np.array([1e300]), c=np.array([0.0]))
        W = hc.Bs(base, L, M, 0.0)['B']; W = 0.5*(W + W.T); ev, V = np.linalg.eigh(W); Wm = V @ np.diag(ev**-0.5) @ V.T
        ln, c, wf = D[nm]; S = np.round(np.arange(1.0, 14.0 + 1e-9, 0.05), 4); rows = []
        for s in S:
            ap = hc.Bs(base, L, M, s); i0, i1 = np.searchsorted(ln, s - L), np.searchsorted(ln, s + L); u = ln[i0:i1]
            cw = c[i0:i1] if wf is None else c[i0:i1]*wf(u)
            B = ap['arch'] + ap['pole'] - prime_mat(u - s, cw*np.exp(-u/2)); Bs = 0.5*(B + B.T)
            rho = float(np.max(np.abs(np.linalg.eigvalsh(Wm @ Bs @ Wm)))); rows.append((float(s), rho))
        r = np.array(rows); viol = r[r[:, 1] > 1]
        out[nm] = dict(lam_min_W=float(ev[0]), max_rho=float(r[:, 1].max()), s_at_max=float(r[np.argmax(r[:, 1]), 0]),
                       n_violating=int(len(viol)), n_scanned=len(r), first_violation=float(viol[0, 0]) if len(viol) else None,
                       rho_table=[dict(s=a, rho=b) for a, b in rows[::4]])
        print(nm, {k: v for k, v in out[nm].items() if k != 'rho_table'}, flush=True)
    json.dump(out, open('hasse_matrix16.json', 'w'), indent=1)
if __name__ == '__main__': main()
