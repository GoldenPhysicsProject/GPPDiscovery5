"""S27 [design, UNTRUSTED float LP]: cutting-plane version of s26_lp_odd.lp_enriched (PL-W exterior, no DC trapezoid).
Solve the LP on a coarse xi set, evaluate m on a fine grid (step FINE), add all local minima with m < c_LP + tol
(plus +-FINE neighbours), re-solve; stop when fine-grid min m >= (1 - rtol) c_LP. Same multiplier family as S26. RH OPEN."""
import numpy as np, json, sys, time, argparse
from scipy.optimize import linprog
from s26_lp_odd import sigA, wPL_columns, choose_Xi

def rows(xi, h, knots, u):
    A0, Ap = wPL_columns(xi, h, knots)
    return np.hstack([np.ones((xi.size, 1)), -A0[:, None]] + [-a[:, None] for a in Ap] + [np.cos(np.outer(xi, u))])

def m_eval(xi, h, knots, u, v, a, chunk=20000):
    out = np.empty_like(xi)
    for s in range(0, xi.size, chunk):
        x = xi[s:s+chunk]; A0, Ap = wPL_columns(x, h, knots)
        out[s:s+chunk] = sigA(x) + v[0]*A0 + sum(vp*A for vp, A in zip(v[1:], Ap)) - np.cos(np.outer(x, u))@a
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h"); ap.add_argument("--U", type=float, default=12.0); ap.add_argument("--D", type=float, default=0.05)
    ap.add_argument("--B", type=float, default=3.0); ap.add_argument("--pl", default="16,6"); ap.add_argument("--pl-knots", default="")
    ap.add_argument("--coarse", type=float, default=0.01); ap.add_argument("--fine", type=float, default=0.0005)
    ap.add_argument("--rtol", type=float, default=0.02); ap.add_argument("--iters", type=int, default=30)
    ap.add_argument("--out", required=True); ap.add_argument("--tag", default="")
    A = ap.parse_args(); t0 = time.time()
    h = np.log(3)/2 if A.h == "hstar" else float(A.h)
    if A.pl_knots: knots = 2*h + np.array([float(x) for x in A.pl_knots.split(",")])
    else:
        M, L = A.pl.split(","); knots = 2*h + np.linspace(0, float(L), int(M)+1)
    K = int(round((A.U - 2*h)/A.D)) + 1; u = 2*h + A.D*np.arange(K); Xi = choose_Xi(A.B)
    nph, bs = 8, 1500; th = np.linspace(0, 2*np.pi, bs, endpoint=False)
    npl = len(knots) - 1
    A2 = np.hstack([np.zeros((bs*nph, 1 + npl)), np.vstack([np.cos(np.outer(th, np.arange(K)) + ph) for ph in 2*np.pi*np.arange(nph)/nph])])
    X = np.arange(0.0, Xi, A.coarse); xf = np.arange(0.0, Xi, A.fine)
    bnds = [(None, None)]*(1 + npl + K); cobj = np.r_[-1.0, np.zeros(npl + K)]
    hist = []
    for it in range(A.iters):
        A1 = rows(X, h, knots, u)
        r = linprog(cobj, A_ub=np.vstack([A1, A2]), b_ub=np.r_[sigA(X), np.full(bs*nph, A.B)], bounds=bnds, method="highs")
        if r.status != 0: print("LP status", r.status, r.message, flush=True); break
        c = -r.fun; v = r.x[1:1+npl]; a = r.x[1+npl:]
        mf = m_eval(xf, h, knots, u, v, a); i = int(mf.argmin()); fmin = float(mf[i])
        hist.append(dict(it=it, n=int(X.size), c_LP=float(c), fine_min=fmin, at=float(xf[i])))
        print(f"it={it} rows={X.size} c_LP={c:.6e} fine_min={fmin:.6e} at {xf[i]:.4f} ({time.time()-t0:.0f}s)", flush=True)
        if fmin >= (1 - A.rtol)*c: break
        thr = c + 0.1*abs(c)
        loc = np.where((mf[1:-1] <= mf[:-2]) & (mf[1:-1] <= mf[2:]) & (mf[1:-1] < thr))[0] + 1
        new = np.unique(np.r_[xf[loc], xf[np.maximum(loc-1, 0)], xf[np.minimum(loc+1, xf.size-1)]])
        X = np.unique(np.r_[X, new])
    out = dict(h=float(h), h_spec="log3/2" if A.h == "hstar" else A.h, U=A.U, D=A.D, B=A.B, Dd_list=[], useDC=False, C=[],
               pl_knots=list(map(float, knots)), xi_step=f"cut coarse={A.coarse} fine={A.fine}", Xi=Xi, K=K, nd=0, npl=npl,
               status=int(r.status), c_LP=float(-r.fun) if r.status == 0 else None,
               v_pl=r.x[1:1+npl].tolist() if r.status == 0 else None, a=r.x[1+npl:].tolist() if r.status == 0 else None,
               u0=float(u[0]), hist=hist, tag=A.tag, time_s=time.time()-t0)
    json.dump(out, open(A.out, "w"), indent=1); print("wrote", A.out, flush=True)
if __name__ == "__main__": main()
