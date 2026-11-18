"""Round 14E2: replaces the withdrawn 14E Gaussian fit.  lambda_min(W_zeta) at 80 digits (mpW.py logs, best = largest M per L),
fits on 0.25 <= L <= 2.0: geometric a+bL, quadratic, a - B e^L, a + bL - B e^L; local slope d log(lambda)/d(e^L).  NUMERICAL."""
import re, json, math, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); best = {}
for f in ['mp80/W.log', 'mp80/W_fit.log']:
    for line in open(os.path.join(HERE, f)):
        m = re.match(r"ZETA (\S+) (\d+) deg 7 low \['([^']+)'", line)
        if m:
            L, M, lam = float(m.group(1)), int(m.group(2)), float(m.group(3))
            if L not in best or M > best[L][0]: best[L] = (M, lam)
Ls = np.array(sorted(best)); y = np.log([best[L][1] for L in Ls]); E = np.exp(Ls)
out = dict(data=[dict(L=float(L), M=best[L][0], lam_min=best[L][1]) for L in Ls], fits={})
for name, X in [('a+bL', [np.ones_like(Ls), Ls]), ('a+bL+cL^2', [np.ones_like(Ls), Ls, Ls**2]), ('a-Be^L', [np.ones_like(Ls), -E]),
                ('a+bL-Be^L', [np.ones_like(Ls), Ls, -E])]:
    X = np.array(X).T; c, *_ = np.linalg.lstsq(X, y, rcond=None); r = np.abs(X @ c - y).max()
    out['fits'][name] = dict(coef=[float(x) for x in c], max_log_resid=float(r))
q = np.polyfit(Ls[Ls <= 1.25], y[Ls <= 1.25], 2)    # the 14E Gaussian, refit on its own 5 points, then extrapolated
out['gaussian_14E_extrapolated'] = {str(L): dict(pred=float(np.polyval(q, L)), meas=float(yy)) for L, yy in zip(Ls, y) if L > 1.25}
out['local_slope_dlog_d_eL'] = [dict(Lmid=float((Ls[i]+Ls[i+1])/2), slope=float((y[i+1]-y[i])/(E[i+1]-E[i]))) for i in range(len(Ls)-1)]
out['prolate_heuristic_B'] = 4*math.pi
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'fits80.json'), 'w'), indent=1)
