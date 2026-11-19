"""Round 15 parts (a)+(c): the one-shift Hasse block.
Z-analogue of the F_q 2x2 block on {Delta, Gamma_F}:  span{g, tau_s g}, g in the width-l window (l < log 2: no prime sees g alone).
  Block PSD for all g in window  <=>  W_l +- Bsym_s >= 0  (M x M),  Bsym_s = (B_s + B_s^T)/2, s >= l (disjoint, ||g +- tau g||^2 = 2).
  Hasse ratio rho(s) = max |eig(W_l^-1/2 Bsym_s W_l^-1/2)|  (needs W_l > 0);  Hasse inequality  <=>  rho(s) <= 1.
F_q dictionary (normalised by sqrt q):  d = W(g) <-> 2g (genus),  b_s <-> -a/sqrt q = (N_1 - q - 1)/sqrt q; the pole part of b_s
  <-> (q+1)/sqrt q = sqrt q + 1/sqrt q (here p^{1/2}+p^{-1/2} times |g^(i/2)|^2), the prime part <-> N_1/sqrt q.
Fixed-g decomposition also recorded for g = ground state of W_l and g = first sine.  Usage: python3 hasse_blocks.py MODEL [l]"""
import sys, json, math, numpy as np, hasse_core as hc
name = sys.argv[1]; l = float(sys.argv[2]) if len(sys.argv) > 2 else 0.6
fmax = 125.0 if name == 'DH' else 45.0
mo = hc.MODELS[name](); M = hc.Mof(l, fmax)
W = hc.Bs(mo, l, M, 0.0)['B']; W = 0.5*(W + W.T); ew, Vw = np.linalg.eigh(W)
g0 = Vw[:, 0]; g1 = np.zeros(M); g1[0] = 1.0
Smax = math.log(mo['n'].max()) - l - 0.05
primes = [int(n) for n in mo['n'] if n <= math.exp(Smax) and all(int(n) % q for q in range(2, int(n**0.5)+1)) and n == int(n)] if name in ('ZETA', 'FTH', 'FTAU', 'BEUR2') else []
if name in ('EPS', 'DH'):
    primes = [p for p in range(2, int(math.exp(Smax))+1) if all(p % q for q in range(2, int(p**0.5)+1))]
grid = list(np.round(np.arange(l, Smax, 0.02), 4))
pl = [math.log(p) for p in primes if math.log(p) >= l]
pl = pl[:60] + pl[60::max(1, len(pl[60:])//60)]   # all small primes, then a thinned sample up to the data cutoff
shifts = sorted(set(grid + pl))
rows = []
Wih = Vw @ np.diag(ew**-0.5) @ Vw.T if ew[0] > 0 else None
for s in shifts:
    D = hc.Bs(mo, l, M, s); Bsy = 0.5*(D['B'] + D['B'].T)
    mp_ = np.linalg.eigvalsh(W + Bsy)[0]; mm = np.linalg.eigvalsh(W - Bsy)[0]
    rho = float(np.abs(np.linalg.eigvalsh(Wih @ Bsy @ Wih)).max()) if Wih is not None else None
    dec = {}
    for tag, g in (('ground', g0), ('sine1', g1)):
        sym = lambda A: 0.5*(A + A.T)
        dec[tag] = dict(d=float(g@W@g), b=float(g@Bsy@g), pole=float(g@sym(D['pole'])@g), prime=float(g@sym(D['prime'])@g), arch=float(g@sym(D['arch'])@g))
    isp = any(abs(s - x) < 1e-12 for x in pl)
    rows.append(dict(s=float(s), prime_shift=isp, p=(round(math.exp(s)) if isp else None), mu_plus=float(mp_), mu_minus=float(mm), rho=rho, fixed=dec))
    if isp or abs(s*50 - round(s*50)) < 1e-9 and round(s*50) % 25 == 0:
        print(name, round(s, 4), 'p' if isp else ' ', '%.3e %.3e rho=%s' % (mp_, mm, rho and round(rho, 4)), flush=True)
neg = [r for r in rows if min(r['mu_plus'], r['mu_minus']) < -1e-7*np.abs(W).max()]
out = dict(model=name, l=l, M=M, fmax=fmax, W_lam=[float(x) for x in ew[:4]], n_shifts=len(rows),
           first_violation_s=(neg[0]['s'] if neg else None), n_violations=len(neg),
           min_mu=float(min(min(r['mu_plus'], r['mu_minus']) for r in rows)),
           max_rho=(float(max(r['rho'] for r in rows)) if Wih is not None else None),
           prime_shifts_violating=[r['p'] for r in neg if r['prime_shift']], rows=rows)
json.dump(out, open(f'blocks_{name}_l{l}.json', 'w'), indent=1)
print('DONE', name, {k: v for k, v in out.items() if k != 'rows'}, flush=True)
