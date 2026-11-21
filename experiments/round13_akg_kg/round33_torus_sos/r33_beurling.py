"""Round 33: Beurling (R32 DMV-type system, primes from round32) — exact step integrals of its own Moebius partial
sums. Reads beur_sqf.bin (log g, mu_B(g)) from r33_beur_dump.c. Outputs per L = kH: dg, q, r, M (as r33_pass)."""
import json, math, numpy as np
H = 0.005; NMAX = 1e8
z = np.fromfile('beur_sqf.bin').reshape(-1, 2); o = np.argsort(z[:, 0], kind='stable'); lg = z[o, 0]; s = z[o, 1]; del z, o
a = s*np.exp(-lg/2); M = np.cumsum(a); dl = np.diff(np.append(lg, math.log(NMAX)))
cq = np.cumsum(M*M*dl); cr = np.cumsum(M*dl); cd = np.cumsum(a*a)
L = np.arange(1, int(math.log(NMAX)/H) + 1)*H; i = np.searchsorted(lg, L, side='right') - 1
tail = L - lg[i]
q = cq[i] - M[i]**2*dl[i] + M[i]**2*tail; r = cq[i]*0 + cr[i] - M[i]*dl[i] + M[i]*tail
json.dump(dict(NMAX=NMAX, H=H, L=L.tolist(), dg=cd[i].tolist(), q=q.tolist(), r=r.tolist(), M=M[i].tolist(), count=int(len(lg))),
          open('r33_beurling.json', 'w'))
print('done', len(lg))
