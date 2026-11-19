"""Round 16 model data at extended cutoffs (integer / own-gamma data only; no zeros).
ZETA-family (ZETA, FTH theta=0.2, FTAU tau=3): von Mangoldt from a segmented sieve up to NZ = e^{19.65}.
Fakes (EPS x^2+5y^2, DH, BEUR2 spliced Beurling): the 14B constructions (harness / models14b) rerun with cutoff NF = e^{14.65}.
Cached in data16_*.npz (not committed if large)."""
import sys, os, math, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); R13 = os.path.join(HERE, '..')
sys.path.insert(0, os.path.join(R13, 'round14b_euler')); sys.path.insert(0, R13); sys.path.insert(0, os.path.join(R13, 'round15_hasse'))
LZ, LF = 19.65, 14.65
NZ, NF = int(math.exp(LZ)), int(math.exp(LF))

def primes_upto(N, seg=1 << 24):
    r = int(math.isqrt(N)) + 1; sv = np.ones(r+1, bool); sv[:2] = False
    for p in range(2, int(r**.5)+1):
        if sv[p]: sv[p*p::p] = False
    small = np.nonzero(sv)[0]; out = [small]
    lo = r + 1
    while lo <= N:
        hi = min(lo + seg, N + 1); s = np.ones(hi - lo, bool)
        for p in small:
            st = max(p*p, ((lo + p - 1)//p)*p)
            if st >= hi: continue
            s[st - lo::p] = False
        out.append(np.nonzero(s)[0] + lo); lo = hi
    return np.concatenate(out)

def zeta_lam(N=NZ):
    f = os.path.join(HERE, 'data16_zeta.npz')
    if os.path.exists(f):
        d = np.load(f); return d['logn'], d['lam']
    P = primes_upto(N); lp = np.log(P.astype(float))
    ln, lam = [lp], [lp]
    for k in range(2, 40):
        m = P[P.astype(float)**k <= N]
        if not len(m): break
        ln.append(k*np.log(m.astype(float))); lam.append(np.log(m.astype(float)))
    ln = np.concatenate(ln); lam = np.concatenate(lam); o = np.argsort(ln, kind='stable')
    np.savez(f, logn=ln[o], lam=lam[o]); return ln[o], lam[o]

def fakes():
    f = os.path.join(HERE, 'data16_fakes.npz')
    if os.path.exists(f):
        d = np.load(f); return {k: (d[k+'_logn'], d[k+'_c']) for k in ('EPS', 'DH', 'BEUR2')}
    import harness as h, models14b as md
    h.NMAX = NF; md.NX = NF
    out = {}
    E = h.epstein_model(); out['EPS'] = (np.log(E['n']), E['c'])
    D = md.dh(); out['DH'] = (np.log(D['n']), D['c'])
    B = md.beur2(); out['BEUR2'] = (np.log(B['n']), B['c'])
    np.savez(f, **{k+'_logn': v[0] for k, v in out.items()}, **{k+'_c': v[1] for k, v in out.items()})
    return out

if __name__ == '__main__':
    import time; t = time.time(); ln, lam = zeta_lam(); print('zeta', len(ln), ln[-1], time.time()-t)
    t = time.time(); F = fakes(); print({k: (len(v[0]), float(v[0][-1])) for k, v in F.items()}, time.time()-t)
