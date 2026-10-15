"""Grok S27 H27: invisible-harmonic lift of the one-stroke envelope.
For supp f in [-L,L], cos(u t) with u >= 2L is invisible in Q. Per prime p, the visible block
B_p(th) = sum_{k log p < 2L} c_{p,k} cos(k th), c_{p,k} = 2 log p / p^{k/2}; add invisible
harmonics h_p(th) = sum_{K_p <= k <= Kmax} a_k cos(k th) (k log p >= 2L) and minimise
kappa_p = max_th (B_p - h_p). Then Psi + sum_p h_p(t log p) >= log(t/2pi) - 1/t - sum_p kappa_p.
Dual (Fejer-Riesz): kappa_p >= max{ sum c_k r_k : 1 + 2 sum_{k<K_p} r_k cos k th >= 0 }.
Float LP design; rigorous kappa bound done in Arb separately.  RH OPEN."""
import numpy as np, math, json, sys
from scipy.optimize import linprog

def blocks(L):
    out = {}
    for p in [q for q in range(2, int(math.exp(2*L))+2) if all(q % r for r in range(2, int(q**.5)+1))]:
        if math.log(p) >= 2*L: continue
        ks = [k for k in range(1, 60) if k*math.log(p) < 2*L]
        out[p] = (ks, [2*math.log(p)/p**(k/2) for k in ks])
    return out

def lift_lp(ks, cs, Kmax, ngrid=4000, amax=None):
    K = max(ks) + 1
    inv = list(range(K, Kmax + 1))
    th = np.linspace(0, math.pi, ngrid)
    B = sum(c*np.cos(k*th) for k, c in zip(ks, cs))
    # variables: a_k (inv), s ; minimise s s.t. B - sum a cos <= s
    nA = len(inv)
    Aub = np.hstack([-np.cos(np.outer(th, inv)), -np.ones((ngrid, 1))])
    bub = -B
    c = np.r_[np.zeros(nA), 1.0]
    bnds = [(-amax, amax) if amax else (None, None)]*nA + [(None, None)]
    r = linprog(c, A_ub=Aub, b_ub=bub, bounds=bnds, method='highs')
    return inv, r.x[:nA], r.x[-1]

def dual_bound(ks, cs):
    # max sum c_k r_k s.t. 1 + 2 sum r_k cos k th >= 0 (grid LP; an upper estimate of the dual)
    K = max(ks)
    th = np.linspace(0, math.pi, 20001)
    Aub = -2*np.cos(np.outer(th, ks)); bub = np.ones(len(th))
    r = linprog(-np.array(cs), A_ub=Aub, b_ub=bub, bounds=[(-1, 1)]*len(ks), method='highs')
    return -r.fun

if __name__ == '__main__':
    L = float(eval(sys.argv[1])); Kmax_mult = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    bl = blocks(L); res = {}
    AL = sum(sum(cs) for ks, cs in bl.values())
    tot = 0.0
    for p, (ks, cs) in bl.items():
        Kmax = (max(ks)+1)*Kmax_mult
        inv, a, kap = lift_lp(ks, cs, Kmax)
        db = dual_bound(ks, cs)
        res[p] = {'ks': ks, 'cs': cs, 'inv': inv, 'a': list(a), 'kappa': kap, 'dual': db, 'raw': sum(cs)}
        tot += kap
        print(f'p={p} visible k={ks} raw={sum(cs):.5f} kappa(Kmax={Kmax})={kap:.5f} dual={db:.5f} max|a|={max(abs(a)):.3f}')
    print(f'L={L} A_L={AL:.5f} A_lift={tot:.5f}  T1=2pi e^A: {2*math.pi*math.exp(AL):.1f} -> {2*math.pi*math.exp(tot):.1f}')
    json.dump({'L': L, 'A_L': AL, 'A_lift': tot, 'primes': res}, open(f'lift_L{L}.json', 'w'), indent=1)
