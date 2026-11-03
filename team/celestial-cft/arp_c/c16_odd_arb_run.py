"""C.16 odd (K2) Arb run at 2L = log 3, T = 40, N = 96, analytic majorants (Thm C.16.A + Prop C.15.B => odd (K2) on (log 2, log 3]).
All inputs enclosed: acb.integral of meromorphic integrands (c16_arb_callbacks), interval Cholesky of A - mu I and s I - M in arb.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok. RH is OPEN."""
import sys, json, time
from multiprocessing import Pool
from flint import arb, acb, ctx
ctx.prec = 128
from c16_arb_callbacks import Callbacks
TWO_L = 'log3'; T = 40; N = 96
def mk(): return Callbacks(arb(3).log(), [(2, arb(2).log()/arb(2).sqrt())], 'odd', T, N)
CB = None
def job(jk):
    global CB
    if CB is None: CB = mk()
    j, k = jk; ctx.prec = 128
    return (j, k, CB.A_entry(j, k).str(40, radius=True), CB.M_entry(j, k).str(40, radius=True))
def epsjob(w):
    ctx.prec = 128; return (w, mk().eps(w).str(40, radius=True))
def chol_pd(Mx):
    """interval Cholesky: returns (True, min pivot lower bound) if every pivot ball is > 0 (proves Mx > 0)."""
    n = len(Mx); Lm = [[arb(0)]*n for _ in range(n)]; minp = None
    for i in range(n):
        for j in range(i + 1):
            s = Mx[i][j] - sum((Lm[i][k]*Lm[j][k] for k in range(j)), arb(0))
            if i == j:
                if not s > 0: return False, s
                minp = s if minp is None or s.mid() < minp.mid() else minp
                Lm[i][i] = s.sqrt()
            else: Lm[i][j] = s/Lm[j][j]
    return True, minp
if __name__ == '__main__':
    t0 = time.time(); pairs = [(j, k) for j in range(N) for k in range(j, N)]
    with Pool(7) as P:
        er = P.map_async(epsjob, ['abs', 'neg'])
        res = P.map(job, pairs, chunksize=8); eps = dict(er.get())
    A = [[None]*N for _ in range(N)]; M = [[None]*N for _ in range(N)]
    for j, k, a, m in res: A[j][k] = A[k][j] = arb(a); M[j][k] = M[k][j] = arb(m)
    json.dump({'T': T, 'N': N, 'twoL': TWO_L, 'A': [[A[i][j].str(30, radius=True) for j in range(N)] for i in range(N)],
               'M': [[M[i][j].str(30, radius=True) for j in range(N)] for i in range(N)], 'eps': eps}, open('c16_odd_arb_entries.json', 'w'))
    cb = mk(); qT = cb.QT; eA = arb(eps['abs']); eN = arb(eps['neg'])
    print(f"2L = log 3, odd, T = {T}, N = {N}, prec 128; entries: {time.time()-t0:.0f}s")
    print("max entry radius A:", max(float(A[i][j].rad()) for i in range(N) for j in range(N)), " M:", max(float(M[i][j].rad()) for i in range(N) for j in range(N)))
    print("q_T =", qT.str(20, radius=True)); print("eps^abs =", eA.str(15, radius=True), " eps^neg =", eN.str(15, radius=True))
    import numpy as np
    Af = np.array([[float(A[i][j].mid()) for j in range(N)] for i in range(N)]); Mf = np.array([[float(M[i][j].mid()) for j in range(N)] for i in range(N)])
    muf = np.linalg.eigvalsh(Af)[0]; sf = np.linalg.eigvalsh(Mf)[-1]; print(f"float mid: mu_A = {muf:.6e}, s_V = {sf:.6f}")
    mu = arb(str(float('%.4g' % (0.995*muf)))); s = arb(str(float('%.4g' % (1.002*sf))))
    okA, pA = chol_pd([[A[i][j] - (mu if i == j else 0) for j in range(N)] for i in range(N)])
    okM, pM = chol_pd([[(s if i == j else 0) - M[i][j] for j in range(N)] for i in range(N)])
    print(f"interval Cholesky A - {mu.str(6)} I > 0: {okA} (min pivot {pA.str(6, radius=True)})")
    print(f"interval Cholesky {s.str(6)} I - M > 0: {okM} (min pivot {pM.str(6, radius=True)})")
    c = qT - eN; B2 = s*eA; lam = (mu + c - ((mu - c)**2 + 4*B2).sqrt())/2
    print("lambda_cert (ball) =", lam.str(15, radius=True), " lower endpoint > 0:", bool(lam > 0) and okA and okM)
