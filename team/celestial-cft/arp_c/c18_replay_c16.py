"""C.18 replay of the C.16 odd p=2 certificate at 2L = log 3 (T = 40, N = 96): recompute ALL entries with the independent
panel-Gauss engine c18_entries.py, compare with the acb.integral entries c16_odd_arb_entries.json (S104 run, 257807 CPU-s),
and redo the Thm C.16.A certificate (interval Cholesky + lambda_cert) from the NEW entries.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok. RH is OPEN."""
import json, time
from c18_entries import *
from c16_odd_arb_run import chol_pd
t0 = time.time(); T, N = 40, 96
E = Engine(arb(3).log(), [(2, arb(2).log())], 'odd', T, N)
NOm = E.gram('Om', range(N), range(N)); B = {2: E.gram(2, range(N), range(N))}
A, M = E.assemble(NOm, B); eA, eN = E.eps(N)
print("engine: T=%d N=%d panels h=1, ng=40, b=0.4, prec %d; time %.0fs" % (T, N, ctx.prec, time.time() - t0))
d = json.load(open('c16_odd_arb_entries.json'))
mx = 0.0; ok = True; mr = 0.0
for X, key in ((A, 'A'), (M, 'M')):
    for i in range(N):
        for j in range(N):
            y = arb(d[key][i][j]); x = X[i][j]
            ok &= bool(x.overlaps(y)); mx = max(mx, float(abs(x.mid() - y.mid()))); mr = max(mr, float(x.rad()))
print("all 2x96^2 entries overlap the c16 JSON balls:", ok, "; max |mid diff| = %.2e ; max new radius = %.2e" % (mx, mr))
for key, x in (('abs', eA), ('neg', eN)):
    y = arb(d['eps'][key]); print("eps^%s new %s  json %s  overlap %s" % (key, x.str(15, radius=True), y.str(15, radius=True), x.overlaps(y)))
print("q_T =", E.QT.str(20, radius=True))
mu = arb('0.0006397'); s = arb('2.914')
okA, pA = chol_pd([[A[i][j] - (mu if i == j else 0) for j in range(N)] for i in range(N)])
okM, pM = chol_pd([[(s if i == j else 0) - M[i][j] for j in range(N)] for i in range(N)])
print("interval Cholesky A - 0.0006397 I > 0:", okA, "; 2.914 I - M > 0:", okM)
c = E.QT - eN; lam = (mu + c - ((mu - c)**2 + 4*s*eA).sqrt())/2
print("lambda_cert (new entries) =", lam.str(15, radius=True), " certified > 0:", bool(lam > 0) and okA and okM)
