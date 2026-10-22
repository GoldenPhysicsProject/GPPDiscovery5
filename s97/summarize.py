"""Grok S97 summary of chain runs -> summary.txt. RH OPEN."""
import json, math, mpmath as mp
mp.mp.dps = 30
out = []
for h in (0.02, 0.01):
    R = {r["name"]: r for r in json.load(open(f"chain2_h{h}.json"))}
    out.append(f"== h={h}, T<=3, dps 60 ==")
    for nm, r in R.items():
        k = [mp.mpf(x) for x in r["kappa"]]
        out.append(f"{nm:12s} PD(T<=3)={r['chol_full']} T_fail={r['T_fail']} max|kappa|={mp.nstr(max(abs(x) for x in k),6)} "
                   f"chol_before={r.get('chol_before_fail')} chol_at={r.get('chol_at_fail')}")
    for p in (2, 3, 5, 7):
        r = R[f"ZETA-drop{p}"]; out.append(f"  drop{p}: T_fail - log p = {r['T_fail'] - math.log(p):.3f}")
    kz = [mp.mpf(x) for x in R["ZETA"]["kappa"]]; kf = [mp.mpf(x) for x in R["FT0.1"]["kappa"]]; kd = [mp.mpf(x) for x in R["DH"]["kappa"]]
    first = lambda a, b: next(((i + 1) * h for i in range(len(a)) if (a[i] > 0) != (b[i] > 0)), None)
    out.append(f"  first kappa sign mismatch: zeta vs FT0.1 at T={first(kz, kf)}, zeta vs DH at T={first(kz, kd)}")
    d2 = lambda k: [k[i + 1] - 2 * k[i] + k[i - 1] for i in range(1, len(k) - 1)]
    for nm, k in (("ZETA", kz), ("DH", kd), ("FT0.1", kf)):
        d = d2(k); top = sorted(range(len(d)), key=lambda i: -abs(d[i]))[:10]
        out.append(f"  {nm} top-10 |d2 kappa| lags: {sorted(round((i + 1) * h, 2) for i in top)}")
out.append("log n (prime powers): " + str([round(math.log(n), 3) for n in (2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19)]))
open("summary.txt", "w").write("\n".join(out) + "\n"); print("\n".join(out))
