"""Zhu's empirical Landau-Widom law  -ln lambda*(L) ~ C N(T*)/ln N(T*),  T* = 2 pi e^{2L}  (arXiv:2608.24827 eq.(5), Conj. 12.1)
checked against (a) Zhu's own certified upper bounds (Table 3), (b) the TEAM's certified two-sided brackets
(GPPDiscovery5 s37/s53/s62/s63/s72-s74/s84), (c) fresh Legendre runs (runs/leg_*.jsonl, this folder).
Everything here is NUMERICAL analysis of certified/numerical inputs; the law itself is a CONJECTURE (fitted constant)."""
import json, math, glob
import mpmath as mp
mp.mp.dps = 30
def Tstar(L): return 2*math.pi*math.exp(2*L)
def Nexact(T): return int(mp.nzeros(T))
def Nsmooth(T): return T/(2*math.pi)*math.log(T/(2*math.pi*math.e)) + 7/8
def lawval(C, N): return C*N/math.log(N)
zhu_upper = {0.7:4.99e-13,0.8:2.27e-17,0.9:5.91e-23,1.0:8.08e-30,1.1:2.78e-38,1.2:9.98e-49,1.4:1.91e-77,1.6:8.61e-121,1.8:7.86e-186,2.0:3.19e-283}
zhu_ref_even = {0.5:9.34e-7,0.6:1.61e-9,0.7:4.18e-13,0.8:1.65e-17,0.9:4.14e-23,1.0:5.88e-30,1.1:2.04e-38,1.2:7.94e-49}   # Table 1, converged (1.2: within ~2x)
zhu_ref_odd  = {0.5:1.94e-4,0.6:5.97e-7,0.7:2.56e-10,0.8:1.57e-14,0.9:7.22e-20,1.0:1.49e-26,1.1:7.68e-35,1.2:4.76e-45}
# team brackets (CERTIFIED lower, CERTIFIED upper); None = not available
team_even = {1.30:(3.634e-62,None),1.35:(5.26499e-70,2.7691e-69),1.40:(7.4196e-78,1.6486e-77),1.42:(1.0711e-81,4.44455e-81),1.7125:(None,2.62874e-149)}
team_odd  = {1.30:(1.4910e-58,None),1.35:(1.5659e-65,3.3829e-65),1.40:(1.2459e-73,2.4959e-73),1.42:(2.4134e-77,7.13734e-77)}
out = {"rows": []}
print("L     T*      N  Nsm   | source               lambda            -ln     C_eff(exactN)  C_eff(smoothN)")
def row(L, lab, lam):
    T = Tstar(L); N = Nexact(T); Ns = Nsmooth(T); y = -math.log(lam)
    ce = y*math.log(N)/N if N > 1 else float('nan'); cs = y*math.log(Ns)/Ns
    print("%-5g %7.2f %3d %6.2f | %-20s %-16.4e %7.2f  %7.3f  %7.3f" % (L, T, N, Ns, lab, lam, y, ce, cs))
    out["rows"].append({"L": L, "Tstar": T, "N": N, "Nsmooth": Ns, "source": lab, "lambda": lam, "minus_ln": y, "C_eff_exact": ce, "C_eff_smooth": cs})
for L, v in sorted(zhu_ref_even.items()):
    if L >= 0.7: row(L, "Zhu T1 even (ref)", v)
for L, v in sorted(zhu_upper.items()):
    if L >= 1.2: row(L, "Zhu T3 even (upper)", v)
for L, (lo, hi) in sorted(team_even.items()):
    if lo: row(L, "team even LOWER", lo)
    if hi: row(L, "team even UPPER", hi)
for f in sorted(glob.glob("runs/leg_*.jsonl")):
    best = None
    for l in open(f):
        d = json.loads(l); best = d
    if best:
        L = float(eval(best["L"].replace("/", "/1.0/") if False else str(float(eval(best["L"])))))
        row(L, "Legendre N=%d even" % best["N_per_parity"], float(best["even"]["lam1"]))
# odd/even ratio law (Zhu: 10^{2.07L+1.32} for L<=1.2)
print("\nodd/even ratio: Zhu law 10^(2.07L+1.32) vs team brackets (lower_odd/upper_even .. upper_odd/lower_even)")
for L in (1.35, 1.40, 1.42):
    lo_o, hi_o = team_odd[L]; lo_e, hi_e = team_even[L]
    print("  L=%.2f  law %.3g   team ratio in [%.3g, %.3g]" % (L, 10**(2.07*L+1.32), lo_o/hi_e, hi_o/lo_e))
# predictions
print("\nPredictions at L = 1.7105 and 1.71:")
pred = {}
for L in (1.71, 1.7105):
    T = Tstar(L); N = Nexact(T); Ns = Nsmooth(T)
    for C in (19.74, 20.13, 20.6):
        for lab, NN in (("exactN", N), ("smoothN", Ns)):
            y = lawval(C, NN); pred[(L, C, lab)] = y
            print("  L=%.4f T*=%.2f N=%d Nsm=%.2f  C=%.2f %-8s  -ln lam = %.1f  -> lam ~ 1e%.1f" % (L, T, N, Ns, C, lab, y, -y/math.log(10)))
# log-linear interpolation of Zhu's certified upper bounds in ln(-ln lam)
a, b = math.log(-math.log(zhu_upper[1.6])), math.log(-math.log(zhu_upper[1.8]))
for L in (1.71, 1.7105):
    y = math.exp(a + (b - a)*(L - 1.6)/0.2)
    print("  interpolated Zhu upper-bound curve at L=%.4f: -ln = %.1f -> 1e%.1f" % (L, y, -y/math.log(10)))
    out["interp_zhu_upper_%s" % L] = y
json.dump(out, open("zhu_law_results.json", "w"), indent=1)
