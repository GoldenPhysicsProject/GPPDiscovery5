"""Grok S96: tables + doubling comparisons from h95a.json -> summary.txt. NUMERICAL. RH OPEN."""
import json, mpmath as mp
mp.mp.dps = 90
R = {r["L"]: r for r in json.load(open("h95a.json"))}
Ls = sorted(R); FORMS = ("ZETA", "DH", "FT0.1", "HT-2", "HT-4")
g = lambda L, F, k: mp.mpf(R[L]["forms"][F]["_num"][k]) * mp.mpf(R[L]["_lmax"]) if k in R[L]["forms"][F]["_num"] else None
s = lambda x, k=3: "—" if x is None else mp.nstr(x, k)
out = []
out.append("Per-L table (rel = / lambda_max(zeta,L)); signs from mp Cholesky / LDL^T at dps 90")
for F in FORMS:
    out.append(f"\n== {F} ==\nL | q_e rel (n-) | q_o rel (n-) | A rel (n-) | 1+mu_min [Q_e] | 1-mu_max [Q_o] | 1-||M|| | side | cos(xext,vQmin) e/o | cos(xext,vS+) | veAve/lamA")
    for L in Ls:
        f = R[L]["forms"][F]; n = lambda c: 0 if c["PD_chol"] else c["n_neg_ldl"]
        if isinstance(f.get("M"), str):
            out.append(f"{L} | {f['qmin_e_rel']} ({n(f['Qe'])}) | {f['qmin_o_rel']} ({n(f['Qo'])}) | {f['Amin_rel']} ({n(f['A'])}) | A not PD")
        else:
            out.append(f"{L} | {f['qmin_e_rel']} ({n(f['Qe'])}) | {f['qmin_o_rel']} ({n(f['Qo'])}) | {f['Amin_rel']} ({n(f['A'])}) | {f['1+mu_min (Q_e)']} | {f['1-mu_max (Q_o)']} | {f['1-||M||']} | {f['extremal_side'][:3]} | {f['cos(x_ext_e, v_Qe_min)']}/{f['cos(x_ext_o, v_Qo_min)']} | {f['cos(x_ext, v_Splus_min)']} | {f['ve\'A ve / lamminA']}")
out.append("\nDoubling: A(L) vs full form at L/2.  qhalf = min(q_e,q_o)(L/2)/2 (abs; e_k has norm^2 2).  n-(A(L)) vs n-(Q_e)+n-(Q_o) at L/2.")
out.append("Also law tests: r1 = (1-||M(L)||) / (qmin(L)/2 / lamA(L));  r2 = log(1-||M(L)||)/log(qhalf/lamA...)")
for F in FORMS:
    out.append(f"\n== {F} ==\nL | lamA(L) abs | qhalf(L/2) abs | ratio lamA/qhalf | n-A(L) | n-Q(L/2) | 1-||M(L)|| | qmin(L)/2/lamA(L) | (1-||M||)/that | log(1-||M||(L))/log(1-||M||(L/2))")
    for L in Ls:
        h = round(L / 2, 3)
        if h not in R: continue
        f = R[L]["forms"][F]; fh = R[h]["forms"][F]
        a = g(L, F, "a"); qh = min(g(h, F, "qe"), g(h, F, "qo")) / 2; qL = min(g(L, F, "qe"), g(L, F, "qo"))
        nA = 0 if f["A"]["PD_chol"] else f["A"]["n_neg_ldl"]; nh = sum(0 if fh[c]["PD_chol"] else fh[c]["n_neg_ldl"] for c in ("Qe", "Qo"))
        om = g(L, F, "onem"); om = om / mp.mpf(R[L]["_lmax"]) if om is not None else None
        omh = g(h, F, "onem"); omh = omh / mp.mpf(R[h]["_lmax"]) if omh is not None else None
        rr = (qL / 2 / a) if a > 0 else None
        lg = (mp.log(om) / mp.log(omh)) if (om and omh and om > 0 and omh > 0 and omh < 1) else None
        out.append(f"{L} | {s(a)} | {s(qh)} | {s(a/qh) if qh else '—'} | {nA} | {nh} | {s(om)} | {s(rr)} | {s(om/rr) if (om is not None and rr) else '—'} | {s(lg)}")
open("summary.txt", "w").write("\n".join(out)); print("\n".join(out))
