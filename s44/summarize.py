import json,sys
rows={}
for M in (24,40,56,72):
    try:
        for l in open(f"scan_M{M}.jsonl"): r=json.loads(l); rows[(M,r['h'])]=r
    except FileNotFoundError: pass
hs=sorted({h for (_,h) in rows},key=float)
print("| L=h | n entered | M | c_- − c_true (odd/even) | c_+ − c_true | width | lower-edge sector | upper-edge sector | λmin(Wζ) even/odd | tower a: [a_-−1, a_+−1] |")
print("|---|---|---|---|---|---|---|---|---|---|")
for h in hs:
    for M in (24,40,56,72):
        r=rows.get((M,h))
        if not r: continue
        e,o=r['even'],r['odd']
        ups={'even':float(e['gap_up']) if e['gap_up']!='inf' else 1e99,'odd':float(o['gap_up']) if o['gap_up']!='inf' else 1e99}
        dns={'even':float(e['gap_down']) if e['gap_down']!='inf' else 1e99,'odd':float(o['gap_down']) if o['gap_down']!='inf' else 1e99}
        us=min(ups,key=ups.get); ds=min(dns,key=dns.get); up,dn=ups[us],dns[ds]
        def tw(x): return x if isinstance(x,str) else ''
        ta=[e.get('tower_2k',{'a_up':'-','a_down':'-','k_in':'-'}),o.get('tower_2k',{'a_up':'-','a_down':'-','k_in':'-'})]
        def num(s):
            try: return float(s)
            except: return None
        aup=[num(t['a_up']) for t in ta if num(t['a_up']) is not None]; adn=[num(t['a_down']) for t in ta if num(t['a_down']) is not None]
        tws=f"[{max(adn):.3e}, +{min(aup):.3e}] k={ta[0]['k_in']}" if aup and adn else ("not run" if ta[0]["a_up"]=="-" else "beyond ±0.5")
        print(f"| {h} | {r['prime_powers']} | {M} | −{dn:.4e} | +{up:.4e} | {up+dn:.4e} | {ds} | {us} | {e['lam_min_Wzeta']}/{o['lam_min_Wzeta']} | {tws} |")
