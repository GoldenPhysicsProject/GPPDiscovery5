import json, sys
def load(f):
    try: return [json.loads(l) for l in open(f)]
    except FileNotFoundError: return []
out={}
for M in (24,40):
    rows=load(f'scan_M{M}.jsonl')
    if not rows: continue
    keys=list(rows[0]['configs'].keys()); tab={}
    for k in keys:
        on=None; cert=None; prev=None; sec=None
        for r in rows:
            c=r['configs'][k]; neg=[s for s in ('even','odd') if float(c[s])<0]
            if neg and on is None:
                on=(prev,r['h']); sec=neg
            for s in ('even','odd'):
                if c.get(s+'_certified_neg') and cert is None: cert=(r['h'],s,c[s+'_cert'])
            prev=r['h']
        lam13=[r['configs'][k] for r in rows if r['h']=='1.3']
        tab[k]={"onset_bracket":on,"binding":sec,"first_cert":cert,"lam_at_1.3":({s:lam13[0][s] for s in ('even','odd')} if lam13 else None),"h_max_scanned":rows[-1]['h']}
    out[f"M{M}"]=tab
json.dump(out,open('s46_results.json','w'),indent=1)
for M,tab in out.items():
    print("==",M)
    for k,v in tab.items(): print(f"{k:20s} onset {v['onset_bracket']} {v['binding']} cert {v['first_cert'][:2] if v['first_cert'] else None} hmax {v['h_max_scanned']}")
