import json
cz=json.load(open('certified_zero_boxes.json'))
sc=json.load(open('scan_T100.json'))
try: cb=json.load(open('certified_counts_T201.json'))
except Exception: cb=[]
L=[]
L.append('| form Q | D | zeros found in [0.52,3]×[1,100] | max Re among them | status |')
L.append('|---|---|---|---|---|')
forms={}
for r in sc:
    forms.setdefault(tuple(r['Q']),[]).extend(r['zeros'])
genus={(1,0,1),(1,0,5),(2,2,3),(1,0,21),(3,0,7),(2,2,11),(5,4,5)}
for Q,z in forms.items():
    a,b,c=Q; D=b*b-4*a*c
    mx=max([x for x,y,m in z],default=None)
    st='CERTIFIED (each in a 1e-8 box)' if Q in genus else 'NUMERICAL'
    if not z: st='none found (NUMERICAL scan)' + ('; see certified counts' if Q in genus else '')
    L.append(f'| {Q} | {D} | {len(z)} | {"%.4f"%mx if mx else "-"} | {st} |')
L.append('')
L.append('Certified zero boxes (centre, half-width 1e-8, exactly one zero each):')
L.append('')
L.append('| form | Re | Im |')
L.append('|---|---|---|')
for r in cz:
    x=(r['box'][0]+r['box'][1])/2; y=(r['box'][2]+r['box'][3])/2
    L.append(f'| {tuple(r["Q"])} | {x:.10f} | {y:.10f} |')
if cb:
    L.append('')
    L.append('Certified zero COUNTS in [0.52,3]×[t0,t1] (Arb argument principle on the whole box boundary):')
    L.append('')
    tot={}
    for r in cb:
        Q=tuple(r['Q']); tot.setdefault(Q,[]).append((r['box'][2],r['box'][3],r['count']))
    L.append('| form | total zeros 0.52<Re<3, 1<Im<201 | per-10 counts |')
    L.append('|---|---|---|')
    for Q,v in tot.items():
        v.sort(); L.append(f'| {Q} | {sum(c for _,_,c in v)} | {" ".join(str(c) for _,_,c in v)} |')
open('ZERO_TABLES.md','w').write('\n'.join(L)+'\n')
print('\n'.join(L[:16]))
