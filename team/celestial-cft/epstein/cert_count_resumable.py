"""Resumable driver for the CERTIFIED full-region count (certify_L.count), region [0.52, 3] x [1, T] in boxes of height 10.
Each finished box is appended to certified_counts_T{T}.jsonl, so a restart skips done boxes. Zero-free for Re s >= 3: no_zeros_right.py.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
import sys, json, os
from multiprocessing import Pool
from certify_L import job
T=float(sys.argv[1]); W=int(sys.argv[2]); out=f'certified_counts_T{int(T)}.jsonl'
forms=[(1,0,1),(1,0,5),(2,2,3),(1,0,21),(3,0,7),(2,2,11),(5,4,5)]
done=set()
if os.path.exists(out):
    for l in open(out): r=json.loads(l); done.add((tuple(r['Q']),tuple(r['box'])))
tasks=[]
for Q in forms:
    y=1.0
    while y<T:
        b=(0.52,3.0,y,min(y+10,T)); y+=10
        if (Q,b) not in done: tasks.append((Q,b))
print(len(done),'done,',len(tasks),'to go',flush=True)
with Pool(W) as p:
    for r in p.imap_unordered(job,tasks):
        r['Q']=list(r['Q']); r['box']=list(r['box'])
        open(out,'a').write(json.dumps(r)+'\n'); print(r,flush=True)
