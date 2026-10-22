"""CERTIFIED zero boxes for D=-23 and D=-56 (class groups C3 and C4; no genus-character decomposition).
The argument-principle driver is certify_L.count, unchanged. Only the evaluator is swapped for eps.Z_Q: the Chowla-Selberg / K-Bessel
expansion of E(tau_Q, s) in Arb, with a rigorous tail (see the eps.py docstring), valid for any positive-definite form.
Each numerically located zero (scan_T100.json) gets a box of half-side r. Count = 1 certifies exactly one zero (with multiplicity) in it.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
import sys, json, time
sys.path.insert(0,'.')
import certify_L, eps
from flint import acb
certify_L.ZQ = lambda Q, s: eps.Z_Q(tuple(Q), s)
from multiprocessing import Pool
FORMS={(1,1,6),(2,1,3),(1,0,14),(2,0,7),(3,2,5)}
def job(a):
    Q,box=a; t=time.time(); n,k,ns=certify_L.count(Q,*box)
    return dict(Q=list(Q),box=list(box),count=n,count_ball=str(k),ball_evals=ns,sec=round(time.time()-t,1))
if __name__=='__main__':
    r=float(sys.argv[1]) if len(sys.argv)>1 else 1e-6; W=int(sys.argv[2]) if len(sys.argv)>2 else 3
    zs=[(tuple(b['Q']),z) for b in json.load(open('scan_T100.json')) if tuple(b['Q']) in FORMS for z in b['zeros']]
    if len(sys.argv)>3: zs=zs[:int(sys.argv[3])]
    tasks=[(Q,(z[0]-r,z[0]+r,z[1]-r,z[1]+r)) for Q,z in zs]
    print(len(tasks),'zero boxes',flush=True)
    with Pool(W) as p:
        res=[]
        for x in p.imap_unordered(job,tasks):
            res.append(x); print(x,flush=True)
    json.dump(res,open('certified_zero_boxes_D23_D56.json','w'),indent=1)
