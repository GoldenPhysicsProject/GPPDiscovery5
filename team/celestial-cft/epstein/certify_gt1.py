"""CERTIFIED: a zero of Z_{x^2+5y^2} with Re s > 1. Two independent certificates:
(a) certify_L.count, the rigorous argument principle with Arb Hurwitz L-values on a box of half-side r, which lies entirely in Re s > 1;
(b) certify_rouche.certify, Rouche with the eps.py evaluator and the lattice-point M bound.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
import sys, json
sys.path.insert(0,'.')
import certify_L, certify_rouche as CR
x,y=float(sys.argv[1]),float(sys.argv[2]); r=float(sys.argv[3]) if len(sys.argv)>3 else 1e-6
Q=(1,0,5)
n,k,ns=certify_L.count(Q,x-r,x+r,y-r,y+r)
a=dict(method='argument principle (zq_L, Arb)',box=[x-r,x+r,y-r,y+r],count=n,count_ball=str(k),ball_evals=ns)
b=CR.certify(Q,CR.setup(Q),x,y)
out=dict(Q=list(Q),zero=[x,y],argument_principle=a,rouche=b,re_lower_bound=x-r)
print(json.dumps(out,indent=1)); json.dump(out,open(f'certified_gt1_{y:.0f}.json','w'),indent=1)
