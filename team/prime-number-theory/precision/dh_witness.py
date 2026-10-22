"""DH negative direction at L: eigenvalue (inverse iteration with shift) and mode content of the witness vector.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok.  usage: dh_witness.py dump.json K"""
import sys, json
import mpmath as mp
D=json.load(open(sys.argv[1])); K=int(sys.argv[2]); mp.mp.prec=400
A=mp.matrix([[mp.mpf(D["W"][i][j][0]) for j in range(K)] for i in range(K)])
E,Q=mp.eigsy(A)
idx=sorted(range(K),key=lambda i:E[i])[:3]
Lf=float(mp.mpf(eval(D["L"].replace('/','*1.0/'))))
out={"form":D["form"],"L":D["L"],"K":K,"lowest_eigs":[mp.nstr(E[i],6) for i in idx]}
v=[float(Q[j,idx[0]]) for j in range(K)]
w=[x*x for x in v]; tot=sum(w)
top=sorted(range(K),key=lambda j:-w[j])[:8]
out["witness_top_modes"]=[(j,round(w[j]/tot,4),round((2*j+1)*3.141592653589793/(2*Lf),2)) for j in top]
cum=0; 
for j in range(K):
    cum+=w[j]
    if cum>=0.5*tot: out["median_mode"]=j; out["median_freq"]=round((2*j+1)*3.141592653589793/(2*Lf),2); break
print(json.dumps(out))
