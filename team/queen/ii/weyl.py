import json,sys
for p in sys.argv[1:]:
    D=json.load(open(p)); W=D["W"]; K=D["K"]; from fractions import Fraction
    L=float(Fraction(D["L"]))
    print(p, {k: k*max(W[i][j][1] for i in range(k) for j in range(k))/L for k in (64,128,200,266,320,400) if k<=K})
