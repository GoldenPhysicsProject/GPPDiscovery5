import sys; sys.path.insert(0,'.'); from s14_exact import *; import numpy as np, time
t=time.time(); A,T=build(14,"0.5",256); print(time.time()-t)
An=np.array([[float(A[i,j].mid()) for j in range(14)] for i in range(14)]); Tn=np.array([[float(T[i,j].mid()) for j in range(14)] for i in range(14)])
print(A[0,0],A[13,13].rad())
import os; os.chdir('..'); sys.argv=['x','14','.3','.3','1']
src=open("experiments/weil_small_support/weil_ritz.py").read().split("cases = [")[0]; ns={"__file__":"experiments/weil_small_support/weil_ritz.py"}; exec(src,ns)
Q,Gm=ns["form"](0.5,"Z",0.0); c=np.sqrt(2)*np.log(2)
print(np.abs((An-c*Tn)-Q).max(), np.abs(Q).max())
