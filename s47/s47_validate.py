"""Validate Arb W_DH (s47_dh) entrywise against S45 float form_dh (independent arch via digamma quadrature)."""
import sys, numpy as np
src=open('../s45/s45_dh_float.py').read().split('rows = []')[0]
sys.argv=['x','12','1.0']; g={}; exec(src,g)
sys.argv=['x']; sys.path.insert(0,'.')
from s47_dh import *
for h in ('0.7','1.0','1.3'):
    sys.argv=['x','12',h]
    Qf,Gf=g['form_dh'](float(h),12)
    W,Wz,nn,b,kap=W_DH(12,h,512)
    Wa=np.array([[float(W[i,j].mid()) if (i+j)%2==0 else 0.0 for j in range(12)] for i in range(12)])
    print(h,'max|diff|',np.max(np.abs(Wa-Qf)),'max|Q|',np.max(np.abs(Qf)),'b2',b[2].str(10),'b3',b[3].str(10))
