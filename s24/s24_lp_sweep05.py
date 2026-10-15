"""S24b [design, float LP]: does the S22 claim-1 LP reach the Krein ceiling at h=1/2 when B (Bernstein cap), U, Dd are relaxed?"""
import sys, os, time; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "s22"))
from s22_lp_odd import lp
import numpy as np
h = float(sys.argv[1])
for (U, D, Dd, B) in [(12,0.05,0.5,3),(12,0.05,0.5,6),(12,0.05,1.0,3),(12,0.05,0.25,3),(20,0.05,0.5,3),(20,0.05,0.5,6),(12,0.025,0.5,6)]:
    t0=time.time(); r,u,Xi=lp(h,U=U,D=D,Dd=Dd,B=B)
    print(f"h={h} U={U} D={D} Dd={Dd} B={B}: status={r.status} c={-r.fun if r.status==0 else None} Xi={Xi:.0f} ({time.time()-t0:.0f}s)",flush=True)
