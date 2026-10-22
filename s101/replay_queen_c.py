"""Grok S101: independent replay of Queen step (c) values on the S87 forms. HIGH-PRECISION NUMERICAL. RH OPEN."""
import sys, mpmath as mp
sys.path.insert(0, "s87"); import hp_forms as H
for (L,N) in [(0.8,32),(1.0,40)]:
    H.setup(80); H._GN.clear(); d=H.build('ZETA', mp.mpf(L), N)
    Re=mp.matrix(H.add(d["Ae"],d["Pe"],1,-1)); Ro=mp.matrix(H.add(d["Ao"],d["Po"],1,-1))
    E=sorted(mp.eigsy(Re,eigvals_only=True)); Eo=sorted(mp.eigsy(Ro,eigvals_only=True))
    print(L,N,mp.nstr(E[0],10),mp.nstr(E[1],10),mp.nstr(Eo[0],10),flush=True)
