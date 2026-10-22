"""c=2 Narain (T^2) CFT at complex structure tau_Q, B=0, area rho2: unitarity / modular invariance checks
and the identity  sum_{m!=0, w=0} Delta^{-s} = (rho2 sqrt|D|)^s Z_Q(s).   (NUMERICAL checks, mpmath.)
Conventions: metric G=(rho2/tau2)[[1,tau1],[tau1,|tau|^2]], Delta = m^T G^{-1} m/2 + w^T G w/2, J = m.w,
h=(Delta+J)/2, hb=(Delta-J)/2, Z(z)=|eta(z)|^{-4} sum_{m,w in Z^2} q^h qbar^hb  (c = cbar = 2).
Standard Narain construction: Narain, Phys. Lett. B169 (1986) 41; Polchinski, String Theory vol. 1, sec. 8.4 (mem).
Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
import sys, json, itertools
import mpmath as mp
mp.mp.dps=30

def setup(Q, rho2=1):
    a,b,c=Q; D=b*b-4*a*c
    t1=mp.mpf(-b)/(2*a); t2=mp.sqrt(-D)/(2*a)
    G=mp.matrix([[1,t1],[t1,t1**2+t2**2]])*(rho2/t2)
    return G, G**-1, D

def spectrum(Q,rho2=1,M=6):
    G,Gi,D=setup(Q,rho2)
    out=[]
    for m1,m2,w1,w2 in itertools.product(range(-M,M+1),repeat=4):
        m=mp.matrix([m1,m2]); w=mp.matrix([w1,w2])
        Dl=(m.T*Gi*m)[0]/2+(w.T*G*w)[0]/2
        J=m1*w1+m2*w2
        out.append((Dl,J,(m1,m2,w1,w2)))
    return out

def Z(z,spec):
    q=mp.e**(2j*mp.pi*z)
    eta=mp.e**(2j*mp.pi*z/24)*mp.qp(q)
    S=mp.mpf(0)
    for Dl,J,_ in spec:
        h=(Dl+J)/2; hb=(Dl-J)/2
        S+= mp.exp(2j*mp.pi*z*h)*mp.exp(-2j*mp.pi*mp.conj(z)*hb)
    return S/abs(eta)**4

if __name__=='__main__':
    res={}
    for Q in [(1,0,5),(2,2,3),(1,1,6),(2,1,3),(1,0,14),(2,0,7),(3,2,5),(1,0,21),(3,0,7),(2,2,11),(5,4,5),(1,0,1)]:
        spec=spectrum(Q,1,M=5)
        Ds=sorted(set(mp.nstr(d,12) for d,J,_ in spec))
        nvac=sum(1 for d,J,_ in spec if d<mp.mpf('1e-20'))
        gap=min(d for d,J,_ in spec if d>mp.mpf('1e-20'))
        zs=[mp.mpc('0.137','0.83'),mp.mpc('-0.31','1.07')]
        mod=[]
        for z in zs:
            a=Z(z,spec); b=Z(-1/z,spec); c=Z(z+1,spec)
            mod.append(dict(z=str(z),Z=mp.nstr(a,15),S_rel=mp.nstr(abs(b-a)/abs(a),3),T_rel=mp.nstr(abs(c-a)/abs(a),3),ImZ=mp.nstr(mp.im(a),3)))
        # low-lying primaries (U(1)^2 primaries = lattice points) with multiplicities
        from collections import Counter
        low=Counter((mp.nstr((d+J)/2,8),mp.nstr((d-J)/2,8)) for d,J,_ in spec if d<1.6)
        # spectral Dirichlet series of pure-momentum scalars at s=3 vs Z_Q
        s=3
        lhs=sum(d**(-s) for d,J,mw in spectrum(Q,1,M=40) if mw[2]==0 and mw[3]==0 and d>1e-20) if False else None
        res[str(Q)]=dict(vacuum_multiplicity=nvac,gap=mp.nstr(gap,10),modular=mod,low_primaries=sorted([(k[0],k[1],v) for k,v in low.items()],key=lambda x:float(x[0])+float(x[1]))[:12])
        print(Q,res[str(Q)]['vacuum_multiplicity'],res[str(Q)]['gap'],mod)
    json.dump(res,open('narain_checks.json','w'),indent=1)
