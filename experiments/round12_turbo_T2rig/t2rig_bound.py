"""(T2-rig) explicit constants, evaluated in Arb (python-flint); plus float checks of the lemmas and the phi-row identity."""
import numpy as np, json
from flint import arb, acb, ctx
from scipy import special
ctx.prec=200
l=arb(1)/4; pi=arb.pi(); NS=[(2,2),(3,3),(4,2),(5,5),(7,7),(8,2),(9,3),(11,11)]
c=[arb(3)/16,arb(-1)/16,arb(-1)/4,arb(3)/32,arb(1)/16,arb(-1)/32]   # P = sum c_k cos(k t/2), exact
Cn={n:sum((c[k] for k in range(6) if k/2<np.log(n)),arb(0)) for n,_ in NS}
kap=arb(5)/4*9/(4*arb(3).sqrt())                    # |A'(t)-1/t| <= kap/t^2
A0=acb(arb(1)/4).digamma().real                      # A(0)=psi(1/4)
def Aup(w): return w.log()-arb(2).log()+kap/w          # A(w) <= log(w/2)+kap/w
def G0(w):  return (Aup(w)-A0)/w                       # G(0)=(A(w)-A(0))/w
def Rfull(w): return pi**2/4+G0(w)+kap/w*(1+w).log()
S=sum((2*abs(c[k])/k for k in range(1,6)),arb(0))
def Rup(w): return c[0]*Rfull(w)+G0(w)*S
Pib=pi*sum((arb(p).log()/arb(n).sqrt()*abs(Cn[n]) for n,p in NS),arb(0))
pref=16/pi*(2/l).sqrt()
bI2=16*(l.sinh())**2; sigc=(2/l).sqrt()*(l/2).sinh()       # |sig(j)|<=sigc/omega
om=lambda j: arb(j)*pi/l
Qd=Rup(om(13))+Pib; Qp=Rup(om(1))+Pib
Cdelta=pref*Qd+bI2**2/8*(l/2).sinh()*sigc/om(13)
Cphi_num=16/pi*(2/l)*Qp   # |w^4 rho_j(phi_i)| <= Cphi_num*(w_i+w_j)/(w_j^2-w_i^2) + sigma term
out={k:str(v) for k,v in dict(kappa=kap,A0=A0,S=S,Pi_bound=Pib,pref=pref,Rup_om13=Rup(om(13)),Rup_4pi=Rup(om(1)),
     Qbound_j13=Qd,Qbound_all=Qp,C_delta=Cdelta,C_delta_over_log_om13=Cdelta/om(13).log(),Cphi_num=Cphi_num,
     phi_row_j13=Cphi_num/(om(13)-om(12)),phi_row_j73=Cphi_num/(om(73)-om(12))).items()}
for k,v in out.items(): print(k,v)
# float lemma checks
A=lambda t: special.digamma(0.25+0.5j*t).real
t=np.logspace(-2,5,20000); h=1e-5*t
Ad=(A(t+h)-A(t-h))/(2*h); out['sup_t2|A1-1/t|']=float(np.max(t*t*np.abs(Ad-1/t)))
s=np.linspace(1e-4,1e6,2000001); B=A(np.sqrt(s)); out['max_second_diff_B(should<=0)']=float(np.max(np.diff(B,2)))
print('sup t^2|A1-1/t|',out['sup_t2|A1-1/t|'],'max 2nd diff',out['max_second_diff_B(should<=0)'])
json.dump(out,open('t2rig_bound_results.json','w'),indent=1)
