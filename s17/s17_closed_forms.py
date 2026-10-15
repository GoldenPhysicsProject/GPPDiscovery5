"""S17 step checks for Codex 2026-10-15 uniform all-Dirichlet small-window Weil gap.
(1) T0/T1 closed forms vs numerically integrated nu0/nu1 tails (mpmath dps 50).
(2) psi(3/4)-psi(1/4)=pi; c0,c1 closed forms.
(3) B(h) [Codex] at h=1/25,1/32 as Arb balls (python-flint, rigorous).
(4) NEW (S17): bathtub/Schur bound E_eps >= T_eps(h)||g||^2 (tail from h, not 2h) => B*(h)=T0(h)-c0-2(sinh h-h).
    Arb balls + roots of B, B*, and the nontrivial-character bounds."""
import mpmath as mp
from flint import arb, ctx
mp.mp.dps = 50
nu = lambda e: (lambda t: 2*mp.e**(-(1+2*e)*t/2)/(1-mp.e**(-2*t)))
T0 = lambda L: 2*(mp.atanh(mp.e**(-L/2))+mp.atan(mp.e**(-L/2)))
T1 = lambda L: 2*(mp.atanh(mp.e**(-L/2))-mp.atan(mp.e**(-L/2)))
print("(1) tail closed forms, L=2h: |quad - closed|")
for h in [mp.mpf(1)/32, mp.mpf(1)/25, mp.mpf('0.1'), mp.mpf('0.2'), mp.mpf('0.3'), mp.log(2)/2-mp.mpf('1e-3')]:
    L = 2*h
    q0 = mp.quad(nu(0), [L, 1, 10, mp.inf]); q1 = mp.quad(nu(1), [L, 1, 10, mp.inf])
    print(f"  h={mp.nstr(h,8):>10}  T0={mp.nstr(T0(L),20)}  err0={mp.nstr(abs(q0-T0(L)),3)}  T1={mp.nstr(T1(L),20)}  err1={mp.nstr(abs(q1-T1(L)),3)}")
print("(2) digamma")
c0 = mp.log(mp.pi)-mp.digamma(mp.mpf(1)/4); c1 = mp.log(mp.pi)-mp.digamma(mp.mpf(3)/4)
print("  psi(3/4)-psi(1/4)-pi =", mp.nstr(mp.digamma(0.75)-mp.digamma(0.25)-mp.pi,5))
print("  c0 =", mp.nstr(c0,25), " c0-(log8pi+gamma+pi/2) =", mp.nstr(c0-(mp.log(8*mp.pi)+mp.euler+mp.pi/2),5))
print("  c1 =", mp.nstr(c1,25), " c1-(log8pi+gamma-pi/2) =", mp.nstr(c1-(mp.log(8*mp.pi)+mp.euler-mp.pi/2),5))
print("  (T1-c1)-(T0-c0)-(pi-4atan r) at h=1/25 :", mp.nstr((T1(0.08)-c1)-(T0(0.08)-c0)-(mp.pi-4*mp.atan(mp.e**-0.04)),5))
# (3),(4) Arb
ctx.prec = 256
def Tarb(L, e):
    r = (-arb(L)/2).exp(); at = ((1+r)/(1-r)).log()/2; an = r.atan()
    return 2*(at+an) if e == 0 else 2*(at-an)
C0 = (8*arb.pi()).log()+arb.const_euler()+arb.pi()/2; C1 = C0-arb.pi()
pole = lambda h: 2*(arb(h).sinh()-arb(h))
B  = lambda h: Tarb(2*arb(h),0)-C0-pole(h)       # Codex B(h): tail from 2h
Bs = lambda h: Tarb(arb(h),0)-C0-pole(h)         # S17 B*(h): tail from h (bathtub+Schur)
print("(3)/(4) Arb balls")
for hq in [(1,32),(1,25),(1,16),(1,12)]:
    h = arb(hq[0])/hq[1]
    print(f"  h={hq[0]}/{hq[1]}: B={B(h).str(12)}  B*={Bs(h).str(12)}  nontriv-even T0(h)-c0={ (Tarb(h,0)-C0).str(12)}  nontriv-odd T1(h)-c1={(Tarb(h,1)-C1).str(12)}")
print("  roots (mpmath):")
Bm  = lambda h: T0(2*h)-c0-2*(mp.sinh(h)-h)
Bsm = lambda h: T0(h)-c0-2*(mp.sinh(h)-h)
for nm,f,x0 in [("B=0.27",lambda h:Bm(h)-0.27,0.03),("B=0.01",lambda h:Bm(h)-0.01,0.04),("B=0",Bm,0.045),
                ("B*=0.27",lambda h:Bsm(h)-0.27,0.06),("B*=0.01",lambda h:Bsm(h)-0.01,0.08),("B*=0",Bsm,0.09),
                ("T0(h)-c0=0 (nontrivial even, f=1 formal)",lambda h:T0(h)-c0,0.09),
                ("T1(h)-c1=0 (odd, f=1 formal)",lambda h:T1(h)-c1,0.2)]:
    print(f"   {nm:45s} h = {mp.nstr(mp.findroot(f,(mp.mpf('1e-4'),mp.mpf('0.3465')),solver='illinois'),12)}")
for f in [3,4,5]:
    for e,(T,c) in enumerate([(T0,c0),(T1,c1)]):
        try: print(f"   conductor {f}, parity {e}: T_e(h)+log f-c_e = 0 at h = {mp.nstr(mp.findroot(lambda h:T(h)+mp.log(f)-c,(mp.mpf('1e-4'),mp.mpf(3)),solver='illinois'),12)}")
        except Exception as ex: print("   root fail",f,e,ex)
