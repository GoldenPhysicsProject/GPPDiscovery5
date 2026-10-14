# Pole-renormalized prime jump measure (reply to Codex non-Levy audit, 2026-10-14).
# (A) M_prime((0,R]) = 2 sum_{n<=e^R} Lambda(n)/sqrt(n) vs pole density int_0^R 2e^{u/2}du = 4(e^{R/2}-1).
# (B) DH (Davenport-Heilbronn, q=5): coefficients b(n) of -F'/F are NOT >= 0 and NOT prime-power supported.
# Float (mpmath / python floats), not certified.
from mpmath import mp, mpf, sqrt, log, exp
import math
mp.dps = 30
def mangoldt_table(N):
    L = [0.0]*(N+1); sieve = bytearray([1])*(N+1)
    for p in range(2, N+1):
        if sieve[p]:
            for m in range(p*p, N+1, p): sieve[m] = 0
            q = p
            while q <= N: L[q] = math.log(p); q *= p
    return L
R_list = [2,4,6,8,10,12,14]
N = int(math.exp(max(R_list)))+1
Lam = mangoldt_table(N)
print("(A) R, M_prime((0,R]), 4(e^{R/2}-1), difference, difference/e^{R/2}")
acc = 0.0; n = 1
for R in R_list:
    X = math.exp(R)
    while n+1 <= X:
        n += 1; acc += 2*Lam[n]/math.sqrt(n)
    pole = 4*(math.exp(R/2)-1)
    print(R, round(acc,7), round(pole,7), round(acc-pole,7), round((acc-pole)/math.exp(R/2),7))
# (B) DH coefficients
s5 = sqrt(5); xi = (sqrt(10-2*s5)-2)/(s5-1)
a = {1: mpf(1), 2: xi, 3: -xi, 4: mpf(-1), 0: mpf(0)}
A = lambda n: a[n % 5]
M = 60; b = [mpf(0)]*(M+1)
for n in range(2, M+1):
    b[n] = A(n)*log(n) - sum(A(d)*b[n//d] for d in range(2, n+1) if n % d == 0 and n//d >= 2)
neg = [n for n in range(2, M+1) if b[n] < -mpf(10)**-20]
nonpp = [n for n in range(2, M+1) if abs(b[n]) > mpf(10)**-20 and len({p for p in range(2,n+1) if n%p==0 and all(p%k for k in range(2,int(p**.5)+1))})>1]
print("(B) DH: n<=60 with b(n)<0:", neg)
print("(B) DH: n<=60 with b(n)!=0 and n not a prime power:", nonpp[:15], "...")
print("(B) b(2..7) =", [mp.nstr(b[k], 8) for k in range(2, 8)])
