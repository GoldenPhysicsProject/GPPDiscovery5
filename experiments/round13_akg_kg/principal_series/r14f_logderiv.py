"""Round 14F: prime side of the Weil form = Dirichlet coefficients c_F(n) of -F'/F (exact integer recursion,
a(n) log n = sum_{d|n} c(d) a(n/d), a(1)=1). Tests: prime-power support, additivity over products
(W_{FG} = W_F + W_G at the coefficient level), Ramanujan bound |c(p)|/log p <= degree. No zeros used."""
import json, math
N = 3000
def chi_m4(a): return [0, 1, 0, -1][a % 4]
def chi_5(a): return [0, 1, -1, -1, 1][a % 5]
def chi_m20(a): return chi_m4(a)*chi_5(a)
def dconv(f, g): return lambda n: sum(f(d)*g(n//d) for d in range(1, n+1) if n % d == 0)
one = lambda n: 1
import functools
def tab(f): return [0]+[f(n) for n in range(1, N+1)]
def logder(a):
    c = [0.0]*(N+1)
    for n in range(2, N+1):
        s = a[n]*math.log(n) - sum(c[d]*a[n//d] for d in range(2, n) if n % d == 0)
        c[n] = s/a[1]
    return c
def is_pp(n):
    for p in range(2, n+1):
        if n % p == 0:
            while n % p == 0: n //= p
            return n == 1
    return False
KAP = (math.sqrt(10-2*math.sqrt(5))-2)/(math.sqrt(5)-1)
dh = {0: 0, 1: 1, 2: KAP, 3: -KAP, 4: -1}
th = 0.2
models = {
 "zeta": tab(one),
 "L(chi_5)": tab(chi_5),
 "zeta_K Q(sqrt5)=zeta*L(chi_5)": tab(dconv(one, chi_5)),
 "Epstein x2+5y2 (r_Q/2)": tab(lambda n: dconv(one, chi_m20)(n) + dconv(chi_m4, chi_5)(n)),
 "Epstein pieces zeta_K(-20)": tab(dconv(one, chi_m20)),
 "Epstein pieces L_psi": tab(dconv(chi_m4, chi_5)),
 "DH": tab(lambda n: dh[n % 5]),
 "F_theta=zeta(s-th)zeta(s+th)": tab(lambda n: sum((d**th)*((n//d)**(-th)) for d in range(1, n+1) if n % d == 0)),
}
for k in models: models[k][1] = 1 if k != "Epstein x2+5y2 (r_Q/2)" else models[k][1]
# Epstein r_Q/2 has a(1)=1 already (b1(1)+b2(1)=2 -> divide by 2)
models["Epstein x2+5y2 (r_Q/2)"] = [x/2 for x in models["Epstein x2+5y2 (r_Q/2)"]]
C = {k: logder(v) for k, v in models.items()}
out = {}
for k, c in C.items():
    nz = [n for n in range(2, N+1) if abs(c[n]) > 1e-9]
    npp = [n for n in nz if not is_pp(n)]
    primes = [p for p in range(2, N+1) if is_pp(p) and all(p % q for q in range(2, int(p**.5)+1))]
    ram = max(abs(c[p])/math.log(p) for p in primes)
    out[k] = {"nonzero_atoms<=N": len(nz), "non_prime_power_atoms": len(npp), "first_non_pp": npp[:6],
              "negative_atoms": sum(1 for n in nz if c[n] < -1e-9), "max_p |c(p)|/log p": round(ram, 4)}
add = max(abs(C["zeta_K Q(sqrt5)=zeta*L(chi_5)"][n] - C["zeta"][n] - C["L(chi_5)"][n]) for n in range(2, N+1))
add_eps_sum = max(abs(C["Epstein x2+5y2 (r_Q/2)"][n] - C["Epstein pieces zeta_K(-20)"][n] - C["Epstein pieces L_psi"][n]) for n in range(2, N+1))
out["additivity_over_product: max|c_{zeta L5} - c_zeta - c_L5|"] = add
out["non-additivity_over_sum: max|c_EPS - c_zetaK - c_Lpsi|"] = add_eps_sum
out["N"] = N
json.dump(out, open("r14f_logderiv.json", "w"), indent=1)
print(json.dumps(out, indent=1))
