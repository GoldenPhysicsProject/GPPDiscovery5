/* Round 34: ratio-class LAYER decomposition of the R33 Fejer square for the real Moebius function (no zeros).
   a_n = mu(n)/sqrt(n) - c[n=1], c = 1/zeta(1/2), L = log X, kernel K(m,n) = L - |log m - log n| (m,n <= X).
   Layer k(m,n) = omega(u)+omega(v) for the reduced ratio v/u = n/m  (= omega of the symmetric difference of m,n).
   E_X(z) = sum_{m,n} a_m a_n z^{k(m,n)} K(m,n) = sum_k E_k z^k;   E(1) = G, E(0) = E_0 = H, E(-1) = G of |mu| (Codex/Turbo translate).
   Exact computation via the gcd/d-decomposition (unique factorization):
     sum_{m,n} mu(m)mu(n)(mn)^{-1/2} z^{k} K = sum_{d squarefree} (1-z^2)^{omega(d)} / d * B_d(z),
     B_d(z) = sum_{m',n' <= X/d, (m'n',d)=1} b_m' b_n' K(m',n'),   b_n = mu(n) z^{omega(n)} n^{-1/2},
   [z^{-2 omega(g)} = sum_{d|g} prod_{p|d}(z^{-2}-1)], plus the exact centering correction at n = 1.
   B_d is a bilinear Fejer form: B = L s^2 - 2 sum_i A_i (s - A_i)(x_{i+1}-x_i), A_i prefix sums over sorted x = log n'.
   Evaluated at NZ roots of unity; E_k by inverse DFT.  argv: X  twist(0/1: flip mu(p) for primes in (X/2,3X/4])  out.json */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#define NZ 48
int main(int argc, char **argv) {
  long X = atol(argv[1]); int twist = atoi(argv[2]);
  double L = log((double)X), c = 1.0 / -1.4603545088095868;
  signed char *mu = malloc(X + 1); unsigned char *om = calloc(X + 1, 1); int *spf = calloc(X + 1, sizeof(int));
  for (long i = 0; i <= X; i++) mu[i] = 1;
  for (long p = 2; p <= X; p++) if (!spf[p]) {
    for (long k = p; k <= X; k += p) { if (!spf[k]) spf[k] = p; mu[k] = -mu[k]; om[k]++; }
    if (p <= X / p) for (long k = p * p; k <= X; k += p * p) mu[k] = 0;
  }
  if (twist) for (long p = X / 2 + 1; p <= 3 * X / 4; p++) if (spf[p] == p) mu[p] = -mu[p];
  double complex z[NZ], zp[NZ][24], E[NZ]; double *x = malloc((X + 1) * sizeof(double)), *w = malloc((X + 1) * sizeof(double));
  int *o = malloc((X + 1) * sizeof(int));
  for (int j = 0; j < NZ; j++) { z[j] = cexp(2 * M_PI * I * j / NZ); zp[j][0] = 1; for (int k = 1; k < 24; k++) zp[j][k] = zp[j][k-1] * z[j]; E[j] = 0; }
  for (long d = 1; d <= X; d++) {
    if (mu[d] == 0) continue;
    long Y = X / d; int pd[16], np = 0; long t = d;
    while (t > 1) { int p = spf[t]; pd[np++] = p; t /= p; }
    long cnt = 0;
    for (long n = 1; n <= Y; n++) {
      if (!mu[n]) continue;
      int ok = 1; for (int i = 0; i < np; i++) if (n % pd[i] == 0) { ok = 0; break; }
      if (!ok) continue;
      x[cnt] = log((double)n); w[cnt] = mu[n] / sqrt((double)n); o[cnt] = om[n]; cnt++;
    }
    for (int j = 0; j < NZ; j++) {
      double complex s = 0, A = 0, acc = 0;
      for (long i = 0; i < cnt; i++) s += w[i] * zp[j][o[i]];
      for (long i = 0; i < cnt; i++) {
        A += w[i] * zp[j][o[i]];
        double dx = (i + 1 < cnt ? x[i+1] : L - log((double)d)) - x[i];
        acc += A * (s - A) * dx;
      }
      double complex Bd = L * s * s - 2 * acc;
      double complex f = 1; for (int i = 0; i < np; i++) f *= (1 - z[j] * z[j]);
      E[j] += f / (double)d * Bd;
    }
  }
  /* centering: a_1 = 1 - c instead of 1  ->  diagonal ((1-c)^2 - 1) L  and cross terms -2c sum_{n>=2} mu(n) n^{-1/2} (L - log n) z^{omega(n)} */
  double cen[24] = {0};
  cen[0] = ((1 - c) * (1 - c) - 1) * L;
  for (long n = 2; n <= X; n++) if (mu[n]) cen[om[n]] += -2 * c * mu[n] / sqrt((double)n) * (L - log((double)n));
  for (int j = 0; j < NZ; j++) for (int k = 0; k < 24; k++) E[j] += cen[k] * zp[j][k];
  FILE *f = fopen(argv[3], "w");
  fprintf(f, "{\"X\": %ld, \"twist\": %d, \"L\": %.17g, \"c\": %.17g, \"Ek\": [", X, twist, L, c);
  for (int k = 0; k < NZ; k++) {
    double complex s = 0; for (int j = 0; j < NZ; j++) s += E[j] * cexp(-2 * M_PI * I * (double)j * k / NZ);
    s /= NZ; fprintf(f, "%s[%.17g, %.6g]", k ? ", " : "", creal(s), cimag(s));
  }
  fprintf(f, "], \"E_at_1\": %.17g, \"E_at_minus1\": %.17g}\n", creal(E[0]), creal(E[NZ/2]));
  fclose(f); return 0;
}
