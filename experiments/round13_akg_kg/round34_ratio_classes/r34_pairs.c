/* Round 34: brute-force ratio-class layers for ANY model, from records (sorted by log n):
   double logn, double a_n (centering already applied at n = 1), int nf, int prime ids (with multiplicity, sorted, <= 24).
   For every pair m <= n: k = |multiset symmetric difference| = Omega(u) + Omega(v) of the reduced ratio n/m;
   E[k] += (m==n ? 1 : 2) a_m a_n (L - (log n - log m)).   argv: records.bin  L  out.json */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef struct { double x, a; int nf; int p[24]; } R;
int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "rb"); fseek(f, 0, SEEK_END); long N = ftell(f) / sizeof(R); rewind(f);
  R *r = malloc(N * sizeof(R)); if (fread(r, sizeof(R), N, f) != (size_t)N) return 1; fclose(f);
  double L = atof(argv[2]); double E[64] = {0};
  for (long i = 0; i < N; i++) {
    E[0] += r[i].a * r[i].a * L;
    for (long j = i + 1; j < N; j++) {
      int a = 0, b = 0, k = 0, *P = r[i].p, *Q = r[j].p, na = r[i].nf, nb = r[j].nf;
      while (a < na && b < nb) { if (P[a] == Q[b]) { a++; b++; } else if (P[a] < Q[b]) { a++; k++; } else { b++; k++; } }
      k += (na - a) + (nb - b);
      E[k] += 2 * r[i].a * r[j].a * (L - (r[j].x - r[i].x));
    }
  }
  f = fopen(argv[3], "w"); fprintf(f, "{\"N\": %ld, \"L\": %.17g, \"Ek\": [", N, L);
  for (int k = 0; k < 64; k++) fprintf(f, "%s%.17g", k ? ", " : "", E[k]);
  fprintf(f, "]}\n"); fclose(f); return 0;
}
