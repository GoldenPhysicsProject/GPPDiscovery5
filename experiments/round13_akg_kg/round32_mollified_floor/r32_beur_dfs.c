/* Round 32: Moebius floor of the discrete Beurling (DMV-type) system. Reads sorted Beurling primes (doubles) from
   argv[1], X = argv[2], H = argv[3]; enumerates every squarefree Beurling integer g <= X with mu_B(g) = (-1)^omega,
   and bins mu_B(g) g^{-1/2} (log g)^j, j = 0,1,2, by ceil(log g / H). Writes 3*nb doubles to argv[4]. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
static double *P, X, H, *A; static long np_, nb, cnt = 0;
static void add(double g, int s) {
  double l = log(g); long j = (long)ceil(l / H - 1e-12); if (j < 0) j = 0;
  double w = s / sqrt(g); A[j] += w; A[nb + j] += w * l; A[2*nb + j] += w * l * l; cnt++;
}
static void dfs(long start, double g, int s) {
  add(g, s);
  for (long i = start; i < np_ && g * P[i] <= X; i++) dfs(i + 1, g * P[i], -s);
}
int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "rb"); fseek(f, 0, SEEK_END); np_ = ftell(f) / 8; rewind(f);
  P = malloc(np_ * 8); if (fread(P, 8, np_, f) != (size_t)np_) return 1; fclose(f);
  X = atof(argv[2]); H = atof(argv[3]); nb = (long)ceil(log(X) / H) + 2;
  A = calloc(3 * nb, 8); dfs(0, 1.0, 1);
  f = fopen(argv[4], "wb"); fwrite(A, 8, 3 * nb, f); fclose(f);
  printf("primes %ld  squarefree %ld  bins %ld\n", np_, cnt, nb); return 0;
}
