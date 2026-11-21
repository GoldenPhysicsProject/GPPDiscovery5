/* Round 33: enumerate every squarefree Beurling integer g <= X (primes from argv[1], R32 construction) and write
   (log g, mu_B(g)) as doubles to argv[3]. Sorting and the exact step integrals are done in r33_beurling.py. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
static double *P, X; static long np_, cnt = 0; static FILE *fo;
static void dfs(long start, double g, int s) {
  double rec[2] = {log(g), (double)s}; fwrite(rec, 8, 2, fo); cnt++;
  for (long i = start; i < np_ && g * P[i] <= X; i++) dfs(i + 1, g * P[i], -s);
}
int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "rb"); fseek(f, 0, SEEK_END); np_ = ftell(f) / 8; rewind(f);
  P = malloc(np_ * 8); if (fread(P, 8, np_, f) != (size_t)np_) return 1; fclose(f);
  X = atof(argv[2]); fo = fopen(argv[3], "wb"); dfs(0, 1.0, 1); fclose(fo);
  printf("primes %ld squarefree %ld\n", np_, cnt); return 0;
}
