#include <stdio.h>
#include "cblas.h"
#include "defs.h"

/* Routine for computing C = A * B + C */

void MY_MMult(int m, int n, int k, double *a, int lda,
              double *b, int ldb,
              double *c, int ldc)
{

  fprintf(stderr, __FILE__ "\n");
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1, a, lda, b, ldb, 1, c, ldc);
}
