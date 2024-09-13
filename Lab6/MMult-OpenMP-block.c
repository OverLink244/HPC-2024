#include <stdio.h>
#include "defs.h"
#include "omp.h"
/* Routine for computing C = A * B + C */

void MY_MMult(int m, int n, int k, double *a, int lda,
              double *b, int ldb,
              double *c, int ldc)
{
    int i, j, p, r, s;
    if (m < 8 || n < 8 || k < 8)
    {
        printf("Matrix can't small than 8*8!\n");
        return ;
    }
    int len = 8;
    // int len = m >= 8 ? 8 : BLOCK_LEN;
    fprintf(stderr, __FILE__ "\n");
    // printf("start\n");
#pragma omp parallel for collapse(2) private(i, j, p) shared(a,b,c) schedule(static)
    for (i = 0; i < m; i += len) /* Loop over the rows of C */
    {
        for (j = 0; j < n; j += len) /* Loop over the columns of C */
        {
            for (p = 0; p < k; p++)
            { /* Update C( i,j ) with the inner product of the ith row of A and the jth column of B */
                C(i + 0, j + 0) = C(i + 0, j + 0) + A(i + 0, p) * B(p, j + 0);
                C(i + 0, j + 1) = C(i + 0, j + 1) + A(i + 0, p) * B(p, j + 1);
                C(i + 0, j + 2) = C(i + 0, j + 2) + A(i + 0, p) * B(p, j + 2);
                C(i + 0, j + 3) = C(i + 0, j + 3) + A(i + 0, p) * B(p, j + 3);
                C(i + 0, j + 4) = C(i + 0, j + 4) + A(i + 0, p) * B(p, j + 4);
                C(i + 0, j + 5) = C(i + 0, j + 5) + A(i + 0, p) * B(p, j + 5);
                C(i + 0, j + 6) = C(i + 0, j + 6) + A(i + 0, p) * B(p, j + 6);
                C(i + 0, j + 7) = C(i + 0, j + 7) + A(i + 0, p) * B(p, j + 7);
                C(i + 1, j + 0) = C(i + 1, j + 0) + A(i + 1, p) * B(p, j + 0);
                C(i + 1, j + 1) = C(i + 1, j + 1) + A(i + 1, p) * B(p, j + 1);
                C(i + 1, j + 2) = C(i + 1, j + 2) + A(i + 1, p) * B(p, j + 2);
                C(i + 1, j + 3) = C(i + 1, j + 3) + A(i + 1, p) * B(p, j + 3);
                C(i + 1, j + 4) = C(i + 1, j + 4) + A(i + 1, p) * B(p, j + 4);
                C(i + 1, j + 5) = C(i + 1, j + 5) + A(i + 1, p) * B(p, j + 5);
                C(i + 1, j + 6) = C(i + 1, j + 6) + A(i + 1, p) * B(p, j + 6);
                C(i + 1, j + 7) = C(i + 1, j + 7) + A(i + 1, p) * B(p, j + 7);
                C(i + 2, j + 0) = C(i + 2, j + 0) + A(i + 2, p) * B(p, j + 0);
                C(i + 2, j + 1) = C(i + 2, j + 1) + A(i + 2, p) * B(p, j + 1);
                C(i + 2, j + 2) = C(i + 2, j + 2) + A(i + 2, p) * B(p, j + 2);
                C(i + 2, j + 3) = C(i + 2, j + 3) + A(i + 2, p) * B(p, j + 3);
                C(i + 2, j + 4) = C(i + 2, j + 4) + A(i + 2, p) * B(p, j + 4);
                C(i + 2, j + 5) = C(i + 2, j + 5) + A(i + 2, p) * B(p, j + 5);
                C(i + 2, j + 6) = C(i + 2, j + 6) + A(i + 2, p) * B(p, j + 6);
                C(i + 2, j + 7) = C(i + 2, j + 7) + A(i + 2, p) * B(p, j + 7);
                C(i + 3, j + 0) = C(i + 3, j + 0) + A(i + 3, p) * B(p, j + 0);
                C(i + 3, j + 1) = C(i + 3, j + 1) + A(i + 3, p) * B(p, j + 1);
                C(i + 3, j + 2) = C(i + 3, j + 2) + A(i + 3, p) * B(p, j + 2);
                C(i + 3, j + 3) = C(i + 3, j + 3) + A(i + 3, p) * B(p, j + 3);
                C(i + 3, j + 4) = C(i + 3, j + 4) + A(i + 3, p) * B(p, j + 4);
                C(i + 3, j + 5) = C(i + 3, j + 5) + A(i + 3, p) * B(p, j + 5);
                C(i + 3, j + 6) = C(i + 3, j + 6) + A(i + 3, p) * B(p, j + 6);
                C(i + 3, j + 7) = C(i + 3, j + 7) + A(i + 3, p) * B(p, j + 7);
                C(i + 4, j + 0) = C(i + 4, j + 0) + A(i + 4, p) * B(p, j + 0);
                C(i + 4, j + 1) = C(i + 4, j + 1) + A(i + 4, p) * B(p, j + 1);
                C(i + 4, j + 2) = C(i + 4, j + 2) + A(i + 4, p) * B(p, j + 2);
                C(i + 4, j + 3) = C(i + 4, j + 3) + A(i + 4, p) * B(p, j + 3);
                C(i + 4, j + 4) = C(i + 4, j + 4) + A(i + 4, p) * B(p, j + 4);
                C(i + 4, j + 5) = C(i + 4, j + 5) + A(i + 4, p) * B(p, j + 5);
                C(i + 4, j + 6) = C(i + 4, j + 6) + A(i + 4, p) * B(p, j + 6);
                C(i + 4, j + 7) = C(i + 4, j + 7) + A(i + 4, p) * B(p, j + 7);
                C(i + 5, j + 0) = C(i + 5, j + 0) + A(i + 5, p) * B(p, j + 0);
                C(i + 5, j + 1) = C(i + 5, j + 1) + A(i + 5, p) * B(p, j + 1);
                C(i + 5, j + 2) = C(i + 5, j + 2) + A(i + 5, p) * B(p, j + 2);
                C(i + 5, j + 3) = C(i + 5, j + 3) + A(i + 5, p) * B(p, j + 3);
                C(i + 5, j + 4) = C(i + 5, j + 4) + A(i + 5, p) * B(p, j + 4);
                C(i + 5, j + 5) = C(i + 5, j + 5) + A(i + 5, p) * B(p, j + 5);
                C(i + 5, j + 6) = C(i + 5, j + 6) + A(i + 5, p) * B(p, j + 6);
                C(i + 5, j + 7) = C(i + 5, j + 7) + A(i + 5, p) * B(p, j + 7);
                C(i + 6, j + 0) = C(i + 6, j + 0) + A(i + 6, p) * B(p, j + 0);
                C(i + 6, j + 1) = C(i + 6, j + 1) + A(i + 6, p) * B(p, j + 1);
                C(i + 6, j + 2) = C(i + 6, j + 2) + A(i + 6, p) * B(p, j + 2);
                C(i + 6, j + 3) = C(i + 6, j + 3) + A(i + 6, p) * B(p, j + 3);
                C(i + 6, j + 4) = C(i + 6, j + 4) + A(i + 6, p) * B(p, j + 4);
                C(i + 6, j + 5) = C(i + 6, j + 5) + A(i + 6, p) * B(p, j + 5);
                C(i + 6, j + 6) = C(i + 6, j + 6) + A(i + 6, p) * B(p, j + 6);
                C(i + 6, j + 7) = C(i + 6, j + 7) + A(i + 6, p) * B(p, j + 7);
                C(i + 7, j + 0) = C(i + 7, j + 0) + A(i + 7, p) * B(p, j + 0);
                C(i + 7, j + 1) = C(i + 7, j + 1) + A(i + 7, p) * B(p, j + 1);
                C(i + 7, j + 2) = C(i + 7, j + 2) + A(i + 7, p) * B(p, j + 2);
                C(i + 7, j + 3) = C(i + 7, j + 3) + A(i + 7, p) * B(p, j + 3);
                C(i + 7, j + 4) = C(i + 7, j + 4) + A(i + 7, p) * B(p, j + 4);
                C(i + 7, j + 5) = C(i + 7, j + 5) + A(i + 7, p) * B(p, j + 5);
                C(i + 7, j + 6) = C(i + 7, j + 6) + A(i + 7, p) * B(p, j + 6);
                C(i + 7, j + 7) = C(i + 7, j + 7) + A(i + 7, p) * B(p, j + 7);
            }
        }
    }
    // printf("my:%lf", C(0, 0));
}
