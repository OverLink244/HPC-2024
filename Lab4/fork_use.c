#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include "time.h"


#define A(i, j) a[(i) * lda + (j)]
#define B(i, j) b[(i) * ldb + (j)]
#define C(i, j) c[(i) * ldc + (j)]
void Native_gemm(int m, int n, int k, double *a, int lda,
                 double *b, int ldb,
                 double *c, int ldc)
{
    int i, j, p;
    // fprintf(stderr, __FILE__ "\n");
    for (i = 0; i < m; i++) /* Loop over the rows of C */
    {
        for (j = 0; j < n; j++) /* Loop over the columns of C */
        {
            for (p = 0; p < k; p++)
            { /* Update C( i,j ) with the inner product of the ith row of A and the jth column of B */
                C(i, j) = C(i, j) + A(i, p) * B(p, j);
            }
        }
    }
}
int main(int argc, char *argv[])
{
    printf("pid:%d\n", (int)getpid());
    int i, m, n, k;
    m = n = k = atoi(argv[1]);
    int sizeofa = m * k;
    int sizeofb = k * n;
    int sizeofc = m * n;
    int lda = m;
    int ldb = k;
    int ldc = m;
    double *a = (double *)malloc(sizeof(double) * sizeofa);
    double *b = (double *)malloc(sizeof(double) * sizeofb);
    double *c = (double *)malloc(sizeof(double) * sizeofc);
    srand((unsigned)time(NULL));

    for (i = 0; i < sizeofa; i++)
    {
        a[i] = i % 3 + 1; // (rand() % 100) / 100.0;
    }

    for (i = 0; i < sizeofb; i++)
    {
        b[i] = i % 3 + 1; //(rand()%100)/10.0;
    }

    for (i = 0; i < sizeofc; i++)
    {
        c[i] = 0.1;
    }
    int rc=fork();
    if(rc<0){
        printf("error");
        exit(0);
    }
    else if(rc==0){
        printf("PID child:%d\n",(int)getpid());
        Native_gemm(m, n, k, a, lda, b, ldb, c, ldc);
    }
    else{
        printf("PID parent:%d\n",(int)getpid());
        Native_gemm(m, n, k, a, lda, b, ldb, c, ldc);
    }
    
    return 0;
}