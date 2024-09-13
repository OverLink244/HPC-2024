#include <stdio.h>
#include <pthread.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"
#include "omp.h"
#include "common_threads.h"
/* Routine for computing C = A * B + C */
void *MY_MMult_4(void *arg)
{
  // pthread_t tid = pthread_self();
  // printf("pid:%d, tid:%u\n", (int)getpid(), (unsigned int)tid);

  myarg_t *args = (myarg_t *)arg;
  double *a = args->a;
  double *b = args->b;
  double *c = args->c;
  int m = args->m;
  int n = args->n;
  int k = args->k;
  int lda = args->lda;
  int ldb = args->ldb;
  int ldc = args->ldc;
  int i, j, p, r, s;
  // int len = (m - args->start_i >= 8) ? 8 : BLOCK_LEN;
  int len = 4;
  // while(1);
  for (i = args->start_i; i < m; i += len) /* Loop over the rows of C */
  {
    for (j = args->start_j; j < n; j += len) /* Loop over the columns of C */
    {
      for (p = 0; p < k; p++)
      { /* Update C( i,j ) with the inner product of the ith row of A and the jth column of B */
        // block compute
        // for (r = 0; r < len; r++)
        // {
        //   for (s = 0; s < len; s++)
        //   {
        //     C(i + r, j + s) = C(i + r, j + s) + A(i + r, p) * B(p, j + s);
        //   }
        // }
        C(i + 0, j + 0) = C(i + 0, j + 0) + A(i + 0, p) * B(p, j + 0);
        C(i + 0, j + 1) = C(i + 0, j + 1) + A(i + 0, p) * B(p, j + 1);
        C(i + 0, j + 2) = C(i + 0, j + 2) + A(i + 0, p) * B(p, j + 2);
        C(i + 0, j + 3) = C(i + 0, j + 3) + A(i + 0, p) * B(p, j + 3);
        C(i + 1, j + 0) = C(i + 1, j + 0) + A(i + 1, p) * B(p, j + 0);
        C(i + 1, j + 1) = C(i + 1, j + 1) + A(i + 1, p) * B(p, j + 1);
        C(i + 1, j + 2) = C(i + 1, j + 2) + A(i + 1, p) * B(p, j + 2);
        C(i + 1, j + 3) = C(i + 1, j + 3) + A(i + 1, p) * B(p, j + 3);
        C(i + 2, j + 0) = C(i + 2, j + 0) + A(i + 2, p) * B(p, j + 0);
        C(i + 2, j + 1) = C(i + 2, j + 1) + A(i + 2, p) * B(p, j + 1);
        C(i + 2, j + 2) = C(i + 2, j + 2) + A(i + 2, p) * B(p, j + 2);
        C(i + 2, j + 3) = C(i + 2, j + 3) + A(i + 2, p) * B(p, j + 3);
        C(i + 3, j + 0) = C(i + 3, j + 0) + A(i + 3, p) * B(p, j + 0);
        C(i + 3, j + 1) = C(i + 3, j + 1) + A(i + 3, p) * B(p, j + 1);
        C(i + 3, j + 2) = C(i + 3, j + 2) + A(i + 3, p) * B(p, j + 2);
        C(i + 3, j + 3) = C(i + 3, j + 3) + A(i + 3, p) * B(p, j + 3);
      }
    }
  }

  return NULL;
}
void *MY_MMult_8(void *arg)
{
  // pthread_t tid = pthread_self();
  // printf("pid:%d, tid:%u\n", (int)getpid(), (unsigned int)tid);

  myarg_t *args = (myarg_t *)arg;
  double *a = args->a;
  double *b = args->b;
  double *c = args->c;
  int m = args->m;
  int n = args->n;
  int k = args->k;
  int lda = args->lda;
  int ldb = args->ldb;
  int ldc = args->ldc;
  int i, j, p, r, s;
  // int len = (m - args->start_i >= 8) ? 8 : BLOCK_LEN;
  int len = 8;
  // while(1);
  for (i = args->start_i; i < m; i += len) /* Loop over the rows of C */
  {
    for (j = args->start_j; j < n; j += len) /* Loop over the columns of C */
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

  return NULL;
}
void *MY_MMult_2(void *arg)
{
  // pthread_t tid = pthread_self();
  // printf("pid:%d, tid:%u\n", (int)getpid(), (unsigned int)tid);

  myarg_t *args = (myarg_t *)arg;
  double *a = args->a;
  double *b = args->b;
  double *c = args->c;
  int m = args->m;
  int n = args->n;
  int k = args->k;
  int lda = args->lda;
  int ldb = args->ldb;
  int ldc = args->ldc;
  int i, j, p, r, s;
  // int len = (m - args->start_i >= 8) ? 8 : BLOCK_LEN;
  int len = 2;
  // while(1);
  for (i = args->start_i; i < m; i += len) /* Loop over the rows of C */
  {
    for (j = args->start_j; j < n; j += len) /* Loop over the columns of C */
    {
      for (p = 0; p < k; p++)
      { /* Update C( i,j ) with the inner product of the ith row of A and the jth column of B */
        C(i + 0, j + 0) = C(i + 0, j + 0) + A(i + 0, p) * B(p, j + 0);
        C(i + 0, j + 1) = C(i + 0, j + 1) + A(i + 0, p) * B(p, j + 1);
        C(i + 1, j + 0) = C(i + 1, j + 0) + A(i + 1, p) * B(p, j + 0);
        C(i + 1, j + 1) = C(i + 1, j + 1) + A(i + 1, p) * B(p, j + 1);
      }
    }
  }

  return NULL;
}

void MY_MMult(int m, int n, int k, double *a, int lda,
              double *b, int ldb,
              double *c, int ldc)
{

  pthread_t *threads;
  myarg_t *args;
  threads=(pthread_t *)malloc(PTHREAD_NUM*sizeof(pthread_t));
  args=(myarg_t *)malloc(PTHREAD_NUM*sizeof(myarg_t));
  // init pthread
  for (int i = 0; i < PTHREAD_NUM; i++)
  {
    pthread_t p=i;
    threads[i] = p;
  }
  int len = (int)sqrt((double)PTHREAD_NUM);
  if (len % 2)
  {
    printf("Thread num must be the square of 4!\n");
    return;
  }
  for (int i = 0; i < len; i++)
  {
    for (int j = 0; j < len; j++)
    {
      myarg_t arg = {a, b, c, (m / len) * (i + 1), (n / len) * (j + 1), k, lda, ldb, ldc, (m / len) * (i + 0), (n / len) * (j + 0)};
      args[i * len + j] = arg;
    }
  }
  // create
  for (int i = 0; i < PTHREAD_NUM; i++)
  {
    if (m / len <= 2)
    {
      Pthread_create(&threads[i], NULL, MY_MMult_2, &args[i]);
    }
    else if (m / len <= 4)
    {
      Pthread_create(&threads[i], NULL, MY_MMult_4, &args[i]);
    }
    else if (m / len <= 8)
    {
      Pthread_create(&threads[i], NULL, MY_MMult_8, &args[i]);
    }
    else
    {
      Pthread_create(&threads[i], NULL, MY_MMult_8, &args[i]);
    }
  }
  // join
  for (int i = 0; i < PTHREAD_NUM; i++)
    Pthread_join(threads[i], NULL);
  free(threads);
  free(args);
}