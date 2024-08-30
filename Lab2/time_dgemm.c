#include "stdio.h"
#include "stdlib.h"
#include "sys/time.h"
#include "time.h"
#include <cblas.h>
#include "omp.h"
// 编译 gcc -o time_dgemm time_dgemm.c –lopenblas
// 运行 ./time_dgemm 1024
void native_dgemm(double alpha, double beta, int M, int N, int K, double *A, double *B, double *C)
{
  double sum0 = 0;
  double sum1 = 0;
  double sum2 = 0;
  double sum3 = 0;
  double sum4 = 0;
  double sum5 = 0;
  double sum6 = 0;
  double sum7 = 0;
  double sum8 = 0;
  double sum9 = 0;
  double sum10 = 0;
  double sum11 = 0;
  double sum12 = 0;
  double sum13 = 0;
  double sum14 = 0;
  double sum15 = 0;
  double sum16 = 0;
  double sum17 = 0;
  double sum18 = 0;
  double sum19 = 0;
  double sum20 = 0;
  double sum21 = 0;
  double sum22 = 0;
  double sum23 = 0;
  double sum24 = 0;
  double sum25 = 0;
  double sum26 = 0;
  double sum27 = 0;
  double sum28 = 0;
  double sum29 = 0;
  double sum30 = 0;
  double sum31 = 0;
  double sum32 = 0;
  double sum33 = 0;
  double sum34 = 0;
  double sum35 = 0;
  double sum36 = 0;
  double sum37 = 0;
  double sum38 = 0;
  double sum39 = 0;
  double sum40 = 0;
  double sum41 = 0;
  double sum42 = 0;
  double sum43 = 0;
  double sum44 = 0;
  double sum45 = 0;
  double sum46 = 0;
  double sum47 = 0;
  double sum48 = 0;
  double sum49 = 0;
  double sum50 = 0;
  double sum51 = 0;
  double sum52 = 0;
  double sum53 = 0;
  double sum54 = 0;
  double sum55 = 0;
  double sum56 = 0;
  double sum57 = 0;
  double sum58 = 0;
  double sum59 = 0;
  double sum60 = 0;
  double sum61 = 0;
  double sum62 = 0;
  double sum63 = 0;
#pragma omp parallel for collapse(2)
  for (int i = 0; i < M; i += 8)
  {
    for (int j = 0; j < N; j += 8)
    {
      for (int k = 0; k < K; k++)
      {
        sum0 += alpha * A[(i + 0) * K + k] * B[k * N + j + 0];
        sum1 += alpha * A[(i + 0) * K + k] * B[k * N + j + 1];
        sum2 += alpha * A[(i + 0) * K + k] * B[k * N + j + 2];
        sum3 += alpha * A[(i + 0) * K + k] * B[k * N + j + 3];
        sum4 += alpha * A[(i + 0) * K + k] * B[k * N + j + 4];
        sum5 += alpha * A[(i + 0) * K + k] * B[k * N + j + 5];
        sum6 += alpha * A[(i + 0) * K + k] * B[k * N + j + 6];
        sum7 += alpha * A[(i + 0) * K + k] * B[k * N + j + 7];
        sum8 += alpha * A[(i + 1) * K + k] * B[k * N + j + 0];
        sum9 += alpha * A[(i + 1) * K + k] * B[k * N + j + 1];
        sum10 += alpha * A[(i + 1) * K + k] * B[k * N + j + 2];
        sum11 += alpha * A[(i + 1) * K + k] * B[k * N + j + 3];
        sum12 += alpha * A[(i + 1) * K + k] * B[k * N + j + 4];
        sum13 += alpha * A[(i + 1) * K + k] * B[k * N + j + 5];
        sum14 += alpha * A[(i + 1) * K + k] * B[k * N + j + 6];
        sum15 += alpha * A[(i + 1) * K + k] * B[k * N + j + 7];
        sum16 += alpha * A[(i + 2) * K + k] * B[k * N + j + 0];
        sum17 += alpha * A[(i + 2) * K + k] * B[k * N + j + 1];
        sum18 += alpha * A[(i + 2) * K + k] * B[k * N + j + 2];
        sum19 += alpha * A[(i + 2) * K + k] * B[k * N + j + 3];
        sum20 += alpha * A[(i + 2) * K + k] * B[k * N + j + 4];
        sum21 += alpha * A[(i + 2) * K + k] * B[k * N + j + 5];
        sum22 += alpha * A[(i + 2) * K + k] * B[k * N + j + 6];
        sum23 += alpha * A[(i + 2) * K + k] * B[k * N + j + 7];
        sum24 += alpha * A[(i + 3) * K + k] * B[k * N + j + 0];
        sum25 += alpha * A[(i + 3) * K + k] * B[k * N + j + 1];
        sum26 += alpha * A[(i + 3) * K + k] * B[k * N + j + 2];
        sum27 += alpha * A[(i + 3) * K + k] * B[k * N + j + 3];
        sum28 += alpha * A[(i + 3) * K + k] * B[k * N + j + 4];
        sum29 += alpha * A[(i + 3) * K + k] * B[k * N + j + 5];
        sum30 += alpha * A[(i + 3) * K + k] * B[k * N + j + 6];
        sum31 += alpha * A[(i + 3) * K + k] * B[k * N + j + 7];
        sum32 += alpha * A[(i + 4) * K + k] * B[k * N + j + 0];
        sum33 += alpha * A[(i + 4) * K + k] * B[k * N + j + 1];
        sum34 += alpha * A[(i + 4) * K + k] * B[k * N + j + 2];
        sum35 += alpha * A[(i + 4) * K + k] * B[k * N + j + 3];
        sum36 += alpha * A[(i + 4) * K + k] * B[k * N + j + 4];
        sum37 += alpha * A[(i + 4) * K + k] * B[k * N + j + 5];
        sum38 += alpha * A[(i + 4) * K + k] * B[k * N + j + 6];
        sum39 += alpha * A[(i + 4) * K + k] * B[k * N + j + 7];
        sum40 += alpha * A[(i + 5) * K + k] * B[k * N + j + 0];
        sum41 += alpha * A[(i + 5) * K + k] * B[k * N + j + 1];
        sum42 += alpha * A[(i + 5) * K + k] * B[k * N + j + 2];
        sum43 += alpha * A[(i + 5) * K + k] * B[k * N + j + 3];
        sum44 += alpha * A[(i + 5) * K + k] * B[k * N + j + 4];
        sum45 += alpha * A[(i + 5) * K + k] * B[k * N + j + 5];
        sum46 += alpha * A[(i + 5) * K + k] * B[k * N + j + 6];
        sum47 += alpha * A[(i + 5) * K + k] * B[k * N + j + 7];
        sum48 += alpha * A[(i + 6) * K + k] * B[k * N + j + 0];
        sum49 += alpha * A[(i + 6) * K + k] * B[k * N + j + 1];
        sum50 += alpha * A[(i + 6) * K + k] * B[k * N + j + 2];
        sum51 += alpha * A[(i + 6) * K + k] * B[k * N + j + 3];
        sum52 += alpha * A[(i + 6) * K + k] * B[k * N + j + 4];
        sum53 += alpha * A[(i + 6) * K + k] * B[k * N + j + 5];
        sum54 += alpha * A[(i + 6) * K + k] * B[k * N + j + 6];
        sum55 += alpha * A[(i + 6) * K + k] * B[k * N + j + 7];
        sum56 += alpha * A[(i + 7) * K + k] * B[k * N + j + 0];
        sum57 += alpha * A[(i + 7) * K + k] * B[k * N + j + 1];
        sum58 += alpha * A[(i + 7) * K + k] * B[k * N + j + 2];
        sum59 += alpha * A[(i + 7) * K + k] * B[k * N + j + 3];
        sum60 += alpha * A[(i + 7) * K + k] * B[k * N + j + 4];
        sum61 += alpha * A[(i + 7) * K + k] * B[k * N + j + 5];
        sum62 += alpha * A[(i + 7) * K + k] * B[k * N + j + 6];
        sum63 += alpha * A[(i + 7) * K + k] * B[k * N + j + 7];
      }
      C[(i + 0) * N + j + 0] = sum0 + C[(i + 0) * N + j + 0] * beta;
      C[(i + 0) * N + j + 1] = sum1 + C[(i + 0) * N + j + 1] * beta;
      C[(i + 0) * N + j + 2] = sum2 + C[(i + 0) * N + j + 2] * beta;
      C[(i + 0) * N + j + 3] = sum3 + C[(i + 0) * N + j + 3] * beta;
      C[(i + 0) * N + j + 4] = sum4 + C[(i + 0) * N + j + 4] * beta;
      C[(i + 0) * N + j + 5] = sum5 + C[(i + 0) * N + j + 5] * beta;
      C[(i + 0) * N + j + 6] = sum6 + C[(i + 0) * N + j + 6] * beta;
      C[(i + 0) * N + j + 7] = sum7 + C[(i + 0) * N + j + 7] * beta;
      C[(i + 1) * N + j + 0] = sum8 + C[(i + 1) * N + j + 0] * beta;
      C[(i + 1) * N + j + 1] = sum9 + C[(i + 1) * N + j + 1] * beta;
      C[(i + 1) * N + j + 2] = sum10 + C[(i + 1) * N + j + 2] * beta;
      C[(i + 1) * N + j + 3] = sum11 + C[(i + 1) * N + j + 3] * beta;
      C[(i + 1) * N + j + 4] = sum12 + C[(i + 1) * N + j + 4] * beta;
      C[(i + 1) * N + j + 5] = sum13 + C[(i + 1) * N + j + 5] * beta;
      C[(i + 1) * N + j + 6] = sum14 + C[(i + 1) * N + j + 6] * beta;
      C[(i + 1) * N + j + 7] = sum15 + C[(i + 1) * N + j + 7] * beta;
      C[(i + 2) * N + j + 0] = sum16 + C[(i + 2) * N + j + 0] * beta;
      C[(i + 2) * N + j + 1] = sum17 + C[(i + 2) * N + j + 1] * beta;
      C[(i + 2) * N + j + 2] = sum18 + C[(i + 2) * N + j + 2] * beta;
      C[(i + 2) * N + j + 3] = sum19 + C[(i + 2) * N + j + 3] * beta;
      C[(i + 2) * N + j + 4] = sum20 + C[(i + 2) * N + j + 4] * beta;
      C[(i + 2) * N + j + 5] = sum21 + C[(i + 2) * N + j + 5] * beta;
      C[(i + 2) * N + j + 6] = sum22 + C[(i + 2) * N + j + 6] * beta;
      C[(i + 2) * N + j + 7] = sum23 + C[(i + 2) * N + j + 7] * beta;
      C[(i + 3) * N + j + 0] = sum24 + C[(i + 3) * N + j + 0] * beta;
      C[(i + 3) * N + j + 1] = sum25 + C[(i + 3) * N + j + 1] * beta;
      C[(i + 3) * N + j + 2] = sum26 + C[(i + 3) * N + j + 2] * beta;
      C[(i + 3) * N + j + 3] = sum27 + C[(i + 3) * N + j + 3] * beta;
      C[(i + 3) * N + j + 4] = sum28 + C[(i + 3) * N + j + 4] * beta;
      C[(i + 3) * N + j + 5] = sum29 + C[(i + 3) * N + j + 5] * beta;
      C[(i + 3) * N + j + 6] = sum30 + C[(i + 3) * N + j + 6] * beta;
      C[(i + 3) * N + j + 7] = sum31 + C[(i + 3) * N + j + 7] * beta;
      C[(i + 4) * N + j + 0] = sum32 + C[(i + 4) * N + j + 0] * beta;
      C[(i + 4) * N + j + 1] = sum33 + C[(i + 4) * N + j + 1] * beta;
      C[(i + 4) * N + j + 2] = sum34 + C[(i + 4) * N + j + 2] * beta;
      C[(i + 4) * N + j + 3] = sum35 + C[(i + 4) * N + j + 3] * beta;
      C[(i + 4) * N + j + 4] = sum36 + C[(i + 4) * N + j + 4] * beta;
      C[(i + 4) * N + j + 5] = sum37 + C[(i + 4) * N + j + 5] * beta;
      C[(i + 4) * N + j + 6] = sum38 + C[(i + 4) * N + j + 6] * beta;
      C[(i + 4) * N + j + 7] = sum39 + C[(i + 4) * N + j + 7] * beta;
      C[(i + 5) * N + j + 0] = sum40 + C[(i + 5) * N + j + 0] * beta;
      C[(i + 5) * N + j + 1] = sum41 + C[(i + 5) * N + j + 1] * beta;
      C[(i + 5) * N + j + 2] = sum42 + C[(i + 5) * N + j + 2] * beta;
      C[(i + 5) * N + j + 3] = sum43 + C[(i + 5) * N + j + 3] * beta;
      C[(i + 5) * N + j + 4] = sum44 + C[(i + 5) * N + j + 4] * beta;
      C[(i + 5) * N + j + 5] = sum45 + C[(i + 5) * N + j + 5] * beta;
      C[(i + 5) * N + j + 6] = sum46 + C[(i + 5) * N + j + 6] * beta;
      C[(i + 5) * N + j + 7] = sum47 + C[(i + 5) * N + j + 7] * beta;
      C[(i + 6) * N + j + 0] = sum48 + C[(i + 6) * N + j + 0] * beta;
      C[(i + 6) * N + j + 1] = sum49 + C[(i + 6) * N + j + 1] * beta;
      C[(i + 6) * N + j + 2] = sum50 + C[(i + 6) * N + j + 2] * beta;
      C[(i + 6) * N + j + 3] = sum51 + C[(i + 6) * N + j + 3] * beta;
      C[(i + 6) * N + j + 4] = sum52 + C[(i + 6) * N + j + 4] * beta;
      C[(i + 6) * N + j + 5] = sum53 + C[(i + 6) * N + j + 5] * beta;
      C[(i + 6) * N + j + 6] = sum54 + C[(i + 6) * N + j + 6] * beta;
      C[(i + 6) * N + j + 7] = sum55 + C[(i + 6) * N + j + 7] * beta;
      C[(i + 7) * N + j + 0] = sum56 + C[(i + 7) * N + j + 0] * beta;
      C[(i + 7) * N + j + 1] = sum57 + C[(i + 7) * N + j + 1] * beta;
      C[(i + 7) * N + j + 2] = sum58 + C[(i + 7) * N + j + 2] * beta;
      C[(i + 7) * N + j + 3] = sum59 + C[(i + 7) * N + j + 3] * beta;
      C[(i + 7) * N + j + 4] = sum60 + C[(i + 7) * N + j + 4] * beta;
      C[(i + 7) * N + j + 5] = sum61 + C[(i + 7) * N + j + 5] * beta;
      C[(i + 7) * N + j + 6] = sum62 + C[(i + 7) * N + j + 6] * beta;
      C[(i + 7) * N + j + 7] = sum63 + C[(i + 7) * N + j + 7] * beta;
      sum0 = 0;
      sum1 = 0;
      sum2 = 0;
      sum3 = 0;
      sum4 = 0;
      sum5 = 0;
      sum6 = 0;
      sum7 = 0;
      sum8 = 0;
      sum9 = 0;
      sum10 = 0;
      sum11 = 0;
      sum12 = 0;
      sum13 = 0;
      sum14 = 0;
      sum15 = 0;
      sum16 = 0;
      sum17 = 0;
      sum18 = 0;
      sum19 = 0;
      sum20 = 0;
      sum21 = 0;
      sum22 = 0;
      sum23 = 0;
      sum24 = 0;
      sum25 = 0;
      sum26 = 0;
      sum27 = 0;
      sum28 = 0;
      sum29 = 0;
      sum30 = 0;
      sum31 = 0;
      sum32 = 0;
      sum33 = 0;
      sum34 = 0;
      sum35 = 0;
      sum36 = 0;
      sum37 = 0;
      sum38 = 0;
      sum39 = 0;
      sum40 = 0;
      sum41 = 0;
      sum42 = 0;
      sum43 = 0;
      sum44 = 0;
      sum45 = 0;
      sum46 = 0;
      sum47 = 0;
      sum48 = 0;
      sum49 = 0;
      sum50 = 0;
      sum51 = 0;
      sum52 = 0;
      sum53 = 0;
      sum54 = 0;
      sum55 = 0;
      sum56 = 0;
      sum57 = 0;
      sum58 = 0;
      sum59 = 0;
      sum60 = 0;
      sum61 = 0;
      sum62 = 0;
      sum63 = 0;
    }
  }
}
void printMatrix(double *A, int m, int n)
{
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      printf("%lf\t", A[i * m + j]);
    }
    printf("\n");
  }
}
int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    printf("Input Error\n");
    return 1;
  }

  printf("test!\n");
  int i, m, n, k;
  m = n = k = atoi(argv[1]);

  int sizeofa = m * k;
  int sizeofb = k * n;
  int sizeofc = m * n;
  int lda = m;
  int ldb = k;
  int ldc = m;

  double alpha = 1.2;
  double beta = 0.001;

  struct timeval start, finish, native_start, native_finish;
  double duration, native_duration;

  double *A = (double *)malloc(sizeof(double) * sizeofa);
  double *B = (double *)malloc(sizeof(double) * sizeofb);
  double *C = (double *)malloc(sizeof(double) * sizeofc);
  double *Native_C = (double *)malloc(sizeof(double) * sizeofc);
  srand((unsigned)time(NULL));

  for (i = 0; i < sizeofa; i++)
  {
    A[i] = i % 3 + 1; // (rand() % 100) / 100.0;
  }

  for (i = 0; i < sizeofb; i++)
  {
    B[i] = i % 3 + 1; //(rand()%100)/10.0;
  }

  for (i = 0; i < sizeofc; i++)
  {
    C[i] = 0.1;
  }

  printf("m=%d,n=%d,k=%d,alpha=%lf,beta=%lf,sizeofc=%d\n", m, n, k, alpha, beta, sizeofc);
  gettimeofday(&start, NULL);
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
  gettimeofday(&finish, NULL);

  // native_dgemm
  gettimeofday(&native_start, NULL);
  native_dgemm(alpha, beta, m, n, k, A, B, Native_C);
  gettimeofday(&native_finish, NULL);

  // 对比
  // printf("C:\n");
  // printMatrix(C,m,n);
  // printf("Native_C:\n");
  // printMatrix(Native_C,m,n);
  // 验证答案
  int flag = 1;
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      // printf("%lf-%lf\t", C[i * m + j], Native_C[i * m + j]);
      if (abs(C[i * m + j] - Native_C[i * m + j]) < __FLT_EPSILON__)
      {
        continue;
      }
      else
      {
        flag = 0;
        break;
      }
    }
    // printf("\n");
    if (!flag)
    {
      printf("Wrong answer!");
      break;
    }
  }
  if (flag)
  {
    printf("Answer true!");
  }
  转成成秒数
  duration = (double)(finish.tv_sec - start.tv_sec) + (double)(finish.tv_usec - start.tv_usec) / 1.0e6;
  native_duration = (double)(native_finish.tv_sec - native_start.tv_sec) + (double)(native_finish.tv_usec - native_start.tv_usec) / 1.0e6;
  double gflops = 4.0 * m * n * k, native_gflops = 4.0 * m * n * k;
  gflops = gflops / duration * 1.0e-9;

  native_gflops = native_gflops / native_duration * 1.0e-9;
  FILE *fp;
  fp = fopen("timeDGEMM.txt", "a"); // 追加写
  fprintf(fp, "Cblas:%dx%dx%d\t%lf s\t%lf GFLOPS\n", m, n, k, duration, gflops);
  fprintf(fp, "Native:%dx%dx%d\t%lf s\t%lf GFLOPS\n", m, n, k, native_duration, native_gflops);
  fclose(fp);

  free(A);
  free(B);
  free(C);
  free(Native_C);
  return 0;
}
