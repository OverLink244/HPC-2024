#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include "defs.h"
#include "mpi.h"

void dgemm(int start_i, int start_j, int m, int n, int k, int lda, int ldb, int ldc, double *a, double *b, double *c);
void read_file(int data_col, int data_row, double *a, double *b, double *cref);
void matrix_print(double *a, char *matrix_name, int m, int n);
void read_matrix(FILE *file, double *matrix, int rows, int cols);
double compare_matrices(int m, int n, double *a, int lda, double *b, int ldb);
int main(int argc, char *argv[])
{
    const int m = M, n = N, k = K;
    int lda = k, ldb = n, ldc = n;
    // double *a, *b, *c, *cref;
    // a = (double *)malloc(m * k * sizeof(double));
    // b = (double *)malloc(k * n * sizeof(double));

    // cref = (double *)malloc(m * n * sizeof(double));
    double a[] = {
        1.900078e-02,
        8.548636e-02,
        -3.314322e-01,
        2.440604e-02,
        -6.461445e-01,
        -7.460251e-01,
        5.706861e-01,
        -6.391164e-02,
        -3.681220e-01,
        3.359572e-01,
        -8.358069e-01,
        1.730207e-01,
        -7.762620e-01,
        -9.918196e-01,
        -1.155304e-01,
        -6.060053e-01,
        4.034051e-01,
        -6.577163e-01,
        -4.593707e-01,
        4.775171e-01,
        7.053617e-01,
        -9.989944e-01,
        -4.241922e-01,
        -1.056461e-01,
        -2.655067e-01,
        6.615583e-01,
        -8.514753e-01,
        7.433776e-01,
        6.223437e-01,
        -6.140962e-01,
        7.067059e-01,
        4.568796e-01,
        7.526247e-01,
        7.238108e-01,
        3.524515e-01,
        1.559346e-01,
        9.097746e-01,
        1.735603e-01,
        5.194394e-02,
        9.321860e-01,
        -2.580107e-01,
        5.384423e-01,
        6.203171e-01,
        9.294072e-01,
        8.085664e-01,
        -2.500817e-01,
        4.568334e-01,
        5.136816e-01,
        5.050480e-02,
        -8.948411e-01,
        -7.157754e-02,
        -3.329414e-01,
        -7.221695e-01,
        -9.464281e-01,
        -4.673714e-01,
        -6.788623e-02,
        -1.443309e-01,
        9.334593e-01,
        1.383701e-01,
        -2.143441e-01,
        3.723747e-01,
        -7.046941e-02,
        5.669680e-01,
        7.379329e-01,
    };
    double b[] = {
        -7.319215e-01,
        8.165651e-01,
        -4.029224e-02,
        8.852518e-01,
        5.361704e-01,
        -2.435351e-01,
        -5.578819e-01,
        -3.161686e-01,
        6.921052e-01,
        5.488294e-01,
        -9.377434e-01,
        -2.239507e-01,
        5.038825e-01,
        -1.625671e-01,
        5.086909e-02,
        -8.637625e-02,
        -9.638247e-01,
        -2.678706e-01,
        2.168604e-01,
        8.610245e-01,
        -8.182718e-01,
        -7.231037e-01,
        -1.094956e-02,
        1.762435e-01,
        8.817096e-01,
        8.885760e-01,
        5.845223e-01,
        8.227684e-01,
        7.079334e-01,
        -2.383066e-01,
        7.650440e-01,
        1.702999e-01,
        3.401642e-01,
        7.696433e-01,
        -6.350064e-01,
        -3.245633e-01,
        -6.507111e-01,
        -4.004399e-01,
        5.476602e-01,
        -8.601459e-01,
        5.999431e-01,
        -1.096578e-01,
        5.024259e-02,
        6.562465e-02,
        -9.604593e-01,
        -9.579758e-01,
        -8.819950e-01,
        6.378629e-01,
        -7.135664e-01,
        -6.739431e-01,
        7.975177e-01,
        4.041545e-01,
        -1.618930e-01,
        -5.223948e-01,
        -7.603972e-01,
        -1.936786e-01,
        2.209206e-01,
        -1.503654e-01,
        -3.824202e-01,
        9.420213e-01,
        -4.072602e-01,
        -9.257813e-01,
        9.560251e-02,
        1.709661e-01,
    };
    double cref[] = {
        -7.024897e-01,
        -6.175928e-01,
        7.138586e-01,
        6.358145e-02,
        1.412361e+00,
        9.497799e-01,
        -1.198863e-01,
        -1.091874e-01,
        5.495471e-01,
        -5.828324e-02,
        2.023827e-01,
        -1.409117e+00,
        2.501541e+00,
        2.480544e+00,
        8.435441e-01,
        -7.664723e-02,
        3.326345e-02,
        1.469986e+00,
        -1.598580e-02,
        -6.368559e-02,
        1.211029e+00,
        1.221015e+00,
        1.691703e+00,
        -1.250213e+00,
        1.568243e+00,
        1.036256e+00,
        -3.969600e-01,
        -3.099133e-02,
        1.298360e+00,
        -5.740635e-02,
        1.148578e+00,
        -9.824463e-01,
        3.303513e-01,
        1.561958e+00,
        -1.425547e+00,
        1.551175e+00,
        -5.565087e-01,
        -2.013690e+00,
        1.271708e-01,
        -7.343229e-01,
        6.956050e-01,
        1.009126e+00,
        -1.748624e-01,
        1.339494e+00,
        -2.857656e-01,
        -1.493150e+00,
        1.240700e+00,
        -5.529886e-01,
        -1.375815e+00,
        -8.533851e-01,
        6.912177e-01,
        -1.710138e-01,
        8.812948e-01,
        1.767113e+00,
        4.605126e-01,
        8.840335e-02,
        2.721833e-01,
        -3.181246e-02,
        -1.034846e+00,
        4.047740e-01,
        -4.389420e-01,
        -1.226528e+00,
        -1.319781e-01,
        -3.960062e-01,
    };
    double *c;
    c = (double *)malloc(m * n * sizeof(double));
    double *localA, *localB, *localC;
    // read_file(MATRIX_ROW, MATRIX_COL, a, b, cref);
    // matrix_print(cref,"cref",m,n);
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int namelen, numprocs, myid;
    double startwtime = 0.0;
    // MPI init
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Get_processor_name(processor_name, &namelen);
    MPI_Status status;
    int split_rows = m / (numprocs - 1);
    // fprintf(stdout, "Process %d of %d is on %s\n", myid, numprocs, processor_name);
    // fflush(stdout);
    if (myid == 0)
        startwtime = MPI_Wtime();
    if (myid == 0)
    {
        // dgemm_2_2(0, 0, m, n, k, lda, ldb, ldc, a, b, c);
        // if ((int)sqrt((double)(numprocs - 1)) % 2 != 0)
        // {
        //     printf("Numprocs wrong!\n");
        //     return -1;
        // }
        int block_len = (int)sqrt((double)(numprocs - 1));
        // send
        for (int i = 0; i < numprocs - 1; i++)
        {
            MPI_Send(&a[split_rows * i * k], split_rows * k, MPI_DOUBLE, i + 1, 1, MPI_COMM_WORLD);
            MPI_Send(b, n * k, MPI_DOUBLE, i + 1, 2, MPI_COMM_WORLD);
            // matrix_print(b,"B",k,n);
        }
        // receive
        for (int i = 0; i < numprocs - 1; i++)
        {
            MPI_Recv(&c[split_rows * i * n], split_rows * n, MPI_DOUBLE, i + 1, 3, MPI_COMM_WORLD, &status);
        }
    }
    else
    {

        localA = (double *)malloc(split_rows * k * sizeof(double));
        localB = (double *)malloc(n * k * sizeof(double));
        localC = (double *)malloc(split_rows * n * sizeof(double));
        MPI_Recv(localA, split_rows * k, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(localB, k * n, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &status);
        dgemm(0, 0, split_rows, n, k, lda, ldb, ldc, localA, localB, localC);
        // matrix_print(localA, "Recv_A", split_rows, k);
        // matrix_print(localB, "Recv_B", k, n);
        // matrix_print(localC, "Recv_C", split_rows, k);

        MPI_Send(localC, split_rows * n, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
        free(localA);
        free(localB);
        free(localC);
    }
    if (myid == 0)
    {
        matrix_print(c, "C", m, n);
        double endtime = MPI_Wtime();
        double diff = compare_matrices(m, n, c, n, cref, n);
        printf("Time:%lf\tDiff:%lf", endtime - startwtime, diff);
        free(c);
    }

    MPI_Finalize();
    return 0;
}