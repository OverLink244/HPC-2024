#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
void dgemm(int start_i, int start_j, int m, int n, int k, int lda, int ldb, int ldc, double *a, double *b, double *c)
{
    int i, j, p;
    // int len = 2;
    for (i = start_i; i < m; i++) /* Loop over the rows of C */
    {
        for (j = start_j; j < n; j++) /* Loop over the columns of C */
        {
            C(i + 0, j + 0) = 0.0;
            for (p = 0; p < k; p++)
            { /* Update C( i,j ) with the inner product of the ith row of A and the jth column of B */
                C(i + 0, j + 0) = C(i + 0, j + 0) + A(i + 0, p) * B(p, j + 0);
                // C(i + 0, j + 1) = C(i + 0, j + 1) + A(i + 0, p) * B(p, j + 1);
                // C(i + 1, j + 0) = C(i + 1, j + 0) + A(i + 1, p) * B(p, j + 0);
                // C(i + 1, j + 1) = C(i + 1, j + 1) + A(i + 1, p) * B(p, j + 1);
            }
        }
    }
}
void read_matrix(FILE *file, double *matrix, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            fscanf(file, "%lf,", &matrix[i * cols + j]);
        }
    }
}
void skip_lines(FILE *file, int lines_to_skip)
{
    char buffer[1024]; // 用来存储每行的字符
    for (int i = 0; i < lines_to_skip; i++)
    {
        fgets(buffer, sizeof(buffer), file); // 跳过当前行
    }
}
void read_file(int data_col, int data_row, double *a, double *b, double *cref)
{
    FILE *file = fopen("mpi_test_data.m", "r"); // 打开数据文件
    if (file == NULL)
    {
        printf("无法打开文件\n");
        return;
    }
    int nums, row, col;
    while (fscanf(file, "data%d: %d*%d\n", &nums, &row, &col) != EOF)
    {
        if (row == data_row && col == data_col)
        {
            // 读取A矩阵
            fscanf(file, "A={\n");
            read_matrix(file, a, row, col);
            fscanf(file, "\n}\n");
            // 读取B矩阵
            fscanf(file, "B={\n");
            read_matrix(file, b, row, col);
            fscanf(file, "\n}\n");
            // 读取C矩阵
            fscanf(file, "C={\n");
            read_matrix(file, cref, row, col);
            break;
        }
        else
        {
            // 如果大小不匹配，跳过矩阵数据
            skip_lines(file, (row + 3) * 3);
        }
    }
    fclose(file); // 关闭文件
}
void matrix_print(double *a, char *matrix_name, int m, int n)
{
    printf("Matrix %s:[\n", matrix_name);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%lf\t", a[i * m + j]);
        }
        printf("\n");
    }
    printf("]\n");
}
double compare_matrices(int m, int n, double *a, int lda, double *b, int ldb)
{
    int i, j;
    double max_diff = 0.0, diff;

    for (j = 0; j < n; j++)
        for (i = 0; i < m; i++)
        {
            diff = abs(A(i, j) - B(i, j));
            max_diff = (diff > max_diff ? diff : max_diff);
        }

    return max_diff;
}
