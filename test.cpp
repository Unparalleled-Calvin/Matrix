#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <ctime>

#include "func.h"

#define MEASURE(__ret_ptr, __func, ...)           \
    ((clock_gettime(CLOCK_MONOTONIC, &start), \
      *(__ret_ptr) = __func(__VA_ARGS__),         \
      clock_gettime(CLOCK_MONOTONIC, &end)),  \
     (end.tv_sec - start.tv_sec) + 1e-9 * (end.tv_nsec - start.tv_nsec))

#define MEASURE_VOID(__func, ...)                 \
    ((clock_gettime(CLOCK_MONOTONIC, &start), \
      __func(__VA_ARGS__),                        \
      clock_gettime(CLOCK_MONOTONIC, &end)),  \
     (end.tv_sec - start.tv_sec) + 1e-9 * (end.tv_nsec - start.tv_nsec))

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define ABS(x) (((x) >= 0.0) ? (x) : -(x))

void RandomFill(struct drand48_data* buf_p, double* d, size_t count)
{
    for (size_t i = 0; i < count; ++i)
    {
        drand48_r(buf_p, d + i);
        d[i] = 2 * d[i] - 1.0;
    }
}

void student_gemm(int m, int n, int k, const double* A, const double* B, double* C, double alpha, double beta, int lda, int ldb, int ldc)
{
    /* TODO */
    double *A_, *B_, *C_;
    LoopRowMajorBlockingPre(m, n, k, A, B, C, &A_, &B_, &C_);
    LoopRowMajorBlockingPro(m, n, k, A_, B_, C_, alpha, beta, k, n, n);
    LoopRowMajorBlockingPost(m, n, C, C_);
    delete[] A_;
    delete[] B_;
    delete[] C_;
}

void naive_gemm(int m, int n, int k, const double* A, const double* B, double* C, double alpha, double beta, int lda, int ldb, int ldc)
{
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            C[i + j * ldc] *= beta;
            for (int p = 0; p < k; p++) {
                C[i + j * ldc] += alpha * A[i + p * lda] * B[p + j * ldb];
            }
        }
    }
    return;
}

void mm_test(int m, int n, int k)
{
    int lda = m;
    int ldb = k;
    int ldc = m;

    double* A = (double*)aligned_alloc(64, m * k * sizeof(double));
    double* B = (double*)aligned_alloc(64, k * n * sizeof(double));
    double* C = (double*)aligned_alloc(64, m * n * sizeof(double));
    double* C_ans = (double*)aligned_alloc(64, m * n * sizeof(double));

    struct drand48_data buffer;
    srand48_r(time(NULL), &buffer);

    RandomFill(&buffer, A, m * k);
    RandomFill(&buffer, B, k * n);
    RandomFill(&buffer, C, m * n);

    memcpy(C_ans, C, sizeof(double) * m * n);

    double alpha, beta;

    drand48_r(&buffer, &alpha);
    drand48_r(&buffer, &beta);
    alpha = 2 * alpha - 1.0;
    beta = 2 * beta - 1.0;

    /* test performance */

    const int TRIAL = 5;
    struct timespec start, end;
    double t_min = __DBL_MAX__;

    for (int i = 0; i < TRIAL; i++) {
        double t = MEASURE_VOID(student_gemm, m, n, k, A, B, C, alpha, beta, lda, ldb, ldc);

        t_min = MIN(t, t_min);
    }

    printf("minimal time spent: %.4f ms\n", t_min * 1000);
    fflush(stdout);

    /* test correctness */

    memcpy(C, C_ans, sizeof(double) * m * n);
    student_gemm(m, n, k, A, B, C, alpha, beta, lda, ldb, ldc);
    naive_gemm(m, n, k, A, B, C_ans, alpha, beta, lda, ldb, ldc);

    double max_err = __DBL_MIN__;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            int idx = i + j * ldc;
            double err = ABS(C_ans[idx] - C[idx]);
            max_err = MAX(err, max_err);
        }
    }
    const double threshold = 1e-7;
    const char* judge_s = (max_err < threshold) ? "correct" : "wrong";

    printf("result: %s (err = %e)\n", judge_s, max_err);

    free(A);
    free(B);
    free(C);
    free(C_ans);
}

int main(int argc, const char* argv[])
{
    if (argc != 4)
    {
        printf("Test usage: ./test m n k\n");
        exit(-1);
    }

    int m = atoi(argv[1]);
    int n = atoi(argv[2]);
    int k = atoi(argv[3]);

    printf("input: %d x %d x %d\n", m, n, k);
    fflush(stdout);

    mm_test(m, n, k);
}
