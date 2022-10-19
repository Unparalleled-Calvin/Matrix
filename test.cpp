#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <ctime>

#include "func.h"
#include "debug.h"

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

int BLOCK_SIZE;
int func_no;
const char* func_names[8] = {
    "LoopRowMajorOrdering",
    "LoopRowMajorBlocking",
    "LoopRowMajorPacking",
    "RecursionRowMajorOrdering",
    "RecursionRowMajorBlocking",
    "RecursionRowMajorPacking",
    "RecursionZmortonOrdering",
    "RecursionZmortonPacking"
};

void (*func_pres[8])(FUNC_PARAM_PRE) = {
    &LoopRowMajorOrderingPre,
    &LoopRowMajorBlockingPre,
    &LoopRowMajorPackingPre,
    &RecursionRowMajorOrderingPre,
    &RecursionRowMajorBlockingPre,
    &RecursionRowMajorPackingPre,
    &RecursionZmortonOrderingPre,
    &RecursionZmortonPackingPre
};

void (*func_pros[8])(FUNC_PARAM_PRO) = {
    &LoopRowMajorOrderingPro,
    &LoopRowMajorBlockingPro,
    &LoopRowMajorPackingPro,
    &RecursionRowMajorOrderingPro,
    &RecursionRowMajorBlockingPro,
    &RecursionRowMajorPackingPro,
    &RecursionZmortonOrderingPro,
    &RecursionZmortonPackingPro
};

void (*func_posts[8])(FUNC_PARAM_POST) = {
    &LoopRowMajorOrderingPost,
    &LoopRowMajorBlockingPost,
    &LoopRowMajorPackingPost,
    &RecursionRowMajorOrderingPost,
    &RecursionRowMajorBlockingPost,
    &RecursionRowMajorPackingPost,
    &RecursionZmortonOrderingPost,
    &RecursionZmortonPackingPost
};

void RandomFill(struct drand48_data* buf_p, double* d, size_t count)
{
    for (size_t i = 0; i < count; ++i)
    {
        drand48_r(buf_p, d + i);
        d[i] = 2 * d[i] - 1.0;
    }
}

double student_gemm(int m, int n, int k, const double* A, const double* B, double* C, double alpha, double beta, int lda, int ldb, int ldc)
{
    /* TODO */
    struct timespec start, end;
    double *A_, *B_, *C_;
    func_pres[func_no](m, n, k, A, B, C, &A_, &B_, &C_, beta);
    double t = MEASURE_VOID(func_pros[func_no], m, n, k, A_, B_, C_, alpha, k, n, n);
    func_posts[func_no](m, n, C, C_);
    delete[] A_;
    delete[] B_;
    delete[] C_;
    return t;
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
    double t_min = __DBL_MAX__;

    for (int i = 0; i < TRIAL; i++) {
        double t = student_gemm(m, n, k, A, B, C, alpha, beta, lda, ldb, ldc);

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
    if (argc != 4 && argc != 5 && argc != 6)
    {
        printf("Test usage: ./test m n k [func_no] [block_size]\n");
        exit(-1);
    }

    int m = atoi(argv[1]);
    int n = atoi(argv[2]);
    int k = atoi(argv[3]);
    BLOCK_SIZE = argc == 6 ? atoi(argv[5]) : 1;
    func_no = argc == 5 ? atoi(argv[4]) : 0;

    printf("input: %d x %d x %d\nmethod: %s\nblock size:%d\n", m, n, k, func_names[func_no], BLOCK_SIZE);
    fflush(stdout);

    if (!isPowerOf2(m) || !isPowerOf2(n) || !isPowerOf2(k)) {
        printf("error!\neach side of the matrix should be a power of 2.\n");
        return -1;
    }
    if (BLOCK_SIZE > m || BLOCK_SIZE > n || BLOCK_SIZE > k) {
        printf("error!\neach side of the matrix should be larger than the block size.\n");
        return -1;
    }

    mm_test(m, n, k);

    return 0;
}
