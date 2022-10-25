#pragma once
#include <cstdio>
#include <cstdlib>
#include <cstring>

#define FUNC_PARAM_PRE int M, int N, int K, const double* A, const double* B, const double* C, double** pA, double **pB, double** pC, double beta
#define FUNC_PARAM_PRO int M, int N, int K, const double* __restrict A, const double* __restrict B, double* __restrict C, double alpha, int lda, int ldb, int ldc
#define FUNC_PARAM_POST int M, int N, double* C, double *C_

#define FUNC_DECLARE(NAME) void NAME##Pre(FUNC_PARAM_PRE); void NAME##Pro(FUNC_PARAM_PRO); void NAME##Post(FUNC_PARAM_POST);

FUNC_DECLARE(LoopRowMajorOrdering)
FUNC_DECLARE(LoopRowMajorBlocking)
FUNC_DECLARE(LoopRowMajorPacking)
FUNC_DECLARE(RecursionRowMajorOrdering)
FUNC_DECLARE(RecursionRowMajorBlocking)
FUNC_DECLARE(RecursionRowMajorPacking)
FUNC_DECLARE(RecursionZmortonOrdering)
FUNC_DECLARE(RecursionZmortonPacking)

extern int BLOCK_SIZE;