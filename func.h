#pragma once

#define FUNC_PARAM_PRE int M, int N, int K, double* A, double* B, double* C
#define FUNC_PARAM_PRO int M, int N, int K, const double* A, const double* B, double* C, double alpha, double beta, int lda, int ldb, int ldc
#define FUNC_PARAM_POST int M, int N, double* C

#define FUNC_DECLARE(NAME) void NAME##Pre(FUNC_PARAM_PRE); void NAME##Pro(FUNC_PARAM_PRO); void NAME##Post(FUNC_PARAM_POST);

FUNC_DECLARE(LoopRowMajorOrderingPre)