#include "func.h"

void ColMajor2RowMajor(int M, int N, double* mat) {
	double* mat_ = new double[M * N];
	for (int i = 0;i < M;i++) {
		for (int j = 0;j < N;j++) {
			mat_[i * N + j] = mat[i + j * M];
		}
	}
	for (int i = 0;i < M * N;i++) {
		mat[i] = mat_[i];
	}
	delete[] mat_;
}

void RowMajor2ColMajor(int M, int N, double* mat) {
	double* mat_ = new double[M * N];
	for (int i = 0;i < M;i++) {
		for (int j = 0;j < N;j++) {
			mat_[i + j * M] = mat[i * N + j];
		}
	}
	for (int i = 0;i < M * N;i++) {
		mat[i] = mat_[i];
	}
	delete[] mat_;
}

void LoopRowMajorOrderingPre(FUNC_PARAM_PRE) {
	ColMajor2RowMajor(M, K, A);
	ColMajor2RowMajor(K, N, B);
	ColMajor2RowMajor(M, N, C);
}

void LoopRowMajorOrderingPro(FUNC_PARAM_PRO) { // 三重循环 + 所有元素按行排序
	for (int i = 0;i < M;i++) {
		for (int j = 0;j < N;j++) {
			C[i * ldc + j] *= beta;
			for (int k = 0;k < K;k++) {
				C[i * ldc + j] += alpha * A[i * lda + k] * B[k * ldb + j];
			}
		}
	}
}

void LoopRowMajorOrderingPost(FUNC_PARAM_POST) {
	RowMajor2ColMajor(M, N, C);
}