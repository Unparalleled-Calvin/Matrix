#include "func.h"

double* ColMajor2RowMajor(int M, int N, const double* mat) {
	double* mat_ = new double[M * N];
	for (int i = 0;i < M;i++) {
		for (int j = 0;j < N;j++) {
			mat_[i * N + j] = mat[i + j * M];
		}
	}
	return mat_;
}

double* RowMajor2ColMajor(int M, int N, double* mat) {
	double* mat_ = new double[M * N];
	for (int i = 0;i < M;i++) {
		for (int j = 0;j < N;j++) {
			mat_[i + j * M] = mat[i * N + j];
		}
	}
	return mat_;
}

// LoopRowMajorOrdering

void LoopRowMajorOrderingPre(FUNC_PARAM_PRE) {
	*pA = ColMajor2RowMajor(M, K, A);
	*pB = ColMajor2RowMajor(K, N, B);
	*pC = ColMajor2RowMajor(M, N, C);
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
	double* C__ = RowMajor2ColMajor(M, N, C_);
	memcpy(C, C__, M * N * sizeof(double));
	delete[] C__;
}

// LoopRowMajorBlocking

#define BLOCK_SIZE 5

void LoopRowMajorBlockingPre(FUNC_PARAM_PRE) {
	*pA = ColMajor2RowMajor(M, K, A);
	*pB = ColMajor2RowMajor(K, N, B);
	*pC = ColMajor2RowMajor(M, N, C);
}

void LoopRowMajorBlockingPro(FUNC_PARAM_PRO) {
	for (int ih = 0;ih < M;ih += BLOCK_SIZE) {
		for (int jh = 0;jh < N;jh += BLOCK_SIZE) {
			for (int kh = 0;kh < K;kh += BLOCK_SIZE) {
				for (int il = 0;il < BLOCK_SIZE && ih + il < M;il++) {
					for (int jl = 0;jl < BLOCK_SIZE && jh + jl < N;jl++) {
						C[(ih + il) * ldc + jh + jl] *= kh ? 1 : beta; // 控制使beta只算一次
						for (int kl = 0;kl < BLOCK_SIZE && kh + kl < K;kl++) {
							C[(ih + il) * ldc + jh + jl] += alpha * A[(ih + il) * lda + kh + kl] * B[(kh + kl) * ldb + jh + jl];
						}
					}
				}
			}
		}
	}
}

void LoopRowMajorBlockingPost(FUNC_PARAM_POST) {
	double* C__ = RowMajor2ColMajor(M, N, C_);
	memcpy(C, C__, M * N * sizeof(double));
	delete[] C__;
}