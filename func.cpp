#include "func.h"

inline int max(int a, int b) {
	return a > b ? a : b;
}

void matScale(int M, int N, double* mat, double scale) {
	for (int i = 0;i < M * N;i++) {
		mat[i] *= scale;
	}
}

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

double* RowMajor2Packing(int M, int N, double* mat) { // 根据BLOCK_SIZE做packing
	double* mat_ = new double[M * N];
	int pos = 0;
	for (int ih = 0;ih < M;ih += BLOCK_SIZE) {
		for (int jh = 0;jh < N;jh += BLOCK_SIZE) {
			for (int il = 0;il < BLOCK_SIZE && ih + il < M;il++) {
				for (int jl = 0;jl < BLOCK_SIZE && jh + jl < N;jl++) {
					mat_[pos++] = mat[(ih + il) * N + jh + jl];
				}
			}
		}
	}
	return mat_;
}

double* Packing2RowMajor(int M, int N, double* mat) {
	double* mat_ = new double[M * N];
	int pos = 0;
	for (int ih = 0;ih < M;ih += BLOCK_SIZE) {
		for (int jh = 0;jh < N;jh += BLOCK_SIZE) {
			for (int il = 0;il < BLOCK_SIZE && ih + il < M;il++) {
				for (int jl = 0;jl < BLOCK_SIZE && jh + jl < N;jl++) {
					mat_[(ih + il) * N + jh + jl] = mat[pos++];
				}
			}
		}
	}
	return mat_;
}

// LoopRowMajorOrdering

void LoopRowMajorOrderingPre(FUNC_PARAM_PRE) {
	*pA = ColMajor2RowMajor(M, K, A);
	*pB = ColMajor2RowMajor(K, N, B);
	*pC = ColMajor2RowMajor(M, N, C);
	matScale(M, N, *pC, beta);
}

#define pos(mat, i, j, ld) (mat)[(i) * (ld) + (j)]

void LoopRowMajorOrderingPro(FUNC_PARAM_PRO) { // 三重循环 + 所有元素按行排序
	for (int i = 0;i < M;i++) {
		for (int j = 0;j < N;j++) {
			for (int k = 0;k < K;k++) {
				pos(C, i, j, ldc) += alpha * pos(A, i, k, lda) * pos(B, k, j, ldb);
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

void LoopRowMajorBlockingPre(FUNC_PARAM_PRE) {
	*pA = ColMajor2RowMajor(M, K, A);
	*pB = ColMajor2RowMajor(K, N, B);
	*pC = ColMajor2RowMajor(M, N, C);
	matScale(M, N, *pC, beta);
}

#define pos(mat, ih, il, jh, jl, ld) (mat)[((ih) + (il)) * (ld) + (jh) + (jl)]

void LoopRowMajorBlockingPro(FUNC_PARAM_PRO) { // 三重循环 + 显式分块
	for (int ih = 0;ih < M;ih += BLOCK_SIZE) {
		for (int jh = 0;jh < N;jh += BLOCK_SIZE) {
			for (int kh = 0;kh < K;kh += BLOCK_SIZE) {
				for (int il = 0;il < BLOCK_SIZE && ih + il < M;il++) {
					for (int jl = 0;jl < BLOCK_SIZE && jh + jl < N;jl++) {
						for (int kl = 0;kl < BLOCK_SIZE && kh + kl < K;kl++) {
							pos(C, ih, il, jh, jl, ldc) += alpha * pos(A, ih, il, kh, kl, lda) * pos(B, kh, kl, jh, jl, ldb);
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

// LoopRowMajorPacking

void LoopRowMajorPackingPre(FUNC_PARAM_PRE) {
	double* A_ = ColMajor2RowMajor(M, K, A);
	double* B_ = ColMajor2RowMajor(K, N, B);
	double* C_ = ColMajor2RowMajor(M, N, C);
	*pA = RowMajor2Packing(M, K, A_);
	*pB = RowMajor2Packing(K, N, B_);
	*pC = RowMajor2Packing(M, N, C_);
	delete[] A_;
	delete[] B_;
	delete[] C_;
	matScale(M, N, *pC, beta);
}

// 仅当MNK均整除BLOCK_SIZE时生效
#define pos(mat, ih, il, jh, jl, ld) (mat)[(ih) * (ld) + (jh) * (BLOCK_SIZE) + (il) * (BLOCK_SIZE) + (jl)]

void LoopRowMajorPackingPro(FUNC_PARAM_PRO) { // 三重循环 + 显式分块 + 小块连续
	for (int ih = 0;ih < M;ih += BLOCK_SIZE) {
		for (int jh = 0;jh < N;jh += BLOCK_SIZE) {
			for (int kh = 0;kh < K;kh += BLOCK_SIZE) {
				for (int il = 0;il < BLOCK_SIZE && ih + il < M;il++) {
					for (int jl = 0;jl < BLOCK_SIZE && jh + jl < N;jl++) {
						for (int kl = 0;kl < BLOCK_SIZE && kh + kl < K;kl++) {
							pos(C, ih, il, jh, jl, ldc) += alpha * pos(A, ih, il, kh, kl, lda) * pos(B, kh, kl, jh, jl, ldb);
						}
					}
				}
			}
		}
	}
}

void LoopRowMajorPackingPost(FUNC_PARAM_POST) {
	double* C__ = Packing2RowMajor(M, N, C_);
	double* C___ = RowMajor2ColMajor(M, N, C__);
	memcpy(C, C___, M * N * sizeof(double));
	delete[] C__;
	delete[] C___;
}

// RecursionRowMajorOrdering

void RecursionRowMajorOrderingPre(FUNC_PARAM_PRE) {
	*pA = ColMajor2RowMajor(M, K, A);
	*pB = ColMajor2RowMajor(K, N, B);
	*pC = ColMajor2RowMajor(M, N, C);
	matScale(M, N, *pC, beta);
}

// 仅当MNK均为2幂次时生效
#define pos(mat, i, j, ld) (mat)[(i) * (ld) + (j)]

void RecursionRowMajorOrderingPro(FUNC_PARAM_PRO) { // 递归结构 + 递归到单个元素
	int L = max(max(M, N), K);
	int X = L / 2;
	if (L == 1) {
		pos(C, 0, 0, ldc) += alpha * pos(A, 0, 0, lda) * pos(B, 0, 0, ldb);
	}
	else if (L == K) {
		RecursionRowMajorOrderingPro(M, N, X, A, B, C, alpha, lda, ldb, ldc);
		RecursionRowMajorOrderingPro(M, N, X, A + X, B + X * ldb, C, alpha, lda, ldb, ldc);
	}
	else if (L == M) {
		RecursionRowMajorOrderingPro(X, N, K, A, B, C, alpha, lda, ldb, ldc);
		RecursionRowMajorOrderingPro(X, N, K, A + X * lda, B, C + X * ldc, alpha, lda, ldb, ldc);
	}
	else if (L == N) {
		RecursionRowMajorOrderingPro(M, X, K, A, B, C, alpha, lda, ldb, ldc);
		RecursionRowMajorOrderingPro(M, X, K, A, B + X, C + X, alpha, lda, ldb, ldc);
	}
}

void RecursionRowMajorOrderingPost(FUNC_PARAM_POST) {
	double* C__ = RowMajor2ColMajor(M, N, C_);
	memcpy(C, C__, M * N * sizeof(double));
	delete[] C__;
}

// RecursionRowMajorOrdering

void RecursionRowMajorBlockingPre(FUNC_PARAM_PRE) {
	*pA = ColMajor2RowMajor(M, K, A);
	*pB = ColMajor2RowMajor(K, N, B);
	*pC = ColMajor2RowMajor(M, N, C);
	matScale(M, N, *pC, beta);
}

// 仅当MNK均为2幂次时生效
#define pos(mat, i, j, ld) (mat)[(i) * (ld) + (j)]

void RecursionRowMajorBlockingPro(FUNC_PARAM_PRO) { // 递归结构 + 分块
	int L = max(max(M, N), K);
	int X = L / 2;
	if (L <= BLOCK_SIZE) {
		for (int i = 0;i < M;i++) {
			for (int j = 0; j < N;j++) {
				for (int k = 0; k < K;k++) {
					pos(C, i, j, ldc) += alpha * pos(A, i, k, lda) * pos(B, k, j, ldb);
				}
			}
		}
	}
	else if (L == K) {
		RecursionRowMajorBlockingPro(M, N, X, A, B, C, alpha, lda, ldb, ldc);
		RecursionRowMajorBlockingPro(M, N, X, A + X, B + X * ldb, C, alpha, lda, ldb, ldc);
	}
	else if (L == M) {
		RecursionRowMajorBlockingPro(X, N, K, A, B, C, alpha, lda, ldb, ldc);
		RecursionRowMajorBlockingPro(X, N, K, A + X * lda, B, C + X * ldc, alpha, lda, ldb, ldc);
	}
	else if (L == N) {
		RecursionRowMajorBlockingPro(M, X, K, A, B, C, alpha, lda, ldb, ldc);
		RecursionRowMajorBlockingPro(M, X, K, A, B + X, C + X, alpha, lda, ldb, ldc);
	}
}

void RecursionRowMajorBlockingPost(FUNC_PARAM_POST) {
	double* C__ = RowMajor2ColMajor(M, N, C_);
	memcpy(C, C__, M * N * sizeof(double));
	delete[] C__;
}