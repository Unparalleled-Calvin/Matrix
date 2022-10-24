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

double* RowMajor2Zmorton(int M, int N, int ld, double* mat, double* mat_, int base) {
	if (mat_ == nullptr) {
		mat_ = new double[M * N];
	}
	int L = max(M, N);
	int X = L / 2;
	if (L <= base) {
		int pos = 0;
		for (int i = 0;i < base; i++) {
			for (int j = 0; j < base;j++) {
				mat_[pos++] = mat[i * ld + j];
			}
		}
	}
	else if(L == M) {
		RowMajor2Zmorton(X, N, ld, mat, mat_, base);
		RowMajor2Zmorton(X, N, ld, mat + X * ld, mat_ + X * N, base);
	}
	else if (L == N) {
		RowMajor2Zmorton(M, X, ld, mat, mat_, base);
		RowMajor2Zmorton(M, X, ld, mat + X, mat_ + M * X, base);
	}
	return mat_;
}

double* Zmorton2RowMajor(int M, int N, int ld, double* mat, double* mat_, int base) {
	if (mat_ == nullptr) {
		mat_ = new double[M * N];
	}
	int L = max(M, N);
	int X = L / 2;
	if (L <= base) {
		int pos = 0;
		for (int i = 0;i < base;i++) {
			for (int j = 0;j < base;j++) {
				mat_[i * ld + j] = mat[pos++];
			}
		}
	}
	else if (L == M) {
		Zmorton2RowMajor(X, N, ld, mat, mat_, base);
		Zmorton2RowMajor(X, N, ld, mat + X * N, mat_ + X * ld, base);
	}
	else if (L == N) {
		Zmorton2RowMajor(M, X, ld, mat, mat_, base);
		Zmorton2RowMajor(M, X, ld, mat + M * X, mat_ + X, base);
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
	int step1 = lda * BLOCK_SIZE, step2 = ldb * BLOCK_SIZE, step3 = ldc * BLOCK_SIZE;
	for (int ih = 0, ih_ldc = 0, ih_lda = 0;ih < M;ih += BLOCK_SIZE, ih_lda += step1, ih_ldc += step3) {
		for (int jh = 0;jh < N;jh += BLOCK_SIZE) {
			for (int kh = 0, kh_ldb = 0;kh < K;kh += BLOCK_SIZE, kh_ldb += step2) {
				for (int il = 0, il_ldc = 0, il_lda = 0;il < BLOCK_SIZE;il++, il_lda += lda, il_ldc += ldc) {
					for (int jl = 0;jl < BLOCK_SIZE;jl++) {
						for (int kl = 0, kl_ldb = 0;kl < BLOCK_SIZE;kl++, kl_ldb += ldb) {
							C[ih_ldc + il_ldc + jh + jl] += alpha * A[ih_lda + il_lda + kh + kl] * B[kh_ldb + kl_ldb + jh + jl];
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
	int step1 = BLOCK_SIZE * lda, step2 = BLOCK_SIZE * ldb, step3 = BLOCK_SIZE * ldc, step4 = BLOCK_SIZE * BLOCK_SIZE;
	for (int ih = 0, ih_ldc = 0, ih_lda = 0;ih < M;ih += BLOCK_SIZE, ih_ldc += step3, ih_lda += step1) {
		for (int jh = 0, jh_block_size = 0;jh < N;jh += BLOCK_SIZE, jh_block_size+=step4) {
			for (int kh = 0, kh_ldb = 0, kh_block_size = 0;kh < K;kh += BLOCK_SIZE, kh_ldb += step2, kh_block_size += step4) {
				for (int il = 0, il_block_size = 0;il < BLOCK_SIZE;il++, il_block_size += BLOCK_SIZE) {
					for (int jl = 0;jl < BLOCK_SIZE;jl++) {
						for (int kl = 0, kl_block_size = 0;kl < BLOCK_SIZE;kl++, kl_block_size += BLOCK_SIZE) {
							C[ih_ldc + jh_block_size + il_block_size + jl] += alpha * A[ih_lda + kh_block_size + il_block_size + kl] * B[kh_ldb + jh_block_size + kl_block_size + jl];
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
	else if (L == M) {
		RecursionRowMajorOrderingPro(X, N, K, A, B, C, alpha, lda, ldb, ldc);
		RecursionRowMajorOrderingPro(X, N, K, A + X * lda, B, C + X * ldc, alpha, lda, ldb, ldc);
	}
	else if (L == K) {
		RecursionRowMajorOrderingPro(M, N, X, A, B, C, alpha, lda, ldb, ldc);
		RecursionRowMajorOrderingPro(M, N, X, A + X, B + X * ldb, C, alpha, lda, ldb, ldc);
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
	else if (L == M) {
		RecursionRowMajorBlockingPro(X, N, K, A, B, C, alpha, lda, ldb, ldc);
		RecursionRowMajorBlockingPro(X, N, K, A + X * lda, B, C + X * ldc, alpha, lda, ldb, ldc);
	}
	else if (L == K) {
		RecursionRowMajorBlockingPro(M, N, X, A, B, C, alpha, lda, ldb, ldc);
		RecursionRowMajorBlockingPro(M, N, X, A + X, B + X * ldb, C, alpha, lda, ldb, ldc);
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

// RecursionRowMajorPacking

void RecursionRowMajorPackingPre(FUNC_PARAM_PRE) {
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

// 仅当MNK均为2幂次时生效
#define pos(mat, i, j) (mat)[(i) * BLOCK_SIZE + (j)]

void RecursionRowMajorPackingPro(FUNC_PARAM_PRO) { // 递归结构 + 分块 + 块与块间连续
	int L = max(max(M, N), K);
	int X = L / 2;
	if (L <= BLOCK_SIZE) {
		for (int i = 0;i < M;i++) {
			for (int j = 0; j < N;j++) {
				for (int k = 0; k < K;k++) {
					pos(C, i, j) += alpha * pos(A, i, k) * pos(B, k, j);
				}
			}
		}
	}
	else if (L == M) {
		RecursionRowMajorPackingPro(X, N, K, A, B, C, alpha, lda, ldb, ldc);
		RecursionRowMajorPackingPro(X, N, K, A + X * lda, B, C + X * ldc, alpha, lda, ldb, ldc);
	}
	else if (L == K) {
		RecursionRowMajorPackingPro(M, N, X, A, B, C, alpha, lda, ldb, ldc);
		RecursionRowMajorPackingPro(M, N, X, A + X * BLOCK_SIZE, B + X * ldb, C, alpha, lda, ldb, ldc);
	}
	else if (L == N) {
		RecursionRowMajorPackingPro(M, X, K, A, B, C, alpha, lda, ldb, ldc);
		RecursionRowMajorPackingPro(M, X, K, A, B + X * BLOCK_SIZE, C + X * BLOCK_SIZE, alpha, lda, ldb, ldc);
	}
}

void RecursionRowMajorPackingPost(FUNC_PARAM_POST) {
	double* C__ = Packing2RowMajor(M, N, C_);
	double* C___ = RowMajor2ColMajor(M, N, C__);
	memcpy(C, C___, M * N * sizeof(double));
	delete[] C__;
	delete[] C___;
}

// RecursionZmortonOrdering

void RecursionZmortonOrderingPre(FUNC_PARAM_PRE) {
	double* A_ = ColMajor2RowMajor(M, K, A);
	double* B_ = ColMajor2RowMajor(K, N, B);
	double* C_ = ColMajor2RowMajor(M, N, C);
	*pA = RowMajor2Zmorton(M, K, K, A_, nullptr, 1);
	*pB = RowMajor2Zmorton(K, N, N, B_, nullptr, 1);
	*pC = RowMajor2Zmorton(M, N, N, C_, nullptr, 1);
	delete[] A_;
	delete[] B_;
	delete[] C_;
	matScale(M, N, *pC, beta);
}

void RecursionZmortonOrderingPro(FUNC_PARAM_PRO) { // 递归结构 + 递归到单个元素 + 所有元素递归序排列
	int L = max(max(M, N), K);
	int X = L / 2;
	if (L <= 1) {
		C[0] += alpha * A[0] * B[0];
	}
	else if (L == M) {
		RecursionZmortonOrderingPro(X, N, K, A, B, C, alpha, lda, ldb, ldc);
		RecursionZmortonOrderingPro(X, N, K, A + X * K, B, C + X * N, alpha, lda, ldb, ldc);
	}
	else if (L == K) {
		RecursionZmortonOrderingPro(M, N, X, A, B, C, alpha, lda, ldb, ldc);
		RecursionZmortonOrderingPro(M, N, X, A + M * X, B + X * N, C, alpha, lda, ldb, ldc);
	}
	else if (L == N) {
		RecursionZmortonOrderingPro(M, X, K, A, B, C, alpha, lda, ldb, ldc);
		RecursionZmortonOrderingPro(M, X, K, A, B + K * X, C + M * X, alpha, lda, ldb, ldc);
	}
}

void RecursionZmortonOrderingPost(FUNC_PARAM_POST) {
	double* C__ = Zmorton2RowMajor(M, N, N, C_, nullptr, 1);
	double* C___ = RowMajor2ColMajor(M, N, C__);
	memcpy(C, C___, M * N * sizeof(double));
	delete[] C__;
	delete[] C___;
}

// RecursionZmortonPacking

void RecursionZmortonPackingPre(FUNC_PARAM_PRE) {
	double* A_ = ColMajor2RowMajor(M, K, A);
	double* B_ = ColMajor2RowMajor(K, N, B);
	double* C_ = ColMajor2RowMajor(M, N, C);
	*pA = RowMajor2Zmorton(M, K, K, A_, nullptr, BLOCK_SIZE);
	*pB = RowMajor2Zmorton(K, N, N, B_, nullptr, BLOCK_SIZE);
	*pC = RowMajor2Zmorton(M, N, N, C_, nullptr, BLOCK_SIZE);
	delete[] A_;
	delete[] B_;
	delete[] C_;
	matScale(M, N, *pC, beta);
}

// 仅当MNK均为2幂次时生效
#define pos(mat, i, j, ld) (mat)[(i) * (ld) + (j)]

void RecursionZmortonPackingPro(FUNC_PARAM_PRO) { // 递归结构 + 递归到单个元素 + 所有元素递归序排列
	int L = max(max(M, N), K);
	int X = L / 2;
	if (L <= BLOCK_SIZE) {
		for (int i = 0;i < BLOCK_SIZE;i++) {
			for (int j = 0;j < BLOCK_SIZE;j++) {
				for (int k = 0;k < BLOCK_SIZE;k++) {
					pos(C, i, j, BLOCK_SIZE) += alpha * pos(A, i, k, BLOCK_SIZE) * pos(B, k, j, BLOCK_SIZE);
				}
			}
		}
	}
	else if (L == M) {
		RecursionZmortonPackingPro(X, N, K, A, B, C, alpha, lda, ldb, ldc);
		RecursionZmortonPackingPro(X, N, K, A + X * K, B, C + X * N, alpha, lda, ldb, ldc);
	}
	else if (L == K) {
		RecursionZmortonPackingPro(M, N, X, A, B, C, alpha, lda, ldb, ldc);
		RecursionZmortonPackingPro(M, N, X, A + M * X, B + X * N, C, alpha, lda, ldb, ldc);
	}
	else if (L == N) {
		RecursionZmortonPackingPro(M, X, K, A, B, C, alpha, lda, ldb, ldc);
		RecursionZmortonPackingPro(M, X, K, A, B + K * X, C + M * X, alpha, lda, ldb, ldc);
	}
}

void RecursionZmortonPackingPost(FUNC_PARAM_POST) {
	double* C__ = Zmorton2RowMajor(M, N, N, C_, nullptr, BLOCK_SIZE);
	double* C___ = RowMajor2ColMajor(M, N, C__);
	memcpy(C, C___, M * N * sizeof(double));
	delete[] C__;
	delete[] C___;
}