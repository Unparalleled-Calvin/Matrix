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
}

#define pos(mat, i, j, ld) (mat)[(i) * (ld) + (j)]

void LoopRowMajorOrderingPro(FUNC_PARAM_PRO) { // 三重循环 + 所有元素按行排序
	for (int i = 0;i < M;i++) {
		for (int j = 0;j < N;j++) {
			pos(C, i, j, ldc) *= beta;
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
}

#define pos(mat, ih, il, jh, jl, ld) (mat)[((ih) + (il)) * (ld) + (jh) + (jl)]

void LoopRowMajorBlockingPro(FUNC_PARAM_PRO) { // 三重循环 + 显式分块
	for (int ih = 0;ih < M;ih += BLOCK_SIZE) {
		for (int jh = 0;jh < N;jh += BLOCK_SIZE) {
			for (int kh = 0;kh < K;kh += BLOCK_SIZE) {
				for (int il = 0;il < BLOCK_SIZE && ih + il < M;il++) {
					for (int jl = 0;jl < BLOCK_SIZE && jh + jl < N;jl++) {
						pos(C, ih, il, jh, jl, ldc) *= kh ? 1 : beta; // 控制使beta只算一次
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
}

// 仅当MNK均整除BLOCK_SIZE时生效
#define pos(mat, ih, il, jh, jl, ld) (mat)[(ih) * (ld) + (jh) * (BLOCK_SIZE) + (il) * (BLOCK_SIZE) + (jl)]

void LoopRowMajorPackingPro(FUNC_PARAM_PRO) { // 三重循环 + 显式分块 + 小块连续
	if (M % BLOCK_SIZE || N % BLOCK_SIZE || K % BLOCK_SIZE) {
		printf("SIZE ERROR!\n");
		exit(-1);
	}
	for (int ih = 0;ih < M;ih += BLOCK_SIZE) {
		for (int jh = 0;jh < N;jh += BLOCK_SIZE) {
			for (int kh = 0;kh < K;kh += BLOCK_SIZE) {
				for (int il = 0;il < BLOCK_SIZE && ih + il < M;il++) {
					for (int jl = 0;jl < BLOCK_SIZE && jh + jl < N;jl++) {
						pos(C, ih, il, jh, jl, ldc) *= kh ? 1 : beta; // 控制使beta只算一次
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