#include "debug.h"

bool isPowerOf2(int x) {
	return (x) == (x & -x);
}

void Print(int M, int N, const double* mat) {
	for (int i = 0;i < M * N;i++) {
		printf("%f ", mat[i]);
	}
	printf("\n");
}