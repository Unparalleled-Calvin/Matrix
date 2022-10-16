#include "debug.h"

void Print(int M, int N, const double* mat) {
	for (int i = 0;i < M * N;i++) {
		printf("%f ", mat[i]);
	}
	printf("\n");
}