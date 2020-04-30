#include <algorithm>
#include "scan.h"

using namespace std;

int* A2;
int* B;
int* F;
int* e1;
int* e2;

void qsort(int* A, int start, int end) {
	if (start >= end-1) return;
	int pivot = A[start];

	cilk_for (int i = start; i < end; i++) A2[i] = A[i];
	
	cilk_for (int i = start; i < end; i++) {
		if (A2[i] < pivot) F[i] = 1; else F[i] = 0;
	}
	scan(F+start, B+start, e1+start, e2+start, end-start);
	
	cilk_for (int i = start+1; i < end; i++) {
		if (F[i]) A[start+B[i]-1] = A2[i];
	}
	
	int x = B[end-1];
	A[start+x] = pivot;

	cilk_for (int i = start+1; i < end; i++) {
		if (A2[i] >= pivot) F[i] = 1; else F[i] = 0;
	}
	scan(F+start, B+start, e1+start, e2+start, end-start);
	
	cilk_for (int i = start+1; i < end; i++) {
		if (F[i]) A[start+x+B[i]] = A2[i];
	}	
	
	cilk_spawn
	qsort(A, start, start+x);
	qsort(A, start+x+1, end);
	cilk_sync;
	
	return;
}