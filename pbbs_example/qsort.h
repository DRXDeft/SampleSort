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

	parallel_for (start, end, [&](int i) {A2[i] = A[i];});
	
	parallel_for (start, end, [&](int i) {if (A2[i] < pivot) F[i] = 1; else F[i] = 0;});
	
	scan(F+start, B+start, e1+start, e2+start, end-start);
	
	parallel_for (start+1, end, [&](int i) {if (F[i]) A[start+B[i]-1] = A2[i];});
	
	int x = B[end-1];
	A[start+x] = pivot;

	parallel_for (start+1, end, [&](int i){
		if (A2[i] >= pivot) F[i] = 1; else F[i] = 0;
	});
	scan(F+start, B+start, e1+start, e2+start, end-start);
	
	parallel_for (start+1, end, [&](int i) {
		if (F[i]) A[start+x+B[i]] = A2[i];
	});
	
	par_do([&]() {qsort(A, start, start+x);}, 
		[&]() {qsort(A, start+x+1, end);});
	
	return;
}