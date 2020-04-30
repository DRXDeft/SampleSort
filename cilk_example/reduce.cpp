#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include "get_time.h"
using namespace std;

int th = 1;

int reduce(int* A, int n) {
	if (n <= th) {
		int ret = 0;
		for (int i = 0; i < n; i++) {
			ret += A[i];
		}
		return ret;
	}
	int L, R;
	L = cilk_spawn reduce(A, n/2);
	R = reduce(A+n/2, n-n/2);
	cilk_sync;
	return L+R;
}

int main(int argc, char** argv) {
	if (argc < 2) {
		cout << "Usage: ./reduce [num_elements] [threshold=1]" << endl;
		return 0;
	}
	int n = atoi(argv[1]);
	if (argc == 3) {
		th = atoi(argv[2]);}
	int* A = new int[n];
	cilk_for (int i = 0; i < n; i++) A[i] = i;
	timer t0; t0.start();
	int ans = 0;
	for (int i = 0; i < n; i++) ans += A[i];
	t0.stop();
	cout << ans << ". sequential running time = " << t0.get_total() << endl;
	timer t; t.start();
	int x = reduce(A, n);
	t.stop();
	cout << "time: " << t.get_total() << endl;
	cout << x << endl;
	
	return 0;
}
