#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include "pbbslib/get_time.h"
#include "pbbslib/utilities.h"

using namespace std;

void reduce(int* A, int n, int& ret) {
	if (n == 1) {
		ret = A[0];
		return;
	}
	int L, R;
	par_do([&] () {reduce(A, n/2, L);}, 
	   		[&] () {reduce(A+n/2, n-n/2, R);});
	ret = L+R;
}

int main(int argc, char** argv) {
	if (argc != 2) {
		cout << "Usage: ./reduce [num_elements]" << endl;
		return 0;
	}
	int n = atoi(argv[1]);
	cout << n << endl;
	int* A = new int[n];
	parallel_for (0, n, [&] (int i) {A[i] = i;}); 
	int x;
	timer t; t.start();
	reduce(A, n, x);
	t.stop();
	cout << "time: " << t.get_total() << endl;
	cout << x << endl;
	
	return 0;
}
