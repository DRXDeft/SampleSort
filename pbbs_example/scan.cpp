#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include "pbbslib/utilities.h"
#include "pbbslib/get_time.h"
#include "scan.h"
using namespace std;

int main(int argc, char** argv) {
	if (argc != 2) {
		cout << "Usage: ./scan [num_elements]" << endl;
		return 0;
	}
	int n = atoi(argv[1]);
	cout << n << endl;
	int* A = new int[n];
	int* B = new int[n];
	int* e1 = new int[n];
	int* e2 = new int[n];
	parallel_for (0, n, [&](int i) {A[i] = i;});
	
	timer t; t.start();
	scan(A, B, e1, e2, n);
	t.stop(); cout << "time: " << t.get_total() << endl;
	
	//for (int i = 0; i < n; i++) cout << B[i] << " ";
	//cout << endl;
	
	return 0;
}