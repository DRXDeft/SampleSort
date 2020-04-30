#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include "get_time.h"
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
	cilk_for (int i = 0; i < n; i++) A[i] = i;
	
	timer t; t.start();
	scan(A, B, e1, e2, n);
	t.stop(); cout << "time: " << t.get_total() << endl;
	
	//for (int i = 0; i < n; i++) cout << B[i] << " ";
	//cout << endl;
	
	return 0;
}