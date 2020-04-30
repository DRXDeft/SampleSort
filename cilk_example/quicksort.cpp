#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include "get_time.h"
#include "qsort.h"
using namespace std;

inline uint32_t hash32(uint32_t a) {
	a = (a+0x7ed55d16) + (a<<12);
	a = (a^0xc761c23c) ^ (a>>19);
	a = (a+0x165667b1) + (a<<5);
	a = (a+0xd3a2646c) ^ (a<<9);
	a = (a+0xfd7046c5) + (a<<3);
	a = (a^0xb55a4f09) ^ (a>>16);
	if (a<0) a = -a;
	return a;
}

int main(int argc, char** argv) {
	if (argc != 2) {
		cout << "Usage: ./qsort [num_elements]" << endl;
		return 0;
	}
	int n = atoi(argv[1]);
	int* A = new int[n];
	cilk_for (int i = 0; i < n; i++) A[i] = (hash32(i)) % (n*2);
	
	A2 = new int[n];
	B = new int[n];
	F = new int[n];
	e1 = new int[n];
	e2 = new int[n];
	
	timer t; t.start();
	qsort(A, 0, n);
	t.stop(); cout << "time: " << t.get_total() << endl;
	
	//for (int i = 0; i < n; i++) cout << A[i] << " ";
	//cout << endl;
	
	return 0;
}