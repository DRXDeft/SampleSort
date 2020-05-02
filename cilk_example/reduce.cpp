#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <cstring>
#include <fstream>
#include "get_time.h"
using namespace std;
ofstream fout;
int n,leaves,bottom,nodes;
void Reduce_Size(int offset , int n, long long* A, long long* PrefixSum ){
	
	PrefixSum[offset] = A[offset];

	if (n==1){	return;	}

	cilk_for (int i = 0; i < n/2; i++)
		A[i+offset+n] = A[offset + 2*i] + A[offset+1 +2*i];
	Reduce_Size(offset + n, n/2, A , PrefixSum);
	cilk_for (int i = 1; i< n; i++){
		if (i%2)
			PrefixSum[offset+i] = PrefixSum[offset+n+i/2];
		else
			PrefixSum[offset+i] = PrefixSum[offset+n+i/2-1]+A[offset+i];
	}
	//----------------------------------------------------------------------------------------
	int j = offset + n, i=offset;
	cilk_for (int lIndex=0; lIndex<n/2; lIndex++){
		A[j] = A[i]+A[i+1];
		//fout <<"A["<<j<<"]=A["<<i<<"]+A["<<i+1<<"] = "<<A[j]<<" "<<A[i]<<" "<<A[i+1]<<endl;
		i+=2;
		j++;
	}
	Reduce_Size(offset + n, n/2, A , PrefixSum);
	j = offset + n;
	int lIndex=1;
	cilk_for (int k = offset+1; k<offset+n; k++){
		if (lIndex%2)
			PrefixSum[k] = PrefixSum[j];
		else
			PrefixSum[k] = PrefixSum[j++]+A[k];
		lIndex++;
	}
	
}
int main(int argc, char** argv){
    if (argc != 2) {
		cout << "Usage: ./scan [num_elements]" << endl;
		return 0;
	}
	n = atoi(argv[1]);
	
	long long* A = new long long[3 * n];
	long long* PrefixSum = new long long[3*n];
	
	cilk_for (int i = 0; i<n; i++)
		A[i] = i;
	//fout.open("prefixsum.txt");
	timer t; 
	t.start();
	Reduce_Size(0,n,A,PrefixSum);
	t.stop();
	//fout.open("prefixsum.txt");
	//for (int i = 0; i < 2*n; i++) cout <<i<<" " <<PrefixSum[i] <<endl;
	//fout.close();
	cout<<PrefixSum[0]<<" "<<PrefixSum[n/2]<<" "<<PrefixSum[n-1]<<endl;
	cout << "time: " <<fixed<<setprecision(4)<< t.get_total() << endl;

    return 0;
}