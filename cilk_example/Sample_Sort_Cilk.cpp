#include "get_time.h"
#include <algorithm>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <random>
#include <string>
#include <vector>
using namespace std;
#define THRESHOLD_OF_TRANSPOSE 1000
#define THRESHOLD_OF_DISTRIBUTION 100
#define THRESHOLD_OF_SCAN 1000
#define THRESHOLD_OF_MERGE 500
int Pick[10000000], Sample[100001], Offset[100000001], InsertPointer[100000001];
inline uint32_t hash32(uint32_t a) {
    a = (a + 0x7ed55d16) + (a << 12);
    a = (a ^ 0xc761c23c) ^ (a >> 19);
    a = (a + 0x165667b1) + (a << 5);
    a = (a + 0xd3a2646c) ^ (a << 9);
    a = (a + 0xfd7046c5) + (a << 3);
    a = (a ^ 0xb55a4f09) ^ (a >> 16);
    if (a < 0)
        a = -a;
    return a;
}
void random(int *A, int n) {
    cilk_for(int i = 0; i < n; i++) A[i] = (hash32(i)) % (n << 2);
}
void almostsort(int *A, int n) {
    cilk_for(int i = 0; i < n / 2; i++) { A[i] = i + 1; }
    cilk_for(int i = n / 2; i < n; i++) { A[i] = i + 2; }
    A[n - 1] = n / 2 + 1;
}
void fewuniq(int *A, int n) {
    cilk_for(int i = 0; i < n; i++) { A[i] = (1 + i / (n / 4)) * (n / 4); }
    for (int i = 0; i < n; i++) {
        int index = (hash32(i)) % n;
        int temp = A[i];
        A[i] = A[index];
        A[index] = temp;
    }
}
void exp(int *A, int n) {
    int nrolls = n;                 // number of experiments
    int nstars = n;                 // maximum number of stars to distribute
    const int nintervals = sqrt(n); // number of intervals

    std::default_random_engine generator;
    std::exponential_distribution<double> distribution(sqrt(n));

    int p[nintervals] = {};

    for (int i = 0; i < nrolls; ++i) {
        double number = distribution(generator);
        if (number < 1.0)
            ++p[int(nintervals * number)];
    }

    int k = 0;
    for (int i = 0; i < nintervals; ++i) {
        for (int j = 0; j < p[i] * nstars / nrolls; j++)
            A[k++] = i;
        /*std::cout << float(i) /// nintervals << "-" << float(i + 1) /
        nintervals
                  << ": ";
        std::cout << std::string(p[i] * nstars / nrolls, '*') << std::endl;*/
    }
}
void normal(int *A, int m) {
    std::default_random_engine e;                   //引擎
    std::normal_distribution<double> n(m / 2, 1.5); //均值, 方差
    std::vector<unsigned> vals(m);
    for (std::size_t i = 0; i != m; ++i) {
        unsigned v = std::lround(n(e)); //取整-最近的整数
        if (v < vals.size())
            ++vals[v];
    }
    int k = 0;
    for (std::size_t j = 0; j != vals.size(); ++j) {
        for (int i = 0; i < vals[j]; i++)
            A[k++] = i;
    }
    // std::cout << j << " : " << vals[j] << std::string(vals[j], '*') <<
    // std::endl; int sum = std::accumulate(vals.begin(), vals.end(), 0); cout
    // << "sum = " << sum << endl;
}

void zipfan(int *List, int n) {
    int R = 2000;
    double A = 1.25;
    double C = 1.0;
    double pf[2000];
    double sum = 0.0;
    for (int i = 0; i < R; i++) {
        sum += C / pow((double)(i + 2),
                       A); //位置为i的频率,一共有r个(即秩), 累加求和
    }
    cilk_for(int i = 0; i < R; i++) {
        if (i == 0)
            pf[i] = C / pow((double)(i + 2), A) / sum;
        else
            pf[i] = pf[i - 1] + C / pow((double)(i + 2), A) / sum;
    }
    srand(time(00));
    //产生n个数
    for (int i = 0; i < n; i++) {
        List[i] = 0;
        double data = (double)rand() / RAND_MAX; //生成一个0~1的数
        while (
            data >
            pf[List[i]]) //找索引,直到找到一个比他小的值,那么对应的index就是随机数了
            List[i]++;
    }
}
int log2_up(int k) {
    int a = 0;
    while (k) {
        k = k >> 1;
        a++;
    }
    return a - 1;
}
bool Verification(int *A, int n) {
    for (int i = 0; i < n - 1; i++)
        if (A[i] > A[i + 1])
            return false;
    return true;
}
int binary_search(int *A, int start, int end, int k) {
    int p;
    if (start == end)
        return start;
    int mid = (start + end) / 2;
    if (A[mid] > k)
        return binary_search(A, start, mid - 1, k);
    else if (A[mid] <= k && A[mid + 1] > k)
        return mid;
    else
        return binary_search(A, mid + 1, end, k);
}
void Merge(int *A, int n, int *C, int start, int end) {
    int i = 0, j = 0;
    while ((i < n) && (j < end - start + 1)) {
        if (A[i] <= Sample[j]) {
            C[j]++;
            i++;
        } else {
            j++;
        }
    }
}

void Scan(int* In, int* Out, int* B, int* C, int n) {
	if (n==0) return;
    if (n < THRESHOLD_OF_SCAN){
        Out[0] = In[0];
        for (int i = 1; i<n;i++)
            Out[i] = Out[i-1] + In[i];
        return;
    }
	
	cilk_for(int i = 0; i < n/2; i++) 
		B[i] = In[2*i] + In[2*i+1];
	
	Scan(B, C, B+n/2, C+n/2, n/2);	
	Out[0] = In[0];
	
	cilk_for(int i = 1; i < n; i++) {
		if (i%2) Out[i] = C[i/2];
		else Out[i] = C[i/2-1] + In[i];
	}
}
void Transpose(int* C, int n, int x, int lx, int y, int ly) {
    if ((lx <= THRESHOLD_OF_TRANSPOSE) && (ly <= THRESHOLD_OF_TRANSPOSE)) {
        for (int i = x; i < x + lx; i++)
            for (int j = y; j < y + ly; j++)
                if (i < j) {
                    int tmp = C[(n * j) + i];
                    C[(n * j) + i] = C[(n * i) + j];
                    C[(n * i) + j] = tmp;
                    // cout <<"Transfer: "<<"x="<<x<<" y="<<y<<endl;
                }
    } else if (lx >= ly) {
        int midx = lx / 2;
        cilk_spawn Transpose(C, n, x, midx, y, ly);
        Transpose(C, n, x + midx, lx - midx, y, ly);
        cilk_sync;
    } else {
        int midy = ly / 2;
        cilk_spawn Transpose(C, n, x, lx, y, midy);
        Transpose(C, n, x, lx, y + midy, ly - midy);
        cilk_sync;
    }
}
void Sample_Sort(int *A, int *B, int *C, int *D, int n) {
    int bucket_quotient = 1;
    int buckets = sqrt(n);
    int bucket_size = n / buckets;

    //--------------------Step 1--------------------
    timer t1;
    t1.start();
    cilk_for(int i = 0; i < buckets; i++) {
        if (i == buckets - 1)
            sort(A + (buckets - 1) * bucket_size, A + n);
        else
            sort(A + i * bucket_size, A + (i + 1) * bucket_size);
    }
    t1.stop();
    cout << "Split:        " << t1.get_total() << endl;
    // for (int i = 0;i<n;i++) cout<< A[i]<<" "; cout<<endl;

    //---------------------Step 2--------------------
    timer t2;
    t2.start();
    int logn = log2_up(n);
    int random_pick = logn * bucket_quotient * buckets;
    //-------Randomly Pick cRootnLogn samples
    cilk_for(int i = 0; i < random_pick; i++) Pick[i] = A[hash32(i) % n];
    // for (int i = 0;i<random_pick;i++) cout<<Pick[i]<<" "; cout<<"\n";
    sort(Pick, Pick + random_pick);
    //-------Randomly Pick every cLogn samples
    // for (int i = 0, j=0; j<buckets-1;  j++,i+=bucket_quotient * logn)
    //  Sample[j] = Pick[i];
    int interval = bucket_quotient * logn;
    for (int i = 0; i < buckets - 1; i++)
        Sample[i] = Pick[i * interval];
    // sort(Sample, Sample+buckets-1);
    Sample[buckets - 1] = INT_MAX; // for (int i = 0; i < buckets; i++) cout<<
                                   // Sample[i]<<" "; cout<<"\n";
    t2.stop();
    cout << "Pivot:        " << t2.get_total() << endl;

    // Step 3
    //-----------------Step 3_1--------------------
    // First, get the count for each subarray in each bucket. I store them in C
    timer t3_1;
    t3_1.start();
    cilk_for(int i = 0; i < buckets; i++) {
        // cout<<"-----------Bucket "<<i/buckets<<"------------------\n"; for
        // (int j = i;j<i+buckets;j++) cout<<A[j]<<" "; cout<<"\n";

        if (i == buckets - 1)
            Merge(A + i * bucket_size, n - i * bucket_size, C + i * buckets, 0,
                  buckets - 1);
        else
            Merge(A + i * bucket_size, bucket_size, C + i * buckets, 0,
                  buckets - 1);
        // for (int j = i;j<i+buckets;j++) cout<<C[j]<<" "; cout<<"\n";
    }
    t3_1.stop();
    cout << "Merge:        " << t3_1.get_total() << endl;

    //-----------------Step 3_2--------------------
    // Then, transpose the array
    // for (int i = 0; i < buckets*buckets; i++) if ((i+1)%buckets == 0) cout<<
    // C[i]<<"\n"; else cout<< C[i]<<" "; cout<<"\n";
    timer t3_2;
    t3_2.start();
    Transpose(C, buckets, 0, buckets, 0, buckets);
    t3_2.stop();
    cout << "Transpose:    " << t3_2.get_total() << endl;
    // for (int i = 0; i < buckets*buckets; i++) if ((i+1)%buckets == 0) cout<<
    // C[i]<<"\n"; else cout<< C[i]<<" "; cout<<"\n";

    //-----------------Step 3_3--------------------
    // scan to compute the offsets
    timer t3_3;
    t3_3.start();
    //cilk_for (int i = 0; i< buckets; i++)
      //  Scan(C,D,InsertPointer,Offset,buckets);
    Scan(C,D,InsertPointer,Offset,buckets*buckets);
    //cilk_for (int i = 0; i < buckets * buckets; i++) D[i] -= C[i];
    cout << "Scan:         " << t3_3.get_total() << endl;
    // for (int i = 0; i<buckets; i++) cout<<InsertPointer[i]<<" "; cout<<"\n";

    //-----------------Step 3_4--------------------
    // Lastly, move each element to the corresponding bucket
    timer t3_4;
    t3_4.start();
    cilk_for (int i = 0; i < buckets; i++){
        int End;
        if (i == buckets - 1)
            End = n;
        else 
            End = bucket_size * (i+1);
        int pivot = 0, j = i * bucket_size;
        while (pivot < buckets && j < End){
            while (A[j] > Sample[pivot]) pivot++;
            int tmp = pivot*buckets+i;
            while (j<End && A[j] <= Sample[pivot]){
                B[D[tmp]++] = A[j++]; 
            }
        }
    }
    t3_4.stop();
    cout << "Distribution: " << t3_4.get_total() << endl;

    //-----------------Step 4--------------------
    timer t4;
    t4.start();
    cilk_for(int i = 0; i < buckets; i++) {
         if (i == buckets-1)
            sort(B+(D[i*buckets]-C[i*buckets]),B+n);
        else
            sort(B+(D[i*buckets]-C[i*buckets]),B+D[(i+1)*buckets]-C[(i+1)*buckets]);
    }
    t4.stop();
    cout << "Sort Buckets: " << t4.get_total() << endl;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        cout << "Usage: ./Sample_Sort_Cilk [num_elements] [distribution]" << endl;
        return 0;
    }
    int n = atoi(argv[1]);
    int *A = new int[n];
    int *B = new int[n];
    int *C = new int[n];
    int *D = new int[n];
    int distribution = atoi(argv[2]);
    if (distribution == 1)
        random(A, n);
    else if (distribution == 2)
        almostsort(A, n);
    else if (distribution == 3)
        zipfan(A, n);
    else if (distribution == 4)
        fewuniq(A, n);
    else if (distribution == 5)
        exp(A, n);
    else if (distribution == 6)
        normal(A, n);
    memset(B, 0, sizeof(B));    
    memset(C, 0, sizeof(C));
    memset(D, 0, sizeof(D));
    // cilk_for (int i = 0; i < n; i++) A[i] = n-i;
    // Verification(A,n);
    cout << "Data Generation is Complete, Start Sorting and Timing!\n";
    timer t;
    t.start();

    Sample_Sort(A, B, C, D, n);
    t.stop();
    cout << "\nTotal Time:   " << t.get_total()
         << "\nStarting Verification Now!\n";
    if (Verification(B, n))
        cout << "The result is correct!" << endl;
    else
        cout << "The result is incorrect!" << endl;

    return 0;
}
