CC = g++
CFLAGS = -O2 -I../c++ -mcx16 -march=native -DCILK -fcilkplus -std=c++11

all:	Sample_Sort

reduce:	reduce.cpp
	$(CC) $(CFLAGS) -DNDEBUG reduce.cpp -o reduce
	
scan:	scan.cpp scan.h
	$(CC) $(CFLAGS) -DNDEBUG scan.cpp -o scan

qsort:	quicksort.cpp qsort.h
	$(CC) $(CFLAGS) -DNDEBUG quicksort.cpp -o qsort
Sample_Sort: Sample_Sort.cpp
	$(CC) $(CFLAGS) -DNDEBUG Sample_Sort.cpp -o Sample_Sort
clean:
	rm -f Sample_Sort
