CC = g++
CPFLAGS = -DHOMEGROWN
LPFLAGS = -pthread
CFLAGS = -mcx16 -O3 -std=c++11 -march=native -Wall $(CPFLAGS) 

all:	reduce scan qsort

reduce:	reduce.cpp
	$(CC) $(CFLAGS) reduce.cpp -o reduce $(LPFLAGS)
	
scan: scan.cpp scan.h
	$(CC) $(CFLAGS) scan.cpp -o scan $(LPFLAGS)
	
qsort:	quicksort.cpp qsort.h
	$(CC) $(CFLAGS) quicksort.cpp -o qsort $(LPFLAGS)

clean:
	rm -f reduce scan qsort
