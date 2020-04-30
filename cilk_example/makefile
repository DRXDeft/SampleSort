CC = g++
CFLAGS = -O2 -I../c++ -mcx16 -march=native -DCILK -fcilkplus -std=c++11

all:	reduce scan qsort

reduce:	reduce.cpp
	$(CC) $(CFLAGS) -DNDEBUG reduce.cpp -o reduce
	
scan:	scan.cpp scan.h
	$(CC) $(CFLAGS) -DNDEBUG scan.cpp -o scan

qsort:	quicksort.cpp qsort.h
	$(CC) $(CFLAGS) -DNDEBUG quicksort.cpp -o qsort

	
clean:
	rm -f reduce scan qsort
