#include ../Makefile.inc

NPROC=8
PSIZE=50000000

CXXFLAGS=-O4 

all: pattern.out

pattern.out: pattern
	mpiexec -n $(NPROC) ./pattern $(PSIZE)

pattern: pattern.o 
	mpicxx $(CXXFLAGS) -o $@ pattern.o 

pattern.o: pattern.cpp
	mpicxx $(CXXFLAGS) -c -o $@ pattern.cpp


clean: 
	rm -f *.out* pattern *.o 
