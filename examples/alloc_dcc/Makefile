include ../Makefile.inc

all: t2_a1_alloc alloc 


alloc: driver.o alloc.o t2_a1_alloc.o a1_alloc.o dcc_alloc.o dcc_mpi.o
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -L$(AMPI_LIBDIR) -o $@ $^ -lm -lAMPI 
driver.o: driver.cpp
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
dcc_alloc.o: dcc_alloc.cpp
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
dcc_mpi.o: dcc_mpi.cpp
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
alloc.o: alloc.cpp
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
a1_alloc.o: a1_alloc.cpp 
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
t2_a1_alloc.o: t2_a1_alloc.cpp 
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
t2_a1_alloc:
	cp alloc.cpp alloc_dcc.cpp
	# no includes
	sed '/#include/g' alloc_dcc.cpp > alloc_dcc_tmp.cpp
	mv alloc_dcc_tmp.cpp alloc_dcc.cpp
	# no using namespace
	sed '/using/g' alloc_dcc.cpp > alloc_dcc_tmp.cpp
	mv alloc_dcc_tmp.cpp alloc_dcc.cpp
	# ready for DCC
	dcc -a -d 1 ./alloc_dcc.cpp
	rm alloc_dcc.cpp
	#cp a1_alloc_dcc.cpp a1_alloc_dcc_out.cpp
	# remove includes again, since we don't do any checkpointing
	sed '/#include/g' a1_alloc_dcc.cpp > a1_alloc_dcc_tmp.cpp
	mv a1_alloc_dcc_tmp.cpp a1_alloc_dcc.cpp
	# make a copy for second order
	cp a1_alloc_dcc.cpp a1_alloc_dcc_t2.cpp
	# add headers
	sed '1i\#include "dcc_alloc.hpp"' a1_alloc_dcc.cpp > a1_alloc_dcc_tmp.cpp
	mv a1_alloc_dcc_tmp.cpp a1_alloc_dcc.cpp
	# remove the forward bmode=2 call
	sed '/alloc(myres/g' a1_alloc_dcc.cpp > a1_alloc_dcc_tmp.cpp
	mv a1_alloc_dcc_tmp.cpp a1_alloc_dcc.cpp
	mv a1_alloc_dcc.cpp a1_alloc.cpp
	# 2nd order
	mv a1_alloc_dcc_t2.cpp a1_alloc_dcc.cpp
	#dcc a1_alloc_dcc.cpp < t2_a1_alloc.conf
	gcc -E -P -I. -C ./a1_alloc_dcc.cpp > ./a1_alloc_dcc_tmp.cpp
	mv a1_alloc_dcc_tmp.cpp a1_alloc_dcc.cpp
	dcc -t -d 2 ./a1_alloc_dcc.cpp
	rm a1_alloc_dcc.cpp
	#dcc alloc_dcc.cpp < t2_alloc.conf
	#cp a1_alloc_dcc.cpp a1_alloc_dcc_out.cpp
	# include the original file
	sed '/t2_alloc/g' t2_a1_alloc_dcc.cpp > t2_a1_alloc_dcc_tmp.cpp
	mv t2_a1_alloc_dcc_tmp.cpp t2_a1_alloc_dcc.cpp
	sed '1i\#include "dcc_alloc.hpp"' t2_a1_alloc_dcc.cpp > t2_a1_alloc_dcc_tmp.cpp
	mv t2_a1_alloc_dcc_tmp.cpp t2_a1_alloc_dcc.cpp
	#sed '1i\#include "alloc.cpp"' t2_a1_alloc_dcc.cpp > t2_a1_alloc_dcc_tmp.cpp
	#mv t2_a1_alloc_dcc_tmp.cpp t2_a1_alloc_dcc.cpp
	mv t2_a1_alloc_dcc.cpp t2_a1_alloc.cpp

	
clean:
	rm -f *.o *.dot *.out *.log *.inc alloc a1_alloc.cpp t2_a1_alloc.cpp 

