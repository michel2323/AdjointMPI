include ../Makefile.inc

all: d2_b1_alloc alloc 


alloc: driver.o alloc.o d2_b1_alloc.o b1_alloc.o dcc_alloc.o dcc_mpi.o
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -L$(AMPI_LIBDIR) -o $@ $^ -lm -lAMPI 
driver.o: driver.cpp
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
dcc_alloc.o: dcc_alloc.cpp
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
dcc_mpi.o: dcc_mpi.cpp
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
alloc.o: alloc.cpp
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
b1_alloc.o: b1_alloc.cpp 
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
d2_b1_alloc.o: d2_b1_alloc.cpp 
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
d2_b1_alloc:
	cp alloc.cpp alloc_dcc.cpp
	# no includes
	sed '/#include/g' alloc_dcc.cpp > alloc_dcc_tmp.cpp
	mv alloc_dcc_tmp.cpp alloc_dcc.cpp
	# no using namespace
	sed '/using/g' alloc_dcc.cpp > alloc_dcc_tmp.cpp
	mv alloc_dcc_tmp.cpp alloc_dcc.cpp
	# ready for DCC
	dcc alloc_dcc.cpp < b1_alloc.conf
	rm alloc_dcc.cpp
	#cp b1_alloc_dcc.cpp b1_alloc_dcc_out.cpp
	# remove includes again, since we don't do any checkpointing
	sed '/#include/g' b1_alloc_dcc.cpp > b1_alloc_dcc_tmp.cpp
	mv b1_alloc_dcc_tmp.cpp b1_alloc_dcc.cpp
	# make a copy for second order
	cp b1_alloc_dcc.cpp b1_alloc_dcc_d2.cpp
	# add headers
	sed '1i\#include "dcc_alloc.hpp"' b1_alloc_dcc.cpp > b1_alloc_dcc_tmp.cpp
	mv b1_alloc_dcc_tmp.cpp b1_alloc_dcc.cpp
	# remove the forward bmode=2 call
	sed '/alloc(myres/g' b1_alloc_dcc.cpp > b1_alloc_dcc_tmp.cpp
	mv b1_alloc_dcc_tmp.cpp b1_alloc_dcc.cpp
	mv b1_alloc_dcc.cpp b1_alloc.cpp
	# 2nd order
	mv b1_alloc_dcc_d2.cpp b1_alloc_dcc.cpp
	dcc b1_alloc_dcc.cpp < d2_b1_alloc.conf
	rm b1_alloc_dcc.cpp
	#dcc alloc_dcc.cpp < d2_alloc.conf
	#cp b1_alloc_dcc.cpp b1_alloc_dcc_out.cpp
	# include the original file
	sed '/d2_alloc/g' d2_b1_alloc_dcc.cpp > d2_b1_alloc_dcc_tmp.cpp
	mv d2_b1_alloc_dcc_tmp.cpp d2_b1_alloc_dcc.cpp
	sed '1i\#include "dcc_alloc.hpp"' d2_b1_alloc_dcc.cpp > d2_b1_alloc_dcc_tmp.cpp
	mv d2_b1_alloc_dcc_tmp.cpp d2_b1_alloc_dcc.cpp
	#sed '1i\#include "alloc.cpp"' d2_b1_alloc_dcc.cpp > d2_b1_alloc_dcc_tmp.cpp
	#mv d2_b1_alloc_dcc_tmp.cpp d2_b1_alloc_dcc.cpp
	mv d2_b1_alloc_dcc.cpp d2_b1_alloc.cpp

	
clean:
	rm -f *.o *.dot *.out *.log *.inc alloc b1_alloc.cpp d2_b1_alloc.cpp 

