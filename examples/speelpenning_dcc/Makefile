include ../Makefile.inc

DCC=~/svn/dcc_with_pointer/src/dcc

all: t2_a1_speelpenning t2_a1_distribute speelpenning 


speelpenning: driver.o distribute.o dcc_mpi.o dcc_alloc.o speelpenning.o t2_a1_speelpenning.o t2_a1_distribute.o a1_speelpenning.o
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -L$(AMPI_LIBDIR) -o $@ $^ -lm -lAMPI 
dcc_mpi.o: dcc_mpi.cpp
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
dcc_alloc.o: dcc_alloc.cpp
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
driver.o: driver.cpp
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
distribute.o: distribute.cpp
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
speelpenning.o: speelpenning.cpp
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
a1_speelpenning.o: a1_speelpenning.cpp 
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
t2_a1_speelpenning.o: t2_a1_speelpenning.cpp 
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
t2_a1_distribute.o: t2_a1_distribute.cpp 
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
t2_a1_speelpenning:
	cp speelpenning.cpp speelpenning_dcc.cpp
	# ass =0 to all external value
	sed '/^extern/s/;/=0;/' speelpenning_dcc.cpp > speelpenning_dcc_tmp.cpp
	mv speelpenning_dcc_tmp.cpp speelpenning_dcc.cpp
	# remove "extern" keyword
	sed 's/extern//g' speelpenning_dcc.cpp > speelpenning_dcc_tmp.cpp
	mv speelpenning_dcc_tmp.cpp speelpenning_dcc.cpp
	# MPI_Request become an integer
	sed 's/MPI_Request/int/g' speelpenning_dcc.cpp > speelpenning_dcc_tmp.cpp
	mv speelpenning_dcc_tmp.cpp speelpenning_dcc.cpp
	# MPI_Status become an integer
	sed 's/MPI_Status/int/g' speelpenning_dcc.cpp > speelpenning_dcc_tmp.cpp
	mv speelpenning_dcc_tmp.cpp speelpenning_dcc.cpp
	# The constants MPI_DOUBLE and MPI_COMM_WORLD become a variable
	sed '1i\int ampi_double = 0;' speelpenning_dcc.cpp > speelpenning_dcc_tmp.cpp
	mv speelpenning_dcc_tmp.cpp speelpenning_dcc.cpp
	sed '1i\int ampi_comm_world = 0;' speelpenning_dcc.cpp > speelpenning_dcc_tmp.cpp
	mv speelpenning_dcc_tmp.cpp speelpenning_dcc.cpp
	sed 's/MPI_DOUBLE/ampi_double/g' speelpenning_dcc.cpp > speelpenning_dcc_tmp.cpp
	mv speelpenning_dcc_tmp.cpp speelpenning_dcc.cpp
	sed 's/MPI_COMM_WORLD/ampi_comm_world/g' speelpenning_dcc.cpp > speelpenning_dcc_tmp.cpp
	mv speelpenning_dcc_tmp.cpp speelpenning_dcc.cpp
	# no includes
	sed '/#include/g' speelpenning_dcc.cpp > speelpenning_dcc_tmp.cpp
	mv speelpenning_dcc_tmp.cpp speelpenning_dcc.cpp
	# no using namespace
	sed '/using/g' speelpenning_dcc.cpp > speelpenning_dcc_tmp.cpp
	mv speelpenning_dcc_tmp.cpp speelpenning_dcc.cpp
	# ready for DCC
	$(DCC) -a speelpenning_dcc.cpp
	rm speelpenning_dcc.cpp
	#cp a1_speelpenning_dcc.cpp a1_speelpenning_dcc_out.cpp
	# remove includes again, since we don't do any checkpointing
	sed '/#include/g' a1_speelpenning_dcc.cpp > a1_speelpenning_dcc_tmp.cpp
	mv a1_speelpenning_dcc_tmp.cpp a1_speelpenning_dcc.cpp
	# make a copy for second order
	cp a1_speelpenning_dcc.cpp a1_speelpenning_dcc_t2.cpp
	# add headers
	sed '1i\#include "dcc_mpi.hpp"' a1_speelpenning_dcc.cpp > a1_speelpenning_dcc_tmp.cpp
	mv a1_speelpenning_dcc_tmp.cpp a1_speelpenning_dcc.cpp
	# add extern
	sed 's/int myid=0/extern int myid/g' a1_speelpenning_dcc.cpp > a1_speelpenning_dcc_tmp.cpp
	mv a1_speelpenning_dcc_tmp.cpp a1_speelpenning_dcc.cpp
	sed 's/int numprocs=0/extern int numprocs/g' a1_speelpenning_dcc.cpp > a1_speelpenning_dcc_tmp.cpp
	mv a1_speelpenning_dcc_tmp.cpp a1_speelpenning_dcc.cpp

	sed 's/int ampi_double=0/extern int ampi_double/g' a1_speelpenning_dcc.cpp > a1_speelpenning_dcc_tmp.cpp
	mv a1_speelpenning_dcc_tmp.cpp a1_speelpenning_dcc.cpp
	sed 's/int ampi_comm_world=0/extern int ampi_comm_world/g' a1_speelpenning_dcc.cpp > a1_speelpenning_dcc_tmp.cpp
	mv a1_speelpenning_dcc_tmp.cpp a1_speelpenning_dcc.cpp
	# remove the forward bmode=2 call
	sed '/speelpenning(myres/g' a1_speelpenning_dcc.cpp > a1_speelpenning_dcc_tmp.cpp
	mv a1_speelpenning_dcc_tmp.cpp a1_speelpenning_dcc.cpp
	mv a1_speelpenning_dcc.cpp a1_speelpenning.cpp
	# 2nd order
	mv a1_speelpenning_dcc_t2.cpp a1_speelpenning_dcc.cpp
	$(DCC) -t -d 2 a1_speelpenning_dcc.cpp
	rm a1_speelpenning_dcc.cpp
	#dcc speelpenning_dcc.cpp < t2_speelpenning.conf
	#cp a1_speelpenning_dcc.cpp a1_speelpenning_dcc_out.cpp
	# include the original file
	sed '/t2_speelpenning/g' t2_a1_speelpenning_dcc.cpp > t2_a1_speelpenning_dcc_tmp.cpp
	mv t2_a1_speelpenning_dcc_tmp.cpp t2_a1_speelpenning_dcc.cpp
	sed '1i\#include "dcc_mpi.hpp"' t2_a1_speelpenning_dcc.cpp > t2_a1_speelpenning_dcc_tmp.cpp
	mv t2_a1_speelpenning_dcc_tmp.cpp t2_a1_speelpenning_dcc.cpp
	# add extern
	sed 's/int myid=0/extern int myid/g' t2_a1_speelpenning_dcc.cpp > t2_a1_speelpenning_dcc_tmp.cpp
	mv t2_a1_speelpenning_dcc_tmp.cpp t2_a1_speelpenning_dcc.cpp
	sed 's/int numprocs=0/extern int numprocs/g' t2_a1_speelpenning_dcc.cpp > t2_a1_speelpenning_dcc_tmp.cpp
	mv t2_a1_speelpenning_dcc_tmp.cpp t2_a1_speelpenning_dcc.cpp

	sed 's/int ampi_double=0/extern int ampi_double/g' t2_a1_speelpenning_dcc.cpp > t2_a1_speelpenning_dcc_tmp.cpp
	mv t2_a1_speelpenning_dcc_tmp.cpp t2_a1_speelpenning_dcc.cpp
	sed 's/int ampi_comm_world=0/extern int ampi_comm_world/g' t2_a1_speelpenning_dcc.cpp > t2_a1_speelpenning_dcc_tmp.cpp
	mv t2_a1_speelpenning_dcc_tmp.cpp t2_a1_speelpenning_dcc.cpp
	#sed '1i\#include "speelpenning.cpp"' t2_a1_speelpenning_dcc.cpp > t2_a1_speelpenning_dcc_tmp.cpp
	#mv t2_a1_speelpenning_dcc_tmp.cpp t2_a1_speelpenning_dcc.cpp
	mv t2_a1_speelpenning_dcc.cpp t2_a1_speelpenning.cpp

t2_a1_distribute:
	cp distribute.cpp distribute_dcc.cpp
	# ass =0 to all external value
	sed '/^extern/s/;/=0;/' distribute_dcc.cpp > distribute_dcc_tmp.cpp
	mv distribute_dcc_tmp.cpp distribute_dcc.cpp
	# remove "extern" keyword
	sed 's/extern//g' distribute_dcc.cpp > distribute_dcc_tmp.cpp
	mv distribute_dcc_tmp.cpp distribute_dcc.cpp
	# MPI_Request become an integer
	sed 's/MPI_Request/int/g' distribute_dcc.cpp > distribute_dcc_tmp.cpp
	mv distribute_dcc_tmp.cpp distribute_dcc.cpp
	# MPI_Status become an integer
	sed 's/MPI_Status/int/g' distribute_dcc.cpp > distribute_dcc_tmp.cpp
	mv distribute_dcc_tmp.cpp distribute_dcc.cpp
	# The constants MPI_DOUBLE and MPI_COMM_WORLD become a variable
	sed '1i\int ampi_double = 0;' distribute_dcc.cpp > distribute_dcc_tmp.cpp
	mv distribute_dcc_tmp.cpp distribute_dcc.cpp
	sed '1i\int ampi_comm_world = 0;' distribute_dcc.cpp > distribute_dcc_tmp.cpp
	mv distribute_dcc_tmp.cpp distribute_dcc.cpp
	sed 's/MPI_DOUBLE/ampi_double/g' distribute_dcc.cpp > distribute_dcc_tmp.cpp
	mv distribute_dcc_tmp.cpp distribute_dcc.cpp
	sed 's/MPI_COMM_WORLD/ampi_comm_world/g' distribute_dcc.cpp > distribute_dcc_tmp.cpp
	mv distribute_dcc_tmp.cpp distribute_dcc.cpp
	# no includes
	sed '/#include/g' distribute_dcc.cpp > distribute_dcc_tmp.cpp
	mv distribute_dcc_tmp.cpp distribute_dcc.cpp
	# no using namespace
	sed '/using/g' distribute_dcc.cpp > distribute_dcc_tmp.cpp
	mv distribute_dcc_tmp.cpp distribute_dcc.cpp
	# ready for DCC
	$(DCC) -a distribute_dcc.cpp
	rm distribute_dcc.cpp
	#cp a1_distribute_dcc.cpp a1_distribute_dcc_out.cpp
	# remove includes again, since we don't do any checkpointing
	sed '/#include/g' a1_distribute_dcc.cpp > a1_distribute_dcc_tmp.cpp
	mv a1_distribute_dcc_tmp.cpp a1_distribute_dcc.cpp
	# make a copy for second order
	cp a1_distribute_dcc.cpp a1_distribute_dcc_t2.cpp
	# add headers
	sed '1i\#include "dcc_mpi.hpp"' a1_distribute_dcc.cpp > a1_distribute_dcc_tmp.cpp
	mv a1_distribute_dcc_tmp.cpp a1_distribute_dcc.cpp
	# add extern
	sed 's/int myid=0/extern int myid/g' a1_distribute_dcc.cpp > a1_distribute_dcc_tmp.cpp
	mv a1_distribute_dcc_tmp.cpp a1_distribute_dcc.cpp
	sed 's/int numprocs=0/extern int numprocs/g' a1_distribute_dcc.cpp > a1_distribute_dcc_tmp.cpp
	mv a1_distribute_dcc_tmp.cpp a1_distribute_dcc.cpp

	sed 's/int ampi_double=0/extern int ampi_double/g' a1_distribute_dcc.cpp > a1_distribute_dcc_tmp.cpp
	mv a1_distribute_dcc_tmp.cpp a1_distribute_dcc.cpp
	sed 's/int ampi_comm_world=0/extern int ampi_comm_world/g' a1_distribute_dcc.cpp > a1_distribute_dcc_tmp.cpp
	mv a1_distribute_dcc_tmp.cpp a1_distribute_dcc.cpp
	# remove the forward bmode=2 call
	sed '/distribute(myres/g' a1_distribute_dcc.cpp > a1_distribute_dcc_tmp.cpp
	mv a1_distribute_dcc_tmp.cpp a1_distribute_dcc.cpp
	mv a1_distribute_dcc.cpp a1_distribute.cpp
	# 2nd order
	mv a1_distribute_dcc_t2.cpp a1_distribute_dcc.cpp
	$(DCC) -t -d 2 a1_distribute_dcc.cpp
	rm a1_distribute_dcc.cpp
	#dcc distribute_dcc.cpp < t2_distribute.conf
	#cp a1_distribute_dcc.cpp a1_distribute_dcc_out.cpp
	# include the original file
	sed '/t2_distribute/g' t2_a1_distribute_dcc.cpp > t2_a1_distribute_dcc_tmp.cpp
	mv t2_a1_distribute_dcc_tmp.cpp t2_a1_distribute_dcc.cpp
	sed '1i\#include "dcc_mpi.hpp"' t2_a1_distribute_dcc.cpp > t2_a1_distribute_dcc_tmp.cpp
	mv t2_a1_distribute_dcc_tmp.cpp t2_a1_distribute_dcc.cpp
	# add extern
	sed 's/int myid=0/extern int myid/g' t2_a1_distribute_dcc.cpp > t2_a1_distribute_dcc_tmp.cpp
	mv t2_a1_distribute_dcc_tmp.cpp t2_a1_distribute_dcc.cpp
	sed 's/int numprocs=0/extern int numprocs/g' t2_a1_distribute_dcc.cpp > t2_a1_distribute_dcc_tmp.cpp
	mv t2_a1_distribute_dcc_tmp.cpp t2_a1_distribute_dcc.cpp

	sed 's/int ampi_double=0/extern int ampi_double/g' t2_a1_distribute_dcc.cpp > t2_a1_distribute_dcc_tmp.cpp
	mv t2_a1_distribute_dcc_tmp.cpp t2_a1_distribute_dcc.cpp
	sed 's/int ampi_comm_world=0/extern int ampi_comm_world/g' t2_a1_distribute_dcc.cpp > t2_a1_distribute_dcc_tmp.cpp
	mv t2_a1_distribute_dcc_tmp.cpp t2_a1_distribute_dcc.cpp
	#sed '1i\#include "distribute.cpp"' t2_a1_distribute_dcc.cpp > t2_a1_distribute_dcc_tmp.cpp
	#mv t2_a1_distribute_dcc_tmp.cpp t2_a1_distribute_dcc.cpp
	mv t2_a1_distribute_dcc.cpp t2_a1_distribute.cpp
	
clean:
	rm -f *.o *.dot *.out *.log *.inc speelpenning a1_speelpenning.cpp t2_a1_speelpenning.cpp 
	rm -f a1_distribute.cpp t2_a1_distribute.cpp 
	rm -f speelpenning_dcc.cpp

