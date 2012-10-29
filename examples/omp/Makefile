include ../../Makefile.inc

all: test_omp test_mpi


test_omp: test_omp.o
	$(CXX) $(CFLAGS) -fopenmp -I$(MPI_INCLUDE) -I../../include -L$(MPI_LIBDIR) -L../../lib -o $@ $^ -lm -lAMPI -lgomp
test_mpi: test_mpi.o
	$(CXX) $(CFLAGS) -I$(MPI_INCLUDE) -I../../include -L$(MPI_LIBDIR) -L../../lib -o $@ $^ -lm -lAMPI -lgomp
test_omp.o: test_omp.cpp
	$(CXX) $(CFLAGS) -fopenmp -I$(MPI_INCLUDE) -I../../include -c $^  
test_mpi.o: test_mpi.cpp
	$(CXX) $(CFLAGS) -I$(MPI_INCLUDE) -I../../include -c $^  
clean:
	rm -f ./*.o test_omp test_mpi 
