#include <mpi.h> 
#include <iostream>




int main (int argc, char *argv[]) {
    MPI_Init(&argc,&argv);
    int num_threads=4;
    int myid, nthreads;
#pragma omp parallel private(th_id)
    {
        MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	std::cout << "Hello World from thread" << th_id << std::endl;
#pragma omp barrier
	if ( th_id == 0 ) {
            MPI_Comm_size(MPI_COMM_WORLD,&num_threads);
	    std::cout << "There are " << nthreads << " threads" << std::endl;
	}
    }
    return 0;
}
