#include <mpi.h> 
#include <iostream>




int main (int argc, char *argv[]) {
    MPI_Init(&argc,&argv);
    int myid, num_threads;
#pragma omp parallel private(myid)
    {
        MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	std::cout << "Hello World from thread" << myid << std::endl;
#pragma omp barrier
	if ( myid == 0 ) {
            MPI_Comm_size(MPI_COMM_WORLD,&num_threads);
	    std::cout << "There are " << num_threads << " threads" << std::endl;
	}
    }
    MPI_Finalize();
    return 0;
}
