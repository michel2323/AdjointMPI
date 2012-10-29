#include "dcc_mpi.hpp"

extern int myid;
extern int numprocs;

void distribute(double * x, double * mpi_x, int &mpi_x_size) 
//$ad indep mpi_x x  
//$ad dep x mpi_x 
{
    double *buf;
    MPI_Status status[1];
    MPI_Request request[1];
    int target=0;
    int tag=0;
    double counter = 0;
    int counter2 = 0;
    int i = 0;
    int j = 0;
    dcc_new(buf,mpi_x_size);
    //request.r.aw = 0;
    if(myid == 0) {
	// set x
	counter = 0;
	// set mpi_x
	for(i = 0 ; i < mpi_x_size ; i=i+1) {
	    mpi_x[i] = x[i];
	}
	// send
	
	for(i = 1 ; i < numprocs ; i=i+1) {
	    for(j = 0 ; j < mpi_x_size ; j = j+1) {
	    counter2 = i*mpi_x_size+j;
		buf[j] = x[counter2];
	    }
	    target = i;

	    //MPI_Isend(&x[mpi_x_size*i], mpi_x_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, request);
	    MPI_Isend(buf, mpi_x_size, MPI_DOUBLE, target, tag, MPI_COMM_WORLD, request);
	    MPI_Wait(request, status);
	}
    }
    else {
	target = 0;
	MPI_Irecv(mpi_x, mpi_x_size, MPI_DOUBLE, target, tag, MPI_COMM_WORLD, request);
	MPI_Wait(request, status);
    }
    dcc_delete(buf);
}
