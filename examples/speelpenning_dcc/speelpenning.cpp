#include "dcc_mpi.hpp"
#include <iostream>


extern int myid;
extern int numprocs;

void speelpenning(double *myres, int &mpi_x_size, double *mpi_x) 
//$ad indep mpi_x 
//$ad dep myres 
{
    int i=0;
    int nelements=1;
    double *buf;
    MPI_Request request[1]; // MPI Request
    MPI_Status status[1]; // MPI Status
    int target=0;
    int tag=0;
    dcc_new(buf,nelements);
    if(myid == 0) {
	myres[0] = 1;
	for(i=0; i < mpi_x_size; i=i+1) {
	    myres[0] = myres[0] * mpi_x[i];
	}
	nelements=1;
	for(i=1; i < numprocs; i=i+1) {
	    buf[0]=0;
	    target=i;
	    MPI_Irecv(buf, nelements, MPI_DOUBLE, target, tag, MPI_COMM_WORLD, request);
	    MPI_Wait(request, status);
	    //buf[0]=20;
	    myres[0] = myres[0] * buf[0];
	}
    }
    else {
	myres[0] = 1;
	for(i = 0; i < mpi_x_size; i=i+1) {
	    myres[0] = myres[0] * mpi_x[i];
	}
	target=0;
	nelements=1;
	MPI_Isend(myres, nelements, MPI_DOUBLE, target, tag, MPI_COMM_WORLD, request);
	MPI_Wait(request, status);
    }
 //   print_num(buf[0]);
    dcc_delete(buf);
}
