#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "dco_tape.hpp"

#define TAPE_MODE 1
#define DERIV_MODE 1
int main(int argc, char** argv){

    // INIT AND STUFF

    int  myid, numprocs;
    int n=4;
    int namelen;
    char processor_name[AMPI_MAX_PROCESSOR_NAME];
    AMPI_Init(&argc,&argv);
    AMPI_Comm_size(AMPI_COMM_WORLD,&numprocs);
    AMPI_Comm_rank(AMPI_COMM_WORLD,&myid);
    AMPI_Get_processor_name(processor_name,&namelen);
    AMPI_Status status, status2;

    int i = 0;
    init_mode(TAPE_MODE, DERIV_MODE);
    active myres;
    active res;
    myres = 1;
    res = 1;
    active x[n/numprocs];

    // FORWARD SPEELPENNING
    AMPI_Request request;
    AMPI_Request request2;
    for(int i=0 ; i < n/numprocs ; i++) {
	x[i] = myid * n/numprocs + i + 2;
        independent(x[i]);
	myres=myres*x[i];
    }
    
    //x = (double) (myid+2);

    request.aw = 0;
    request2.aw = 0;
    //if(myid == 0) {
	////myres = myres * x;
	////printf("myres.v: %f, myres.va: %d\n", myres.v, myres.va);
	////AMPI_Send(&myres, 1, MPI_DOUBLE, myid+1, 0, MPI_COMM_WORLD);
	//AMPI_Isend(&myres, 1, MPI_DOUBLE, myid+1, 0, MPI_COMM_WORLD, &request);
	//AMPI_Wait(&request, &status);
    //}	
    //if(myid > 0 && myid < numprocs-1) {
	////AMPI_Recv(&res, 1, MPI_DOUBLE, myid-1, 0, MPI_COMM_WORLD, &status);
	//AMPI_Irecv(&myres, 1, MPI_DOUBLE, myid-1, 0, MPI_COMM_WORLD, &request);
	//AMPI_Wait(&request, &status);
	//myres = myres * res;
	////AMPI_Send(&myres, 1, MPI_DOUBLE, myid+1, 0, MPI_COMM_WORLD);
	//AMPI_Isend(&myres, 1, MPI_DOUBLE, myid+1, 0, MPI_COMM_WORLD, &request);
	//AMPI_Wait(&request, &status2);
    //}
    //if(myid == numprocs - 1) {	
	////AMPI_Recv(&res, 1, MPI_DOUBLE, myid-1, 0, MPI_COMM_WORLD, &status);
	//AMPI_Irecv(&res, 1, MPI_DOUBLE, myid-1, 0, MPI_COMM_WORLD, &request);
	//AMPI_Wait(&request, &status);
	////printf("myres: (%f,%d)\n", myres.v, myres.va);
	////printf("res: (%f,%d)\n", res.v, res.va);
	//res = myres * res;
	////printf("myres: (%f,%d) x: (%f,%d)\n", myres.v, myres.va);
	////printf("res: (%f,%d)\n", res.v, res.va);
    //}
    //AMPI_Reduce(&myres, &res, 1, MPI_DOUBLE, MPI_PROD, numprocs-1, MPI_COMM_WORLD);
    AMPI_Allreduce(&myres, &res, 1, MPI_DOUBLE, MPI_PROD, MPI_COMM_WORLD);
    if(myid == numprocs-1){
	dependent(res);
	seed_dependent(res);
	printf("Res[%d]: %f\n", myid, res.v);
    }

    // PRINT TAPE AFTER FORWARD MODE


    	for(i = 0 ; i < numprocs ; i++) {
	if(myid == i)
	    //print_tape();
	AMPI_Barrier(AMPI_COMM_WORLD);
	}
		
    // REVERSE
 //   printf("finished forward\n");

//    interpret_tape(AMPI_COMM_WORLD, &status);
	reverse_tape_interpreter();

   // printf("finished backward\n");
    // PRINT TAPE AFTER REVERSE MODE

    for(i = 0 ; i < numprocs ; i++) {
	if(myid == i){
	    //print_tape();
	    print_tape_indeps();
	}
	AMPI_Barrier(AMPI_COMM_WORLD);
    }
    AMPI_Barrier(AMPI_COMM_WORLD);
    AMPI_Finalize();
} 

