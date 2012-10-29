#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <tape.hpp>
#include <iostream>

#define TAPE_MODE 1
#define DERIV_MODE 1

using namespace std;

int main(int argc, char** argv){
    // INIT AND STUFF
    int  myid, numprocs;
    int namelen;
    char processor_name[AMPI_MAX_PROCESSOR_NAME];
    AMPI_Init(&argc,&argv);
    AMPI_Comm_size(AMPI_COMM_WORLD,&numprocs);
    AMPI_Comm_rank(AMPI_COMM_WORLD,&myid);
    AMPI_Get_processor_name(processor_name,&namelen);
    AMPI_Status status;
    int n = 2;
    AMPI_Status * status_array = new AMPI_Status[n];
    AMPI_dco_Request * request = new AMPI_dco_Request[n];
    double time0 = MPI_Wtime();
    request[0].r.aw = 0;
    request[1].r.aw = 0;
    init_mode(TAPE_MODE, DERIV_MODE);
    active ** problem = new active*[n];
    active ** result = new active*[n];
    active ** buf;
    active ** res = new active*[n];
    for(int i = 0 ; i < n ; i++) {
	res[i] = new active[n];
	problem[i] = new active[n];
	result[i] = new active[n];
    }
    for(int i = 0 ; i < n ; i++) {
	for(int j = 0 ; j < n ; j++) {
	    //problem[i][j] = 0.1;
	    problem[i][j] = 2;
	    independent(problem[i][j]);
	    //	printf("problem[i]: %f\n", problem[i].v);
	    //		result[i] = 1;
	}
    }
    if(myid == 0) {
	request[0].r.aw = 0;
	request[1].r.aw = 0;
//	AMPI_Awaitall(n, request, status_array);
	//	printf("AWAITALL: %d\n", request[0].aw);
	buf = new active*[n];
	for(int i = 0 ; i < n ; i++) {
	    buf[i] = new active[n];
	    for(int j = 0 ; j < n ; j++)
		buf[i][j] = 0;
	    AMPI_Irecv(buf[i], n, AMPI_DOUBLE, 1, 0, AMPI_COMM_WORLD, &request[i]);
	}
	//	AMPI_Wait(&request[0], &status);
	for(int i = 0 ; i < n ; i++) {
	    for(int j = 0 ; j < n ; j++) {
		//result[i][j] = sin(problem[i][j])* problem[i][j];
		result[i][j] = problem[i][j]* problem[i][j];
	    }
	}
	//		AMPI_Recv(buf, n, AMPI_DOUBLE, 1, 0, AMPI_COMM_WORLD, &status);
	AMPI_Wait(&request[0], &status);
	AMPI_Wait(&request[1], &status);
//	AMPI_Waitall(n, request, status_array); 
    }
    else {
	buf = new active*[n];
	for(int i = 0 ; i < n ; i++) {
	    for(int j = 0 ; j < n ; j++) {
		//result[i][j] = cos(problem[i][j])* problem[i][j];
		result[i][j] = problem[i][j]* problem[i][j];
	    }
	    buf[i] = result[i];
	    AMPI_Isend(buf[i], n, AMPI_DOUBLE, 0, 0, AMPI_COMM_WORLD, &request[i]);
	}
	//AMPI_Send(buf, n, AMPI_DOUBLE, 0, 0, AMPI_COMM_WORLD);
	AMPI_Wait(&request[0], &status);
	AMPI_Wait(&request[1], &status);
//	AMPI_Waitall(n, request, status_array); 
    }
    if(myid == 0) {
	for(int i = 0 ; i < n ; i++){
	    for(int j = 0 ; j < n ; j++) {
		res[i][j] = buf[i][j]*result[i][j];
	    }
	}
	for(int i = 0 ; i < n ; i++){
	    for(int j = 0 ; j < n ; j++) {
		dependent(res[i][j]);
		seed_dependent(res[i][j]);
	    }
	}
    }
    double time1 = MPI_Wtime();
    AMPI_Barrier(AMPI_COMM_WORLD);
    reverse_tape_interpreter();
    for(int i = 0 ; i < numprocs ; i++) {
	if(myid == i)
	    //   print_indeps();
	    print_tape();
	AMPI_Barrier(AMPI_COMM_WORLD);
    }
    AMPI_Barrier(AMPI_COMM_WORLD);
    double time2 = MPI_Wtime();
    if(myid == 0){
	printf("Forward run: %f s\n", time1-time0);
	printf("Backward run: %f s\n", time2-time1);
    }
    AMPI_Finalize();
}
