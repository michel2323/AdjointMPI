#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cassert>
using namespace std;

//#include "dco_tape.hpp"
#include "f.hpp"
#define TAPE_MODE 1
#define DERIV_MODE 1

int myid, numprocs;
/*
 * Euclidean norm of an n-vector
 */
double norm(int n, double *x) {
    double res=0;
    for (int i=0;i<n;i++) res+=x[i]*x[i];
    return sqrt(res);
}

int main(int argc, char* argv[]) {
    /* AMPI Initialization */
    double tmp=0;
    AMPI_Init(&argc, &argv);
    AMPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    AMPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Status status;
    MPI_Status status2;
    //dco::a1s::global_tape->reset();
    dco::a1s::global_tape=dco::a1s::tape::create(1e8);
    double time1,time2;
    time1 = MPI_Wtime();
    int nit=0;
    // for stability of Euler's methode
    // we require delta_t <= (1/nx)^2/(2*c) 
    // =1/(nx^2*2*c)
    // int nx = 50; // default
    int nx=50;
    // int nt = 100; // default
    int nt=100;
    active cost, cost_d, delta_t=1./nt, c=1e-3;
    // max time steps on tape
    int stride=10;
    // boundary condition
    active a=2, b=0;
    // initial condition
    // observations
    active* temp_obs=new active[nx+1];
    temp_obs[0]=a;
   //  for (int i=1;i<nx;i++) temp_obs[i]=((1.*i)/nx-.5)*((1.*i)/nx-.5);
 //   for (int i=1;i<nx;i++) temp_obs[i]=2./i;
    for(int i = 1; i < nx; i++) temp_obs[i] = a - (double) i*a/nx ;
    // for (int i=1;i<nx;i++) temp_obs[i]=nx-i;
    // for (int i=1;i<nx;i++) temp_obs[i]=1;
    temp_obs[nx]=b;

    active* temp=new active[nx+1];
    temp[0]=a;
    for (int i=1;i<nx;i++) temp[i]=0;
    temp[nx]=b;

    active b1_cost, b1_delta_t=0, b1_c=0;
    active b1_a=0, b1_b=0;
    active* b1_temp=new active[nx+1];
    active* b1_temp_obs=new active[nx+1];
    //double* b1_temp_aux= new double[nx-1];
    double* b1_temp_aux= new double[nx+1];

    for (int i=0;i<=nx;i++) b1_temp_obs[i]=0;

    double* temp_aux=new double[nx+1];
    //for (int j=0;j<=nx;j++) temp_aux[j]=temp[j];
    for (int j=0;j<=nx;j++) dco::a1s::get(temp[j],temp_aux[j]);
    double buf[2];
    double alpha=1;
    double ng;
    const double eps=1e-3;
    f(nx,nt,stride,delta_t,a,b,c,temp,temp_obs,cost);
    cout << "start cost=" << cost << endl;

    int mpi_j = 0;
    int mpi_nx = 0;

    for (int i=1;i<nx;i++) temp[i]=0;
    do {
	//reset_tape();
	dco::a1s::global_tape->reset();
	AMPI_Reset_Tape();
	a = 2.0;
	b = 0.0;
	temp[0]=a;
	for (int j=1;j<nx;j++) temp[j]=temp_aux[j];
	temp[nx]=b;
	delta_t=1./nt;
	alpha=1;
	c=1e-3;
	temp_obs[0]=a;
	for (int i=1;i<nx;i++){ 
	  //temp_obs[i]=2.0 - (double) i / 100.0;
            temp_obs[i] = a - (double) i*a/nx ;
	    assert(temp_obs[i] == temp_obs[i]);
	}
	temp_obs[nx]=b;
	cost = 0;
	for (int j=0;j<=nx;j++) dco::a1s::global_tape->register_variable(temp[j]);
	f(nx,nt,stride,delta_t,a,b,c,temp,temp_obs,cost);
	//if(myid == 0)
	cost_d = cost;
	dco::a1s::set(cost_d,1.0,-1);
//	interpret_tape(MPI_COMM_WORLD, &status);
	dco::a1s::global_tape->interpret_adjoint();
	//cout << "cost: " << cost.v << "temp: " << temp[0].d << " " << tape_entry::indeps[0] << endl;
	//cout << "PRINTING TAPE\n";
	//MPI_Barrier(MPI_COMM_WORLD);
	//for(int i = 0 ; i < numprocs ; i++){
	    //if(myid == i)
		//print_tape_indeps();
	    //MPI_Barrier(MPI_COMM_WORLD);
	//}

	mpi_j = (myid * (nx/numprocs))+1;
	mpi_nx = ((myid+1) * (nx/numprocs))+1;
	if(myid == numprocs -1)
	    mpi_nx--;

	for (int i=mpi_j;i<mpi_nx;i++) {
	    dco::a1s::get(temp[i],tmp,-1);
	    temp_aux[i] = temp_aux[i] - (alpha*tmp);
	    //temp_aux[i] = temp_aux[i] - (alpha*tape[tape_entry::indeps[i]].d);
	}

	if(myid != 0){
	    buf[0]=temp_aux[mpi_j];
	    MPI_Send(&buf[0], 1, MPI_DOUBLE, myid-1, 0, AMPI_COMM_WORLD);
	}
	// Recieve from right
	if(myid != numprocs - 1){
	    MPI_Recv(&buf[1], 1, MPI_DOUBLE, myid+1, 0, AMPI_COMM_WORLD, &status);
	    temp_aux[mpi_nx] = buf[1];
	}
	if(myid != 0) {
	}
	// Send to right
	if(myid != numprocs-1){
	    buf[0]=temp_aux[mpi_nx-1];
	    MPI_Send(&buf[0], 1, MPI_DOUBLE, myid+1, 0, AMPI_COMM_WORLD);
	}
	// Recieve from left
	if(myid != 0){
	    MPI_Recv(&buf[1], 1, MPI_DOUBLE, myid-1, 0, AMPI_COMM_WORLD, &status2);
	    temp_aux[mpi_j-1] = buf[1];
	}
	if(myid != numprocs-1) {
	}

	if(myid == 0)
	    mpi_j--;	    
	nit++;

	for (int i=0;i<nx-1;i++){
	    dco::a1s::get(temp[i+1],tmp,-1);
	    b1_temp_aux[i]=tmp;
	    //b1_temp_aux[i]=tape[tape_entry::indeps[i+1]].d;
	}
	if(numprocs > 1){ 
	    if(myid != 0) {
		MPI_Send(&b1_temp_aux[mpi_j], mpi_nx - mpi_j, MPI_DOUBLE, 0, 0, AMPI_COMM_WORLD);
	    }
	    else {
		for(int i = 1 ; i < numprocs ; i++) {
		    mpi_j = (i * (nx/numprocs))+1;
		    mpi_nx = ((i+1) * (nx/numprocs))+1;
		    if(i == numprocs - 1)
			mpi_nx--;
		    MPI_Recv(&b1_temp_aux[mpi_j], mpi_nx - mpi_j, MPI_DOUBLE, i, 0, AMPI_COMM_WORLD, &status);
		}
	    }
	}
	ng=norm(nx-1,b1_temp_aux);
	MPI_Bcast(&ng, 1, MPI_DOUBLE, 0, AMPI_COMM_WORLD );
	//cout << "id: "<< myid <<" ng: " << ng << endl;
	//cout << "id: "<< myid <<" eps: " << eps << endl;
		if ((!(nit%10)) && (myid == 0)) cout << "ng: " << ng << endl;
    } while (ng>eps);
    if(numprocs > 1){ 
	if(myid != 0) {
	    MPI_Send(&temp_aux[mpi_j], mpi_nx - mpi_j, MPI_DOUBLE, 0, 0, AMPI_COMM_WORLD);
	}
	else {
	    for(int i = 1 ; i < numprocs ; i++) {
		mpi_j = (i * (nx/numprocs))+1;
		mpi_nx = ((i+1) * (nx/numprocs))+1;
		if(i == numprocs - 1)
		    mpi_nx--;
		MPI_Recv(&temp_aux[mpi_j], mpi_nx - mpi_j, MPI_DOUBLE, i, 0, AMPI_COMM_WORLD, &status);
	    }
	}
    }
    if(myid == 0) {
	temp[0] = a;
	for (int j=1;j<nx;j++) temp[j]=temp_aux[j];
	temp[nx] = b;
	cout << "solution ..." << endl;
	sf(nx,nt,stride,delta_t,a,b,c,temp,temp_obs,cost);
	//cout << "optimized cost=" << cost.v << endl;
	cout << "optimized cost=" << cost << endl;
	for (int i=0;i<=nx;i++) cout << "temp[" << i << "]=" << temp_aux[i] << endl;
	//for (int i=0;i<=nx;i++) cout << "temp_obs[" << i << "]=" << temp_obs[i].v << endl;
	for (int i=0;i<=nx;i++) cout << "temp_obs[" << i << "]=" << temp_obs[i] << endl;
	//for (int i=0;i<=nx;i++) cout << "temp[" << i << "]=" << temp[i].v << endl;
	for (int i=0;i<=nx;i++) cout << "temp[" << i << "]=" << temp[i] << endl;
	cout << "#iterations=" << nit << endl;
    }

       delete [] b1_temp;
       delete [] b1_temp_aux;
       delete [] b1_temp_obs;
       delete [] temp_obs;
       delete [] temp;
    time2 = MPI_Wtime();
    if(myid == 0)
	cout << "Time: " << time2-time1 << endl;
    AMPI_Finalize();

    return 0;
}

