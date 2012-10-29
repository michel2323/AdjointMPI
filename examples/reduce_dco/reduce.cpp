#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include "dco_tape.hpp"

#define TAPE_MODE 1
#define DERIV_MODE 1

using namespace std;


int main(int argc, char **argv)
{
    int nx = 2;
    int myid, numprocs;
    init_mode(TAPE_MODE, DERIV_MODE);
    active a = 2;
    active b = 0;
    /* MPI Initialization */
    AMPI_Init(&argc, &argv);
    AMPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    AMPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Status status;
    active * temp_obs = new active[nx+1];
    active * temp = new active[nx+1];
    temp_obs[0]=a;
    for (int i=1;i<nx;i++) temp_obs[i]=2.0/(double) i;
    temp_obs[nx]=b;
    temp[0]=a;
    for (int j=1;j<nx;j++) temp[j]=0;
    temp[nx]=b;

    for (int j=0;j<=nx;j++) independent(temp[j]);
    active cost;
    cost = 0;
    active mpi_cost;
    mpi_cost = 0;
    int mpi_j = 0;
    int mpi_nx = 0;
    mpi_j = (myid * (nx/numprocs))+1;
    mpi_nx = ((myid+1) * (nx/numprocs))+1;
    if(myid == numprocs -1)
	mpi_nx--;
    while (mpi_j<mpi_nx) {
	mpi_cost=mpi_cost+0.5*(temp[mpi_j]-temp_obs[mpi_j])*(temp[mpi_j]-temp_obs[mpi_j]);
	mpi_j=mpi_j+1;
    }
    cout << "mpi_cost.v: " << mpi_cost.v << " mpi_cost.d: " << mpi_cost.d << endl;
    AMPI_Reduce(&mpi_cost, &cost, 1, AMPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    cout << "cost: " << cost.v << " cost.d: " << cost.d << endl;
    dependent(cost);
    seed_dependent(cost);
    reverse_tape_interpreter();
    print_tape();
    for(int i = 0 ; i < numprocs ; i++) {
	if(myid == i)
	    print_tape_indeps();
	MPI_Barrier(MPI_COMM_WORLD);
    }
    AMPI_Finalize();
    return 0;
}

