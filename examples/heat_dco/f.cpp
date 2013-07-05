// +++++++++++++++++++++++++++++++++++++++++
// data assimilation with the 
// one-dimensional heat equation 
//
// d temp /d t = c * d^2 temp / d x^2
//
// author: Uwe Naumann
//
// heating stick on both ends with scalars 
// a and b
//
// discrete initial condition 
// temp(0,x)=temp in nx+1 vector
//
// nt time steps
//
// discrete observations temp(T,x)=temp_obs 
// in nx+1 vector
//
// objective as least squares fit of 
// prediction to observations
// as function of initial conditions
//
// +++++++++++++++++++++++++++++++++++++++++

// single time step
//#include "dco_tape.hpp"
#include "f.hpp"
#include <iostream>
#include <cassert>
using namespace std;
extern int myid, numprocs;
int counter=0;
void ts(int& nx, 
	active& delta_t, 
	active& c,
	active* temp) 
{
    //printf("counter: %d\n",counter);
    counter=counter+1;
    int mpi_j;
    int mpi_nx;
    active buf[4];
    MPI_Request request[4];
    MPI_Status status[4];
    for(int i = 0 ; i < 4 ; i++) {
	buf[i] = 0;
    }
    mpi_j = (myid * (nx/numprocs))+1;
    mpi_nx = ((myid+1) * (nx/numprocs))+1;
    if(myid == numprocs -1)
	mpi_nx--;
    for(int i = mpi_j ; i < mpi_nx ; i++) {
	temp[i]=temp[i]+c*nx* nx*delta_t*(temp[i-1]-2*temp[i]+temp[i+1]);
	assert(temp[i] == temp[i]);
    }
    // Send to left
    if(myid != 0){
	buf[0]=temp[mpi_j];
	AMPI_Isend(&buf[0], 1, AMPI_DOUBLE, myid-1, 0, AMPI_COMM_WORLD,&request[0]);
    }
    // Recieve from right
    if(myid != numprocs - 1){
	AMPI_Irecv(&buf[1], 1, AMPI_DOUBLE, myid+1, 0, AMPI_COMM_WORLD, &request[1]);
    }
    // Send to right
    if(myid != numprocs-1){
	buf[2]=temp[mpi_nx-1];
	AMPI_Isend(&buf[2], 1, MPI_DOUBLE, myid+1, 0, AMPI_COMM_WORLD,&request[2]);
    }
    // Recieve from left
    if(myid != 0){
	AMPI_Irecv(&buf[3], 1, MPI_DOUBLE, myid-1, 0, AMPI_COMM_WORLD, &request[3]);
    }
    if(myid != numprocs-1) {
	AMPI_Wait(&request[1],&status[1]);
	AMPI_Wait(&request[2],&status[2]);
	temp[mpi_nx] = buf[1];
    }
    if(myid != 0) {
	AMPI_Wait(&request[0],&status[0]);
	AMPI_Wait(&request[3],&status[3]);
	temp[mpi_j-1] = buf[3];
    }
}

void sts(int& nx, 
	active& delta_t, 
	active& c,
	active* temp) 
//$ad indep temp temp temp
//$ad dep temp
{
    //  AMPI_Status status;
    int j=0;
    int jp1=0;
    int jm1=0;
    j=1;
    while (j<nx) {
	jp1=j+1;
	jm1=j-1;
	temp[j]=temp[j]+c*nx*nx*delta_t*(temp[jp1]-2*temp[j]+temp[jm1]);
	j=j+1;
    }
}

// time stepping scheme
void tss(int& nx,
	int& from,
	int& to,
	int& stride,
	active& delta_t,
	active& a,
	active& b,
	active& c,
	active* temp,
	active* temp_obs,
	active& cost)
//$ad indep temp
//$ad dep temp
{
    int i=0;
    int nt=0;
    int lto=0;
    nt=to-from;
    if (nt>stride) {
	lto=from+nt/2;
	tss(nx,from,lto,stride,delta_t,a,b,c,temp,temp_obs,cost);
	lto=from+nt/2+1;
	tss(nx,lto,to,stride,delta_t,a,b,c,temp,temp_obs,cost);
    } else {
	// time integration
	i=from;
	while (i<=to) {
	    ts(nx,delta_t,c,temp);
	    i=i+1;
	}
    }
}

void stss(int& nx,
	int& from,
	int& to,
	int& stride,
	active& delta_t,
	active& a,
	active& b,
	active& c,
	active* temp,
	active* temp_obs,
	active& cost)
//$ad indep temp
//$ad dep temp
{
    int i=0;
    int nt=0;
    int lto=0;
    nt=to-from;
    if (nt>stride) {
	lto=from+nt/2;
	stss(nx,from,lto,stride,delta_t,a,b,c,temp,temp_obs,cost);
	lto=from+nt/2+1;
	stss(nx,lto,to,stride,delta_t,a,b,c,temp,temp_obs,cost);
    } else {
	// time integration
	i=from;
	while (i<=to) {
	    sts(nx,delta_t,c,temp);
	    i=i+1;
	}
    }
}

void af(int& nx, 
	int& nt, 
	int& stride, 
	active& delta_t, 
	active& a,
	active& b,
	active& c,
	active* temp,
	active* temp_obs,
	active& cost) 
//$ad indep temp
//$ad dep cost
{
    int zero=0;
    int j=0;
    tss(nx,zero,nt,stride,delta_t,a,b,c,temp,temp_obs,cost);
    // objective at final time
    j=1;
    cost=0;
    while (j<nx) {
	cost=cost+0.5*(temp[j]-temp_obs[j])*(temp[j]-temp_obs[j])/nx;
	j=j+1;
    }
}
void f(int& nx, 
	int& nt, 
	int& stride, 
	active& delta_t, 
	active& a,
	active& b,
	active& c,
	active* temp,
	active* temp_obs,
	active& cost) 
//$ad indep temp
//$ad dep cost
{
    int zero=0;
    tss(nx,zero,nt,stride,delta_t,a,b,c,temp,temp_obs,cost);
    // objective at final time
    cost=0;
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
    AMPI_Reduce(&mpi_cost, &cost, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}
void sf(int& nx, 
	int& nt, 
	int& stride, 
	active& delta_t, 
	active& a,
	active& b,
	active& c,
	active* temp,
	active* temp_obs,
	active& cost) 
//$ad indep temp
//$ad dep cost
{
    int zero=0;
    int j=0;
    stss(nx,zero,nt,stride,delta_t,a,b,c,temp,temp_obs,cost);
    // objective at final time
    j=1;
    cost=0;
    while (j<nx) {
	cost=cost+0.5*(temp[j]-temp_obs[j])*(temp[j]-temp_obs[j]);
	j=j+1;
    }
}
