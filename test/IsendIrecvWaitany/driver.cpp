#include <iostream>
#include <cmath>
#include <cstdlib>
#include <mpi.h>

using namespace std;

void comp(double *x, double &y, int &n) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  double *buf=new double[n];
  y=0;
  int indx=-1;
  MPI_Request request[2];
  if(rank==0) {
    for(int i=0;i<n;i++) x[i]=x[i]*x[i];
    MPI_Irecv(buf,n,MPI_DOUBLE,1,0,MPI_COMM_WORLD,&request[0]);
    MPI_Isend(x,n,MPI_DOUBLE,1,0,MPI_COMM_WORLD,&request[1]);
    MPI_Wait(&request[0],MPI_STATUS_IGNORE);
    for(int i=0;i<n;i++) {
      y+=buf[i];
    }
  }
  if(rank==1) {
    MPI_Irecv(buf,n,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&request[0]);
    MPI_Wait(&request[0],MPI_STATUS_IGNORE);
    for(int i=0;i<n;i++) { 
      buf[i]=sin(buf[i]);
    }
    MPI_Isend(buf,n,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&request[1]);
  }
  MPI_Waitany(2,request,&indx,MPI_STATUS_IGNORE);
  MPI_Waitany(2,request,&indx,MPI_STATUS_IGNORE);
  delete [] buf;
}


int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  double h=1e-6;
  if(argc<2) {
    cout << "Not enough arguments. Missing problem size." << endl;
    MPI_Finalize();
    return 0;
  }
  int numprocs;
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  if(numprocs!=2) {
    cout << "Must be run with number of processes set to 2." << endl;
    MPI_Finalize();
    return 0;
  }
  int n=atoi(argv[1]);
  cout << "Problem size: " << n << endl;
  double *x=new double[n];
  double y=0;
  for(int i=0;i<n;i++) x[i]=(double) i;
  comp(x,y,n);
  cout << "Result:" << y << endl;
  MPI_Finalize();
  delete [] x;
  return 0;
}
