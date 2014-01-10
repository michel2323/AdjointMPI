#include <iostream>
#include <cmath>
#include <cstdlib>
#include <mpi.h>

using namespace std;

void comp(double *x, double &y, int n) {
  int numprocs, rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  int n_local=n/numprocs;
  double *buf=new double[n_local];
  y=0;
  for(int i=0;i<n;i++) x[i]=x[i]*x[i];
  MPI_Scatter(x,n_local,MPI_DOUBLE,buf,n_local,MPI_DOUBLE,0,MPI_COMM_WORLD);
  for(int i=0;i<n_local;i++) { 
    buf[i]=sin(buf[i]);
  }
  MPI_Gather(buf,n_local,MPI_DOUBLE,x,n_local,MPI_DOUBLE,0,MPI_COMM_WORLD);
  if(rank==0) {
    for(int i=0;i<n;i++) {
      y+=x[i];
    }
  }
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
  int numprocs,rank;
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  int n=atoi(argv[1]);
  cout << "Problem size: " << n << endl;
  double *x=new double[n*numprocs];
  double y=0;
  if(rank==0) {
    for(int i=0;i<numprocs;i++) {
      for(int j=0;j<n;j++) x[i*n+j]=(double) j;
    }
  }
  else {
    for(int i=0;i<n*numprocs;i++) x[i]=0;
  }
  comp(x,y,n*numprocs);
  cout << "Result:" << y << endl;
  MPI_Finalize();
  delete [] x;
  return 0;
}
