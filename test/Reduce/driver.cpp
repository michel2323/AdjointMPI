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
  if(rank==0) {
    for(int i=0;i<n;i++) x[i]=x[i]*x[i];
    for(int i=0;i<n;i++) buf[i]=x[i];
    MPI_Bcast(x,n,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Reduce(buf,x,n,MPI_DOUBLE,MPI_PROD,0,MPI_COMM_WORLD);
    MPI_Bcast(x,n,MPI_DOUBLE,0,MPI_COMM_WORLD);
    for(int i=0;i<n;i++) { 
      x[i]=sin(x[i]);
    }
    MPI_Reduce(x,buf,n,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    for(int i=0;i<n;i++) {
      y+=buf[i];
    }
  }
  else {
    MPI_Bcast(buf,n,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Reduce(buf,x,n,MPI_DOUBLE,MPI_PROD,0,MPI_COMM_WORLD);
    MPI_Bcast(buf,n,MPI_DOUBLE,0,MPI_COMM_WORLD);
    for(int i=0;i<n;i++) { 
      buf[i]=sin(buf[i]);
    }
    MPI_Reduce(buf,x,n,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
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
  int numprocs;
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  int n=atoi(argv[1]);
  cout << "Problem size: " << n << endl;
  double *x=new double[n];
  double y=0;
  for(int i=0;i<n;i++) x[i]=(double) i+0.1;
  comp(x,y,n);
  cout << "Result:" << y << endl;
  MPI_Finalize();
  delete [] x;
  return 0;
}
