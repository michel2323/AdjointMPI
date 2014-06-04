#include <iostream>
#include <cmath>
#include <cstdlib>
#include <mpi.h>

int rank;

using namespace std;

void passive_pattern(double *x, int &n) {
  if(rank==0)
    MPI_Send(x,n,MPI_DOUBLE,1,0,MPI_COMM_WORLD);
  if(rank==1)
    MPI_Recv(x,n,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}

void adjoint_forward_pattern(double *x, int &n) {
  if(rank==0)
    MPI_Send(x,n,MPI_DOUBLE,1,0,MPI_COMM_WORLD);
  if(rank==1)
    MPI_Recv(x,n,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}

void adjoint_reverse_pattern(double *x, double *z, int &n) {
  if(rank==1) {
    MPI_Send(z,n,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
    for(int i=0;i<n;i++) x[i]+=z[i];
  }
  if(rank==0)
    MPI_Recv(x,n,MPI_DOUBLE,1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  int n=atoi(argv[1]);
  double *x=new double[n]; 
  double t0=MPI_Wtime();
  passive_pattern(x,n);
  double t1=MPI_Wtime();
  cout << "Passive pattern: " << t1-t0 << " s." << endl;
  delete [] x; 

  x=new double[n]; 
  double t2=MPI_Wtime();
  adjoint_forward_pattern(x,n);
  double t3=MPI_Wtime();
  cout << "Active forward pattern: " << t3-t2 << " s." << endl;
  delete [] x; 

  x=new double[n]; 
  double *z=new double[n]; 
  double t4=MPI_Wtime();
  adjoint_reverse_pattern(x,z,n);
  double t5=MPI_Wtime();
  cout << "Active reverse pattern: " << t5-t4 << " s." << endl;
  delete [] x; delete [] z; 

  cout << "Ratio:" << (t5-t4+t3-t2)/(t1-t0) << endl; 
  MPI_Finalize();
  return 0;
}
