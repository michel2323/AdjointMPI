#include <iostream>
#include <cmath>
#include <cstdlib>
#include <mpi.h>

int rank;
MPI_Win win;

using namespace std;

void passive_pattern(double *x, int &n) {
  if(rank==0) {
    MPI_Win_fence(0,win);
    MPI_Accumulate(x,n,MPI_DOUBLE,1,0,n,MPI_DOUBLE,MPI_SUM,win);
    MPI_Win_fence(0,win);
  }
  if(rank==1) {
    MPI_Win_fence(0,win);
    MPI_Win_fence(0,win);
  }
}

void adjoint_forward_pattern(double *x, int &n) {
  if(rank==0) {
    MPI_Win_fence(0,win);
    MPI_Accumulate(x,n,MPI_DOUBLE,1,0,n,MPI_DOUBLE,MPI_SUM,win);
    MPI_Win_fence(0,win);
  }
  if(rank==1) {
    MPI_Win_fence(0,win);
    MPI_Win_fence(0,win);
  }
}

void adjoint_reverse_pattern(double *x, double *z, int &n) {
  if(rank==0) {
    MPI_Win_fence(0,win);
    MPI_Get(z,n,MPI_DOUBLE,1,0,n,MPI_DOUBLE,win);
    MPI_Win_fence(0,win);
    for(int i=0;i<n;i++) x[i]+=z[i];
    MPI_Win_fence(0,win);
  }
  if(rank==1) {
    MPI_Win_fence(0,win);
    MPI_Win_fence(0,win);
    MPI_Win_fence(0,win);
  }
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  int n=atoi(argv[1]);

  double *x=new double[n]; 
  if(rank==0) {
    MPI_Win_create(NULL,0,1,MPI_INFO_NULL,MPI_COMM_WORLD,&win);
  }
  if(rank==1) {
    MPI_Win_create(x,n*sizeof(double), sizeof(double),MPI_INFO_NULL,MPI_COMM_WORLD,&win);
  }
  double t0=MPI_Wtime();
  passive_pattern(x,n);
  double t1=MPI_Wtime();
  cout << "Passive pattern: " << t1-t0 << " s." << endl;
  MPI_Win_free(&win);
  delete [] x; 

  x=new double[n]; 
  if(rank==0) {
    MPI_Win_create(NULL,0,1,MPI_INFO_NULL,MPI_COMM_WORLD,&win);
  }
  else {
    MPI_Win_create(x,n*sizeof(double), sizeof(double),MPI_INFO_NULL,MPI_COMM_WORLD,&win);
  }
  double t2=MPI_Wtime();
  adjoint_forward_pattern(x,n);
  double t3=MPI_Wtime();
  cout << "Active forward pattern: " << t3-t2 << " s." << endl;
  MPI_Win_free(&win);
  delete [] x; 

  x=new double[n]; 
  double *z=new double[n]; 
  if(rank==0) {
    MPI_Win_create(NULL,0,1,MPI_INFO_NULL,MPI_COMM_WORLD,&win);
  }
  else {
    MPI_Win_create(z,n*sizeof(double), sizeof(double),MPI_INFO_NULL,MPI_COMM_WORLD,&win);
  }
  x=new double[n]; 
  double t4=MPI_Wtime();
  adjoint_reverse_pattern(x,z,n);
  double t5=MPI_Wtime();
  cout << "Active reverse pattern: " << t5-t4 << " s." << endl;
  MPI_Win_free(&win);
  delete [] x; delete [] z; 

  cout << "Ratio:" << (t5-t4+t3-t2)/(t1-t0) << endl; 
  MPI_Finalize();
  return 0;
}
