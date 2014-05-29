#include <iostream>
#include <cmath>
#include <cstdlib>
#include <mpi.h>

int rank;
int numprocs;

using namespace std;

void passive_pattern(double *sendbuf, double *recvbuf, int &n) {
  MPI_Reduce(sendbuf,recvbuf,n,MPI_DOUBLE,MPI_PROD,0,MPI_COMM_WORLD);
}

void adjoint_forward_pattern(double *sendbuf, double *recvbuf, int &n) {
  MPI_Allreduce(sendbuf,recvbuf,n,MPI_DOUBLE,MPI_PROD,MPI_COMM_WORLD);
  MPI_Alltoall(sendbuf,n,MPI_DOUBLE,recvbuf,n,MPI_DOUBLE,MPI_COMM_WORLD);
}

void adjoint_reverse_pattern(double *sendbuf, double *recvbuf, int &n) {
  MPI_Allreduce(sendbuf,recvbuf,n,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  int n=atoi(argv[1]);
  double *sendbuf=new double[n]; 
  double *recvbuf=new double[n*numprocs]; 
  double t0=MPI_Wtime();
  passive_pattern(sendbuf,recvbuf,n);
  double t1=MPI_Wtime();
  cout << "Passive pattern: " << t1-t0 << " s." << endl;
  delete [] sendbuf; delete [] recvbuf; 

  sendbuf=new double[n]; 
  recvbuf=new double[n*numprocs]; 
  double t2=MPI_Wtime();
  adjoint_forward_pattern(sendbuf,recvbuf,n);
  double t3=MPI_Wtime();
  cout << "Active forward pattern: " << t3-t2 << " s." << endl;
  delete [] sendbuf; delete [] recvbuf; 

  sendbuf=new double[n]; 
  recvbuf=new double[n]; 
  double t4=MPI_Wtime();
  adjoint_reverse_pattern(sendbuf,recvbuf,n);
  double t5=MPI_Wtime();
  cout << "Active reverse pattern: " << t5-t4 << " s." << endl;
  delete [] sendbuf; delete [] recvbuf; 

  cout << "Ratio:" << (t5-t4+t3-t2)/(t1-t0) << endl; 
  MPI_Finalize();
  return 0;
}
