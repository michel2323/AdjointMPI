#include <stdio.h>
#include <mpi.h>
#include <math.h>

using namespace std;

int head(double &x, double &y) { 
  MPI_Request r; 
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_rank == 0) {
    x=x*2;
    MPI_Irecv(&y, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &r);
    MPI_Send(&x, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
    MPI_Wait(&r,MPI_STATUS_IGNORE);
    y=y*3;
  } else if (world_rank == 1) {
    double local;
    MPI_Irecv(&local, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &r);
    MPI_Wait(&r,MPI_STATUS_IGNORE);
    local=sin(local);
    MPI_Send(&local, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  } 
}

int main(int argc, char** argv) {
  MPI_Init(&argc,&argv);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  double x,y;
  double xv,h,yp;
  xv=3.5; 
  h=1.0e-8;
  if (world_rank == 0) {
    x=xv;
    head(x,y);
    cout << __FILE__  << ": process 0 got number" << y << endl;
  } else if (world_rank == 1) {
    head(x,y);
  } 
  if (world_rank == 0) {
    x=xv+h;
    head(x,yp);
    cout << __FILE__ << ": process 0 got fd value" <<  (yp-y)/h << endl;
  } else if (world_rank == 1) {
    head(x,y);
  } 
  MPI_Finalize();
  return 0;
}
