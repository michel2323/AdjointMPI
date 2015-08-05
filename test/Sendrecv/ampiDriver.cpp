#include <iostream>
#include <cmath>
#include <cstdlib>
#include "dco.hpp"
#include "ampi_tape.hpp"

using namespace std;

typedef dco::ga1s<double>::type active;



void comp(active *x, active &y, int &n) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  active *buf=new active[n];
  y=0;
  if(rank==0) {
    for(int i=0;i<n;i++) x[i]=x[i]*x[i];
    AMPI_Sendrecv(x,n,MPI_DOUBLE,1,0,buf,n,MPI_DOUBLE,1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    for(int i=0;i<n;i++) {
      y+=buf[i];
    }
  }
  if(rank==1) {
    AMPI_Recv(buf,n,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    for(int i=0;i<n;i++) { 
      buf[i]=sin(buf[i]);
    }
    AMPI_Send(buf,n,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
  }
  delete [] buf;
}


int main(int argc, char *argv[]) {
  AMPI_Init(&argc, &argv);
  dco::ga1s<double>::global_tape = dco::ga1s<double>::tape_t::create(1e4);
  active h=1e-6;
  if(argc<2) {
    cout << "Not enough arguments. Missing problem size." << endl;
    AMPI_Finalize();
    return 0;
  }
  int numprocs,rank;
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if(numprocs!=2) {
    cout << "Must be run with number of processes set to 2." << endl;
    AMPI_Finalize();
    return 0;
  }
  int n=atoi(argv[1]);
  cout << "Problem size: " << n << endl;
  active *x=new active[n];
  active y=0;
  for(int i=0;i<n;i++) x[i]=(active) i;
  if(rank==0) {
    for(int i=0;i<n;i++)
      dco::ga1s<double>::global_tape->register_variable(x[i]);
  }
  active *x_saved=new active[n];
  for(int i=0;i<n;i++) x_saved[i]=x[i];
  comp(x,y,n);
  cout << "Result:" << y << endl;
  cout << "Derivatives:" << endl;
  if(rank == 0) dco::ga1s<double>::set(y, 1., -1);
  if(rank == 1) dco::ga1s<double>::set(y, 0., -1);
  dco::ga1s<double>::global_tape->interpret_adjoint();
  //if(rank == 0) {
    double g=0;
    for(int i=0;i<n;i++) {
      dco::ga1s<double>::get(x_saved[i], g, -1);
      cout << g << endl;
    }
    //}
    
  AMPI_Finalize();
  delete [] x; delete [] x_saved;
  return 0;
}
