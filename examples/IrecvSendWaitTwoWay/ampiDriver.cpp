#include <stdio.h>
#include <math.h>
#include "ampi_tape.hpp"
#include "dco.hpp"

using namespace std;

typedef dco::a1s::type active;

int head(active &x, active &y) { 
  MPI_Request r; 
  int world_rank;
  AMPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_rank == 0) {
    x=x*2;
    AMPI_Irecv (&y, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &r);
    AMPI_Send(&x, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
    AMPI_Wait(&r,MPI_STATUS_IGNORE);
    y=y*3;
  } else if (world_rank == 1) {
    active local;
    AMPI_Irecv (&local, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &r);
    AMPI_Wait(&r,MPI_STATUS_IGNORE);
    local=sin(local);
    AMPI_Send (&local, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  } 
}

int main(int argc, char** argv) {
  AMPI_Init(&argc,&argv);
  dco::a1s::global_tape = dco::a1s::tape::create(1e4);
  int world_rank;
  AMPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  active x,y;
  active x_save=0;
  double xp,yp,w, g ;
  if (world_rank == 0) {
    xp=3.5;
    x=xp;
    dco::a1s::global_tape->register_variable(x);
    x_save=x;
    head(x,y);
    cout << __FILE__ << ":process 0 got number " << y << endl;
  } 
  else if (world_rank == 1) {
    head(x,y);
  } 
  if (world_rank == 0) {
    dco::a1s::set(y, 1., -1);
    dco::a1s::global_tape->interpret_adjoint();
    dco::a1s::get(x_save, g, -1);
    cout << __FILE__ << ":process 0 got gradient " << g << endl;
    
  } 
  else if (world_rank == 1) {
    dco::a1s::set(y, 0., -1);
    dco::a1s::global_tape->interpret_adjoint();
    dco::a1s::get(x, g, -1);
    cout << __FILE__ << ":process 1 got gradient " << g << endl;
  }   
  AMPI_Finalize();
  return 0;
}
