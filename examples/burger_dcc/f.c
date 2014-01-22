// Viscous Burger's equation
// E. Kalnay: Atmospheric modeling, data assimilation and predictability
// Cambridge University Press

int myid=0;
int numprocs=0;

void f(int& nx, // number of grid points
       int& n, // number of time steps
       double& cost, // cost function
       double** uob, // observations
       double** ub, // basic states
       double** us, // step i-1 saved states
       double** u, // model solutions
       double* ui, // initial conditions
       double& dx, // space increment
       double& dt, // time increment
       double& r, // Reynolds number
       double& dtdx,
       double& c0,
       double& c1,
       double* buf // MPI buffer
)
#pragma ad indep ui
#pragma ad dep cost

{
    // MPI
  int mpi_i=0;
  int mpi_j=0;
  int mpi_k=0;
  int mpi_n=0;
  int mpi_x=0;
  int mpi_nx=0;
  int nelements = 0;
  int target = 0;
  int tag = 0;
  int ampi_double = 0;
  int ampi_comm = 0;
  int request[1]; // MPI Request
  int request2[1]; // MPI Request
  int status[1]; // MPI Request

  double dtdxsq=0;
  double cost_mpi=0;
  int i=0;
  int j=0;
  int k=0;
  int aux1=0;
  int aux2=0;
  int aux3=0;
  int aux4=0;
  int aux5=0;
  r=1000.0;
  dx=1.0;
  dt=0.1;
  dtdx=dt/dx;
  dtdxsq=dt/(dx*dx);
  c1=(2.0/r)*dtdxsq;
  c0=1.0/(1.0+c1);
  cost=0.0;
  cost_mpi=0.0;
  request[0]=0;
  status[0]=0;

  // initial conditions
  i=0;
  while (i<nx) {
    u[i][0]=ui[i];
    us[i][0]=u[i][0];
    i=i+1;
  }

  // boundary conditions
  j=0;
  while (j<n) {
    u[0][j]=0.0;
    us[0][j]=0.0;
    aux1=nx-1;
    u[aux1][j]=0.0;
    us[aux1][j]=0.0;
    j=j+1;
  }

  // numerical integration
  
  // FTCS for first time step
  i=1;
  k=nx-1;
  while (i<k) {
    aux1=i+1;
    aux2=i-1;
    u[i][1]=u[i][0]-0.5*dtdx*u[i][0]*(u[aux1][0]-u[aux2][0])+0.5*c1*(u[aux1][0]-2.0*u[i][0]+u[aux2][0]);
    us[i][1]=u[i][1];
    i=i+1;
  }

  // Leap Frog / DuPont-Frankel afterwards 

//  mpi_i = (myid) * (nx/numprocs);
//  mpi_nx = (myid+1) * (nx/numprocs);

j=2;
while (j<n) {
    //  while (i<k) {
    i=1;
    mpi_i = i;
    mpi_nx=nx;
    if(mpi_i == 0) {
	mpi_i = 1;
    }
    if(mpi_nx == nx) {
	mpi_nx = nx-1;
    }
    i=mpi_i;
    k=mpi_nx;
    while (i<k) {
	aux1=i+1;
	aux2=i-1;
	aux3=j-1;
	aux4=j-2;
	us[i][j]=u[i][j];
	u[i][j]=c0*(u[i][aux4]+c1*(us[aux1][aux3]-u[i][aux4]+u[aux2][aux3])-dtdx*u[i][aux3]*(us[aux1][aux3]-u[aux2][aux3]));
	i=i+1;
    }

    j=j+1;
}

  // cost function
  i=0;
  mpi_i = (myid) * (nx/numprocs);
  mpi_nx = (myid+1) * (nx/numprocs);
  cost_mpi=0;
  while (mpi_i < mpi_nx) {
    j=0;
    while (j<n) {
      cost_mpi=cost_mpi+0.5*(u[mpi_i][j]-uob[mpi_i][j])*(u[mpi_i][j]-uob[mpi_i][j]);
      j=j+1;
    }
  mpi_i=mpi_i+1;
  cost = cost_mpi;
  }
  for(i=1; i < numprocs ; i=i+1) {
      if(myid == 0) {
	  target=i;
	  nelements=1;
	  buf[0] = 0;
	  MPI_Irecv(buf, nelements, ampi_double, target, tag, ampi_comm, request);
	  MPI_Wait(request,status);
//	  print_num(buf[0]);
	  cost = cost + buf[0];
      } else {
	  if(myid == i) {
	      target=0;
	      nelements=1;
	      buf[0]=cost_mpi;
	      MPI_Isend(buf, nelements, ampi_double, target, tag, ampi_comm, request);
	      MPI_Wait(request,status);
	      cost_mpi = 0;
	  }
      }
  }


  // save nonlinear solutions to the basic field
  i=0;
  while (i<nx) {
    j=0;
    while (j<n) {
      ub[i][j]=u[i][j];
      j=j+1;
    }
    i=i+1;
  }
}
