#include <mpi.h>
#include<iostream>
#include<fstream>
#include<cmath>
#include "dcc_mpi.hpp"
using namespace std;

int myid;
int numprocs;

#include "t2_a1_f.c"


int main(int argc, char* argv[]) {
  int namelen;
  int mpi_x_size;
  char processor_name[AMPI_MAX_PROCESSOR_NAME];
  MPI_Init(&argc,&argv);
  MPI_Comm_size(AMPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(AMPI_COMM_WORLD,&myid);
  MPI_Get_processor_name(processor_name,&namelen);
  MPI_Status status;
  int nx=10; 
  int n=100; 

  double cost, dx, dt, r, dtdx, c0, c1;
  double **uob, **ub, **u, **us, *ui, *buf;
  uob=new double*[nx]; ub=new double*[nx]; u=new double*[nx]; us=new double*[nx]; 
  for (int i=0;i<nx;i++) {
    uob[i]=new double[n];
    ub[i]=new double[n];
    us[i]=new double[n];
    u[i]=new double[n];
  }
  ui=new double[nx];
  buf=new double[1];

  // tangent-linear model
  double d1_cost, d1_dx, d1_dt, d1_r, d1_dtdx, d1_c0, d1_c1;
  double **d1_uob, **d1_ub, **d1_u, **d1_us, *d1_ui, *d1_buf;
  d1_uob=new double*[nx]; d1_ub=new double*[nx]; d1_u=new double*[nx]; d1_us=new double*[nx]; 
  for (int i=0;i<nx;i++) {
    d1_uob[i]=new double[n];
    d1_ub[i]=new double[n];
    d1_us[i]=new double[n];
    d1_u[i]=new double[n];
  }
  d1_ui=new double[nx];
  d1_buf=new double[1];
  // adjoint model
  double b0_cost, b0_dx, b0_dt, b0_r, b0_dtdx, b0_c0, b0_c1;
  double **b0_uob, **b0_ub, **b0_u, **b0_us, *b0_ui, *b0_buf;
  b0_uob=new double*[nx]; b0_ub=new double*[nx]; b0_u=new double*[nx]; b0_us=new double*[nx]; 
  for (int i=0;i<nx;i++) {
    b0_uob[i]=new double[n];
    b0_ub[i]=new double[n];
    b0_us[i]=new double[n];
    b0_u[i]=new double[n];
  }
  b0_ui=new double[nx];
  b0_buf=new double[1];

  // Hessian
  double ** H= new double*[n];
  for (int i=0;i<n;i++) H[i]=new double[n];



  // second-order adjoint model
  double d1_b0_cost, d1_b0_dx, d1_b0_dt, d1_b0_r, d1_b0_dtdx, d1_b0_c0, d1_b0_c1;
  double **d1_b0_uob, **d1_b0_ub, **d1_b0_u, **d1_b0_us, *d1_b0_ui, *d1_b0_buf;
  d1_b0_uob=new double*[nx]; d1_b0_ub=new double*[nx]; d1_b0_u=new double*[nx]; d1_b0_us=new double*[nx]; 
  for (int i=0;i<nx;i++) {
    d1_b0_uob[i]=new double[n];
    d1_b0_ub[i]=new double[n];
    d1_b0_u[i]=new double[n];
    d1_b0_us[i]=new double[n];
  }
  for(int i=0;i<nx;i++) 
    for(int j=0;j<n;j++) { 
      d1_b0_us[i][j]=0;
      b0_us[i][j]=0;
    }
  d1_b0_ui=new double[nx];
  d1_b0_buf=new double[1];
  ofstream soadm_out;
    if(myid == 0)
    	soadm_out.open("soadm0.out");
    if(myid == 1)
    	soadm_out.open("soadm1.out");
  int bmode0=1;
  for (int j=0;j<nx;j++) {
    for (int i=0;i<nx;i++) { 
      ui[i]=1.0/nx; 
      b0_ui[i]=0; 
      d1_ui[i]=0; 
      d1_b0_ui[i]=0; 
      for (int j2=0;j2<n;j2++) { 
        uob[i][j2]=1.0/(i*j2+1); 
        b0_uob[i][j2]=0; 
        d1_uob[i][j2]=0; 
        d1_b0_uob[i][j2]=0; 
        u[i][j2]=0; 
        b0_u[i][j2]=0; 
        d1_u[i][j2]=0; 
        d1_b0_u[i][j2]=0; 
        ub[i][j2]=0; 
        b0_ub[i][j2]=0; 
        d1_ub[i][j2]=0; 
        d1_b0_ub[i][j2]=0; 
      }
    }
    d1_dx=0;
    b0_dx=0;
    d1_b0_dx=0;
    d1_dt=0;
    b0_dt=0;
    d1_b0_dt=0;
    d1_r=0;
    b0_r=0;
    d1_b0_r=0;
    d1_dtdx=0;
    b0_dtdx=0;
    d1_b0_dtdx=0;
    d1_c0=0;
    b0_c0=0;
    d1_b0_c0=0;
    d1_c1=0;
    b0_c1=0;
    d1_b0_c1=0;
    d1_cost=0;
    d1_b0_cost=0;
    if(myid==0) {
	b0_cost=1;
	d1_ui[j]=1;
    }
    else {
	b0_cost=0;
	d1_ui[j]=1;
    }
    t2_a1_f(bmode0,nx, n, 
              cost, d1_cost, b0_cost, d1_b0_cost,
              uob, d1_uob, b0_uob, d1_b0_uob,
              ub, d1_ub, b0_ub, d1_b0_ub,
              us, d1_us, b0_us, d1_b0_us,
              u, d1_u, b0_u, d1_b0_u,
              ui, d1_ui, b0_ui, d1_b0_ui,
              dx, d1_dx, b0_dx, d1_b0_dx,
              dt, d1_dt, b0_dt, d1_b0_dt,
              r, d1_r, b0_r, d1_b0_r,
              dtdx, d1_dtdx, b0_dtdx, d1_b0_dtdx,
              c0, d1_c0, b0_c0, d1_b0_c0,
              c1, d1_c1, b0_c1, d1_b0_c1,
	      buf, d1_buf, b0_buf, d1_b0_buf);
    for (int i=0;i<nx;i++) H[j][i]=d1_b0_ui[i];
  }
  double *res=new double[nx];
  for (int j=0;j<nx;j++) res[j]=0; 
  for (int j=0;j<nx;j++) { 
    MPI_Reduce(H[j],res,nx,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    if(myid==0)
      for (int i=0;i<nx;i++) 
	soadm_out << "H[" << j << "," << i << "]=" << res[i] <<endl;
  }
  delete [] res;

  // deallocate
 
  for (int i=0;i<nx;i++) {
    delete [] H[i];
    delete [] b0_uob[i];
    delete [] d1_uob[i];
    delete [] uob[i];
    delete [] b0_ub[i];
    delete [] d1_ub[i];
    delete [] ub[i];
    delete [] b0_u[i];
    delete [] d1_u[i];
    delete [] u[i];
    delete [] b0_us[i];
    delete [] d1_us[i];
    delete [] us[i];
  }
  delete [] H;
  delete [] b0_ui; delete [] b0_buf; delete [] b0_u; delete [] b0_ub; delete [] b0_uob;
  delete [] d1_ui; delete [] d1_buf; delete [] d1_u; delete [] d1_ub; delete [] d1_uob;
  delete [] ui; delete [] u; delete [] ub; delete [] uob; delete [] us;
  MPI_Barrier(AMPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}
