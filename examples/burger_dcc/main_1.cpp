#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include "dcc_mpi.hpp"

int myid;
int numprocs;

using namespace std;

//#include "d1_f.c"
#include "a1_f.c"

int main(int argc, char* argv[]) {
    // init

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
    ofstream adm_out;
    int bmode0;
    if(myid == 0)
    	adm_out.open("adm0.out");
    if(myid == 1)
    	adm_out.open("adm1.out");
    for (int i=0;i<nx;i++) { 
	ui[i]=1.0/nx; 
	b0_ui[i]=0; 
    }
    for (int j1=0;j1<nx;j1++) { 
	for (int j2=0;j2<n;j2++) { 
	    uob[j1][j2]=1.0/(j1*j2+1); 
	    b0_uob[j1][j2]=0; 
	    u[j1][j2]=0; 
	    b0_u[j1][j2]=0; 
	    ub[j1][j2]=0; 
	    b0_ub[j1][j2]=0; 
	    b0_us[j1][j2]=0;
	}
    }
    if(myid == 0)
	b0_cost=1;
    else
	b0_cost=0;
    bmode0=1;
    a1_f(bmode0, nx, n, cost, b0_cost, uob, b0_uob, ub, b0_ub, us, b0_us, u, b0_u, 
	    ui, b0_ui, dx, b0_dx, dt, b0_dt, r, b0_r, dtdx, b0_dtdx, 
	    c0, b0_c0, c1, b0_c1, buf, b0_buf);
    double *res=new double[nx];
    MPI_Reduce(b0_ui, res, nx, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
    if(myid==0)
      for (int i=0;i<nx;i++) adm_out << "g[" << i << "]=" << res[i] <<endl;
    delete [] res;

    for (int i=0;i<nx;i++) {
	delete [] b0_uob[i];
	delete [] uob[i];
	delete [] b0_ub[i];
	delete [] ub[i];
	delete [] b0_us[i];
	delete [] us[i];
	delete [] b0_u[i];
	delete [] u[i];
    }
    delete [] b0_ui; delete [] b0_u; delete [] b0_ub; delete [] b0_us; delete [] b0_uob;
    delete [] ui; delete [] u; delete [] ub; delete [] us; delete [] uob;
    MPI_Barrier(AMPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
