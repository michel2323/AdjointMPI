#include <mpi.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ampi.h>
using namespace std;

int  myid;
int  numprocs;
int ampi_comm_world=0; 
int ampi_double=0; 
extern void distribute(double *, double *, int &);
extern void speelpenning(double *, int &, double *);
extern void a1_speelpenning(int &, double *, double *, int &, double *, double *);
extern void t2_a1_speelpenning(int &, double *, double *, double *, double *, int &, double *, double *, double *, double*);
extern void t2_a1_distribute(int &, double *, double *, double *, double *, double *, double *, double* , double*, int &);

int main(int argc, char** argv){

    // init

    int namelen;
    char processor_name[AMPI_MAX_PROCESSOR_NAME];
    MPI_Init(&argc,&argv);
    MPI_Comm_size(AMPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(AMPI_COMM_WORLD,&myid);
    MPI_Get_processor_name(processor_name,&namelen);
    MPI_Status status;
    int n=4;
    int bmode=1;
    int mpi_x_size = n/numprocs;
    // mpi_x
    double * mpi_x = new double[mpi_x_size];
    double * b1_mpi_x = new double[mpi_x_size];
    double * d2_mpi_x = new double[mpi_x_size];
    double * d2_b1_mpi_x = new double[mpi_x_size];
    // buffer for distribution ()
    //buf_mpi = new double[mpi_x_size]; 
    //d2_buf_mpi = new double[mpi_x_size];
    //b1_buf_mpi = new double[mpi_x_size];
    //d2_b1_buf_mpi = new double[mpi_x_size];
    //double * cp_buf_mpi = new double[mpi_x_size];
    // myres
    double myres = 1;
    double b1_myres = 1;
    double d2_myres = 0;
    double d2_b1_myres = 0;
    // x
    double * x = new double[n];
    double * b1_x = new double[n];
    double * d2_x = new double[n];
    double * d2_b1_x = new double[n];
    double * cp_x = new double[n];
    // buf for speelpenning ()


    // UNDIFFERENTIATED CODE

    if(myid==0) {
	cout << "Original code" << endl;
	cout << "-----------------------" << endl;
    }
    if(myid == 0) {
	for(int i = 0 ; i < n ; i++)
	    x[i] = i + 2;
    }
    else {
	for(int i = 0 ; i < n ; i++)
	    x[i] = 0;
    }


    // distribute x

    // Checkpointing
    for(int i = 0 ; i < n; i++) {
	cp_x[i] = x[i]; 
    }
    distribute(x, mpi_x, mpi_x_size);    
    //distribute(x, mpi_x, mpi_x_size, buf_mpi);    

    // speelpenning

    //speelpenning(&myres, mpi_x_size, mpi_x, buf_mpi);
    speelpenning(&myres, mpi_x_size, mpi_x);

    if(myid==0) {
	cout << "Result: " <<myres << endl;
    }

    // Seeding 
    //
    myres = 1;
    b1_myres = 1;
    d2_b1_myres = 0;
    d2_myres = 0;
    bmode=1;
    for(int i=0 ; i < mpi_x_size ; i++) {
	d2_mpi_x[i]=0;
	b1_mpi_x[i]=0;
	d2_b1_mpi_x[i]=0;
    }
    if(myid==0)
	d2_mpi_x[0]=1;

    //Speelpenning

    t2_a1_speelpenning(bmode, &myres, &d2_myres, &b1_myres, &d2_b1_myres, mpi_x_size, mpi_x, d2_mpi_x, b1_mpi_x, d2_b1_mpi_x);
    //d2_b1_speelpenning(bmode, &myres, &d2_myres, &b1_myres, &d2_b1_myres, mpi_x_size, mpi_x, d2_mpi_x, b1_mpi_x, d2_b1_mpi_x, buf_mpi, d2_buf_mpi, b1_buf_mpi, d2_b1_buf_mpi);

    // Output
    //cout << "First order adjoints:" << endl;
    //cout << "-----------------------" << endl;
    //for(int i=0 ; i < mpi_x_size ; i++) {
	//cout << b1_mpi_x[i] << " " ;
    //}
    //cout << endl;
    //cout << "Forward over Reverse wrt. x[0]: " << endl;
    //cout << "-----------------------" << endl;
    //for(int i=0 ; i < mpi_x_size ; i++) {
	//cout << d2_b1_mpi_x[i] << " " ;
    //}
    //cout << endl;
    //Restore Checkpoint for distribute ()
    for(int i = 0 ; i < mpi_x_size; i++) {
	//buf_mpi[i] = 0; 
	//d2_b1_buf_mpi[i] = 0; 
	//b1_buf_mpi[i] = 0; 
	//d2_buf_mpi[i] = 0; 
	mpi_x[i] = 0; 
	d2_mpi_x[i] = 0; 
    }
    bmode = 1;
    for(int i=0 ; i < n ; i++) {
	x[i]=cp_x[i];
	d2_x[i]=0;
	b1_x[i]=0;
	d2_b1_x[i]=0;
    }

    //d2_b1_distribute(bmode, x, d2_x, b1_x, d2_b1_x, mpi_x, d2_mpi_x, b1_mpi_x, d2_b1_mpi_x,mpi_x_size, buf_mpi, d2_buf_mpi, b1_buf_mpi, d2_b1_buf_mpi);    
    t2_a1_distribute(bmode, x, d2_x, b1_x, d2_b1_x, mpi_x, d2_mpi_x, b1_mpi_x, d2_b1_mpi_x,mpi_x_size);    
    // Output
    if(myid == 0) {
	cout << "First order adjoints:" << endl;
	cout << "-----------------------" << endl;
	for(int i=0 ; i < n ; i++) {
	    cout << b1_x[i] << " " ;
	}
	cout << endl;
	cout << "Forward over Reverse wrt. x[0]: " << endl;
	cout << "-----------------------" << endl;
	for(int i=0 ; i < n ; i++) {
	    cout << d2_b1_x[i] << " " ;
	}
	cout << endl;
    }
    delete [] mpi_x;
    delete [] b1_mpi_x;
    delete [] d2_mpi_x;
    delete [] d2_b1_mpi_x;

    delete [] x;
    delete [] b1_x;
    delete [] d2_x;
    delete [] d2_b1_x;
    delete [] cp_x;

    MPI_Barrier(AMPI_COMM_WORLD);
    MPI_Finalize();
} 

