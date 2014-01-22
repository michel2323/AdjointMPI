#include "dcc_mpi.hpp"
#include <cassert>
#include <iostream>

using namespace std;

#define STACK_SIZE 200 

// Request stacks

MPI_Request mpi_rst[STACK_SIZE];
int mpi_rstc = 0;

AMPI_Request ampi_rst[STACK_SIZE];
int ampi_rstc = 0;

double *buf_ampi_rst[STACK_SIZE];
double **ptr_buf_ampi_rst[STACK_SIZE];

// DUMMY MPI CALLS FOR DCC ADAPTATION

// Communication

void MPI_Irecv(double * buf, int &nelements, int &ampi_double, int &target, int &tag, int &ampi_comm, int *request) {
    *request = mpi_rstc;
    MPI_Irecv( buf, nelements, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &mpi_rst[mpi_rstc++]);
}

void MPI_Isend(double * buf, int &nelements, int &ampi_double, int &target, int &tag, int &ampi_comm, int *request) {
    *request = mpi_rstc;
    MPI_Isend(buf, nelements, MPI_DOUBLE, target, 0 , MPI_COMM_WORLD, &mpi_rst[mpi_rstc++]);
}

void MPI_Wait(int *request, int *status) {
    MPI_Status mpi_status;
    MPI_Wait(&mpi_rst[*request], &mpi_status);
}

void MPI_Reduce(double *sendbuf, double *recvbuf, int &nelements, int &ampi_double, int &op, int &root, int &ampi_comm) {
    MPI_Reduce(sendbuf, recvbuf, nelements, MPI_DOUBLE, MPI_SUM , root, MPI_COMM_WORLD);
}

// Buffer allocation

void print_num(double &num) {
    //cout << "num: " << num << endl;
}

// FIRST ORDER ADJOINTS

// Communication

void a1_MPI_Irecv(int bmode, double * buf, double *&a1_buf, int &nelements, int &ampi_double, int &target, int &tag, int &ampi_comm, int *request) {
    int ierr=0;
    if(bmode == 2) {
        *request = ampi_rstc;
	ampi_rst[ampi_rstc].request=(MPI_Request*) malloc(sizeof(MPI_Request));
	ierr = AMPI_Irecv_f( buf, nelements, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &ampi_rst[ampi_rstc]);
	ptr_buf_ampi_rst[ampi_rstc] = &a1_buf;
        ampi_rstc++;
	assert(ampi_rstc<STACK_SIZE);
    }
    if(bmode == 1) {
	ierr = AMPI_Irecv_b( a1_buf, nelements, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &ampi_rst[*request]);
	free(ampi_rst[*request].request);
	*request = *request - 1;
    }
}

void a1_MPI_Isend(int bmode, double * buf, double *&a1_buf, int &nelements, int &ampi_double, int &target, int &tag, int &ampi_comm, int *request) {
    int ierr=0;
    if(bmode == 2) {
        *request = ampi_rstc;
	ampi_rst[ampi_rstc].request=(MPI_Request*) malloc(sizeof(MPI_Request));
	ierr = AMPI_Isend_f( buf, nelements, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &ampi_rst[ampi_rstc]);
	ptr_buf_ampi_rst[ampi_rstc] = &a1_buf;
        ampi_rstc++;
	assert(ampi_rstc<STACK_SIZE);
    }
    if(bmode == 1) {
	ierr = AMPI_Isend_b( a1_buf, nelements, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &ampi_rst[*request]);
	free(ampi_rst[*request].request);
	*request = *request - 1;
    }
}

void a1_MPI_Wait(int bmode, int *request, int *status) {
    MPI_Status mpi_status;
    int ierr = 0;
    if(bmode == 2) {
	ierr = AMPI_Wait_f(&ampi_rst[*request], &mpi_status);
    }
    if(bmode == 1) {
	ampi_rst[*request].a = *ptr_buf_ampi_rst[*request];
	ierr = AMPI_Wait_b(&ampi_rst[*request], &mpi_status);
    }
}

void a1_MPI_Reduce(int bmode, double *sendbuf, double *a1_sendbuf, double *recvbuf, double *a1_recvbuf, int &nelements, int &ampi_double, int &op, int &root, int &ampi_comm) {
    if(bmode == 2) {
	if(op == AMPI_SUM)
	    AMPI_Reduce_f(sendbuf, recvbuf, nelements, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
	if(op == AMPI_PROD)
	    AMPI_Reduce_f(sendbuf, recvbuf, nelements, MPI_DOUBLE, MPI_PROD, root, MPI_COMM_WORLD);
    }
    if(bmode == 1) {
	if(op == AMPI_SUM)
	    AMPI_Reduce_b(a1_sendbuf, a1_recvbuf, nelements, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
	if(op == AMPI_PROD)
	    AMPI_Reduce_b(a1_sendbuf, a1_recvbuf, nelements, MPI_DOUBLE, MPI_PROD, root, MPI_COMM_WORLD);
    }
}

// FORWARD OVER REVERSE

// Communication

void t2_a1_MPI_Irecv(int bmode, double * buf, double * t2_buf, double *&a1_buf, double *&t2_a1_buf, int &nelements, int &ampi_double, int &target, int &tag, int &ampi_comm, int *request) {
    int ierr=0;
    if(bmode == 2) {
        *request = ampi_rstc;
	ampi_rst[ampi_rstc].request=(MPI_Request*) malloc(sizeof(MPI_Request));
	ierr = AMPI_Irecv_f( buf, nelements, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &ampi_rst[ampi_rstc]);
	ptr_buf_ampi_rst[ampi_rstc] = &a1_buf;
        ampi_rstc++;
	assert(ampi_rstc<STACK_SIZE);
	ampi_rst[ampi_rstc].request=(MPI_Request*) malloc(sizeof(MPI_Request));
	ierr = AMPI_Irecv_f( t2_buf, nelements, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &ampi_rst[ampi_rstc]);
	ptr_buf_ampi_rst[ampi_rstc] = &t2_a1_buf;
	ampi_rstc++;
	assert(ampi_rstc<STACK_SIZE);
    }
    if(bmode == 1) {
	ierr = AMPI_Irecv_b( t2_a1_buf, nelements, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &ampi_rst[*request+1]);
	free(ampi_rst[*request+1].request);
	ierr = AMPI_Irecv_b( a1_buf, nelements, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &ampi_rst[*request]);
	free(ampi_rst[*request].request);
    }
}

void t2_a1_MPI_Isend(int bmode, double * buf, double * t2_buf, double *&a1_buf, double *&t2_a1_buf, int &nelements, int &ampi_double, int &target, int &tag, int &ampi_comm, int *request) {
    int ierr=0;
    if(bmode == 2) {
	*request=ampi_rstc;
	ampi_rst[ampi_rstc].request=(MPI_Request*) malloc(sizeof(MPI_Request));
	ierr = AMPI_Isend_f( buf, nelements, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &ampi_rst[ampi_rstc]);
	ptr_buf_ampi_rst[ampi_rstc] = &a1_buf;
        ampi_rstc++;
	assert(ampi_rstc<STACK_SIZE);
	ampi_rst[ampi_rstc].request=(MPI_Request*) malloc(sizeof(MPI_Request));
	ierr = AMPI_Isend_f( t2_buf, nelements, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &ampi_rst[ampi_rstc]);
	ptr_buf_ampi_rst[ampi_rstc] = &t2_a1_buf;
	ampi_rstc++;
	assert(ampi_rstc<STACK_SIZE);
    }
    if(bmode == 1) {
	ierr = AMPI_Isend_b( t2_a1_buf, nelements, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &ampi_rst[*request+1]);
	free(ampi_rst[*request+1].request);
	ierr = AMPI_Isend_b( a1_buf, nelements, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &ampi_rst[*request]);
	free(ampi_rst[*request].request);
    }
}

void t2_a1_MPI_Wait(int bmode, int *request, int *status) {
    MPI_Status mpi_status;
    int ierr = 0;
    if(bmode == 2) {
	ierr = AMPI_Wait_f(&ampi_rst[*request+1], &mpi_status);
	ierr = AMPI_Wait_f(&ampi_rst[*request], &mpi_status);
    }
    if(bmode == 1) {
	ampi_rst[*request+1].a = *ptr_buf_ampi_rst[*request+1];
	ierr = AMPI_Wait_b(&ampi_rst[*request+1], &mpi_status);
	ampi_rst[*request].a = *ptr_buf_ampi_rst[*request];
	ierr = AMPI_Wait_b(&ampi_rst[*request], &mpi_status);
    }
}

void t2_a1_MPI_Reduce(int bmode, double *sendbuf, double *t2_sendbuf, double *a1_sendbuf, double *t2_a1_sendbuf, double *recvbuf, double *t2_recvbuf, double *a1_recvbuf, double *t2_a1_recvbuf, int &nelements, int &ampi_double, int &op, int &root, int &ampi_comm) {
    if(bmode == 2) {
	AMPI_Reduce_f(sendbuf, recvbuf, nelements, MPI_DOUBLE, MPI_SUM , root, MPI_COMM_WORLD);
	AMPI_Reduce_f(t2_sendbuf, t2_recvbuf, nelements, MPI_DOUBLE, MPI_SUM , root, MPI_COMM_WORLD);
    }
    if(bmode == 1) {
	AMPI_Reduce_b(t2_a1_sendbuf, t2_a1_recvbuf, nelements, MPI_DOUBLE, MPI_SUM , root, MPI_COMM_WORLD);
	AMPI_Reduce_b(a1_sendbuf, a1_recvbuf, nelements, MPI_DOUBLE, MPI_SUM , root, MPI_COMM_WORLD);
    }
}

void a1_print_num(int bmode, double &num, double &a1_num) {
    if(bmode == 2)
	cout << "num: " << num << endl;
    if(bmode == 1)
	cout << "a1_num: " << a1_num << endl;
}

void t2_a1_print_num(int bmode, double &num, double &t2_num, double &a1_num, double &t2_a1_num) {
    if(bmode == 2) {
	cout << "num: " << num << endl;
	cout << "t2_num: " << t2_num << endl;
    }
    if(bmode == 1) {
	cout << "a1_num: " << a1_num << endl;
	cout << "t2_a1_num: " << t2_a1_num << endl;
    }
}
