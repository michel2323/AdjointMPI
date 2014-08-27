#include "dcc_mpi.hpp"

using namespace std;

// Request stacks

MPI_Request mpi_rst[200];
int mpi_rstc = 0;

AMPI_Request ampi_rst[200];
int ampi_rstc = 0;

double *buf_ampi_rst[200];
double **ptr_buf_ampi_rst[200];

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

void a1_MPI_Irecv(int bmode, double * buf, double *&b1_buf, int &nelements, int &ampi_double, int &target, int &tag, int &ampi_comm, int *request) {
    int ierr=0;
    if(bmode == 2) {
        *request = ampi_rstc;
	ierr = AMPI_Irecv_f( buf, nelements, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &ampi_rst[ampi_rstc]);
	ptr_buf_ampi_rst[ampi_rstc] = &b1_buf;
        ampi_rstc++;
    }
    if(bmode == 1) {
	ierr = AMPI_Irecv_b( b1_buf, nelements, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &ampi_rst[*request]);
	*request = *request - 1;
    }
}

void a1_MPI_Isend(int bmode, double * buf, double *&b1_buf, int &nelements, int &ampi_double, int &target, int &tag, int &ampi_comm, int *request) {
    int ierr=0;
    if(bmode == 2) {
        *request = ampi_rstc;
	ierr = AMPI_Isend_f( buf, nelements, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &ampi_rst[ampi_rstc]);
	ptr_buf_ampi_rst[ampi_rstc] = &b1_buf;
        ampi_rstc++;
    }
    if(bmode == 1) {
	ierr = AMPI_Isend_b( b1_buf, nelements, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &ampi_rst[*request]);
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

void a1_MPI_Reduce(int bmode, double *sendbuf, double *b1_sendbuf, double *recvbuf, double *b1_recvbuf, int &nelements, int &ampi_double, int &op, int &root, int &ampi_comm) {
    if(bmode == 2) {
	if(op == AMPI_SUM)
	    AMPI_Reduce_f(sendbuf, recvbuf, nelements, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
	if(op == AMPI_PROD)
	    AMPI_Reduce_f(sendbuf, recvbuf, nelements, MPI_DOUBLE, MPI_PROD, root, MPI_COMM_WORLD);
    }
    if(bmode == 1) {
	if(op == AMPI_SUM)
	    AMPI_Reduce_b(b1_sendbuf, b1_recvbuf, nelements, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
	if(op == AMPI_PROD)
	    AMPI_Reduce_b(b1_sendbuf, b1_recvbuf, nelements, MPI_DOUBLE, MPI_PROD, root, MPI_COMM_WORLD);
    }
}

// FORWARD OVER REVERSE

// Communication

void t2_a1_MPI_Irecv(int bmode, double * buf, double * d2_buf, double *&b1_buf, double *&d2_b1_buf, int &nelements, int &ampi_double, int &target, int &tag, int &ampi_comm, int *request) {
    int ierr=0;
    if(bmode == 2) {
        *request = ampi_rstc;
	ierr = AMPI_Irecv_f( buf, nelements, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &ampi_rst[ampi_rstc]);
	ptr_buf_ampi_rst[ampi_rstc] = &b1_buf;
        ampi_rstc++;
	ierr = AMPI_Irecv_f( d2_buf, nelements, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &ampi_rst[ampi_rstc]);
	ptr_buf_ampi_rst[ampi_rstc] = &d2_b1_buf;
	ampi_rstc++;
    }
    if(bmode == 1) {
	ierr = AMPI_Irecv_b( d2_b1_buf, nelements, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &ampi_rst[*request+1]);
	ierr = AMPI_Irecv_b( b1_buf, nelements, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &ampi_rst[*request]);
    }
}

void t2_a1_MPI_Isend(int bmode, double * buf, double * d2_buf, double *&b1_buf, double *&d2_b1_buf, int &nelements, int &ampi_double, int &target, int &tag, int &ampi_comm, int *request) {
    int ierr=0;
    if(bmode == 2) {
	*request=ampi_rstc;
	ierr = AMPI_Isend_f( buf, nelements, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &ampi_rst[ampi_rstc]);
	ptr_buf_ampi_rst[ampi_rstc] = &b1_buf;
        ampi_rstc++;
	ierr = AMPI_Isend_f( d2_buf, nelements, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &ampi_rst[ampi_rstc]);
	ptr_buf_ampi_rst[ampi_rstc] = &d2_b1_buf;
	ampi_rstc++;
    }
    if(bmode == 1) {
	ierr = AMPI_Isend_b( d2_b1_buf, nelements, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &ampi_rst[*request+1]);
	ierr = AMPI_Isend_b( b1_buf, nelements, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &ampi_rst[*request]);
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

void t2_a1_MPI_Reduce(int bmode, double *sendbuf, double *d2_sendbuf, double *b1_sendbuf, double *d2_b1_sendbuf, double *recvbuf, double *d2_recvbuf, double *b1_recvbuf, double *d2_b1_recvbuf, int &nelements, int &ampi_double, int &op, int &root, int &ampi_comm) {
    if(bmode == 2) {
	AMPI_Reduce_f(sendbuf, recvbuf, nelements, MPI_DOUBLE, MPI_SUM , root, MPI_COMM_WORLD);
	AMPI_Reduce_f(d2_sendbuf, d2_recvbuf, nelements, MPI_DOUBLE, MPI_SUM , root, MPI_COMM_WORLD);
    }
    if(bmode == 1) {
	AMPI_Reduce_b(d2_b1_sendbuf, d2_b1_recvbuf, nelements, MPI_DOUBLE, MPI_SUM , root, MPI_COMM_WORLD);
	AMPI_Reduce_b(b1_sendbuf, b1_recvbuf, nelements, MPI_DOUBLE, MPI_SUM , root, MPI_COMM_WORLD);
    }
}

void a1_print_num(int bmode, double &num, double &b1_num) {
    if(bmode == 2)
	cout << "num: " << num << endl;
    if(bmode == 1)
	cout << "b1_num: " << b1_num << endl;
}

void t2_a1_print_num(int bmode, double &num, double &d2_num, double &b1_num, double &d2_b1_num) {
    if(bmode == 2) {
	cout << "num: " << num << endl;
	cout << "d2_num: " << d2_num << endl;
    }
    if(bmode == 1) {
	cout << "b1_num: " << b1_num << endl;
	cout << "d2_b1_num: " << d2_b1_num << endl;
    }
}
