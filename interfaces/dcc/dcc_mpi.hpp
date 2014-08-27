#ifndef DCC_MPI_H
#define DCC_MPI_H
#include <mpi.h>
#include <math.h>
#include <cassert>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <stack>
#include "dcc_alloc.hpp"
#define AMPI
extern "C" {
#include "ampi.h"
}
using namespace std;


// Communciation

void MPI_Irecv(double *buf, int &nelements, int &ampi_double, int &target, int &tag, int &ampi_comm, int *request);
void MPI_Isend(double *buf, int &nelements, int &ampi_double, int &target, int &tag, int &ampi_comm, int *request);
void MPI_Wait(double *buf, int *request, int *status);
void MPI_Reduce(double *sendbuf, double *recvbuf, int &nelements, int &ampi_double, int &op, int &root, int &ampi_comm);

// First order

// Communciation

void a1_MPI_Irecv(int bmode, double *buf, double *&b1_buf, int &nelements, int &ampi_double, int &target, int &tag, int &ampi_comm, int *request);
void a1_MPI_Isend(int bmode, double *buf, double *&b1_buf, int &nelements, int &ampi_double, int &target, int &tag, int &ampi_comm, int *request);
void a1_MPI_Wait(int bmode, int *request, int *status);
void a1_MPI_Reduce(int bmode, double *sendbuf, double *b1_sendbuf, double *recvbuf, double *b1_recvbuf, int &nelements, int &ampi_double, int &op, int &root, int &ampi_comm);


// Forward over Reverse

// Communciation

void t2_a1_MPI_Irecv(int bmode, double * buf, double * d2_buf, double *&b1_buf, double *&d2_b1_buf, int &nelements, int &ampi_double, int &target, int &tag, int &ampi_comm, int *request);
void t2_a1_MPI_Isend(int bmode, double * buf, double * d2_buf, double *&b1_buf, double *&d2_b1_buf, int &nelements, int &ampi_double, int &target, int &tag, int &ampi_comm, int *request);
void t2_a1_MPI_Wait(int bmode, int *request, int *status);
void t2_a1_MPI_Reduce(int bmode, double *sendbuf, double *d2_sendbuf, double *b1_sendbuf, double *d2_b1_sendbuf, double *recvbuf, double *d2_recvbuf, double *b1_recvbuf, double *d2_b1_recvbuf, int &nelements, int &ampi_double, int &op, int &root, int &ampi_comm);

// Buffer allocation

void t2_a1_dcc_malloc(int bmode, double *&buf, double *&d2_buf, double *&b1_buf, double *&d2_b1_buf, int &size);
void t2_a1_dcc_free(int bmode, double *&buf, double *&d2_buf, double *&b1_buf, double *&d2_b1_buf);

void t2_a1_print_num(int bmode, double &num, double &d2_num, double &b1_num, double &d2_b1_num);
#endif
