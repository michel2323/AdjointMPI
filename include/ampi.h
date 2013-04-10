#ifndef AMPI_H
#define AMPI_H
/* Basic AMPI C library used for overloading and source transformation. All MPI routines
 * are subdivided into their forward _f and backward _b counterpart. The forward routines
 * are called during the forward/taping run. The backward routines are called during the
 * reverse/interpretation run.
 */
#include <mpi.h>
#include <ampi_stack.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#define AMPI_IS 1
#define AMPI_IR 2

#define AMPI_COMM_WORLD MPI_COMM_WORLD
#define AMPI_MAX_PROCESSOR_NAME MPI_MAX_PROCESSOR_NAME
#define AMPI_Status MPI_Status
#define AMPI_DOUBLE MPI_DOUBLE
/*#define INT64 double**/
#define INT64 int

/* Reduce operations */

#define AMPI_SUM   1
#define AMPI_PROD  2
#define AMPI_MIN   3
#define AMPI_MAX   4


/* AMPI request, replacing the MPI request */
typedef struct AMPI_Request {
    MPI_Request *request;
    MPI_Status status;
    MPI_Comm comm;

    /*special for tape*/
    void *buf;

    int tag;
    double *v;
    double *a;
    int va;
    int oc;
    int dest;
    int size;
    int aw;
} AMPI_Request;

typedef struct AMPI_Tupel {
    double v;
    int j;
} AMPI_Tupel;

/* Non communication routines (Init, Finalize...)*/

int AMPI_Init_f(int*, char***);
int AMPI_Init_b(int*, char***);
int AMPI_Comm_size(MPI_Comm, int*);
int AMPI_Comm_rank(MPI_Comm, int*);
int AMPI_Get_processor_name(char*, int*);
int AMPI_Barrier(MPI_Comm);
int AMPI_Finalize_f();
int AMPI_Finalize_b();

/* Blocking communication */

int AMPI_Send_f(double *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
int AMPI_Recv_f(double *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Status *status);
int AMPI_Send_b(double *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
int AMPI_Recv_b(double *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Status *status);

/* Non blocking communication */

int AMPI_Isend_f(double *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, AMPI_Request *request);
int AMPI_Irecv_f(double *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, AMPI_Request *request);
int AMPI_Wait_f(AMPI_Request *request , MPI_Status *status);
int AMPI_Isend_b(double *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, AMPI_Request *request);
int AMPI_Irecv_b(double *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, AMPI_Request *request);
int AMPI_Wait_b(AMPI_Request *request, MPI_Status *status);
int AMPI_Waitall_f(int count, AMPI_Request *requests, MPI_Status *status);
int AMPI_Waitall_b(int count, AMPI_Request *requests);

/* Anti Waitall. Experimental */

int AMPI_Awaitall_f(int count, AMPI_Request *requests, MPI_Status *status);
int AMPI_Awaitall_b(int count, AMPI_Request *request);

/* Collective */

int AMPI_Bcast_f(double *buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm);
int AMPI_Bcast_b(double *buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm);
int AMPI_Reduce_f(double *sendbuf, double *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
int AMPI_Reduce_b(double *sendbuf, double *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
int AMPI_Allreduce_f(double *sendbuf, double *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int AMPI_Allreduce_b(double *sendbuf, double *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int AMPI_Gather_f(void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm);
int AMPI_Gather_b(void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm);
int AMPI_Scatter_f(void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm);
int AMPI_Scatter_b(void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm);

/* Special comms */

int AMPI_Sendrecv_replace_f(double *buf, int count, MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag, MPI_Comm comm, MPI_Status *status);
int AMPI_Sendrecv_replace_b(double *buf, int count, MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag, MPI_Comm comm, MPI_Status *status);

/* MPI derived datatypes for AMPI */

int AMPI_Reduc_Tupel();
void AMPI_Tupel_Max(void *invec_in, void *outvec_in, int *len, MPI_Datatype *datatype);
void AMPI_Tupel_Min(void *invec_in, void *outvec_in, int *len, MPI_Datatype *datatype);

int AMPI_get_stack_counter();
void AMPI_set_stack_counter(int counter);
#endif
