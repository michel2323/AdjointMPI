#ifndef AMPI_H
#define AMPI_H
/** \file
 * \brief Header file for source transformation tools.
 *
 * This file constitutes the interface for source transformation tools. The AMPI
 * tape is dependent on these routines. All MPI routines
 * are subdivided into their forward _f and backward _b counterpart. The forward routines
 * are called during the forward/taping run. The backward routines are called during the
 * reverse/interpretation run.
 */
#include <mpi.h>
#include <ampi_stack.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
/**
 * \def INT64
 * Defines the type of the tape index
 */
#define INT64 long int

/**
 * @{
 * \name AMPI_Request Send/Receive
 * Defines tor distinguishing send and receive in an AMPI_Request
 */
#define AMPI_IS 1
#define AMPI_IR 2
/**@}*/


/**
 * @{
 * \name Legacy defines with active predefined constants
 */
#define AMPI_COMM_WORLD MPI_COMM_WORLD
#define AMPI_MAX_PROCESSOR_NAME MPI_MAX_PROCESSOR_NAME
#define AMPI_Status MPI_Status
#define AMPI_DOUBLE MPI_DOUBLE
/**@}*/


/**
 * @{
 * \name Active Reduce Operations
 */
#define AMPI_SUM   1
#define AMPI_PROD  2
#define AMPI_MIN   3
#define AMPI_MAX   4
/**@}*/

/**
 * AMPI request, replacing the MPI request 
 */
typedef struct AMPI_Request {
    MPI_Request *request; /** Original request */
    MPI_Status status; /** Original status */
    MPI_Comm comm; /** Original communicator */

    /*special for tape*/
    void *buf;            /** Only used in the tape. Here we store the pointer to the active buffer that is conveyed to the Wait. */

    int tag;
    double *v;            /** Incoming or outgoing values of the mapped active buffer. */
    double *a;            /** Incoming or outgoing adjoints of the mapped tape. */
    long int va;          /** Tape index */
    int oc;
    int dest;             /** Destination or source. */
    int size;             /** Size of the buffer. */
    long int aw;          /** Anti wait flag*/
} AMPI_Request;


/**
 * Tupel used in the MPI_MAX and MPI_MIN reduction to save the rank of
 * the process with the maximum or minimum
 */
typedef struct AMPI_Tupel {
    double v; /** value */
    int j;    /** rank */
} AMPI_Tupel;


/**
 * Active forward MPI_Init.
 * First chunk of the AMPI tape is allocated.
 *
 * @param int* argc
 * @param char*** argv
 *
 * @return error code 
 */
int AMPI_Init_f(int*, char***);

/**
 * Active reverse MPI Init.
 * AMPI data structures are destroyed and MPI_Finalize() is called.
 *
 * @param int* Dummy argc. Not used.
 * @param char*** Dummy argv. Not used.
 *
 * @return error code
 */
int AMPI_Init_b(int*, char***);

/**
 * Legacy active variant of MPI_Comm_size. Only a wrapper of MPI_Comm_size.
 *
 * @param MPI_Comm Communicator.
 * @param int* Number of processes.
 *
 * @return error code 
 */
int AMPI_Comm_size(MPI_Comm, int*);

/**
 * Active variant of MPI_Comm_rank. The rank is saved in a global variable
 * to avoid repeated calls to MPI_Comm_rank.
 *
 * @param MPI_Comm Communicator.
 * @param int* Rank.
 *
 * @return error code 
 */
int AMPI_Comm_rank(MPI_Comm, int*);

/**
 * Legacy. Wrapper to MPI_Get_processor_name 
 *
 * @param char* Processor name.
 * @param int* Processor name length. 
 *
 * @return error code
 */
int AMPI_Get_processor_name(char*, int*);

/**
 * Active barrier. Amounts to a wrapper of MPI_Barrier. Does 
 * not need to be traced.
 *
 * @param MPI_Comm Communicator.
 *
 * @return error code
 */
int AMPI_Barrier(MPI_Comm);

/**
 * Active forward MPI_Finalize(). Nothing is done here. 
 *
 * @return error code 
 */
int AMPI_Finalize_f();
/**
 * Active reverse MPI_Finalize(). Nothing is done here. 
 *
 * @return error code 
 */
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
