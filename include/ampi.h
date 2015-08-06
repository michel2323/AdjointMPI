#ifndef AMPI_H
#define AMPI_H
/** \file
 * \brief Header file for source transformation tools.
 * This file constitutes the interface for source transformation tools. The AMPI
 * tape is dependent on these routines. All MPI routines
 * are subdivided into their forward _f and backward _b counterpart. The forward routines
 * are called during the forward/taping run. The backward routines are called during the
 * reverse/interpretation run.
 */
#include <mpi.h>
#include <ampi_interface.h>
#include <ampi_stack.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

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

struct ampi_tape_entry;

/**
 * AMPI request, replacing the MPI request 
 */
typedef struct AMPI_Request {
    MPI_Request *mpiRequest; /**< Original request */
    MPI_Status status; /**< Original status */
    MPI_Comm comm; /**< Original communicator */

    /* AMPI tape only */

    void *buf;            /**< Only used in the tape. Here we store the pointer to the active buffer that is conveyed to the Wait. */
    int tag;              /**< MPI tag */
    double *v;            /**< Incoming or outgoing values of the mapped active buffer. */
    double *a;            /**< Incoming or outgoing adjoints of the mapped tape. */
    struct ampi_tape_entry* va;  /**< Tape index */
    int oc;               /**< Operation code */
    int dest;             /**< Destination or source. */
    int size;             /**< Size of the buffer. */
} AMPI_Request;


/**
 * Tupel used in the MPI_MAX and MPI_MIN reduction to save the rank of
 * the process with the maximum or minimum
 */
typedef struct AMPI_Tupel {
    double v; /**< value */
    int j;    /**< rank */
} AMPI_Tupel;


/**
 * Active forward MPI_Init.
 * First chunk of the AMPI tape is allocated.
 *
 * @param argc argc of passive code
 * @param argv argv of passive code
 *
 * @return error code 
 */
int AMPI_Init_f(int* argc, char*** argv);

/**
 * Active reverse MPI Init.
 * AMPI data structures are destroyed and MPI_Finalize() is called.
 *
 * @param argv Dummy argc. Not used
 * @param argc Dummy argv. Not used
 *
 * @return error code
 */
int AMPI_Init_b(int* argc, char*** argv);

/**
 * Legacy active variant of MPI_Comm_size. Only a wrapper of MPI_Comm_size.
 *
 * @param comm Communicator
 * @param numprocs Number of processes
 *
 * @return error code 
 */
int AMPI_Comm_size(MPI_Comm comm, int* numprocs);

/**
 * Active variant of MPI_Comm_rank. The rank is saved in a global variable
 * to avoid repeated calls to MPI_Comm_rank.
 *
 * @param comm Communicator
 * @param rank Rank of calling process
 *
 * @return error code 
 */
int AMPI_Comm_rank(MPI_Comm comm, int* rank);

/**
 * Legacy. Wrapper to MPI_Get_processor_name 
 *
 * @param name Processor name
 * @param namelength Processor name length. 
 *
 * @return error code
 */
int AMPI_Get_processor_name(char* name, int* namelength);

/**
 * Active barrier. Amounts to a wrapper of MPI_Barrier. Does 
 * not need to be traced.
 *
 * @param comm Communicator with the processes that execute a barrier
 *
 * @return error code
 */
int AMPI_Barrier(MPI_Comm comm);

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

/**
 * Active forward send. 
 *
 * The forward send amounts to a wrapper of the MPI_Send.
 *
 * @param buf Buffer with values that are to be sent
 * @param count Number of buffer elements
 * @param datatype MPI data type of the buffer elements
 * @param dest Rank of destination process
 * @param tag Message tag
 * @param comm MPI communicator
 *
 * @return error code 
 */
int AMPI_Send_f(double *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);

/**
 * Active forward buffered send.
 *
 * The forward send amounts to a wrapper of the MPI_Send.
 *
 * @param buf Buffer with values that are to be sent
 * @param count Number of buffer elements
 * @param datatype MPI data type of the buffer elements
 * @param dest Rank of destination process
 * @param tag Message tag
 * @param comm MPI communicator
 *
 * @return error code
 */
int AMPI_Bsend_f(double *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);

/**
 * Active reverse send. 
 *
 * The active reverse send amounts to an MPI_Receive with MPI_STATUS_IGNORE.
 *
 * @param buf Buffer with adjoints that are received
 * @param count Number of buffer elements
 * @param datatype MPI data type of the buffer elements
 * @param src Rank of source process
 * @param tag Message tag
 * @param comm MPI communicator
 *
 * @return error code 
 */
int AMPI_Send_b(double *buf, int count, MPI_Datatype datatype, int src, int tag, MPI_Comm comm);

/**
 * Active forward receive.
 *
 * The active forward receive is a wrapper of MPI_Recv.
 *
 * @param buf Buffer with values that are received
 * @param count Number of buffer elements
 * @param datatype MPI data type of the buffer elements
 * @param src Rank of source process
 * @param tag Message tag
 * @param comm MPI communicator
 * @param status MPI status of the received value message
 *
 * @return error code 
 */
int AMPI_Recv_f(double *buf, int count, MPI_Datatype datatype, int src, int tag, MPI_Comm comm, MPI_Status *status);

/**
 * Active reverse receive.
 *
 * The active reverse receive is a wrapper of MPI_Send. The status is ignored.
 *
 * @param buf Buffer with adjoints that are sent
 * @param count Number of buffer elements
 * @param datatype MPI data type of the buffer elements
 * @param dest Rank of destination process
 * @param tag Message tag
 * @param comm MPI communicator
 * @param status Ignored. Only present for consistency with the MPI signatures
 *
 * @return error code
 */
int AMPI_Recv_b(double *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Status *status);


/**
 * Documentation TODO for Michel
 */
int AMPI_Brecv_b(double *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Status *status);

/* Non blocking communication */

/**
 * Active forward non blocking send
 *
 * @param buf Buffer with values that are sent
 * @param count Number of buffer elements
 * @param datatype MPI data type of the buffer elements
 * @param dest Rank of destination process
 * @param tag Message tag
 * @param comm MPI communicator
 * @param request Active MPI_Request. Additional information (buf, count, dest,
 * oc, comm and the original MPI_Request) is saved in here for later use in
 * AMPI_Wait_f or AMPI_Waitall_f.
 *
 * @return error code 
 */
int AMPI_Isend_f(double *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, AMPI_Request *request);

/**
 * Active reverse non blocking send. This amounts to an MPI_Wait. The received
 * adjoints are copied to the adjoint buffer that is conveyed together with the
 * MPI_Request through the
 * AMPI_Request.
 *
 * @param buf Points to the buffer address for the received adjoints in the request
 * @param count Dummy argument
 * @param datatype Dummy argument
 * @param src Dummy argument
 * @param tag Message tag
 * @param comm MPI Dummy argument
 * @param request Active MPI_Request. It conveys the buffer address for the
 * adjoints from the reverse wait to the reverse send (=receive).
 *
 * @return error code 
 */
int AMPI_Isend_b(double *buf, int count, MPI_Datatype datatype, int src, int tag, MPI_Comm comm, AMPI_Request *request);

/**
 * Active forward non blocking receive
 *
 * @param buf Buffer with values that are to be received
 * @param count Number of buffer elements
 * @param datatype MPI data type of the buffer elements
 * @param dest Rank of source process
 * @param tag Message tag
 * @param comm MPI communicator
 * @param request Active MPI_Request. Additional information (buf, count, dest,
 * oc, comm and the original MPI_Request) is saved in here for later use in
 * AMPI_Wait_f or AMPI_Waitall_f.
 *
 * @return error code 
 */
int AMPI_Irecv_f(double *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, AMPI_Request *request);

/**
 * Active reverse non blocking receive. This amounts to an MPI_Wait. For
 * sementical reasons, the sent adjoints are copied to the adjoint buffer buf
 * that is conveyed together with the MPI_Request through the AMPI_Request.
 *
 * @param buf Points to the buffer address of the sent adjoints in the request
 * @param count Dummy argument
 * @param datatype Dummy argument
 * @param dest Dummy argument
 * @param tag Message tag
 * @param comm MPI Dummy argument
 * @param request Active MPI_Request. It conveys the buffer address for the
 * adjoints from the reverse wait to the reverse send (=receive).
 *
 * @return error code 
 */
int AMPI_Irecv_b(double *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, AMPI_Request *request);


/**
 * Active forward wait. An MPI_Wait is executed. In the request, additional information
 * is provided (destination, source, value buffer address, communicator, count,
 * type of communication).
 *
 * @param request Active MPI request
 * @param status The original status
 *
 * @return error code 
 */
int AMPI_Wait_f(AMPI_Request *request , MPI_Status *status);

/**
 * Active reverse wait. Depending on the associated communication a non
 * blocking send (for a receive) or a receive (for a send) is executed with the
 * adjoint buffer address. The information for the reverse communication is
 * assumed to be stored in the active request. 
 *
 * @param request Active MPI request with all the information needed for the
 * reverse communication
 * @param status Dummy argument. Status is saved in the active requests
 *
 * @return error code 
 */
int AMPI_Wait_b(AMPI_Request *request, MPI_Status *status);

/**
 * Active forward waitall. An MPI_Waitall is executed. The original
 * MPI_Request requests are saved into the active requests. 
 *
 * @param count Number of requests
 * @param requests Active requests
 * @param status Original statuses
 *
 * @return error code 
 */
int AMPI_Waitall_f(int count, AMPI_Request *requests, MPI_Status *status);

/**
 * Forward waitany. 
 * 
 * @param count Number of requests
 * @param requests MPI requests
 * @param status Original statuses
 *
 * @return error code 
 */
int AMPI_Waitany_f(int count, MPI_Request array_of_request[], int *index, MPI_Status *status);

/**
 * Active reverse waitall. This marks a performance loss in the active
 * case. For each forward wait, a send or receive communication has to be
 * started. Therefore we call an active reverse wait for each active request. 
 *
 * @param count Number of requests
 * @param requests Active requests
 * @param status Dummy argument. Status is saved in the active requests
 *
 * @return error code 
 */
int AMPI_Waitall_b(int count, AMPI_Request *requests, MPI_Status *status);

/* Anti Waitall. Experimental */

/**
 * Experimental implementation of the anti-wait. See corresponding paper
 * for additional information. In the forward routine, the anti-wait flag in the
 * active requests are set.
 *
 * @param count Number of requests
 * @param requests Active requests
 * @param status Original statuses. Dummy argument
 *
 * @return error code 
 */
int AMPI_Awaitall_f(int count, AMPI_Request *requests, MPI_Status *status);

/**
 * Experimental implementation of the anti-wait. See corresponding paper
 * for additional information. In the reverse routine, all the wait operations
 * for the adjoint communications are executed here, as opposed to single waits
 * in the non blocking sends and receives.
 *
 * @param count Number of requests
 * @param requests Active requests
 * @param status Dummy argument. Status is inside active requests
 *
 * @return error code 
 */
int AMPI_Awaitall_b(int count, AMPI_Request *requests, MPI_Status *status);

/* Collective */

/**
 * Forward active broadcast. An MPI_Broadcast is executed.
 *
 * @param buf Buffer of broadcast value
 * @param count Number of buffer elements
 * @param datatype MPI data type of the buffer elements
 * @param root Process with the buffer that is broadcast
 * @param comm MPI communicator
 *
 * @return error code 
 */
int AMPI_Bcast_f(double *buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm);

/**
 * Reverse active broadcast. An MPI_Reduce with the MPI_SUM operation is
 * executed on the adjoints.
 *
 * @param buf Send buffer with adjoints that are to be reduced. After
 * AMPI_Bcast_f is called this is the received buffer of reduced adjoint on the
 * root process.
 * @param count Number of buffer elements
 * @param datatype MPI data type of the buffer elements
 * @param root Process with the buffer that the adjoints are reduced to
 * @param comm MPI communicator
 *
 * @return error code
 */
int AMPI_Bcast_b(double *buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm);

/**
 * Forward active reduce. See corresponding papers for additional
 * information. Depending on the operation op, an MPI communication is executed
 * and the information necessary for the adjoint operation is saved.
 *
 * MPI_SUM: An MPI_Reduce is executed.
 * MPI_PROD: An MPI_Allreduce is executed, because all the processes need to
 * result of the multiplication.
 * MPI_MAX and MPI_MIN: The rank of the process which holds the maximum oder
 * minimum is saved. Hence, we execute an MPI_MAXLOC or MPI_MINLOC after copying
 * the data accordingly.
 *
 * @param sendbuf The sent value buffer
 * @param recvbuf The received value buffer. In case of the MPI_PROD all the
 * processes have the result.
 * @param count Number of buffer elements
 * @param datatype MPI data type of the buffer elements
 * @param op Reduction operation (MPI_SUM, MPI_PROD, MPI_MAX or MPI_MIN)
 * @param root Process with the buffer that the values are reduced to
 * Irreleveant in case of MPI_PROD.
 * @param comm MPI communicator
 *
 * @return error code 
 */
int AMPI_Reduce_f(double *sendbuf, double *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm, ampi_stack** stack);

/**
 * Reverse active reduction. The adjoint of a reduction is a broadcast
 * with an additional operation on the incoming distributed adjoint. There is no
 * MPI routine that implements this. So we have to decompose the communication
 * into the broadcast and a local operation on each process. In case of
 * MPI_PROD, we have currently implemented the short solution, where the local
 * adjoint is computed by dividing the global result with the local result.
 *
 * @param sendbuf The incoming adjoint buffer address
 * @param recvbuf The broadcast adjoint buffer address
 * @param count Number of buffer elements
 * @param datatype MPI data type of the buffer elements
 * @param op Reduction operation (MPI_SUM, MPI_PROD, MPI_MAX or MPI_MIN)
 * @param root Process with the buffer that the adjoints are broadcast from
 * @param comm MPI communicator
 *
 * @return error code 
 */
int AMPI_Reduce_b(double *sendbuf, double *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm, ampi_stack* stack);

/**
 * Forward active allreduce. See corresponding papers for additional
 * information. The information necessary for the adjoint operation is saved.
 *
 * MPI_SUM and MPI_PROD: The result and the input send buffer are saved.
 * MPI_MAX and MPI_MIN: The rank of the process which holds the maximum oder
 * minimum is saved. Hence, we execute an MPI_MAXLOC or MPI_MINLOC after copying
 * the data accordingly.
 *
 * @param sendbuf The sent value buffer
 * @param recvbuf The received value buffer. In case of the MPI_PROD all the
 * processes have the result.
 * @param count Number of buffer elements
 * @param datatype MPI data type of the buffer elements
 * @param op Reduction operation (MPI_SUM, MPI_PROD, MPI_MAX or MPI_MIN)
 * @param comm MPI communicator
 *
 * @return error code 
 */
int AMPI_Allreduce_f(double *sendbuf, double *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, ampi_stack** stack);

/**
 * Reverse active allreduce. All the adjoints of all the processes have
 * to be summed up and sent to all processes, essentially amounting to an allreduce with an MPI_SUM operation. We
 * then apply the adjoint computation like in the common reduction. 
 *
 * @param sendbuf The incoming adjoint buffer address
 * @param recvbuf The broadcast adjoint buffer address
 * @param count Number of buffer elements
 * @param datatype MPI data type of the buffer elements
 * @param op Reduction operation (MPI_SUM, MPI_PROD, MPI_MAX or MPI_MIN)
 * @param comm MPI communicator
 *
 * @return error code 
 */
int AMPI_Allreduce_b(double *sendbuf, double *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, ampi_stack* stack);

/**
 * Active forward gather. MPI_Gather wrapper. 
 *
 */
int AMPI_Gather_f(void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm);

/**
 * Active reverse gather. MPI_Scatter wrapper. sendbuf are the received
 * adjoints while recvbuf are the sent adjoints. 
 *
 */
int AMPI_Gather_b(void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm);

/**
 * Active forward scatter. MPI_Scatter wrapper. 
 *
 */
int AMPI_Scatter_f(void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm);

/** Active reverse scatter. MPI_Gather wrapper. sendbuf are the received
 * adjoints while recvbuf are the sent adjoints. 
 *
 */
int AMPI_Scatter_b(void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm);

/* Special comms */

/** Active forward combined send and receive where the same buffer is used.
 * MPI_Sendrecv_replace wrapper. 
 *
 */
int AMPI_Sendrecv_replace_f(double *buf, int count, MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag, MPI_Comm comm, MPI_Status *status);


/**
 * Active reverse combined send and receive where the same buffer is used.
 * MPI_Sendrecv_replace wrapper. The reversal amounts to switching destination
 * and source 
 *
 */
int AMPI_Sendrecv_replace_b(double *buf, int count, MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag, MPI_Comm comm, MPI_Status *status);

int AMPI_Sendrecv_f(void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf, int recvcount, MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status *status);

int AMPI_Sendrecv_b(void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf, int recvcount, MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status *status);
#endif
