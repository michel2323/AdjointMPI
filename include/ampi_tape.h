#ifndef AMPI_TAPE_INCLUDE
#define AMPI_TAPE_INCLUDE

/* Generic AMPI C tape. Don't touch this. The code is always mirrored with the AMPI repo.
 * Changing code here will result in merge conflicts. 
 */

/* Tape is static at the moment. Should become dynamic based on a chunk scheme */
#define AMPI_CHUNK_SIZE 500000
#define ASSERT 

#define REDUCE_ADD 1
#define REDUCE_MUL 2
#define REDUCE_MIN 3
#define REDUCE_MAX 4

#define SEND 1
#define RECV 2
#define ISEND 3
#define IRECV 4
#define WAIT 5
#define WAITALL 6
#define AWAITALL 7
#define BCAST 8
#define REDUCE 9
#define ALLREDUCE 10
#define MPI_DUMMY 11
#define MPI_ADUMMY 12
#define SENDRECVREPLACE 13
#define SCATTER 14
#define GATHER 15
#define SCATTERV 16
#define GATHERV 17
#define SEND_INIT 18
#define RECV_INIT 19
#define START 20
#define STARTALL 21

#include <stdlib.h>
#include <assert.h>
#include <ampi.h>
#include <uthash.h>

/*int ampi_vac=0;*/

typedef struct {
    void *key;
    AMPI_Request request;
    UT_hash_handle hh; 
} AMPI_ht_el;


typedef struct ampi_tape_entry {
    int oc;
    int *arg;
    /*int arg1;*/
    /*int arg2;*/
    /*double v;*/
    /*double d;*/
    INT64 idx;
    AMPI_Request *request;
    MPI_Comm comm;
    int tag;
}ampi_tape_entry;

/* AMPI taping routines which have an MPI counterpart that is adjoined.*/

int AMPI_Reset_Tape();
int AMPI_Init(int*, char***);
int AMPI_Finalize();

int AMPI_Send(void*, int, MPI_Datatype, int, int, MPI_Comm);
int AMPI_Recv(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status*);

int AMPI_Isend(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
int AMPI_Irecv(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);

int AMPI_Wait(MPI_Request *, MPI_Status *);
int AMPI_Waitall(int , MPI_Request *, MPI_Status *);
int AMPI_Waitany(int count, MPI_Request array_of_requests[], int *index, MPI_Status *status);
int AMPI_Awaitall(int , AMPI_Request *, MPI_Status *);

int AMPI_Bcast(void *buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm);
int AMPI_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
int AMPI_Allreduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int AMPI_Scatter(void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm);
int AMPI_Gather(void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm);

int AMPI_Send_init(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
int AMPI_Recv_init(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request);
int AMPI_Start(MPI_Request *request);
int AMPI_Startall(int count, MPI_Request array_of_requests[]);

int AMPI_Sendrecv_replace(void *buf, int count, MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag, MPI_Comm comm, MPI_Status *status);
/* AMPI routines that have to be called by the external tape */

/* Call AMPI tape interpreter from external tape */
void ampi_interpret_tape(void);

/* Call AMPI tape printer from external tape printer */
void ampi_print_tape(void);
/* Call AMPI tape printer from external tape printer for one tape entry */
void ampi_print_tape_entry(int *j);
void ampi_check_tape_size(int size);

/* AMPI routines which are defined as external. These routines need to be implemented by
 * the external tape library. They should implement the data flow between the external
 * tape and the AMPI tape.
 */

/* Get a value v from a specific tape entry with index idx */
extern void ampi_get_val(void* buf, int* i, double* v);

/* Get a value v from a specific tape variable buf[i] */
extern void ampi_set_val(void* buf, int* i, double* v);

/* Get an adjoint a from a specific tape entry with index idx */
extern void ampi_get_adj(INT64* idx, double* a);

/* Set an adjoint a from a specific tape entry with index idx */
extern void ampi_set_adj(INT64*, double*);

/* Get a tape index from a variable buf[i] */
extern void ampi_get_idx(void* buf, int* i, INT64* idx);

/* Create a tape entry in the external tape, indicating an external AMPI call */
extern void ampi_create_tape_entry(int* i);

/* Create size tape entries to store the values of buf. Refer to receive buffer without
 * initialization */
extern void ampi_create_dummies(void* buf, int *size);
#endif
