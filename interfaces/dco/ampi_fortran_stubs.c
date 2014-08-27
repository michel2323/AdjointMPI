
/* Fortran wrappers for C AMPI routines
 * ----------------------------------------
 * The signatures of the MPI Fortran routines have no return value. Additionally, requests,
 * comms, datatypes... are of type integer. This interface maps these signatures to the
 * AMPI signatures of the AMPI C library.
 */


#include "ampi_fortran_stubs.h"


int ampi_request_chunks=1;

/* Active datatype of COMPAD. Needed to do cast from void* */


/* AMPI routines have an MPI counterpart that is adjoined */

void ampi_reset_tape_fort(int * ierr){
    double temp;
    *ierr = AMPI_Reset_Tape();
}

void ampi_init_fort(int *ierr) {
    crequest = 0;
    request_idx=malloc(REQUEST_CHUNK_SIZE*sizeof(AMPI_Request));
    char *** argv = NULL;
    int * argc=NULL; 
    *ierr =  AMPI_Init(argc, argv);
}

void ampi_comm_size_fort(int * comm, int * numprocs, int * ierr) {
    *ierr = AMPI_Comm_size(MPI_COMM_WORLD, numprocs);
}


void ampi_comm_rank_fort(int * comm, int * myid, int * ierr) {
    *ierr = AMPI_Comm_rank(MPI_COMM_WORLD, myid);
}


void ampi_finalize_fort(int * ierr){
    double temp;
    *ierr = AMPI_Finalize();
}

void ampi_send_fort(compad_type *buf, int *count, int *datatype, int *dest, int *tag, int *comm, int *ierr) {
    int i = 0;
    *ierr = AMPI_Send(buf, *count, MPI_DOUBLE, *dest, *tag, MPI_COMM_WORLD);
}


void ampi_recv_fort(compad_type *buf, int *count, int *datatype, int *dest, int *tag, int *comm, int *status, int *ierr) {
    MPI_Status tmp_status;
    *ierr = AMPI_Recv(buf, *count, MPI_DOUBLE, *dest, *tag, MPI_COMM_WORLD, &tmp_status);
}

void ampi_isend_fort(compad_type *buf, int *count, int *datatype, int *dest, int *tag, int *comm, int *request, int *ierr) {
    /*#ifdef ASSERT*/
    /*assert(crequest < REQUEST_IDX_SIZE);*/
    ampi_check_request_size();
    /*#endif*/
    request_idx[crequest].v = buf;
    request_idx[crequest].dest= dest;
    request_idx[crequest].oc = AMPI_IS;
    request_idx[crequest].size = count;
    request_idx[crequest].tag = tag;
    request_idx[crequest].comm = comm;
    *request = crequest;
    *ierr = AMPI_Isend(buf, *count, MPI_DOUBLE, *dest, *tag, MPI_COMM_WORLD, &request_idx[crequest]);
    crequest++;
}

void ampi_irecv_fort(compad_type *buf, int *count, int *datatype, int *dest, int *tag, int *comm, int *request, int *ierr) {
    /*#ifdef ASSERT*/
    /*assert(crequest < REQUEST_IDX_SIZE);*/
    ampi_check_request_size();
    /*#endif*/
    request_idx[crequest].v = buf;
    request_idx[crequest].dest= dest;
    request_idx[crequest].oc = AMPI_IR;
    request_idx[crequest].size = count;
    request_idx[crequest].tag = tag;
    request_idx[crequest].comm = comm;
    *request = crequest;
    *ierr = AMPI_Irecv(buf, *count, MPI_DOUBLE, *dest, *tag, MPI_COMM_WORLD, &request_idx[crequest]);
    crequest++;
}

void ampi_wait_fort(int *request, int *status, int * ierr) {
    *ierr = AMPI_Wait(&request_idx[*request], &request_idx[*request].status);
}

void ampi_waitall_fort(int *count, int *request, int *status, int * ierr) {
    int i=0;
    for(i=0;i<*count;i=i+1) {
    	*ierr = AMPI_Wait(&request_idx[request[i]], &request_idx[request[i]].status);
    }
}

void ampi_reduce_fort(compad_type *sendbuf, compad_type *recvbuf, int *count, int *datatype, int *op, int *root, int *comm, int *ierr) {
    int i = 0;
    MPI_Op tmp_op=-1;
    if(*op == AMPI_SUM)
	tmp_op = MPI_SUM;
    if(*op == AMPI_PROD)
	tmp_op = MPI_PROD;
    if(*op == AMPI_MIN)
	tmp_op = MPI_MIN;
    if(*op == AMPI_MAX)
	tmp_op = MPI_MAX;
    if(tmp_op == -1) {
	printf("Unsupported reduction operation in AMPI\n");
	*ierr=-1;
	return;
    }
    *ierr = AMPI_Reduce(sendbuf, recvbuf, *count, MPI_DOUBLE, tmp_op, *root, MPI_COMM_WORLD);
}

void ampi_allreduce_fort(compad_type *sendbuf, compad_type *recvbuf, int *count, int *datatype, int *op, int *comm, int *ierr) {
    int i = 0;
    MPI_Op tmp_op = -1;
    if(*op == AMPI_SUM)
	tmp_op = MPI_SUM;
    if(*op == AMPI_PROD)
	tmp_op = MPI_PROD;
    if(*op == AMPI_MIN)
	tmp_op = MPI_MIN;
    if(*op == AMPI_MAX)
	tmp_op = MPI_MAX;
    if(tmp_op == -1) {
	printf("Unsupported reduction operation in AMPI\n");
	*ierr=-1;
	return;
    }
    *ierr = AMPI_Allreduce(sendbuf, recvbuf, *count, MPI_DOUBLE, tmp_op, MPI_COMM_WORLD);
}

/* The fundamental AMPI tape routines that need to be implemented in Fortran. An
 * interface is needed to cast compad_type to void
 */

void ampi_get_val(void *buf, int *i, double *v) {
   compad_type *buf_compad = buf;
   ampi_get_val_fort(buf_compad, i, v); 
}

void ampi_set_val(void *buf, int *i, double *v) {
   compad_type *buf_compad = buf;
   ampi_set_val_fort(buf_compad, i, v); 
}

void ampi_get_adj(INT64 *idx, double *v) {
   ampi_get_adj_fort(idx, v); 
}

void ampi_set_adj(INT64 *idx, double *v) {
   ampi_set_adj_fort(idx, v); 
}

void ampi_get_idx(void *buf, int *i, INT64 *idx) {
   compad_type *buf_compad = buf;
   ampi_get_idx_fort(buf_compad, i, idx); 
}

void ampi_create_tape_entry(int *i) {
    ampi_create_tape_entry_fort(i);
}

void ampi_create_dummies(void *buf, int *size) {
    compad_type *buf_compad = buf;
    ampi_create_dummies_fort(buf_compad, size);
}

void ampi_check_request_size() {
    if(crequest>=REQUEST_CHUNK_SIZE) {
	AMPI_Request *tmp;
        tmp=realloc(request_idx,(ampi_request_chunks+1)*REQUEST_CHUNK_SIZE*sizeof(AMPI_Request));
	if(tmp != NULL) {
	    request_idx=tmp;
	    ampi_request_chunks=ampi_request_chunks+1;
	}
	else {
	    printf("AMPI Fortran request allocation error.\n");
	}
    }
}

/*void ampi_print_tape_entry_(int *i) {*/
/*ampi_print_tape_entry(i);*/
/*}*/
