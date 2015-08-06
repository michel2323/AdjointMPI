/* Generic AMPI C tape. Don't touch this. The code is always mirrored with the AMPI repo.
 * Changing code here will result in merge conflicts.
 *
 * See header for more information. 
 */
/*#define DEBUG*/
#define NO_COMM_WORLD

#include <stddef.h>
#include <ampi_tape.h>
#ifdef __cplusplus
#include <ampi_interface.hpp> 
#endif
#ifdef AMPI_COUNT_COMMS
int ampi_comm_count=0;
#endif

AMPI_ht_el *AMPI_ht=NULL;
ampi_tape_entry *ampi_tape;

int AMPI_Init(int* argc, char*** argv) {
    return AMPI_Init_f(argc, argv);
}

int AMPI_Finalize() {
#ifdef AMPI_COUNT_COMMS
    printf("AMPI comunications executed: %d\n", ampi_comm_count);
#endif
    return AMPI_Init_b(NULL, NULL);

}

int AMPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
#ifdef AMPI_COUNT_COMMS
    ampi_comm_count=ampi_comm_count+1;
#endif
    if(datatype!=AMPI_DOUBLE) {
      return MPI_Send(buf, count, datatype, dest, tag, comm);
    }
    int i=0;
    double *primalValues = (double*) malloc(sizeof(double)*count);

    for(i=0;i<count;i=i+1) {
        ampi_get_val(buf,&i,&primalValues[i]);
    }

    if (ampi_is_tape_active()){

      ampi_tape_entry* ampi_tape = ampi_create_tape(count+1);

      ampi_tape->arg=(int*) malloc(sizeof(int)*2);

      ampi_create_tape_entry((void*)ampi_tape);

      for(i=0;i<count;i=i+1) {
        ampi_get_idx(buf, &i, &ampi_tape->idx[i]);
      }

      /*create actual MPI entry*/

      ampi_tape->oc = SEND;
      ampi_tape->arg[0] = count;
      ampi_tape->arg[1] = dest;
      ampi_tape->comm = comm;
      ampi_tape->tag = tag;
    }

    int exitCode = AMPI_Send_f(primalValues, count, datatype, dest, tag, comm);

    free(primalValues);

    return exitCode;
}

int AMPI_Bsend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
#ifdef AMPI_COUNT_COMMS
    ampi_comm_count=ampi_comm_count+1;
#endif
    if(datatype!=AMPI_DOUBLE) {
      return MPI_Bsend(buf, count, datatype, dest, tag, comm);
    }
    int i=0;
    /*Allocate one more element and hint that it is a Bsend*/
    double *primalValues = (double*) malloc(sizeof(double)*(count+1));
    primalValues[count]=1;

    for(i=0;i<count;i=i+1) {
        ampi_get_val(buf,&i,&primalValues[i]);
    }

    if (ampi_is_tape_active()){

      ampi_tape_entry* ampi_tape = ampi_create_tape(count+1);

      ampi_tape->arg=(int*) malloc(sizeof(int)*2);

      ampi_create_tape_entry((void*)ampi_tape);

      for(i=0;i<count;i=i+1) {
        ampi_get_idx(buf, &i, &ampi_tape->idx[i]);
      }

      /*create actual MPI entry*/

      /*BSEND is adjoined like a SEND. The corresponding RECV is adjoined as a BSEND!*/
      ampi_tape->oc = SEND;
      ampi_tape->arg[0] = count;
      ampi_tape->arg[1] = dest;
      ampi_tape->comm = comm;
      ampi_tape->tag = tag;
    }

    int exitCode = AMPI_Bsend_f(primalValues, count+1, datatype, dest, tag, comm);

    free(primalValues);

    return exitCode;
}

int AMPI_Recv(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Status *status) {
#ifdef AMPI_COUNT_COMMS
    ampi_comm_count=ampi_comm_count+1;
#endif
    if(datatype!=AMPI_DOUBLE) {
      return MPI_Recv(buf, count, datatype, dest, tag, comm, status);
    }

    int i=0;
    /*One element longer to distinguish between Bsend and Send of the received message*/
    double *primalValues = (double*) malloc(sizeof(double)*(count+1));
    /* set to 0. if this is a 1, the received message was a Bsend */
    primalValues[count]=0;

    int exitCode = AMPI_Recv_f(primalValues, count+1, datatype, dest, tag, comm, status);
    /*if 1, it is a BSEND on the other end*/
    if (ampi_is_tape_active()){
      ampi_tape_entry* ampi_tape = ampi_create_tape(count+1);
      ampi_tape->arg=(int*) malloc(sizeof(int)*2);

      ampi_create_dummies(buf, &count);

      ampi_create_tape_entry((void*)ampi_tape);

      ampi_tape->arg[0] = count;
      ampi_tape->arg[1] = dest;
      ampi_tape->comm = comm;
      if(status==MPI_STATUS_IGNORE) {
        ampi_tape->tag = tag;
      } else {
        ampi_tape->tag = status->MPI_TAG;
      }

      for(i=0;i<count;i=i+1) {
        ampi_get_idx(buf, &i, &ampi_tape->idx[i]);
      }
      if(primalValues[count]==1) {
        ampi_tape->oc=BRECV;
      }
      else {
        ampi_tape->oc=RECV;
      }
    }

    for(i=0;i<count;i=i+1) {
        ampi_set_val(buf, &i, &primalValues[i]);
    }
    free(primalValues);
    return exitCode;
}

int AMPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *mpi_request) {
#ifdef AMPI_COUNT_COMMS
    ampi_comm_count=ampi_comm_count+1;
#endif
    if(datatype!=AMPI_DOUBLE) {
      return MPI_Isend(buf, count, datatype, dest, tag, comm, mpi_request);
    }
    int i=0, temp;
    double *primalValues = (double*) malloc(sizeof(double)*count);

    for(i=0 ; i<count ; i=i+1) {
      ampi_get_val(buf,&i,&primalValues[i]);
    }

    /*new AMPI primalRequest*/
    AMPI_Request primalRequest;
    /*new hash element*/
    AMPI_ht_el *ht_el=(AMPI_ht_el*) malloc(sizeof(AMPI_ht_el));
    ht_el->key=mpi_request;

    primalRequest.buf = buf;
    primalRequest.mpiRequest =mpi_request;

    if (ampi_is_tape_active()){

      ampi_tape_entry* ampi_tape = ampi_create_tape(count+1);
      ampi_tape->arg=(int*) malloc(sizeof(int)*2);

      /*create dummy of each element*/

      for(i=0 ; i<count ; i=i+1) {
        ampi_get_idx(buf, &i, &ampi_tape->idx[i]);
      }
      ampi_tape->oc = ISEND;
      ampi_tape->arg[0] = count;
      ampi_tape->arg[1] = dest;
      ampi_tape->comm = comm;
      ampi_tape->tag = tag;

      ampi_create_tape_entry((void*)ampi_tape);
      ampi_tape->request = (AMPI_Request*) malloc(sizeof(AMPI_Request));
      primalRequest.va = ampi_tape;

    }

    /*point current primalRequest index to this tape entry*/
    primalRequest.v = primalValues;
    temp = AMPI_Isend_f(primalValues, count, datatype, dest, tag, comm, &primalRequest);
    ht_el->request= primalRequest;
    HASH_ADD_PTR(AMPI_ht,key,ht_el);

    return temp;
}

int AMPI_Irecv(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *mpi_request) {
#ifdef AMPI_COUNT_COMMS
    ampi_comm_count=ampi_comm_count+1;
#endif
    if(datatype!=AMPI_DOUBLE) {
      return MPI_Irecv(buf, count, datatype, dest, tag, comm, mpi_request);
    }
    int i=0;    
    /* One more element. Could be a Bsend on the other side. This is only
     * important to avoid deadlocks in the blocking case. Irrelevant in this
     * case, however the element needs to be there */
    double * tmp = (double*) malloc(sizeof(double)*(count+1));
    tmp[count]=0;

    AMPI_Request primalRequest;
    /*INT64 tmp_int64 = 0;*/
    AMPI_ht_el *ht_el=(AMPI_ht_el*) malloc(sizeof(AMPI_ht_el));
    ht_el->key=mpi_request;
    primalRequest.buf = buf;
    primalRequest.mpiRequest =mpi_request;

    if (ampi_is_tape_active()){
      ampi_tape_entry* ampi_tape = ampi_create_tape(count+1);
      ampi_tape->arg=(int*) malloc(sizeof(int)*2);

      ampi_create_tape_entry((void*)ampi_tape);

      ampi_create_dummies(buf, &count);

      ampi_tape->arg[0] = count;
      ampi_tape->arg[1] = dest;
      ampi_tape->oc = IRECV;
      ampi_tape->comm = comm;
      ampi_tape->tag = tag;

      for(i=0 ; i<count ; i=i+1) {
        ampi_get_idx(buf, &i, &ampi_tape->idx[i]);
      }
      ampi_tape->request = (AMPI_Request*) malloc(sizeof(AMPI_Request));
      primalRequest.va = ampi_tape;
    }

    int temp = AMPI_Irecv_f(tmp, count+1, datatype, dest, tag, comm, &primalRequest);

    /*point current primalRequest index to this tape entry*/
    ht_el->request= primalRequest;
    HASH_ADD_PTR(AMPI_ht,key,ht_el);

    return temp;
}

int AMPI_Wait(MPI_Request *mpi_request, MPI_Status *status) {
#ifdef AMPI_COUNT_COMMS
    ampi_comm_count=ampi_comm_count+1;
#endif
  int i = 0;
  int ret = 0;
  AMPI_ht_el *ht_req = NULL;
  void *addr = mpi_request;
  HASH_FIND_PTR(AMPI_ht, &addr, ht_req);
  if (ht_req) {
    AMPI_Request primalRequest;
    /*AMPI_Request *primalRequest=(AMPI_Request*) malloc(sizeof(AMPI_Request));*/
    primalRequest = ht_req->request;
    double *primalValues = (double *) primalRequest.v;
    /*get the corresponding isend or irecv tape entry*/

    ret = AMPI_Wait_f(&primalRequest, status);
    /*finally copy the primalRequest to the tape*/
    if (status == MPI_STATUS_IGNORE) {
      primalRequest.tag = primalRequest.tag;
    } else {
      primalRequest.tag = status->MPI_TAG;
    }

    if (primalRequest.oc == AMPI_IR) {
      for (i = 0; i < primalRequest.size; i = i + 1) {
        ampi_set_val(primalRequest.buf, &i, &primalValues[i]);
      }
    }

    if (ampi_is_tape_active()){
      ampi_tape_entry *ampi_tape = ampi_create_tape(1);
      ampi_tape->oc = WAIT;
      ampi_tape->tag=primalRequest.tag;
      ampi_tape->request = primalRequest.va->request;
      *primalRequest.va->request = primalRequest;
      /*ampi_tape->primalRequest->a = &ampi_tape[primalRequest->va].d;*/
      ampi_tape->request->size = primalRequest.va->arg[0];
      ampi_tape->request->va = primalRequest.va; /* link here the information backward */
      ampi_create_tape_entry((void *) ampi_tape);
    }
    HASH_DEL(AMPI_ht, ht_req);
    free(ht_req);
    free(primalValues);

    return ret;
  } else {
    return MPI_Wait(mpi_request, status);
  }
}

int AMPI_Waitall(int count, MPI_Request *mpi_request, MPI_Status *status) {
#ifdef AMPI_COUNT_COMMS
    ampi_comm_count=ampi_comm_count+1;
#endif
    int i=0;
    int exitCode=-1;
    if (status != MPI_STATUSES_IGNORE) {
        for(i=0; i<count; i=i+1) {
            exitCode=AMPI_Wait(&mpi_request[i],&status[i]);
        }
    } else {
        for(i=0; i<count; i=i+1) {
            exitCode=AMPI_Wait(&mpi_request[i],status);
        }
    }
    return exitCode;
}
int AMPI_Waitany(int count, MPI_Request array_of_requests[], int *index, MPI_Status *status) {
#ifdef AMPI_COUNT_COMMS
    ampi_comm_count=ampi_comm_count+1;
#endif
    int exitCode=AMPI_Waitany_f(count,array_of_requests,index,status);
    if(*index < count && *index >=0) {
      int i=0;
      AMPI_ht_el *ht_req = NULL;
      void *addr = &array_of_requests[*index];
      HASH_FIND_PTR(AMPI_ht, &addr, ht_req);
      if (ht_req) {
        AMPI_Request primalRequest;
        /*AMPI_Request *primalRequest=(AMPI_Request*) malloc(sizeof(AMPI_Request));*/
        primalRequest = ht_req->request;
        double *primalValues = (double *) primalRequest.v;
        /*get the corresponding isend or irecv tape entry*/
        /*finally copy the primalRequest to the tape*/
        if (status == MPI_STATUS_IGNORE) {
          primalRequest.tag = primalRequest.tag;
        } else {
          primalRequest.tag = status->MPI_TAG;
        }

        if (primalRequest.oc == AMPI_IR) {
          for (i = 0; i < primalRequest.size; i = i + 1) {
            ampi_set_val(primalRequest.buf, &i, &primalValues[i]);
          }
        }

        if (ampi_is_tape_active()){
          ampi_tape_entry *ampi_tape = ampi_create_tape(1);
          ampi_tape->oc = WAIT;
          ampi_tape->tag=primalRequest.tag;
          ampi_tape->request = primalRequest.va->request;
          *primalRequest.va->request = primalRequest;
          /*ampi_tape->primalRequest->a = &ampi_tape[primalRequest->va].d;*/
          ampi_tape->request->size = primalRequest.va->arg[0];
          ampi_tape->request->va = primalRequest.va; /* link here the information backward */
          ampi_create_tape_entry((void *) ampi_tape);
        }
        HASH_DEL(AMPI_ht, ht_req);
        free(ht_req);
        free(primalValues);
      }
    }
    return exitCode; 
}
int AMPI_Bcast(void *buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm) {
#ifdef AMPI_COUNT_COMMS
    ampi_comm_count=ampi_comm_count+1;
#endif
    if(datatype!=AMPI_DOUBLE) {
      return MPI_Bcast(buf, count, datatype, root, comm);
    }
    int rank=0;
    int i=0;
    double *primalValues = (double*) malloc(sizeof(double)*count);
    for(i=0;i<count;i++) primalValues[i]=0;
    MPI_Comm_rank(comm,&rank);

    if(rank==root) {
      for(i = 0 ; i < count ; i=i+1) {
        ampi_get_val(buf,&i,&primalValues[i]);
      }
    }

    if (ampi_is_tape_active()){

      ampi_tape_entry* ampi_tape = ampi_create_tape(count+1);

      if(rank==root) {
        for(i = 0 ; i < count ; i=i+1) {
          ampi_get_idx(buf, &i, &ampi_tape->idx[i]);
        }
      }

      /*create actual MPI entry*/
      if(rank!=root) {
        ampi_create_dummies(buf, &count);
      }

      ampi_tape->arg=(int*) malloc(sizeof(int)*2);
      ampi_tape->oc = BCAST;
      ampi_tape->arg[0] = count;
      ampi_tape->arg[1] = root;
      ampi_tape->comm = comm;

      if(rank!=root) {
        for(i=0;i<count;i=i+1) {
          ampi_get_idx(buf, &i, &ampi_tape->idx[i]);
        }
      }
      ampi_create_tape_entry((void*)ampi_tape);
    }

    int temp = AMPI_Bcast_f(primalValues, count, datatype, root, comm);

    if(rank!=root) {
      for(i=0;i<count;i=i+1) {
          ampi_set_val(buf, &i, &primalValues[i]);
        }
      }

    free(primalValues);

    return temp;
}

int AMPI_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {
#ifdef AMPI_COUNT_COMMS
    ampi_comm_count=ampi_comm_count+1;
#endif
    if(datatype!=AMPI_DOUBLE) {
      return MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
    }
    int i=0;
    double * tmp_send = (double*) malloc(sizeof(double)*count);
    double * tmp_recv = (double*) malloc(sizeof(double)*count);

    for(i=0 ; i<count ; i=i+1) {
        ampi_get_val(sendbuf,&i,&tmp_send[i]);
        ampi_get_val(recvbuf,&i,&tmp_recv[i]);
    }


    if (ampi_is_tape_active()){
      ampi_tape_entry* ampi_tape = ampi_create_tape(2*count+1);


      ampi_tape->arg=(int*) malloc(sizeof(int)*3);

      ampi_create_dummies(recvbuf, &count);

      ampi_create_tape_entry((void*)ampi_tape);

      /*sendbuf dummies*/

      for(i=0 ; i<count ; i=i+1) {
          ampi_get_idx(sendbuf, &i, &ampi_tape->idx[i]);
      }

      /*actual reduce entry*/

      ampi_tape->oc = REDUCE;
      ampi_tape->arg[0] = count;
      ampi_tape->arg[1] = root;
      if(op == MPI_SUM){
          ampi_tape->arg[2] = AMPI_REDUCE_ADD;
      }
      if(op == MPI_PROD){
          ampi_tape->arg[2] = AMPI_REDUCE_MUL;
      }
      if(op == MPI_MIN){
          ampi_tape->arg[2] = AMPI_REDUCE_MIN;
      }
      if(op == MPI_MAX){
          ampi_tape->arg[2] = AMPI_REDUCE_MAX;
      }
      ampi_tape->comm = comm;

      AMPI_Reduce_f(tmp_send, tmp_recv, count, datatype, op, root, comm, &ampi_tape->stack);

      /*recvbuf entry*/

      for(i=0 ; i<count ; i=i+1) {
          ampi_get_idx(recvbuf, &i, &ampi_tape->idx[count + i]);
      }

    }else{
       MPI_Reduce(tmp_send, tmp_recv, count, datatype, op, root, comm);
    }


    for(i=0 ; i<count ; i=i+1) {
        ampi_set_val(recvbuf, &i, &tmp_recv[i]);
    }

    free(tmp_send);
    free(tmp_recv);
    return 0;
}

int AMPI_Allreduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
  int ierr=0;
#ifdef AMPI_COUNT_COMMS
  ampi_comm_count=ampi_comm_count+1;
#endif
  if(datatype!=AMPI_DOUBLE) {
    return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
  }

  int i=0;
  double * tmp_send = (double*) malloc(sizeof(double)*count);
  double * tmp_recv = (double*) malloc(sizeof(double)*count);

  for(i=0 ; i<count ; i=i+1) {
    ampi_get_val(sendbuf,&i,&tmp_send[i]);
  }

  if(ampi_is_tape_active()) {

    ampi_tape_entry* ampi_tape = ampi_create_tape(2*count+1);

    ampi_tape->arg=(int*) malloc(sizeof(int)*3);

    ampi_create_dummies(recvbuf, &count);
    ampi_create_tape_entry((void*)ampi_tape);;
    /*actual reduce entry*/

    /*sendbuf dummies*/

     for(i=0 ; i<count ; i=i+1) {
       ampi_get_idx(sendbuf, &i, &ampi_tape->idx[i]);
     }

    ampi_tape->oc = ALLREDUCE;
    ampi_tape->arg[0] = count;
    ampi_tape->comm = comm;
    /*ampi_tape->arg[1] = root;*/
    if(op == MPI_SUM){
      ampi_tape->arg[2] = AMPI_REDUCE_ADD;
    }
    if(op == MPI_PROD){
      ampi_tape->arg[2] = AMPI_REDUCE_MUL;
    }
    if(op == MPI_MIN){
      ampi_tape->arg[2] = AMPI_REDUCE_MIN;
    }
    if(op == MPI_MAX){
      ampi_tape->arg[2] = AMPI_REDUCE_MAX;
    }

    ierr=AMPI_Allreduce_f(tmp_send, tmp_recv, count, datatype, op, comm, &ampi_tape->stack);

    for(i=0 ; i<count ; i=i+1) {
      ampi_get_idx(recvbuf, &i, &ampi_tape->idx[count + i]);
    }

  }else{
    ierr = MPI_Allreduce(tmp_send, tmp_recv, count, datatype, op, comm);
  }


  for(i=0 ; i<count ; i=i+1) {
    ampi_set_val(recvbuf, &i, &tmp_recv[i]);
  }


  free(tmp_send);
  free(tmp_recv);
  return ierr;
}

int AMPI_Scatter(void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm) {
#ifdef AMPI_COUNT_COMMS
    ampi_comm_count=ampi_comm_count+1;
#endif
  if(sendtype !=AMPI_DOUBLE || recvtype != AMPI_DOUBLE) {
    return MPI_Scatter(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm);
  }
  int i=0;
  int size=0;
  MPI_Comm_size(comm,&size);
  double * tmp_send = (double*) malloc(sizeof(double)*sendcnt*size);
  double * tmp_recv = (double*) malloc(sizeof(double)*recvcnt);

  for(i=0 ; i<sendcnt*size ; i=i+1) {
    ampi_get_val(sendbuf,&i,&tmp_send[i]);
  }

  if (ampi_is_tape_active()){

    ampi_tape_entry* ampi_tape = ampi_create_tape(recvcnt+sendcnt*size+1);
    ampi_tape->arg=(int*) malloc(sizeof(int)*3);
    ampi_tape->comm = comm;
    ampi_create_dummies(recvbuf, &recvcnt);

    /*sendbuf dummies*/

    for(i=0 ; i<sendcnt*size ; i=i+1) {
      ampi_get_idx(sendbuf, &i, &ampi_tape->idx[i]);
    }

    ampi_tape->oc = SCATTER;
    ampi_tape->arg[2] = sendcnt;
    ampi_tape->arg[1] = recvcnt;
    ampi_tape->arg[0] = root;
    ampi_tape->comm = comm;
    ampi_create_tape_entry((void*)ampi_tape);;
    /*recvbuf entry*/

    for(i=0 ; i<recvcnt ; i=i+1) {
      ampi_get_idx(recvbuf, &i, &ampi_tape->idx[sendcnt*size+i]);
    }
  }

  AMPI_Scatter_f(tmp_send, sendcnt, sendtype, tmp_recv, recvcnt, recvtype, root, comm);

  for(i=0 ; i<recvcnt ; i=i+1) {
    ampi_set_val(recvbuf, &i, &tmp_recv[i]);
  }

  free(tmp_send);
  free(tmp_recv);
  return 0;

}

int AMPI_Scatterv(void *sendbuf, int *sendcnts, int *displs, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm) {
#ifdef AMPI_COUNT_COMMS
    ampi_comm_count=ampi_comm_count+1;
#endif
  if(sendtype !=AMPI_DOUBLE || recvtype != AMPI_DOUBLE) {
    return MPI_Scatterv(sendbuf, sendcnts, displs, sendtype, recvbuf, recvcnt, recvtype, root, comm);
  }
    int i=0;
    int size=0;
    int rank=0;
    int max_size=0;
    MPI_Comm_size(comm,&size);
    MPI_Comm_size(comm,&rank);
    /* Allocate maximum size of sendbuf */
    for(i=0;i<size;i=i+1) {
      if (sendcnts[i] + displs[i] > max_size) {
        max_size = sendcnts[i] + displs[i];
      }
    }
    double * tmp_send = (double*) malloc(sizeof(double)*max_size);
    double * tmp_recv = (double*) malloc(sizeof(double)*recvcnt);


    for(i=0 ; i<max_size ; i=i+1) {
      ampi_get_val(sendbuf,&i,&tmp_send[i]);
    }

    if (ampi_is_tape_active()){
      ampi_tape_entry* ampi_tape = ampi_create_tape(max_size + recvcnt + 1);
      ampi_tape->arg=(int*) malloc(sizeof(int)*(3+2*size));
      ampi_tape->comm = comm;

      /*sendbuf dummies*/

      for(i=0 ; i<max_size ; i=i+1) {
        ampi_get_val(sendbuf,&i,&tmp_send[i]);
      }

      for(i=0 ; i<size*2 ; i=i+2) {
        ampi_tape->arg[i+3] = sendcnts[i/2];
        ampi_tape->arg[i+4] = displs[i/2];
      }
      ampi_tape->oc = SCATTERV;
      ampi_tape->arg[2] = max_size;
      ampi_tape->arg[1] = recvcnt;
      ampi_tape->arg[0] = root;
      ampi_tape->comm = comm;

      /*recvbuf entry*/

      for(i=0 ; i<recvcnt ; i=i+1) {
          ampi_get_idx(recvbuf, &i, &ampi_tape->idx[max_size + i]);
      }
    }

    /*AMPI_Scatterv_f(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm);*/

    for(i=0 ; i<recvcnt ; i=i+1) {
        ampi_set_val(recvbuf, &i, &tmp_recv[i]);
    }
    free(tmp_send);
    free(tmp_recv);
    return 0;

}

int AMPI_Gather(void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm) {
#ifdef AMPI_COUNT_COMMS
    ampi_comm_count=ampi_comm_count+1;
#endif
  if(sendtype !=AMPI_DOUBLE || recvtype != AMPI_DOUBLE) {
    return MPI_Gather(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm);
  }
  int i=0;
  int size=0;
  int rank=0;
  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);
  double * tmp_send = (double*)malloc(sizeof(double)*sendcnt);
  double * tmp_recv = 0;
  if(rank == root) {
    tmp_recv = (double*)malloc(sizeof(double)*recvcnt*size);
  }

  /* check for size */
  int recvSize=recvcnt*size;
  int bufferSize = sendcnt+1;
  if(rank == root) {
          bufferSize = recvcnt*size+sendcnt+1;
  }

  for(i=0 ; i<sendcnt ; i=i+1) {
    ampi_get_val(sendbuf,&i,&tmp_send[i]);
  }

  if (ampi_is_tape_active()){
    ampi_tape_entry* ampi_tape = ampi_create_tape(bufferSize);

    ampi_tape->arg=(int*)malloc(sizeof(int)*3);
    ampi_tape->comm = comm;


    if(rank == root) {
      ampi_create_dummies(recvbuf, &recvSize);
    }

    /*sendbuf dummies*/

    for(i=0 ; i<sendcnt ; i=i+1) {
      ampi_get_idx(sendbuf, &i, &ampi_tape->idx[i]);
    }

    ampi_tape->oc = GATHER;
    ampi_tape->arg[2] = sendcnt;
    ampi_tape->arg[1] = recvcnt;
    ampi_tape->arg[0] = root;
    ampi_tape->comm = comm;
    ampi_create_tape_entry((void*)ampi_tape);;

    /*recvbuf entry*/

    if(rank == root) {
      for(i=0 ; i<recvcnt*size ; i=i+1) {
        ampi_get_idx(recvbuf, &i, &ampi_tape->idx[sendcnt + i]);
      }
    }
  }

  AMPI_Gather_f(tmp_send, sendcnt, sendtype, tmp_recv, recvcnt, recvtype, root, comm);

  if(rank == root) {
    for(i=0 ; i<recvcnt*size ; i=i+1) {
      ampi_set_val(recvbuf, &i, &tmp_recv[i]);
    }
  }
  free(tmp_send);
  if( rank == root ) {
    free(tmp_recv);
  }
  return 0;

}

int AMPI_Gatherv(void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int *recvcnts, int *displs, MPI_Datatype recvtype, int root, MPI_Comm comm) {
#ifdef AMPI_COUNT_COMMS
    ampi_comm_count=ampi_comm_count+1;
#endif
  if(sendtype !=AMPI_DOUBLE || recvtype != AMPI_DOUBLE) {
    return MPI_Gatherv(sendbuf, sendcnt, sendtype, recvbuf, recvcnts, displs, recvtype, root, comm);
  }

  int i=0;
  int size=0;
  int rank=0;
  int max_size=0;
  MPI_Comm_size(comm,&size);
  MPI_Comm_size(comm,&rank);
  /* Determine and allocate maximum size of recvbuf */
  for(i=0;i<size;i=i+1) {
    if (recvcnts[i] + displs[i] > max_size) {
      max_size = recvcnts[i] + displs[i];
    }
  }
  double * tmp_send = (double*) malloc(sizeof(double)*sendcnt);
  double * tmp_recv = (double*) malloc(sizeof(double)*max_size);

  for(i=0 ; i<sendcnt ; i=i+1) {
    ampi_get_val(sendbuf,&i,&tmp_send[i]);
  }

  if (ampi_is_tape_active()){
    ampi_tape_entry* ampi_tape = ampi_create_tape(max_size*size+sendcnt+1);
    ampi_tape->arg=(int*) malloc(sizeof(int)*3);
    ampi_tape->comm = comm;

    /*sendbuf dummies*/

    for(i=0 ; i<sendcnt ; i=i+1) {
      ampi_get_idx(sendbuf, &i, &ampi_tape->idx[i]);
    }
    ampi_tape->oc = GATHERV;
    ampi_tape->arg[2] = sendcnt;
    ampi_tape->arg[1] = max_size;
    ampi_tape->arg[0] = root;
    ampi_tape->comm = comm;


    /*recvbuf entry*/

    for(i=0 ; i<max_size ; i=i+1) {
        ampi_get_idx(recvbuf, &i, &ampi_tape->idx[sendcnt+i]);
    }
  }

  /* TODO: implement primal evaluation */
  /*AMPI_Gather(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm);*/

  for(i=0 ; i<max_size ; i=i+1) {
      ampi_set_val(recvbuf, &i, &tmp_recv[i]);
  }
  free(tmp_send);
  free(tmp_recv);
  return 0;
}

int AMPI_Recv_init(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *mpi_request) {
#ifdef AMPI_COUNT_COMMS
    ampi_comm_count=ampi_comm_count+1;
#endif
  if(datatype!=AMPI_DOUBLE) {
    return MPI_Recv_init(buf, count, datatype, source, tag, comm, mpi_request);
  }
    AMPI_Request request;
    AMPI_ht_el *ht_el=(AMPI_ht_el*) malloc(sizeof(AMPI_ht_el));
    ht_el->key=mpi_request;
    request.buf = buf;
    ht_el->request=request;

    if (ampi_is_tape_active()){
      ampi_tape_entry* ampi_tape = ampi_create_tape(1);
      ampi_tape->arg=(int*) malloc(sizeof(int)*2);
      request.va=ampi_tape;

      ampi_create_tape_entry((void*)ampi_tape);;
      ampi_tape->oc = RECV_INIT;
      ampi_tape->arg[0] = count;
      ampi_tape->arg[1] = source;
      ampi_tape->comm = comm;
      ampi_tape->tag = tag;
    }

    HASH_ADD_PTR(AMPI_ht,key,ht_el);
    return 0;
}
int AMPI_Send_init(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *mpi_request) {
#ifdef AMPI_COUNT_COMMS
    ampi_comm_count=ampi_comm_count+1;
#endif
  if(datatype!=AMPI_DOUBLE) {
    return MPI_Send_init(buf, count, datatype, dest, tag, comm, mpi_request);
  }
    AMPI_Request request;
    AMPI_ht_el *ht_el=(AMPI_ht_el*) malloc(sizeof(AMPI_ht_el));
    ht_el->key=mpi_request;
    request.buf = buf;
    ht_el->request=request;

    if (ampi_is_tape_active()){
      ampi_tape_entry* ampi_tape = ampi_create_tape(1);
      ampi_tape->arg=(int*) malloc(sizeof(int)*2);
      request.va=ampi_tape;

      ampi_create_tape_entry((void*)ampi_tape);;
      ampi_tape->oc = SEND_INIT;
      ampi_tape->arg[0] = count;
      ampi_tape->arg[1] = dest;
      ampi_tape->comm = comm;
      ampi_tape->tag = tag;
    }
    HASH_ADD_PTR(AMPI_ht,key,ht_el);
    return 0;
}
int AMPI_Sendrecv_replace(void *buf, int count, MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag, MPI_Comm comm, MPI_Status *status) {
#ifdef AMPI_COUNT_COMMS
    ampi_comm_count=ampi_comm_count+1;
#endif
  if(datatype!=AMPI_DOUBLE) {
    return MPI_Sendrecv_replace(buf, count, datatype, dest, sendtag, source, recvtag, comm, status);
  }
    int i=0;
    double * tmp = (double*) malloc(sizeof(double)*count);
    for(i=0;i<count;i=i+1) {
        ampi_get_val(buf,&i,&tmp[i]);
    }

    if (ampi_is_tape_active()){
      /*create actual MPI entry*/
      ampi_tape_entry* ampi_tape = ampi_create_tape(2*count+2);
      ampi_tape->arg=(int*) malloc(sizeof(int)*4);

      for(i=0;i<count;i=i+1) {
          ampi_get_idx(buf, &i, &ampi_tape->idx[i]);
      }
      ampi_tape->oc = SENDRECVREPLACE;
      ampi_tape->arg[0] = count;
      ampi_tape->arg[1] = dest;
      ampi_tape->arg[2] = source;

      ampi_tape->comm = comm;
      ampi_tape->tag = sendtag;

      ampi_create_dummies(buf, &count);
      ampi_create_tape_entry((void*)ampi_tape);;


      if(status!=MPI_STATUS_IGNORE) {
        ampi_tape->arg[3] = status->MPI_TAG;
      } else {
        ampi_tape->arg[3] = recvtag;
      }

      for(i=0;i<count;i=i+1) {
          ampi_get_idx(buf, &i, &ampi_tape->idx[count+i]);
      }
    }

    int temp=AMPI_Sendrecv_replace_f(tmp, count, datatype, dest, sendtag, source, recvtag, comm, status);

    for(i=0;i<count;i=i+1) {
        ampi_set_val(buf, &i, &tmp[i]);
    }
    free(tmp);
    return temp;
}

int AMPI_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf, int recvcount, MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status *status){
#ifdef AMPI_COUNT_COMMS
    ampi_comm_count=ampi_comm_count+1;
#endif
  if(recvtype!=AMPI_DOUBLE) {
    return MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf,recvcount,recvtype,source,recvtag,comm,status);
  }

    double * sendtmp = (double*) malloc(sizeof(double)*sendcount);
    double * recvtmp = (double*) malloc(sizeof(double)*recvcount);

    int i=0;

    for(i=0;i<sendcount;i=i+1) {
        ampi_get_val(sendbuf,&i,&sendtmp[i]);
    }

    if (ampi_is_tape_active()){
      /*create actual MPI entry*/
      ampi_tape_entry* ampi_tape = ampi_create_tape(sendcount+recvcount+2);
      ampi_tape->arg = (int*) malloc(sizeof(int)*6);


      for(i=0;i<sendcount;i=i+1) {
          ampi_get_idx(sendbuf, &i, &ampi_tape->idx[i]);
      }

      ampi_tape->oc = SENDRECV;
      ampi_tape->arg[0] = sendcount;
      ampi_tape->arg[1] = dest;
      ampi_tape->arg[2] = sendtag;
      ampi_tape->arg[3] = recvcount;
      ampi_tape->arg[4] = source;
      ampi_tape->arg[5] = recvtag;
      ampi_tape->comm   = comm;

      ampi_create_dummies(recvbuf, &recvcount);
      ampi_create_tape_entry((void*)ampi_tape);

      for(i=0;i<recvcount;i=i+1) {
        ampi_get_idx(recvbuf, &i, &ampi_tape->idx[sendcount+i]);
      }
    }

    int temp=AMPI_Sendrecv_f(sendtmp, sendcount, sendtype, dest, sendtag, recvtmp,recvcount,recvtype,source,recvtag,comm,status);

    for(i=0;i<recvcount;i=i+1) {
      ampi_set_val(recvbuf, &i, &recvtmp[i]);
    }

    free(recvtmp);
    free(sendtmp);
    return temp;
}

int AMPI_Start(MPI_Request *mpi_request) {
#ifdef AMPI_COUNT_COMMS
    ampi_comm_count=ampi_comm_count+1;
#endif
    AMPI_ht_el *ht_req=NULL;
    int ret=0;
    AMPI_Request request;
    /*AMPI_Request *request=malloc(sizeof(AMPI_Request));*/
    void *addr=mpi_request;
    HASH_FIND_PTR(AMPI_ht,&addr,ht_req);
    if(ht_req) {
        request=ht_req->request;
        if (ampi_is_tape_active()){
          ampi_tape_entry* ampi_tape=0;
          ampi_tape=request.va;
          if(ampi_tape->oc==SEND_INIT) {
              ret=AMPI_Isend(request.buf,ampi_tape->arg[0],MPI_DOUBLE,ampi_tape->arg[1],ampi_tape->tag,ampi_tape->comm,mpi_request);
              if(addr!=mpi_request) printf("changed\n");
              HASH_FIND_PTR(AMPI_ht,&addr,ht_req);
              return ret;
          }
          else if(ampi_tape->oc==RECV_INIT) {
              ret=AMPI_Irecv(request.buf,ampi_tape->arg[0],MPI_DOUBLE,ampi_tape->arg[1],ampi_tape->tag,ampi_tape->comm,mpi_request);
              if(addr!=mpi_request) printf("changed\n");
              HASH_FIND_PTR(AMPI_ht,&addr,ht_req);
              return ret;
          }
          else {
              printf("Active start No opcode: %d\n",ampi_tape->oc);
          }
        }
    } else {
        return MPI_Start(mpi_request);
    }
    return -1;
}

int AMPI_Startall(int count, MPI_Request array_of_requests[]) {
#ifdef AMPI_COUNT_COMMS
    ampi_comm_count=ampi_comm_count+1;
#endif
    int i=0;
    for(i=0;i<count;i=i+1) AMPI_Start(&array_of_requests[i]);
    return 0;
}

void ampi_interpret_tape(void* handle){
    ampi_tape_entry* ampi_tape = (ampi_tape_entry*)handle;
    int j=0;
    double *tmp_d;
    double *tmp_d_recv;
    double *tmp_d_send;
    MPI_Comm comm;
    MPI_Op op = MPI_SUM;
    comm = MPI_COMM_WORLD;
    MPI_Status status;
    /*ampi_vac=ampi_vac-1;*/
    /*while(ampi_tape[ampi_vac].oc == MPI_DUMMY)*/
    /*ampi_vac=ampi_vac-1;*/
    /*i=ampi_vac;*/
#ifdef DEBUG
    printf("AMPI_TAPE Interpreter OC: %d\n", ampi_tape[i].oc);
    printf("--------------------------------START\n");
#endif
    switch(ampi_tape->oc){
  case SEND : {
      tmp_d = (double*) malloc(sizeof(double)*ampi_tape->arg[0]);
#ifdef NO_COMM_WORLD
      comm=ampi_tape->comm;
#endif
      AMPI_Send_b(tmp_d, ampi_tape->arg[0], MPI_DOUBLE, ampi_tape->arg[1], ampi_tape->tag, comm);
      for(j=0;j<ampi_tape->arg[0];j=j+1) {
          ampi_set_adj(&ampi_tape->idx[j], &tmp_d[j]);
#ifdef DEBUG
          printf("SEND: %d %e\n", tmp_entry->idx, tmp_d[j]);
#endif
      }
      free(tmp_d);
      break;
        }
  case RECV : {
      tmp_d = (double*) malloc(sizeof(double)*ampi_tape->arg[0]);
      for(j=0;j<ampi_tape->arg[0];j=j+1) {
          ampi_get_adj(&ampi_tape->idx[j], &tmp_d[j]);
#ifdef DEBUG
          printf("RECV: %d %e\n", tmp_entry->idx, tmp_d[j]);
#endif
      } 
#ifdef NO_COMM_WORLD
      comm=ampi_tape->comm;
#endif
      AMPI_Recv_b(tmp_d, ampi_tape->arg[0], MPI_DOUBLE, ampi_tape->arg[1], ampi_tape->tag, comm,&status);
      free(tmp_d);
      break;
        }
  case BRECV : {
      tmp_d = (double*) malloc(sizeof(double)*ampi_tape->arg[0]);
      for(j=0;j<ampi_tape->arg[0];j=j+1) {
          ampi_get_adj(&ampi_tape->idx[j], &tmp_d[j]);
#ifdef DEBUG
          printf("BRECV: %d %e\n", tmp_entry->idx, tmp_d[j]);
#endif
      } 
#ifdef NO_COMM_WORLD
      comm=ampi_tape->comm;
#endif
      AMPI_Brecv_b(tmp_d, ampi_tape->arg[0], MPI_DOUBLE, ampi_tape->arg[1], ampi_tape->tag, comm,&status);
      free(tmp_d);
      break;
        }
  case ISEND : {
       /*tmp_d = malloc(sizeof(double)*ampi_tape[i].arg[0]);*/
       /*if(!tape[i].request->r.aw) {*/
       tmp_d = (double*) ampi_tape->request->a;
#ifdef NO_COMM_WORLD
      comm=ampi_tape->comm;
#endif
       AMPI_Isend_b(tmp_d, ampi_tape->arg[0], MPI_DOUBLE, ampi_tape->arg[1], ampi_tape->tag, comm,ampi_tape->request);
       for(j = 0 ; j < ampi_tape->arg[0] ; j++) {
           ampi_set_adj(&ampi_tape->idx[j], &tmp_d[j]);
         }
       free(tmp_d);
       free(ampi_tape->request->mpiRequest);
       break;
         }
  case IRECV : {
       /*tmp_d = malloc(sizeof(double)*ampi_tape->arg[0]);*/
       tmp_d = (double*) ampi_tape->request->a;
#ifdef NO_COMM_WORLD
      comm=ampi_tape->comm;
#endif
       /*if(tape->request->r.aw) {*/

       /*}*/
       /*else {*/
           AMPI_Irecv_b(tmp_d, ampi_tape->arg[0], MPI_DOUBLE, ampi_tape->arg[1], ampi_tape->tag, comm,ampi_tape->request);
           /*}*/
       free(tmp_d);
       free(ampi_tape->request->mpiRequest);
       break;
         }
  case WAIT : {
      tmp_d = (double*) malloc(sizeof(double)*ampi_tape->request->va->arg[0]);
      if(ampi_tape->request->oc == AMPI_IR) {
          for(j = 0 ; j < ampi_tape->request->va->arg[0] ; j++) {
            ampi_get_adj(&ampi_tape->request->va->idx[j], &tmp_d[j]);
          }
#ifdef DEBUG
          printf("AMPI_Wait_interpret: ");
          printf("%d ", ampi_tape[ampi_tape->arg[0]].arg[0]);
          for(j = 0 ; j < ampi_tape[ampi_tape->arg[0]].arg[0] ; j++) {
        printf("%e ", tmp_d[j]);
          }
          printf("\n");
#endif
      }
      ampi_tape->request->a = tmp_d;
      ampi_tape->request->mpiRequest =(MPI_Request*) malloc(sizeof(MPI_Request));
      ampi_tape->request->tag=ampi_tape->request->va->tag;
#ifndef NO_COMM_WORLD 
      ampi_tape->request->comm=MPI_COMM_WORLD;
#endif
      AMPI_Wait_b(ampi_tape->request,&status);
      ampi_tape->request->va->request = ampi_tape->request;
      break;
        }
  case BCAST : {
           int rank=0;
           int root=ampi_tape->arg[1];
#ifdef NO_COMM_WORLD
      comm=ampi_tape->comm;
#endif
           MPI_Comm_rank(comm,&rank);
           
      tmp_d = (double*) malloc(sizeof(double)*ampi_tape->arg[0]);
      if(rank!=root) {
        for(j=0;j<ampi_tape->arg[0];j=j+1) {
          ampi_get_adj(&ampi_tape->idx[j], &tmp_d[j]);
        } 
      }
      else {
        for(j=0;j<ampi_tape->arg[0];j=j+1) tmp_d[j]=0;
      }
       AMPI_Bcast_b(tmp_d, ampi_tape->arg[0], MPI_DOUBLE, root, comm);
      if(rank==root) {
         for(j=0 ; j<ampi_tape->arg[0] ; j++) {
           ampi_set_adj(&ampi_tape->idx[j],&tmp_d[j]);
         }
      }
       free(tmp_d);
       break;
         }
  case REDUCE : {
        if(ampi_tape->arg[2] == AMPI_REDUCE_ADD)
            op = MPI_SUM;
        if(ampi_tape->arg[2] == AMPI_REDUCE_MUL)
            op = MPI_PROD;
        if(ampi_tape->arg[2] == AMPI_REDUCE_MIN)
            op = MPI_MIN;
        if(ampi_tape->arg[2] == AMPI_REDUCE_MAX)
            op = MPI_MAX;
        tmp_d_send = (double*) malloc(sizeof(double)*ampi_tape->arg[0]);
        tmp_d_recv = (double*) malloc(sizeof(double)*ampi_tape->arg[0]);
        for(j=0;j<ampi_tape->arg[0];j=j+1) {
          ampi_get_adj(&ampi_tape->idx[ampi_tape->arg[0] + j], &tmp_d_recv[j]);
        }
        for(j=0;j<ampi_tape->arg[0];j=j+1) {
          tmp_d_send[j]=0;
        }

#ifdef NO_COMM_WORLD
      comm=ampi_tape->comm;
#endif
        AMPI_Reduce_b(tmp_d_send, tmp_d_recv, ampi_tape->arg[0], MPI_DOUBLE, op, ampi_tape->arg[1], comm, ampi_tape->stack);
        for(j=0 ; j<ampi_tape->arg[0] ; j++) {
            ampi_set_adj(&ampi_tape->idx[j],&tmp_d_send[j]);
        }
        free(tmp_d_send);
        free(tmp_d_recv);
        break;
          }
  case ALLREDUCE : {
        if(ampi_tape->arg[2] == AMPI_REDUCE_ADD)
            op = MPI_SUM;
        if(ampi_tape->arg[2] == AMPI_REDUCE_MUL)
            op = MPI_PROD;
        if(ampi_tape->arg[2] == AMPI_REDUCE_MIN)
            op = MPI_MIN;
        if(ampi_tape->arg[2] == AMPI_REDUCE_MAX)
            op = MPI_MAX;
        tmp_d_send = (double*) malloc(sizeof(double)*ampi_tape->arg[0]);
        tmp_d_recv = (double*) malloc(sizeof(double)*ampi_tape->arg[0]);
        for(j=0;j<ampi_tape->arg[0];j=j+1)
            ampi_get_adj(&ampi_tape->idx[ampi_tape->arg[0] + j], &tmp_d_recv[j]);
#ifdef NO_COMM_WORLD
      comm=ampi_tape->comm;
#endif

#ifdef DEBUG
        printf("AMPI_Allreduce tmp_d_recv: ");
        for(j=0 ; j<ampi_tape->arg[0] ; j++) {
            /*if(tmp_d_recv[j]!=tmp_d_recv[j]) tmp_d_recv[j]=0;*/
            printf("%e ",tmp_d_recv[j]);
        }
        printf("\n");
#endif
        AMPI_Allreduce_b(tmp_d_send, tmp_d_recv, ampi_tape->arg[0], MPI_DOUBLE, op, comm, ampi_tape->stack);
#ifdef DEBUG
        printf("AMPI_Allreduce tmp_d_send: ");
        for(j=0 ; j<ampi_tape->arg[0] ; j++) {
            printf("%e ",tmp_d_send[j]);
        }
        printf("\n");
#endif
        for(j=0 ; j<ampi_tape->arg[0] ; j++) {
            ampi_set_adj(&ampi_tape->idx[j],&tmp_d_send[j]);
        }
        free(tmp_d_send);
        free(tmp_d_recv);
        break;
          }
  case SCATTER : {
#ifdef NO_COMM_WORLD
       comm=ampi_tape->comm;
#endif
       int size=0;
       int root=ampi_tape->arg[0];
       int crecv=ampi_tape->arg[1];
       int csend=ampi_tape->arg[2];
       MPI_Comm_size(comm,&size);
       double *sendbuf=(double*) malloc(sizeof(double)*csend*size);
       double *recvbuf=(double*) malloc(sizeof(double)*crecv);
       for(j=0;j<crecv;j++) {
         ampi_get_adj(&ampi_tape->idx[csend*size+j],&recvbuf[j]);
       }
       AMPI_Scatter_b(sendbuf,csend,MPI_DOUBLE,recvbuf,crecv,MPI_DOUBLE,root,comm);
       for(j=0;j<csend*size;j++) {
         ampi_set_adj(&ampi_tape->idx[j],&sendbuf[j]);
       }
       free(sendbuf);
       free(recvbuf);
       break;
           }
  case SCATTERV : {
         break;
           }
  case GATHER : {
#ifdef NO_COMM_WORLD
       comm=ampi_tape->comm;
#endif
       int size=0;
       int rank=0;
       int root=ampi_tape->arg[0];
       int crecv=ampi_tape->arg[1];
       int csend=ampi_tape->arg[2];
       MPI_Comm_size(comm,&size);
       MPI_Comm_rank(comm,&rank);
       double *sendbuf=(double*)malloc(sizeof(double)*csend);
       double *recvbuf=NULL;
       if(rank == root) {
         recvbuf = (double*)malloc(sizeof(double)*crecv*size);
         for(j=0;j<crecv*size;j++) {
           ampi_get_adj(&ampi_tape->idx[csend + j],&recvbuf[j]);
         }
       }

       AMPI_Gather_b(sendbuf,csend,MPI_DOUBLE,recvbuf,crecv,MPI_DOUBLE,root,comm);
       for(j=0;j<csend;j++) {
         ampi_set_adj(&ampi_tape->idx[j],&sendbuf[j]);
       }
       free(sendbuf);
       if(rank == root) {
         free(recvbuf);
       }
         break;
           }
  case GATHERV : {
         break;
           }
  case SENDRECVREPLACE : {
    int count = ampi_tape->arg[0];
    tmp_d = (double*) malloc(sizeof(double)* count);
      /*take adjoints out of the tape in the send buffer*/
      for(j=0;j<count;j=j+1) {
          ampi_get_adj(&ampi_tape->idx[count + j], &tmp_d[j]);
#ifdef DEBUG
          printf("SENDRECV_B: %d %e\n", tmp_entry->idx, tmp_d[j]);
#endif
      } 
#ifdef NO_COMM_WORLD
      comm=ampi_tape->comm;
#endif
      AMPI_Sendrecv_replace_b(tmp_d, count, MPI_DOUBLE, ampi_tape->arg[2], ampi_tape->arg[3], ampi_tape->arg[1], ampi_tape->tag, comm, &status);

                  /*two entries for sendrecvreplace, decrease tape index*/
      for(j=0;j<count;j=j+1) {
          ampi_set_adj(&ampi_tape->idx[j], &tmp_d[j]);
#ifdef DEBUG
          printf("SEND: %d %e\n", tmp_entry->idx, tmp_d[j]);
#endif
      }
      free(tmp_d);
      break;
       }
      case SENDRECV: {
          int recvcount = ampi_tape->arg[3];
          int sendcount = ampi_tape->arg[0];
          double *sendbuf = (double*) malloc(sizeof(double)*sendcount);
          double *recvbuf = (double*) malloc(sizeof(double)*recvcount);
          for (j=0; j < recvcount; j=j+1){
              ampi_get_adj(&ampi_tape->idx[sendcount+j],&recvbuf[j]);
            }
          comm = ampi_tape->comm;
          AMPI_Sendrecv_b(sendbuf, sendcount, MPI_DOUBLE, ampi_tape->arg[1], ampi_tape->arg[2], recvbuf, recvcount, MPI_DOUBLE, ampi_tape->arg[4], ampi_tape->arg[5],comm,&status);

          for (j=0;j<sendcount;j=j+1){
              ampi_set_adj(&ampi_tape->idx[j], &sendbuf[j]);
            }
          free(sendbuf);
          free(recvbuf);
          break;
        }
  case SEND_INIT : {
      break;
       }
  case RECV_INIT : {
      break;
       }
  default: {
         printf("Warning: Missing opcode in the AMPI tape interpreter for %d.\n", ampi_tape->oc);
         break;
          }

    }
}
void ampi_reset_entry(void* handle){
  ampi_release_tape((ampi_tape_entry*)handle);
}

void ampi_release_tape(ampi_tape_entry* ampi_tape) {
  if(ampi_tape->arg != NULL) {
    free(ampi_tape->arg);
    ampi_tape->arg = NULL;
  }
  if(ampi_tape->idx != NULL) {
    free(ampi_tape->idx);
    ampi_tape->idx = NULL;
  }
  if(ampi_tape->stack != NULL) {
    AMPI_stack_delete(ampi_tape->stack);
    ampi_tape->stack = NULL;
  }
  if(ampi_tape->oc == WAIT) {
    /* wait does not need to delete its request */
    ampi_tape->request = NULL;
  } else {
    if(ampi_tape->request != NULL) {
      free(ampi_tape->request);
      ampi_tape->request = NULL;
    }
  }

  free(ampi_tape);
}

ampi_tape_entry* ampi_create_tape(long int size) {
  ampi_tape_entry* ampi_tape = (ampi_tape_entry*) calloc(1, sizeof(ampi_tape_entry));
  ampi_tape->arg = NULL;
  ampi_tape->idx = NULL;
  ampi_tape->stack = NULL;
  ampi_tape->request = NULL;
  if(0 != size) {
    ampi_tape->idx = (INT64*) malloc(sizeof(INT64) * size);
  }

  return ampi_tape;
}
