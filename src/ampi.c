/* Basic AMPI C library used for overloading and source transformation. All MPI routines
 * are subdivided into their forward _f and backward _b counterpart. The forward routines
 * are called during the forward/taping run. The backward routines are called during the
 * reverse/interpretation run.
 */

/*#define DEBUG*/
#include <ampi.h>

int * AMPI_myid;
/* Stack to compute the adjoints of the MPI reduce operation */
ampi_stack reduce_stack;

/* Non communication routines (Init, Finalize...)*/

int AMPI_Init_f(int *argc, char ***argv) {
    int ierr;
    stack_init(&reduce_stack);
    ierr=MPI_Init(argc, argv);
    return ierr;
}

int AMPI_Init_b(int *argc, char ***argv) {
    /*destroy(&reduce_stack);*/
    return MPI_Finalize();
}

int AMPI_Comm_size(MPI_Comm comm, int * numprocs) {
    return MPI_Comm_size(comm, numprocs);
}


int AMPI_Comm_rank(MPI_Comm comm, int * myid) {
    AMPI_myid = myid;
    return MPI_Comm_rank(comm, myid);
}

int AMPI_Get_processor_name(char *processor_name, int *namelen) {
    return MPI_Get_processor_name(processor_name, namelen);
}

int AMPI_Barrier(MPI_Comm comm) {
    return MPI_Barrier(comm);
}

int AMPI_Finalize_f(){
    return 1;
}

int AMPI_Finalize_b(){
    return 1;
}

/* Blocking communication */

int AMPI_Send_f(double *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
#ifdef DEBUG
    int i = 0;
    printf("AMPI_Send_f: ");
    for(i = 0 ; i < count ; i++) {
	printf("%f ", buf[i]);
    }
    printf("\n");
#endif
    return MPI_Send(buf, count, datatype, dest, tag, comm);
}

int AMPI_Recv_f(double *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Status *status) {
    int tmp = 0;
#ifdef DEBUG
    int i = 0;
#endif
    tmp = MPI_Recv(buf, count, datatype, dest, tag, comm, status);
#ifdef DEBUG
    printf("AMPI_Recv_f: ");
    for(i = 0 ; i < count ; i++) {
	printf("%f ", buf[i]);
    }
    printf("\n");
#endif
    return tmp;

}

int AMPI_Send_b(double *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
    MPI_Status status;
    return MPI_Recv(buf, count, datatype, dest, tag, comm, &status);
    
}

int AMPI_Recv_b(double *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Status *status) {
    int ierr=0;
    ierr=MPI_Send(buf, count, datatype, dest, tag, comm);
    return ierr;
}

/* Non blocking communication */

int AMPI_Isend_f(double *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, AMPI_Request *request) {
#ifdef DEBUG
    int i = 0;
#endif
    request->v = buf;
    request->dest= dest;
    request->oc = AMPI_IS;
    request->size = count;
    /*request->tag = tag;*/
    request->comm = comm;
#ifdef DEBUG
    printf("AMPI_Isend_f: ");
    for(i = 0 ; i < count ; i++) {
	printf("%f ", buf[i]);
    }
    printf("\n");
#endif
    return MPI_Isend(buf, count, datatype, dest, tag, comm, request->request);
}

int AMPI_Irecv_f(double *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, AMPI_Request *request) {
    request->v = buf;
    request->dest= dest;
    request->oc = AMPI_IR;
    request->size = count;
    /*request->tag = tag;*/
    request->comm = comm;
    return MPI_Irecv(buf, count, datatype, dest, tag, comm, request->request);
}

int AMPI_Wait_f(AMPI_Request *request, MPI_Status *status) {
    int ierr=0;
    ierr=MPI_Wait(request->request, status);
    return ierr;
}

int AMPI_Isend_b(double *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, AMPI_Request *request) {
    int i = 0;
    if(!request->aw) {
	MPI_Wait(request->request, &request->status);
	for(i = 0 ; i < request->size ; i++) {
	    buf[i] = request->a[i];
	}
#ifdef DEBUG
	printf("AMPI_Isend_b: ");
	for(i = 0 ; i < request->size ; i++) {
	    printf("%f ", buf[i]);
	}
	printf("\n");
#endif
	return MPI_SUCCESS;
    } else {
	return MPI_SUCCESS;
    }
}

int AMPI_Irecv_b(double *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, AMPI_Request *request) {
    int i = 0;
    if(!request->aw) {
	MPI_Wait(request->request, &request->status);
	for(i = 0 ; i < request->size ; i++) {
	    buf[i] = request->a[i];
	}
#ifdef DEBUG
	printf("AMPI_Irecv_b: ");
	for(i = 0 ; i < request->size ; i++) {
	    printf("%f ", buf[i]);
	}
	printf("\n");
#endif
	return MPI_SUCCESS;
    } else {
	return MPI_SUCCESS;
    }
}

int AMPI_Wait_b(AMPI_Request *request, MPI_Status * status) {
#ifdef DEBUG
    int i=0;
#endif
    if(request->oc == AMPI_IS) {
	return MPI_Irecv(request->a, request->size, MPI_DOUBLE, request->dest, request->tag, request->comm, request->request);
    }
    else {
	if(request->oc == AMPI_IR) {
#ifdef DEBUG
	    printf("AMPI_Wait_recv: ");
	    for(i=0;i<request->size;i=i+1) {
		printf("%f ",request->a[i]);
	    }
	    printf("\n");
#endif
	    return MPI_Isend(request->a, request->size, MPI_DOUBLE, request->dest, request->tag, request->comm, request->request);
	} else { 
	    printf("Error: OC invalid\n");
	}
    }
    return MPI_SUCCESS;
}

int AMPI_Waitall_f(int count, AMPI_Request *requests, MPI_Status *status) {
    int i = 0;
    int ierr;
    MPI_Request * reqs = malloc(count * sizeof(MPI_Request*)); 

    for(i = 0 ; i < count ; i++) {
	reqs[i] = *requests[i].request;
    }

    ierr = MPI_Waitall(count, reqs, status);

    for(i = 0 ; i < count ; i++) {
	*requests[i].request = reqs[i];
    }
    free(reqs);
    return ierr;
}

int AMPI_Waitall_b(int count, AMPI_Request *requests) {
    MPI_Status status[count];
    int i = 0;
    for(i = 0 ; i < count ; i++) {
	AMPI_Wait_b(&requests[i], &status[i]);
    }
    return MPI_SUCCESS;
}

/* Anti Waitall. Experimental */

int AMPI_Awaitall_f(int count, AMPI_Request *requests, MPI_Status *status) {
    int i = 0;
    for(i = 0 ; i < count ; i++) {
	requests[i].aw = 1;
    }
    return MPI_SUCCESS;
}

int AMPI_Awaitall_b(int count, AMPI_Request *requests) {
    int i = 0;
    for(i = 0 ; i < count ; i++) {
	if(requests[i].aw) {
	    MPI_Wait(requests[i].request, &requests[i].status);
	}
    }
    return MPI_SUCCESS;
}

/* Collective Communication */

int AMPI_Bcast_f(double *buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm) {
    int ierr=0;
#ifdef DEBUG
    printf("AMPI_Bcast_f before: ");
    for(i = 0 ; i < count ; i++) {
	printf("%f ", buf[i]);
    }
    printf("\n");
#endif
    ierr=MPI_Bcast(buf, count, datatype, root, comm);
    return ierr;
}

int AMPI_Bcast_b(double *buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm){
    int ierr=0;
    int i=0;
    int rank=0;
    MPI_Comm_rank(comm,&rank);

    double *sendbuf=malloc(sizeof(double)*count);
    double *recvbuf=malloc(sizeof(double)*count);
    if(rank!=root) {
      for(i=0;i<count;i=i+1) sendbuf[i]=buf[i];
      for(i=0;i<count;i=i+1) recvbuf[i]=0;
    }
    else {
      for(i=0;i<count;i=i+1) {
	sendbuf[i]=0;
      }
      for(i=0;i<count;i=i+1) recvbuf[i]=0;
    }
    ierr=MPI_Reduce(sendbuf,recvbuf,1,MPI_DOUBLE,MPI_SUM,root,comm);
    if(rank==root) {
      for(i=0;i<count;i=i+1) buf[i]=recvbuf[i];
    }
    free(sendbuf);
    free(recvbuf);
    return ierr;
}

int AMPI_Reduce_f(double *sendbuf, double *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {
    int myid=0;
    double *recvbuf_tmp=0;
    AMPI_Tupel *sendtupel=0;
    AMPI_Tupel *recvtupel=0;
    int i=0;
    if(op == MPI_PROD || op == MPI_SUM) {
	recvbuf_tmp=malloc(sizeof(double)*count);
	for(i=0; i<count ; i=i+1) {
	    push(&reduce_stack,sendbuf[i]);
	}
	MPI_Allreduce(sendbuf, recvbuf_tmp, count, datatype, op, comm);
	for(i=0; i<count; i=i+1)
	    push(&reduce_stack,recvbuf_tmp[i]);
	/*MPI_Comm_rank(MPI_COMM_WORLD, &myid);*/
	if(myid==root) {
	    for(i=0; i<count; i=i+1) {
		recvbuf[i]=recvbuf_tmp[i];
	    }
	}
	free(recvbuf_tmp);
	return MPI_SUCCESS;
    }
    if(op == MPI_MIN || op == MPI_MAX) {
	sendtupel=malloc(sizeof(AMPI_Tupel)*count);
	recvtupel=malloc(sizeof(AMPI_Tupel)*count);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	for(i=0; i<count; i=i+1) {
	    sendtupel[i].v=sendbuf[i];
	    sendtupel[i].j=myid;
	}
	if(op == MPI_MAX)
	    MPI_Allreduce(sendtupel, recvtupel, count, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
	else
	    MPI_Allreduce(sendtupel, recvtupel, count, MPI_DOUBLE_INT, MPI_MINLOC, comm);
	for(i=0; i<count; i=i+1) {
	    push(&reduce_stack,(double) recvtupel[i].j);
	}
	/*push(&reduce_stack,(double) count);*/
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	if(myid==root) {
	    for(i=0; i<count; i=i+1) {
		recvbuf[i]=recvtupel[i].v;
	    }
	}
	free(sendtupel);
	free(recvtupel);
	return MPI_SUCCESS;
    }
    return -1;
}

int AMPI_Reduce_b(double *sendbuf, double *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {
    int i=0;
    double *s = 0;
    double *s_d = 0;
    double *r = 0;
    int idx=0;
    int myid=0;
    if(op == MPI_PROD) {
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	    s = malloc(sizeof(double)*count);
	    s_d = malloc(sizeof(double)*count);
	    r = malloc(sizeof(double)*count);
	    MPI_Bcast(recvbuf, count, datatype, root, comm);
	    for(i=count-1 ; i>=0 ; i=i-1) {
		r[i] = pop(&reduce_stack);
	    }
	    for(i=count-1 ; i>=0 ; i=i-1) {
		s[i] = pop(&reduce_stack);
	    }
	    for(i=0; i<count; i=i+1) {
		s_d[i] = recvbuf[i];
	    }
	    for(i=0 ; i<count ; i=i+1) {
		if(r[i]*s_d[i] == 0.0){
		    sendbuf[i]=0;
		    printf("--------------------------------------\n");
		    printf("Arithmetic Exeption in adjoint reduce!\n");
		    printf("Result set to 0. Hoping for the best.\n");
		    printf("--------------------------------------\n");
		    return 1;
		}
		else {
		    sendbuf[i] = (r[i]/s[i])*s_d[i];
		}
	    }
	    free(s);
	    free(s_d);
	    free(r);
	    return MPI_SUCCESS;
    }
    if(op == MPI_SUM) {
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	s = malloc(sizeof(double)*count);
	s_d = malloc(sizeof(double)*count);
	r = malloc(sizeof(double)*count);
	MPI_Bcast(recvbuf, count, datatype, root, comm);
	for(i=count-1 ; i>=0 ; i=i-1) {
	    r[i] = pop(&reduce_stack);
	}
	for(i=count-1 ; i>=0 ; i=i-1) {
	    s[i] = pop(&reduce_stack);
	}
	for(i=0; i<count; i=i+1)
	    s_d[i] = recvbuf[i];
	for(i=0 ; i<count ; i=i+1) {
	    sendbuf[i] = s_d[i];
	}
	free(s);
	free(s_d);
	free(r);
	return MPI_SUCCESS;
    }
    if(op == MPI_MIN || op == MPI_MAX) {
	MPI_Bcast(recvbuf, count, datatype, root, comm);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	for(i=count-1;i>=0;i--) {
	    idx=(int)pop(&reduce_stack);
	    if(myid!=idx)
		sendbuf[i]=0;
	    else {
		sendbuf[i]=recvbuf[i];
	    }
	}
	return MPI_SUCCESS;
    }
    return -1;
}

int AMPI_Allreduce_f(double *sendbuf, double *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
    int i=0;
    int myid=0;
    AMPI_Tupel *sendtupel=0;
    AMPI_Tupel *recvtupel=0;
    if(op == MPI_PROD || op == MPI_SUM) {
	for(i=0; i<count ; i=i+1) {
	    push(&reduce_stack,sendbuf[i]);
	}

	MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
	for(i=0; i<count; i=i+1)
	    push(&reduce_stack,recvbuf[i]);
	return MPI_SUCCESS;
    }
    if(op == MPI_MIN || op == MPI_MAX) {
	sendtupel=malloc(sizeof(AMPI_Tupel)*count);
	recvtupel=malloc(sizeof(AMPI_Tupel)*count);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	for(i=0; i<count; i=i+1) {
	    sendtupel[i].v=sendbuf[i];
	    sendtupel[i].j=myid;
	}
	if(op == MPI_MAX)
	    MPI_Allreduce(sendtupel, recvtupel, count, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
	else
	    MPI_Allreduce(sendtupel, recvtupel, count, MPI_DOUBLE_INT, MPI_MINLOC, comm);
	for(i=0; i<count; i=i+1) {
	    push(&reduce_stack,(double) recvtupel[i].j);
	}

	/*push(&reduce_stack,(double) count);*/
	for(i=0; i<count; i=i+1) {
	    recvbuf[i]=recvtupel[i].v;
	}
	free(sendtupel);
	free(recvtupel);
	return MPI_SUCCESS;
    }
    return -1;
}

int AMPI_Allreduce_b(double *sendbuf, double *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
    int i=0;
    int j=0;
    int myid=0;
    int numprocs=0;
    int idx=0;
    double **recvbuf_tmp = 0;
    double *minmaxbuf_tmp = 0;
    MPI_Request *requests = 0;
    MPI_Status *status = 0;
    double *s = 0;
    double *s_d = 0; 
    double *r = 0;
    /*double tmp=0;*/
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
#ifdef DEBUG
    printf("AMPI_Allreduce_b count: %d\n",count);
    for(i=0 ; i<count ; i++) {
	printf("AMPI_Allreduce_b recvbuf: %e\n", recvbuf[i]);
    }
#endif
    if(op == MPI_PROD || op == MPI_SUM) {        
    if(myid==0) {
        requests = malloc(sizeof(MPI_Request)*(numprocs-1));
        status = malloc(sizeof(MPI_Status)*(numprocs-1));
	recvbuf_tmp = malloc(sizeof(double*)*numprocs);
	recvbuf_tmp[0] = recvbuf;
        for(i=1;i<numprocs;i++) {
	    recvbuf_tmp[i] = malloc(sizeof(double)*count);
	    MPI_Irecv(recvbuf_tmp[i], count, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &requests[i-1]);
        }
	MPI_Waitall(numprocs-1,requests,status);
	free(requests);
	for(j=0;j<count;j++) {
	    for(i=1;i<numprocs;i++) {
		recvbuf[j]=recvbuf[j]+recvbuf_tmp[i][j];
	    }
	}
	free(recvbuf_tmp);
    }
    else {
	i=0;
	MPI_Send(recvbuf, count, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
    }
    MPI_Bcast(recvbuf, count, datatype, 0, MPI_COMM_WORLD);
/*MPI_Bcast(recvbuf, count, datatype, root, comm);*/
	s = malloc(sizeof(double)*count);
	s_d = malloc(sizeof(double)*count);
	r = malloc(sizeof(double)*count);
	for(i=count-1 ; i>=0 ; i=i-1) {
	    r[i] = pop(&reduce_stack);
	}
	for(i=count-1 ; i>=0 ; i=i-1) {
	    s[i] = pop(&reduce_stack);
	}
	for(i=0; i<count; i=i+1) {
	    s_d[i] = recvbuf[i];
	}
	if(op == MPI_PROD) {
#ifdef DEBUG
	    printf("AMPI_Allreduce_b MPI_PROD\n");
#endif
	    for(i=0 ; i<count ; i=i+1) {
		if(r[i]*s_d[i] == 0.0){
		    sendbuf[i]=0;
		    printf("--------------------------------------\n");
		    printf("Arithmetic Exeption in adjoint reduce!\n");
		    printf("Result set to 0. Hoping for the best.\n");
		    printf("--------------------------------------\n");
		}
		else {
		    sendbuf[i] = (r[i]/s[i])*s_d[i];
		}
	    }
	}
	if(op == MPI_SUM) {
#ifdef DEBUG
	    printf("AMPI_Allreduce_b MPI_SUM\n");
#endif
	    for(i=0 ; i<count ; i=i+1) {
		sendbuf[i] = s_d[i];
	    }
	}
	free(s);
	free(s_d);
	free(r);
    }
    if(op == MPI_MIN || op == MPI_MAX) {
	minmaxbuf_tmp = malloc(sizeof(double)*count);
#ifdef DEBUG
	printf("AMPI_Allreduce_b MPI_MIN or MPI_MAX\n");
#endif
	MPI_Allreduce(recvbuf, minmaxbuf_tmp, count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	for(i=count-1;i>=0;i--) {
	    idx=(int)pop(&reduce_stack);
	    if(myid==idx) {
		sendbuf[i]=minmaxbuf_tmp[i];
	    }
	    else {
		sendbuf[i]=0;
	    }
            /*MPI_Bcast(&tmp, 1, datatype, idx,
	     * MPI_COMM_WORLD);*/
	    /*sendbuf[i]=tmp;*/
	    
	}
	free(minmaxbuf_tmp);
    }
#ifdef DEBUG
    printf("AMPI_Allreduce_b result: ");
    for(i=0;i<count;i++) { 
	/*sendbuf[i]=1;*/
	/*recvbuf[i]=1;*/
	printf("%f ",sendbuf[i]);
	printf("%f ",recvbuf[i]);
    }
    printf("\n");
#endif
    return MPI_SUCCESS;
}

int AMPI_Sendrecv_replace_f(double *buf, int count, MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag, MPI_Comm comm, MPI_Status *status) {
	return MPI_Sendrecv_replace(buf, count, datatype, dest, sendtag, source, recvtag, comm, status);
}

int AMPI_Sendrecv_replace_b(double *buf, int count, MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag, MPI_Comm comm, MPI_Status *status) {
    /*switch dest and source*/
	return MPI_Sendrecv_replace(buf, count, datatype, source, sendtag, dest, recvtag, comm, status);
}

int AMPI_Gather_f(void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm) {
    return MPI_Gather(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm);
}

int AMPI_Gather_b(void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm) {
    return MPI_Scatter(recvbuf, recvcnt, recvtype, sendbuf, sendcnt, sendtype, root, comm);
}

int AMPI_Scatter_f(void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm) {
    return MPI_Scatter(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm);
}

int AMPI_Scatter_b(void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm) {
    return MPI_Gather(recvbuf, recvcnt, recvtype, sendbuf, sendcnt, sendtype, root, comm);
}

