#include "dco_tape.hpp"

// definitions of static elements of class tape_entry:
int                       tape_entry::vac = 0;
vector<int>               tape_entry::indeps;
vector<int>               tape_entry::deps;
bool TAPE_MODE = false;
bool DERIV_MODE = false;
//int TAPE_SIZE = 100000;
 tape_entry* tape = new tape_entry [TAPE_SIZE];
//tape_entry* tape=0;

#if (DEBUG_ACTIVE == 0)
 int a_ctr=0;
 int ai_ctr=0;
 int al_ctr=0;
 int ad_ctr=0;
 int ald_ctr=0;
 int aa_ctr=0;
 int aeqa_ctr=0;
 int aeqp_ctr=0;
 int nsign_ctr=0;
//  int aeqp_ctr=0;
 int ainc_ctr=0;
 int adec_ctr=0;
 int atimesa_ctr=0;
 int ptimesa_ctr=0;
 int atimesp_ctr=0;
 int adiva_ctr=0;
 int pdiva_ctr=0;
//  int adivp_ctr=0;
 int aplusa_ctr=0;
 int pplusa_ctr=0;
 int aplusp_ctr=0;
 int aminusa_ctr=0;
 int pminusa_ctr=0;
 int aminusp_ctr=0;
 int atimeseqa_ctr=0;
 int atimeseqp_ctr=0;
 int adiveqa_ctr=0;
 int adiveqp_ctr=0;
 int apluseqa_ctr=0;
 int apluseqp_ctr=0;
 int aminuseqa_ctr=0;
 int aminuseqp_ctr=0;
 int sin_ctr=0;
 int asin_ctr=0;
 int sinh_ctr=0;
 int asinh_ctr=0;
 int cos_ctr=0;
 int acos_ctr=0;
 int cosh_ctr=0;
 int acosh_ctr=0;
 int tan_ctr=0;
 int atan_ctr=0;
 int tanh_ctr=0;
 int atanh_ctr=0;
 int aatan2a_ctr=0;
 int aatan2p_ctr=0;
 int patan2a_ctr=0;
 int afmaxa_ctr=0;
 int afmaxp_ctr=0;
 int pfmaxa_ctr=0;
 int amaxa_ctr=0;
 int amaxp_ctr=0;
 int pmaxa_ctr=0;
 int afmina_ctr=0;
 int afminp_ctr=0;
 int pfmina_ctr=0;
 int amina_ctr=0;
 int aminp_ctr=0;
 int pmina_ctr=0;
 int exp_ctr=0;
 int apowa_ctr=0;
 int apowp_ctr=0;
 int ppowa_ctr=0;
 int sqrt_ctr=0;
 int ahypota_ctr=0;
 int ahypotp_ctr=0;
 int phypota_ctr=0;
 int log_ctr=0;
 int log10_ctr=0;
 int log1p_ctr=0;
 int fabs_ctr=0;
 int abs_ctr=0;
 int aldexpa_ctr=0;
 int aldexpp_ctr=0;
 int pldexpa_ctr=0;
 int afrexpi_ctr=0;
 int floor_ctr=0;
 int ceil_ctr=0;
 int setv_ctr=0;
#endif

// Adjoint MPI AMPI

void ampi_get_val(void *buf, int *i, double *x) {
    active *tmp;
    active *tmp_x;
    tmp=(active*) buf;
    tmp_x = &tmp[*i];
    *x = tmp_x->v;
}
void ampi_set_val(void*x, int *i, double *v) {
    active *tmp;
    active *tmp_x;
    tmp = (active*) x;
    tmp_x = &tmp[*i];
    tape[tmp_x->va].v = *v;
    tmp_x->v=*v;
    //printf("SETVAL: %f, %d, %f\n", tmp->v, tmp->va, *v);
}
void ampi_get_adj(INT64 *idx, double *x) {
    *x = tape[*idx].d;
}
void ampi_set_adj(INT64 *idx, double *x) {
    tape[*idx].d = *x;
}
void ampi_get_idx(void *x, int *i, INT64 *idx) {
    active *tmp;
    active *tmp_x;
    tmp = (active*) x;
    tmp_x = &tmp[*i];
    *idx = tmp_x->va;
    //printf("tmp_x.v, tmp_x.va: %f, %d\n", tmp_x->v, tmp_x->va);
}

void ampi_create_tape_entry(int *i) {
    tape[tape_entry::vac].oc = AMPI;
    tape[tape_entry::vac].arg1 = *i;
    tape[tape_entry::vac].arg2 = 0;
    tape[tape_entry::vac].v = 0;
    tape[tape_entry::vac].d = 0;
    tape_entry::vac++;
}

void ampi_create_dummies(void *buf, int *size) {
    active *buf_tmp = (active *) buf;
    for(int i=0;i<*size;i++) {
	tape[tape_entry::vac].oc = ASG;
	tape[tape_entry::vac].arg1 = buf_tmp[i].va;
	tape[tape_entry::vac].arg2 = 0;
	tape[tape_entry::vac].v = buf_tmp[i].v;
	tape[tape_entry::vac].d = 0;
	buf_tmp[i].va=tape_entry::vac;
	tape_entry::vac++;
    }
}

//void ampi_get_elem(void* buf, int i, void* x) {
    //active *tmp;
    //active *tmp_x;
    //tmp = (active*) buf;
    //tmp_x = &tmp[i];
    //x = tmp_x;
//}
    
//int AMPI_Send(active *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {

    ////double * tmp = new double[count];

    ////// create dummy of each element

    ////for(int i = 0 ; i < count ; i++) {
	////tape[tape_entry::vac].oc = MPI_DUMMY;
	////tape[tape_entry::vac+i].v = buf[i].v;
	////tmp[i] = buf[i].v;
	////tape[tape_entry::vac+i].d = 0;
	////tape[tape_entry::vac+i].arg1 = buf[i].va;
	////tape[tape_entry::vac].arg2=-1;
	////buf[i].va = tape_entry::vac+i;
    ////}

    ////// create actual MPI entry

    ////tape[tape_entry::vac+count].oc = SEND;
    ////tape[tape_entry::vac+count].arg1 = count;
    ////tape[tape_entry::vac+count].arg2 = dest;
    ////tape[tape_entry::vac+count].v = 0;
    ////tape[tape_entry::vac+count].d = 0;

    ////int temp = AMPI_Send_f(tmp, count, datatype, dest, tag, comm);

    ////tape_entry::vac+=count+1;
    ////delete [] tmp;
    ////return temp;
////}

//int AMPI_Recv(active *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Status *status) {

    ////double * tmp = new double[count];

    ////int temp = AMPI_Recv_f(tmp, count, datatype, dest, tag, comm, status);

    ////tape[tape_entry::vac].oc = RECV;
    ////tape[tape_entry::vac].arg1 = count;
    ////tape[tape_entry::vac].arg2 = dest;
    ////tape[tape_entry::vac].v = 0;
    ////tape[tape_entry::vac].d = 0;

    ////for(int i = 0 ; i < count ; i++) {
	////tape[tape_entry::vac+i+1].oc = MPI_DUMMY;
	////buf[i].v = tmp[i];
	////tape[tape_entry::vac+i+1].v = tmp[i];
	////tape[tape_entry::vac+i+1].d = 0;
	////tape[tape_entry::vac+i+1].arg1 = buf[i].va;
	////tape[tape_entry::vac].arg2=-1;
	////buf[i].va = tape_entry::vac+i+1;
    ////} 
    ////tape_entry::vac+=count+1;
    ////delete [] tmp;
    ////return temp;
////}

//int AMPI_Isend(active *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, AMPI_dco_Request *request) {
    ////request->buf = buf;
    ////double * tmp = new double[count];

    ////// create dummy of each element

    ////for(int i = 0 ; i < count ; i++) {
	////tape[tape_entry::vac+i].oc = MPI_DUMMY;
	////tape[tape_entry::vac+i].v = buf[i].v;
	////tmp[i] = buf[i].v;
	////tape[tape_entry::vac+i].d = 0;
	////tape[tape_entry::vac+i].arg1 = buf[i].va;
	////tape[tape_entry::vac+i].arg2 = -1;
	////buf[i].va = tape_entry::vac+i;
    ////}
    ////tape[tape_entry::vac+count].oc = ISEND;
    ////tape[tape_entry::vac+count].arg1 = count;
    ////tape[tape_entry::vac+count].arg2 = dest;
    ////tape[tape_entry::vac+count].v = 0;
    ////tape[tape_entry::vac+count].d = 0;

    ////if(request->r.aw) { // if antiwait, memory was already allocated
	////tape[tape_entry::vac+count].request = new AMPI_dco_Request;
	////tape[request->va].request = tape[tape_entry::vac+count].request;
	////tape[request->va].arg1 = tape_entry::vac+count;
    ////}
    ////else {
	////tape[tape_entry::vac+count].request = new AMPI_dco_Request;
    ////}
    ////request->va = tape_entry::vac+count; // point current request index to this tape entry
    ////request->r.v = tmp;
    ////int temp = AMPI_Isend_f(tmp, count, datatype, dest, tag, comm, &request->r);
    ////tape_entry::vac+=count+1;
    ////return temp;
////}

//int AMPI_Irecv(active *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, AMPI_dco_Request *request) {
    ////request->buf = buf;
    ////double * tmp = new double[count];
    ////int temp = AMPI_Irecv_f(tmp, count, datatype, dest, tag, comm, &request->r);

    ////tape[tape_entry::vac].arg1 = count;
    ////tape[tape_entry::vac].arg2 = dest;
    ////tape[tape_entry::vac].oc = IRECV;
    ////tape[tape_entry::vac].v = 0;
    ////tape[tape_entry::vac].d = 0;

    ////for(int i = 0 ; i < count ; i++) {
	////tape[tape_entry::vac+i+1].arg1 = buf[i].va;
	////tape[tape_entry::vac+i+1].arg2 = -1;
	////tape[tape_entry::vac+i+1].oc = MPI_DUMMY;
	////tape[tape_entry::vac+i+1].v = 0;
	////tape[tape_entry::vac+i+1].d = 0;
	////buf[i].va = tape_entry::vac+i+1;


    ////} 
    ////if(request->r.aw) { // if antiwait, memory was already allocated
	////tape[tape_entry::vac].request = new AMPI_dco_Request;
	////tape[request->va].request = tape[tape_entry::vac].request;
	////tape[request->va].arg1 = tape_entry::vac;
    ////}
    ////else {
	////tape[tape_entry::vac].request =  new AMPI_dco_Request;
    ////}
    //////request->r.v = tmp;
    ////request->va = tape_entry::vac; // point current request index to this tape entry
    ////tape_entry::vac+=count+1;
    ////return temp;
////}

//int AMPI_Wait(AMPI_dco_Request *request, MPI_Status *status) {
    ////double * tmp = (double*) request->r.v;
    ////tape[tape_entry::vac].arg1 = request->va; // get the corresponding isend or irecv tape entry
    ////tape[tape_entry::vac].arg2=-1;
    ////tape[tape_entry::vac].oc = WAIT;
    ////tape[tape_entry::vac].v = 0;
    ////tape[tape_entry::vac].d = 0;
    ////AMPI_Wait_f(&request->r, status);
    ////*tape[request->va].request = *request;			// finally copy the request to the tape
    ////tape[tape_entry::vac].request = tape[request->va].request;
    ////if(request->r.oc == IR) {
	////for(int i = 0 ; i < tape[request->va].arg1 ; i++) {
	    ////request->buf[i].v = tmp[i];
	    ////tape[request->va+i+1].v = tmp[i];
	////}	 
    ////}
    ////tape[tape_entry::vac].request->r.a = &tape[request->va].d;
    ////tape[tape_entry::vac].request->r.size = tape[request->va].arg1;
    ////tape_entry::vac++;
    ////delete [] tmp;
    ////return 1;
////}

//int AMPI_Waitall(int count, AMPI_dco_Request *request, MPI_Status *status) {
    ////double ** tmp = new double*[count];
    ////AMPI_Request * tmp_r = new AMPI_Request[count];
    ////for(int i = 0 ; i < count ; i++) {
	////tmp[i] = (double*) request[i].r.v;
	////tape[tape_entry::vac+i].arg1 = request[i].va;
	////tape[tape_entry::vac+i].arg2 = -1;
	////tape[tape_entry::vac+i].oc = MPI_DUMMY;
	////tape[tape_entry::vac+i].v = 0;
	////tape[tape_entry::vac+i].d = 0;
    ////}
    ////tape[tape_entry::vac+count].arg1 = count;
    ////tape[tape_entry::vac+count].arg2 = -1;
    ////tape[tape_entry::vac+count].oc = WAITALL;
    ////tape[tape_entry::vac+count].v = 0;
    ////tape[tape_entry::vac+count].d = 0;
    ////for(int i = 0 ; i < count ; i++) {
	////tmp_r[i] = request[i].r;
    ////}	
    ////AMPI_Waitall_f(count, tmp_r, status);
    ////for(int i = 0 ; i < count ; i++) {
	////request[i].r = tmp_r[i];
    ////}	

    ////for(int i = 0 ; i < count ; i++) {
	////*tape[request[i].va].request = request[i];
	////tape[tape_entry::vac+i].request = tape[request[i].va].request;
    ////}
    ////for(int j = 0 ; j < count ; j++) {
	////if(request[j].r.oc == IR) {
	    ////for(int i = 0 ; i < tape[request[j].va].arg1 ; i++) {
		////request[j].buf[i].v = tmp[j][i];
		////tape[request[j].va+i+1].v = tmp[j][i];
	    ////}	 
	////}
	////tape[tape_entry::vac+j].request->r.a = &tape[request[j].va].d;
	////tape[tape_entry::vac+j].request->r.size = tape[request[j].va].arg1;
	////tape[tape_entry::vac+j].arg2 = tape[request[j].va].arg1;
    ////}
    ////tape_entry::vac+=count+1;
    ////delete [] tmp_r;
    ////return 1;
////}

//int AMPI_Awaitall(int count, AMPI_dco_Request *request, MPI_Status *status) {
    ////for(int i = 0 ; i < count ; i++) {
	////tape[tape_entry::vac+i].arg1 = -1;
	////tape[tape_entry::vac+i].arg2 = -1;
	////tape[tape_entry::vac+i].oc = MPI_ADUMMY;
	////tape[tape_entry::vac+i].v = 0;
	////tape[tape_entry::vac+i].d = 0;
	////request[i].va = tape_entry::vac+i;
    ////}
    ////tape[tape_entry::vac+count].arg1 = count;
    ////tape[tape_entry::vac+count].arg2 = -1;
    ////tape[tape_entry::vac+count].oc = AWAITALL;
    ////tape[tape_entry::vac+count].v = 0;
    ////tape[tape_entry::vac+count].d = 0;
    ////tape_entry::vac+=count+1;
    ////return 1;
////}

//int AMPI_Bcast(active *buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm) {

    ////double * tmp = new double[count];
    ////// create dummy of each element

    ////for(int i = 0 ; i < count ; i++) {
	////tape[tape_entry::vac].oc = MPI_DUMMY;
	////tape[tape_entry::vac+i].v = buf[i].v;
	////tmp[i] = buf[i].v;
	////tape[tape_entry::vac+i].d = 0;
	////tape[tape_entry::vac+i].arg1 = buf[i].va;
	////buf[i].va = tape_entry::vac+i;
    ////}

    ////// create actual MPI entry

    ////tape[tape_entry::vac+count].oc = BCAST;
    ////tape[tape_entry::vac+count].arg1 = count;
    ////tape[tape_entry::vac+count].arg2 = root;
    ////tape[tape_entry::vac+count].v = 0;
    ////tape[tape_entry::vac+count].d = 0;

    ////int temp = AMPI_Bcast_f(tmp, count, datatype, root, comm);

    ////tape_entry::vac+=count+1;
    ////delete [] tmp;
    ////return temp;
////}

//int AMPI_Reduce(active *sendbuf, active *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {

    ////double * tmp = new double[count];

    ////// sendbuf dummies

    ////for(int i = 0 ; i < count ; i++) {
	////tape[tape_entry::vac+i].oc = MPI_DUMMY;
	////tape[tape_entry::vac+i].v = sendbuf[i].v;
	////tmp[i] = sendbuf[i].v;
	////tape[tape_entry::vac+i].d = 0;
	////tape[tape_entry::vac+i].arg1 = sendbuf[i].va;
	////sendbuf[i].va = tape_entry::vac+i;
    ////}

    ////// actual reduce entry

    ////tape[tape_entry::vac+count].oc = REDUCE;
    ////tape[tape_entry::vac+count].arg1 = count;
    ////tape[tape_entry::vac+count].arg2 = root;
    ////tape[tape_entry::vac+count].v = 0;
    ////tape[tape_entry::vac+count].d = 0;

    ////// recvbuf entry

    ////tape[tape_entry::vac+count+1].oc = MPI_DUMMY;
    ////tape[tape_entry::vac+count+1].d = 0;
    ////tape[tape_entry::vac+count+1].v = 0;
    ////tape[tape_entry::vac+count+1].arg1 = -1;
    ////recvbuf->va = tape_entry::vac+count+1;


    ////if(op == MPI_SUM)
	////tape[tape_entry::vac+count+1].arg2 = ADD;
    ////if(op == MPI_PROD)
	////tape[tape_entry::vac+count+1].arg2 = MUL;
    ////AMPI_Reduce_f(tmp, &tape[tape_entry::vac+count+1].v, count, datatype, op, root, comm);
    ////recvbuf->v = tape[tape_entry::vac+count+1].v;
    ////tape_entry::vac+=count+2;
    ////delete [] tmp;
    ////return 1;
////}

void forward_tape_interpreter(int i)
{
  switch (tape[i].oc) {
    case INDEP : {} break;
    case DEP : {
      #if (DEBUG_ACTIVE == 0)
      cout << "WARN > tape_entry " << i << "is marked as dependent" << endl;
      assert(false);
      #endif
    } break;
    case CONST : {
      tape[i].d=0;
    } break;
    case COPY : {
      tape[i].d = tape[tape[i].arg1].d;
      // separation of passive rhs
//       if (tape[i].arg1 >= 0)
//         tape[i].d = tape[tape[i].arg1].d;
//       else tape[i].d = 0;
    } break;
    case ASG : {
      tape[i].d = tape[tape[i].arg1].d; 
      // separation of passive rhs
/*      if (tape[i].arg1 >= 0)
        tape[i].d = tape[tape[i].arg1].d;
      else tape[i].d = 0;*/
    } break;
    case PLUSEQ : {
      tape[i].d=tape[tape[i].arg1].d + tape[tape[i].arg2].d;
    } break;
    case MINUSEQ : {
      tape[i].d=tape[tape[i].arg1].d - tape[tape[i].arg2].d;
    } break;
    case MULTEQ : {
      tape[i].d= tape[tape[i].arg2].v*tape[tape[i].arg1].d + tape[tape[i].arg1].v*tape[tape[i].arg2].d;
    } break;
    case DIVEQ : {
      if (tape[i].arg2 == UNDEF || tape[tape[i].arg2].v == 0)  tape[i].d = 0;
      else  tape[i].d=(tape[tape[i].arg1].d/tape[tape[i].arg2].v) - 
                      (tape[tape[i].arg2].d*(tape[tape[i].arg1].v/tape[tape[i].arg2].v*tape[tape[i].arg2].v));
    } break;
    case ADD : {
      if (tape[i].arg1 != UNDEF && tape[i].arg2 != UNDEF ) tape[i].d=tape[tape[i].arg1].d + tape[tape[i].arg2].d;
      else tape[i].d = (tape[i].arg2== UNDEF) ? tape[tape[i].arg1].d : tape[tape[i].arg2].d;
    } break;
    case SUB : {
      if (tape[i].arg1 != UNDEF && tape[i].arg2 != UNDEF ) tape[i].d=tape[tape[i].arg1].d - tape[tape[i].arg2].d;
      else tape[i].d = (tape[i].arg2== UNDEF) ? tape[tape[i].arg1].d : -tape[tape[i].arg2].d;
    } break;
    case NEGSIGN : {
      tape[i].d=-tape[tape[i].arg1].d;
    } break;
    case MUL :  {
      if(tape[i].arg1 != UNDEF && tape[i].arg2 != UNDEF ) 
        tape[i].d=tape[tape[i].arg2].v*tape[tape[i].arg1].d + tape[tape[i].arg1].v*tape[tape[i].arg2].d;
      else if (tape[i].arg1 != UNDEF && tape[tape[i].arg1].v) 
        tape[i].d = (tape[i].v/tape[tape[i].arg1].v)*tape[tape[i].arg1].d;
      else if (tape[i].arg2 != UNDEF && tape[tape[i].arg2].v) 
        tape[i].d = (tape[i].v/tape[tape[i].arg2].v)*tape[tape[i].arg2].d;
      else tape[i].d = 0;
    } break;
    case DIV : {
      if (tape[i].arg2 == UNDEF || tape[tape[i].arg2].v == 0) tape[i].d = 0;
      else if (tape[i].arg1 != UNDEF && tape[i].arg2 != UNDEF )
        tape[i].d=tape[tape[i].arg1].d/tape[tape[i].arg2].v - 
        tape[tape[i].arg2].d * (tape[tape[i].arg1].v/(tape[tape[i].arg2].v * tape[tape[i].arg2].v));
//       else if (tape[i].arg1 != UNDEF && tape[tape[i].arg1].v) 
//         tape[i].d = (tape[i].v/tape[tape[i].arg1].v)*tape[tape[i].arg1].d;
      else if (tape[i].arg2 != UNDEF && tape[tape[i].arg2].v) 
        tape[i].d = (-tape[i].v/tape[tape[i].arg2].v)*tape[tape[i].arg2].d;
      else tape[i].d = 0;
    } break;
    case SIN : {
      tape[i].d=cos(tape[tape[i].arg1].v)*tape[tape[i].arg1].d;
    } break;
    case ASIN : {
      double tmp = 1.0-tape[tape[i].arg1].v*tape[tape[i].arg1].v;
      if (tmp <= 0) tape[i].d=0;
      else tape[i].d = tape[tape[i].arg1].d/sqrt(tmp);
    } break;
    case SINH : {
      tape[i].d=cosh(tape[tape[i].arg1].v)*tape[tape[i].arg1].d;
    } break;
    case ASINH : {
      double tmp = tape[tape[i].arg1].v * tape[tape[i].arg1].v+1.0;
      if (tmp <= 0)  tape[i].d=0;
      else tape[i].d=tape[tape[i].arg1].d/sqrt(tmp);
    } break;
    case COS : {
      tape[i].d=-sin(tape[tape[i].arg1].v)*tape[tape[i].arg1].d;
    } break;
    case ACOS : {
      double tmp = 1.0-tape[tape[i].arg1].v*tape[tape[i].arg1].v;
      if (tmp <= 0)  tape[i].d=0;
      else tape[i].d=-tape[tape[i].arg1].d/sqrt(tmp);
    } break;
    case COSH : {
      tape[i].d=sinh(tape[tape[i].arg1].v)*tape[tape[i].arg1].d;
    } break;
    case ACOSH : {
      double tmp = tape[tape[i].arg1].v*tape[tape[i].arg1].v-1.0;
      if (tmp <= 0)  tape[i].d=0;
      else tape[i].d=tape[tape[i].arg1].d/sqrt(tmp);
    } break;
    case TAN : {
      if (!cos(tape[tape[i].arg1].v))  tape[i].d=0;
      else tape[i].d=tape[tape[i].arg1].d/(cos(tape[tape[i].arg1].v)*cos(tape[tape[i].arg1].v));
    } break;
    case ATAN : {
      double tmp = 1.0+tape[tape[i].arg1].v*tape[tape[i].arg1].v;
      if (tmp == 0)  tape[i].d=0;
      else tape[i].d=tape[tape[i].arg1].d/tmp;
    } break;
    case TANH : {
      double tmp = tape[tape[i].arg1].v*tape[tape[i].arg1].v;
      if (tmp == 0)  tape[i].d=0;
      else tape[i].d=tape[tape[i].arg1].d/tmp;
    } break;
    case ATANH : {
      double tmp = 1.0-tape[tape[i].arg1].v*tape[tape[i].arg1].v;
      if (tmp == 0)  tape[i].d=0;
      else tape[i].d=tape[tape[i].arg1].d/tmp;
    } break;
    case EXP : {
      tape[i].d=tape[i].v*tape[tape[i].arg1].d;
    } break;
    case LOG : {
      if (tape[tape[i].arg1].v <= 0)  tape[i].d=0;
      else tape[i].d = tape[tape[i].arg1].d/tape[tape[i].arg1].v;
    } break;
    case LOG1P : {
      double tmp = tape[tape[i].arg1].v+1.0;
      if (tmp == 0)  tape[i].d=0;
      else tape[i].d=tape[tape[i].arg1].d/tmp;
    } break;
    case LOG10 : {
      double tmp = log(10.0)*tape[tape[i].arg1].v;
      if (tmp == 0)  tape[i].d=0;
      else tape[i].d=tape[tape[i].arg1].d/tmp;
    } break;
    case POW : {
      if (tape[i].arg1 != UNDEF && tape[i].arg2 != UNDEF)
        tape[i].d=tape[tape[i].arg2].v*pow(tape[tape[i].arg1].v,tape[tape[i].arg2].v-1)*tape[tape[i].arg1].d +
                  log(tape[tape[i].arg1].v)*tape[i].v*tape[tape[i].arg2].d;
      else if (tape[i].arg1 != UNDEF && log(tape[tape[i].arg1].v)) {
        double a = log(tape[i].v)/log(tape[tape[i].arg1].v);
        tape[i].d=a*pow(tape[tape[i].arg1].v, a-1)*tape[tape[i].arg1].d;
      }
      else if (tape[i].arg2 != UNDEF)  tape[i].d=tape[i].v*log(tape[tape[i].arg2].v)*tape[tape[i].arg2].d;
      else tape[i].d=0;
    } break;
    case SQRT : {
      double tmp = 2.0*tape[i].v;
      if (tmp == 0)  tape[i].d=0;
      else tape[i].d=tape[tape[i].arg1].d/tmp;
    } break;
    case FABS : {
      if (tape[tape[i].arg1].v > 0)  tape[i].d=tape[tape[i].arg1].d;
      else if (tape[tape[i].arg1].v < 0)  tape[i].d=-tape[tape[i].arg1].d;
      else tape[i].d=0;
    } break;
    case MAX : {
      if (tape[i].arg1 != UNDEF && tape[i].arg2 != UNDEF) {
        if (tape[tape[i].arg1].v >= tape[tape[i].arg2].v)  tape[i].d=tape[tape[i].arg1].d;
        else tape[i].d=tape[tape[i].arg2].d;
      }
      else if (tape[i].arg1 != UNDEF && tape[i].v == tape[tape[i].arg1].v)  tape[i].d=tape[tape[i].arg1].d;
      else if (tape[i].arg2 != UNDEF && tape[i].v == tape[tape[i].arg2].v)  tape[i].d=tape[tape[i].arg2].d;
      else tape[i].d=0;
    } break;
    case MIN : {
      if (tape[i].arg1 != UNDEF && tape[i].arg2 != UNDEF) {
        if (tape[tape[i].arg1].v < tape[tape[i].arg2].v)  tape[i].d=tape[tape[i].arg1].d;
        else  tape[i].d=tape[tape[i].arg2].d;
      }
      else if (tape[i].arg1 != UNDEF && tape[i].v == tape[tape[i].arg1].v)  tape[i].d=tape[tape[i].arg1].d;
      else if (tape[i].arg2 != UNDEF && tape[i].v == tape[tape[i].arg2].v)  tape[i].d=tape[tape[i].arg2].d;
      else  tape[i].d=0;
    } break;
    case HYPOT : {
      if (tape[i].v == 0) {
        tape[i].d=0;
      } 
      else if (tape[i].arg1 != UNDEF && tape[i].arg2 != UNDEF)
        tape[i].d=tape[tape[i].arg1].v/tape[i].v*tape[tape[i].arg1].d +  
                  tape[tape[i].arg2].v/tape[i].v*tape[tape[i].arg2].d;
      else if (tape[i].arg1 != UNDEF)  tape[i].d=tape[tape[i].arg1].v/tape[i].v*tape[tape[i].arg1].d;
      else if (tape[i].arg2 != UNDEF)  tape[i].d=tape[tape[i].arg2].v/tape[i].v*tape[tape[i].arg2].d;
      else tape[i].d=0;
    } break;
    case FLOOR : {
      tape[i].d=0;
    } break;
    case CEIL : {
      tape[i].d=0;
    } break;
    default : {
    }
  }
}
/**
 * reverse interpretation of tape with index idx
 */
void reverse_tape_interpreter(int i) {
  switch (tape[i].oc) {
    case COPY : {
      if (tape[i].arg1 != UNDEF)  tape[tape[i].arg1].d+=tape[i].d;
      tape[i].d=0.;
    } break;
    case ASG : {
      if (tape[i].arg1 != UNDEF)  tape[tape[i].arg1].d+=tape[i].d;
      tape[i].d=0.;
    } break;
    case PLUSEQ : {
      tape[tape[i].arg1].d=tape[i].d;
      tape[tape[i].arg2].d+=tape[i].d;
      tape[i].d=0.;
    } break;
    case MINUSEQ : {
      tape[tape[i].arg1].d=tape[i].d;
      tape[tape[i].arg2].d-=tape[i].d;
      tape[i].d=0.;
    } break;
    case MULTEQ : {
      tape[tape[i].arg1].d=tape[tape[i].arg2].v*tape[i].d;
      tape[tape[i].arg2].d+=tape[tape[i].arg1].v*tape[i].d;
      tape[i].d=0.;
    } break;
    case DIVEQ : {
      #if (DEBUG_ACTIVE == 4)
      if (tape[i].arg2 == UNDEF || tape[tape[i].arg2].v == 0) {print_tape(i,"WARN");}
      #endif
      tape[tape[i].arg1].d=tape[i].d/tape[tape[i].arg2].v;
      tape[tape[i].arg2].d-=tape[i].d*(tape[tape[i].arg1].v/(tape[tape[i].arg2].v*tape[tape[i].arg2].v));
      tape[i].d=0.;
    } break;
    case ADD : {
      if (tape[i].arg1 != UNDEF) tape[tape[i].arg1].d+=tape[i].d;
      if (tape[i].arg2 != UNDEF) tape[tape[i].arg2].d+=tape[i].d;
      tape[i].d=0.;
    } break;
    case SUB : {
      if(tape[i].arg1 != UNDEF) tape[tape[i].arg1].d+=tape[i].d;
      if(tape[i].arg2 != UNDEF) tape[tape[i].arg2].d-=tape[i].d;
      tape[i].d=0.;
    } break;
    case NEGSIGN : {
      tape[tape[i].arg1].d-=tape[i].d;
      tape[i].d=0.;
      break;
    }
    case MUL : {
      if(tape[i].arg1 != UNDEF && tape[i].arg2 != UNDEF ) {
        tape[tape[i].arg1].d+=tape[tape[i].arg2].v*tape[i].d;
        tape[tape[i].arg2].d+=tape[tape[i].arg1].v*tape[i].d;
        tape[i].d=0.;
      }
      else if(tape[i].arg1 != UNDEF && tape[tape[i].arg1].v) tape[tape[i].arg1].d+=(tape[i].v/tape[tape[i].arg1].v)*tape[i].d;
      else if(tape[i].arg2 != UNDEF && tape[tape[i].arg2].v) tape[tape[i].arg2].d+=(tape[i].v/tape[tape[i].arg2].v)*tape[i].d;
      tape[i].d=0.;
    } break;
    case DIV : {
      #if (DEBUG_ACTIVE == 4)
      if (tape[i].arg2 == UNDEF || tape[tape[i].arg2].v == 0) {print_tape(i,"WARN");}
      #endif
      if(tape[i].arg1 != UNDEF && tape[i].arg2 != UNDEF) {
        tape[tape[i].arg1].d+= tape[i].d/tape[tape[i].arg2].v;
        tape[tape[i].arg2].d-=tape[i].d*(tape[tape[i].arg1].v/(tape[tape[i].arg2].v * tape[tape[i].arg2].v));
      }
      else if(tape[i].arg2 != UNDEF  && tape[tape[i].arg2].v)
        tape[tape[i].arg2].d-= tape[i].d*(tape[i].v/tape[tape[i].arg2].v);
      tape[i].d=0; 
    }  break;
    case SIN : {
      tape[tape[i].arg1].d+=cos(tape[tape[i].arg1].v)*tape[i].d;
      tape[i].d=0.;
    }  break;
    case ASIN : {
      double tmp = 1.0-tape[tape[i].arg1].v*tape[tape[i].arg1].v;
      #if (DEBUG_ACTIVE == 4)
      if (tmp <= 0) {print_tape(i,"WARN");}
      #endif
      tape[tape[i].arg1].d+=tape[i].d/sqrt(tmp);
      tape[i].d=0;
    }  break;
    case SINH : {
      tape[tape[i].arg1].d+=cosh(tape[tape[i].arg1].v)*tape[i].d;
      tape[i].d=0;
    }  break;
    case ASINH : {
      double tmp = tape[tape[i].arg1].v * tape[tape[i].arg1].v+1.0;
      if (tmp <= 0) {
        #if (DEBUG_ACTIVE == 4)
        print_tape(i,"WARN");
        #endif
      }
      else tape[tape[i].arg1].d+=tape[i].d/sqrt(tmp);
      tape[i].d=0;
    }  break;
    case COS : {
      tape[tape[i].arg1].d-=sin(tape[tape[i].arg1].v)*tape[i].d;
      tape[i].d=0.;
    }  break;
    case ACOS : {
      double tmp = 1.0-tape[tape[i].arg1].v*tape[tape[i].arg1].v;
      #if (DEBUG_ACTIVE == 4)
      if (tmp <= 0) {print_tape(i,"WARN");}
      #endif
      tape[tape[i].arg1].d-=tape[i].d/sqrt(tmp);
      tape[i].d=0;
    }  break;
    case COSH : {
      tape[tape[i].arg1].d+=sinh(tape[tape[i].arg1].v)*tape[i].d;
      tape[i].d=0;
    }  break;
    case ACOSH : {
      double tmp = tape[tape[i].arg1].v*tape[tape[i].arg1].v-1.0;
      #if (DEBUG_ACTIVE == 4)
      if (tmp <= 0) {print_tape(i,"WARN");}
      #endif
      tape[tape[i].arg1].d+=tape[i].d/sqrt(tmp);
      tape[i].d=0;
    }  break;
    case TAN : {
      #if (DEBUG_ACTIVE == 4)
      if (!cos(tape[tape[i].arg1].v)) {print_tape(i,"WARN");}
      #endif
      tape[tape[i].arg1].d+=tape[i].d/(cos(tape[tape[i].arg1].v)*cos(tape[tape[i].arg1].v));
      tape[i].d=0;
    }  break;
    case ATAN : {
      double tmp = 1.0+tape[tape[i].arg1].v*tape[tape[i].arg1].v;
      #if (DEBUG_ACTIVE == 4)
      if (tmp == 0) {print_tape(i,"WARN");}
      #endif
      tape[tape[i].arg1].d+=tape[i].d/tmp;
      tape[i].d=0;
    }  break;
    case TANH : {
      double tmp = tape[tape[i].arg1].v*tape[tape[i].arg1].v;
      #if (DEBUG_ACTIVE == 4)
      if (tmp == 0) {print_tape(i,"WARN");}
      #endif
      tape[tape[i].arg1].d+=tape[i].d/tmp;
      tape[i].d=0;
    }  break;
    case ATANH : {
      double tmp = 1.0-tape[tape[i].arg1].v*tape[tape[i].arg1].v;
      #if (DEBUG_ACTIVE == 4)
      if (tmp == 0) {print_tape(i,"WARN");}
      #endif
      tape[tape[i].arg1].d+=tape[i].d/tmp;
      tape[i].d=0;
    }  break;
    case EXP : {
      tape[tape[i].arg1].d+=tape[i].v*tape[i].d;
      tape[i].d=0;
    }  break;
    case LOG : {
      #if (DEBUG_ACTIVE == 4)
      if (tape[tape[i].arg1].v <= 0) {print_tape(i,"WARN");}
      #endif
      tape[tape[i].arg1].d+=tape[i].d/tape[tape[i].arg1].v;
      tape[i].d=0;
    }  break;
    case LOG1P : {
      double tmp = tape[tape[i].arg1].v+1.0;
      #if (DEBUG_ACTIVE == 4)
      if (tmp == 0) {print_tape(i,"WARN");}
      #endif
      tape[tape[i].arg1].d+=tape[i].d/tmp;
      tape[i].d=0;
    }  break;
    case LOG10 : {
      double tmp = log(10.0)*tape[tape[i].arg1].v;
      #if (DEBUG_ACTIVE == 4)
      if (tmp == 0) {print_tape(i,"WARN");}
      #endif
      tape[tape[i].arg1].d+=tape[i].d/tmp;
      tape[i].d=0;
    }  break;
    case POW : {
      if(tape[i].arg1 != UNDEF && tape[i].arg2 != UNDEF ) {
        tape[tape[i].arg1].d+=tape[tape[i].arg2].v*pow(tape[tape[i].arg1].v, tape[tape[i].arg2].v-1)*tape[i].d;
        tape[tape[i].arg2].d+=log(tape[tape[i].arg2].v)*tape[i].v*tape[i].d;  
      }
      else if(tape[i].arg1 != UNDEF && log(tape[tape[i].arg1].v)) {
        double a = log(tape[i].v)/log(tape[tape[i].arg1].v); 
        tape[tape[i].arg1].d+=a*pow(tape[tape[i].arg1].v, a-1);
      }
      else if(tape[i].arg2 != UNDEF) tape[tape[i].arg2].d+=tape[i].v*log(tape[tape[i].arg2].v);
      tape[i].d=0;
    }  break;
    case SQRT : {
      double tmp = 2.0*tape[i].v;
      #if (DEBUG_ACTIVE == 4)
      if (tmp == 0) {print_tape(i,"WARN");}
      #endif
      tape[tape[i].arg1].d+=tape[i].d/tmp;
      tape[i].d=0;
    }  break;
    case FABS : {
      if (tape[tape[i].arg1].v > 0)  tape[tape[i].arg1].d+=tape[i].d;
      else if (tape[tape[i].arg1].v < 0)  tape[tape[i].arg1].d-=tape[i].d;
      else  tape[i].d=0;
    }  break;
    case MAX : {
      if (tape[i].arg1 != UNDEF && tape[i].arg2 != UNDEF) {
        if (tape[tape[i].arg1].v >= tape[tape[i].arg2].v)  tape[tape[i].arg1].d+=tape[i].d;
        else  tape[tape[i].arg2].d+=tape[i].d;
      }
      else if (tape[i].arg1 != UNDEF && tape[i].v == tape[tape[i].arg1].v)  tape[tape[i].arg1].d+=tape[i].d;
      else if (tape[i].arg2 != UNDEF && tape[i].v == tape[tape[i].arg2].v)  tape[tape[i].arg2].d+=tape[i].d;
      tape[i].d=0;
    }  break;
    case MIN : {
      if (tape[i].arg1 != UNDEF && tape[i].arg2 != UNDEF) {
        if (tape[tape[i].arg1].v < tape[tape[i].arg2].v)  tape[tape[i].arg1].d+=tape[i].d;
        else  tape[tape[i].arg2].d+=tape[i].d;
      }
      else if (tape[i].arg1 != UNDEF && tape[i].v == tape[tape[i].arg1].v)  tape[tape[i].arg1].d+=tape[i].d;
      else if (tape[i].arg2 != UNDEF && tape[i].v == tape[tape[i].arg2].v)  tape[tape[i].arg2].d+=tape[i].d;
      tape[i].d=0;
    }  break;
    case HYPOT : { 
      #if (DEBUG_ACTIVE == 4)
      if (tape[i].v == 0) {print_tape(i,"WARN");}
      #endif
      if (tape[i].arg1 != UNDEF && tape[i].arg2 != UNDEF) {
        tape[tape[i].arg1].d+=tape[tape[i].arg1].v/tape[i].v*tape[i].d;
        tape[tape[i].arg2].d+=tape[tape[i].arg2].v/tape[i].v*tape[i].d;
      }
      else if (tape[i].arg1 != UNDEF)  tape[tape[i].arg1].d+=tape[tape[i].arg1].v/tape[i].v*tape[i].d;
      else if (tape[i].arg2 != UNDEF)  tape[tape[i].arg2].d+=tape[tape[i].arg2].v/tape[i].v*tape[i].d;
      tape[i].d=0;
    }  break;
    case FLOOR : {
      tape[i].d=0;
      break;
    }
    case CEIL : {
      tape[i].d=0;
      break;
    }
    case AMPI : {
		    //printf("AMPI Entry\n");
		    ampi_interpret_tape();
		}
		/* AMPI */
		//case MPI_DUMMY : {
				     ////if(tape[i].arg1 != -1) // recv buffer of reduce
				          ////tape[tape[i].arg1].d+=tape[i].d;
				     ////break;
				 ////}
		//case SEND : {
				////tmp_d = new double[tape[i].arg1];
				////AMPI_Recv_b(tmp_d, tape[i].arg1, MPI_DOUBLE, tape[i].arg2, 0, comm, &status);
				////for(int j = 0 ; j < tape[i].arg1 ; j++) {
				    ////tape[i - tape[i].arg1 + j].d += tmp_d[j];
				////}
				////delete [] tmp_d;
				////break;
			    ////}
		//case RECV : {
				////tmp_d = new double[tape[i].arg1];
				////for(int j = 0 ; j < tape[i].arg1 ; j++) {
				    ////tmp_d[j] = tape[i+j+1].d;
				////} 
				////AMPI_Send_b(tmp_d, tape[i].arg1, MPI_DOUBLE, tape[i].arg2, 0, comm);
				////delete [] tmp_d;
				////break;
			    ////}
		//case ISEND : {
				 ////tmp_d = new double[tape[i].arg1];
				 ////if(!tape[i].request->r.aw) {
				     ////tmp_d = (double*) tape[i].request->r.a;
				     //////AMPI_Isend_b(&tape[i].request->r);
				     ////AMPI_Isend_b(tmp_d, tape[i].arg1, MPI_DOUBLE, tape[i].arg2, 0, comm,&tape[i].request->r);
				     ////for(int j = 0 ; j < tape[i].arg1 ; j++) {
					 ////tape[i-j-1].d = tmp_d[j];
				     ////}
				     //////				     delete [] tmp_d;
				     //////				     delete tape[i].request;
				 ////}
				 ////else {
//////				     tape[i].request->va = i;
				 ////}
				////delete [] tmp_d;
				 ////break;
			     ////}
		//case IRECV : {
				////tmp_d = new double[tape[i].arg1];
				 ////if(tape[i].request->r.aw) {
				     
				 ////}
				 ////else {
				     ////AMPI_Irecv_b(tmp_d, tape[i].arg1, MPI_DOUBLE, tape[i].arg2, 0, comm,&tape[i].request->r);
				     //////AMPI_Irecv_b(&tape[i].request->r);
				     //////				      delete [] tmp_d;
				     //////				      delete tape[i].request;
				 ////}
				////delete [] tmp_d;
				 ////break;
			     ////}
		//case WAIT : {
				//tmp_d = new double[tape[tape[i].arg1].arg1];
				//if(tape[i].request->r.oc == IR) {
				    //for(int j = 0 ; j < tape[tape[i].arg1].arg1 ; j++) {
					//tmp_d[j] = tape[tape[i].arg1+j+1].d;
				    //}
				//}
				//tape[i].request->r.a = tmp_d;
				//AMPI_Wait_b(&tape[i].request->r,&status);
				//tape[tape[i].arg1].request = tape[i].request;
				//break;
			    //}

		//case WAITALL : {
				   ////count = tape[i].arg1;
				   ////reqs = new AMPI_Request[count];
				   ////tmp_dd = new double*[count];
				   ////for(int k = 0 ; k < count ; k++) {
				       ////tmp_dd[k] = new double[tape[tape[i-count+k].arg1].arg1];
				       ////if(tape[i-count+k].request->r.oc == IR) {
					   ////for(int j = 0 ; j < tape[tape[i-count+k].arg1].arg1 ; j++) {
					       ////tmp_dd[k][j] = tape[tape[i-count+k].arg1+j+1].d;
					   ////}
				       ////}
				       ////reqs[k] = tape[i-count+k].request->r;
				       ////reqs[k].a = tmp_dd[k];
				   ////}
				   ////AMPI_Waitall_b(count , reqs);
				   ////for(int k = 0 ; k < count ; k++) {
				       ////tape[tape[i-count+k].arg1].request->r = reqs[k];
				   ////}
				   ////break;
			       ////}
		//case AWAITALL : {
				   ////count = tape[i].arg1;
				   ////reqs = new AMPI_Request[count];
				   ////for(int k = 0 ; k < count ; k++) {
				       ////reqs[k] = tape[i-count+k].request->r;
				   ////}
				   ////AMPI_Awaitall_b(count, reqs);
				   ////for(int k = 0 ; k < count ; k++) {
				     ////tmp_d = (double*) reqs[k].a;
				     ////if(reqs[k].oc == IS) {
					////for(int j = tape[i-count+k].arg1 ; j < (tape[i-count+k].arg1 + tape[tape[i-count+k].arg1].arg1) ; j++) {
					////tape[tape[j].arg1].d = tmp_d[j];
					////}
				     ////}
				   ////}


				    ////break;

				////}
		//case BCAST : {
				 ////tmp_d = new double[tape[i].arg1];
				 ////AMPI_Bcast_b(tmp_d, tape[i].arg1, MPI_DOUBLE, tape[i].arg2, comm);
				 ////if(*AMPI_myid == tape[i].arg2) { 
				     ////for(int j = 0 ; j < tape[i].arg1 ; j++) {
					 ////tape[i-tape[i].arg1+j].d += tmp_d[j];
				     ////}
				 ////}
				 ////delete [] tmp_d;
				 ////break;
			     ////}
		//case REDUCE : {
				  ////if(tape[i+1].arg2 == ADD)
				      ////op = MPI_SUM;
				  ////if(tape[i+1].arg2 ==MUL)
				      ////op = MPI_PROD;
				  ////tmp_d = new double[tape[i].arg1];
				  ////AMPI_Reduce_b(tmp_d, &tape[i+1].d, tape[i].arg1, MPI_DOUBLE, op, tape[i].arg2, comm);
				  ////switch(tape[i+1].arg2) {
				      ////case ADD : {
						     ////for(int j = 0 ; j < tape[i].arg1 ; j++) {
							 //////tape[i-tape[i].arg1+j].d += tape[i+1].d;
							 ////tape[i-tape[i].arg1+j].d += tmp_d[j];
						     ////}
						     ////break;
						 ////}
				      ////case MUL : {
						     ////for(int j = 0 ; j < tape[i].arg1 ; j++) {
							 //////tape[i-tape[i].arg1+j].d = (tape[i+1].v/tape[i-tape[i].arg1+j].v)*tape[i+1].d;
							 ////tape[i-tape[i].arg1+j].d = tmp_d[j];
						     ////}
						     ////break;
						 ////}
				  ////}
				  ////delete [] tmp_d;
				  ////break;
			      ////}
  }
}
#if (DEBUG_ACTIVE == 3)
void compute_jacobian_pattern(unsigned int**& P, int& m) {
  cout << "Entering compute_jacobian_pattern" << endl;
  assert(tape);

  int n = int(tape_entry::indeps.size());
  m = int(tape_entry::deps.size());
  P = new unsigned int* [m];
  bool** bJ = new bool* [m];
  for (int i=0; i<m; i++)  bJ[i]=new bool [n];
  int* nz_ctr = new int [m];
  memset(nz_ctr, 0, m*sizeof(int));

  for (int i=0; i<n; i++)  tape[tape_entry::indeps[i]].ds=0;

  for (int i=0; i<n; i++) {

    tape[tape_entry::indeps[i]].ds=1;

    for (int j=0; j<tape_entry::vac; j++) {
      int arg1 = tape[j].arg1;
      int arg2 = tape[j].arg2;
      if (arg1 != UNDEF && arg2 != UNDEF) tape[j].ds = tape[arg1].ds || tape[arg2].ds;
      else if (arg1 != UNDEF) tape[j].ds = tape[arg1].ds;
      else if (arg2 != UNDEF) tape[j].ds = tape[arg2].ds; 
//       else tape[j].ds=0;  // ERROR SINCE INDEPENDENTS ARE RESET TO ZERO
    }
    tape[tape_entry::indeps[i]].ds=0;
    for (int j=0; j<m; j++) {
//       cout << tape[tape_entry::deps[j]].ds << " ";
      bJ[j][i] = tape[tape_entry::deps[j]].ds;
      if (bJ[j][i])  nz_ctr[j]++;
    }
  }

  for (int j=0; j<m; j++) {
    P[j] = new unsigned int [nz_ctr[j]+1];
    P[j][0] = nz_ctr[j];
    int next=1;
    for (int i=0; i<n; i++)
      if (bJ[j][i])  P[j][next++] = i;
  }
  cout << "Leaving compute_jacobian_pattern" << endl;  
}
#endif
/*
 * See active.h
 */
void print_tape(int i, string kind) {
  cout << kind << " > " << i << " th tape entry : [ " 
       << tape[i].oc << ", "
       << tape[i].arg1 << ", "
       << tape[i].arg2 << ", "
       << tape[i].v << ", "
       << tape[i].d << " ]" << endl;
  if (kind == "WARN") {
    if (tape[i].arg1 != UNDEF)  print_tape(tape[i].arg1);
    if (tape[i].arg2 != UNDEF)  print_tape(tape[i].arg2);
  }
}
/*
 * See active.h
 */
void print_tape() {
  cout << "THE TAPE : " << endl;
  for (int i=0; i<tape_entry::vac; i++)  print_tape(i);
  ampi_print_tape();
}
/*
 * See active.h
 */
void print_tape_indeps() {
  cout << "INDEPENDET TAPES:" << endl;
  for (unsigned int i=0; i<tape_entry::indeps.size(); i++)
    print_tape(tape_entry::indeps[i]);
}
/*
 * See active.h
 */
void print_tape_deps() {
  cout << "DEPENDET TAPES:" << endl;
  for (unsigned int i=0; i<tape_entry::deps.size(); i++)
    print_tape(tape_entry::deps[i]);
}
/*
 * See active.h
 */
void print_tape_extremes() {
  print_tape_indeps();
  print_tape_deps();
}
/*
 * See active.h
 */
void print_deriv_of_indeps(int dep_indx) {
  for (unsigned int i=0; i<tape_entry::indeps.size(); i++)
//     if (tape[tape_entry::indeps[i]].d)
      cout << dep_indx << " " << i << " " << tape[tape_entry::indeps[i]].d << endl;

  cout << endl;
}
void print_deriv_of_deps(int indep_indx) {
  for (unsigned int i=0; i<tape_entry::deps.size(); i++)
//     if (tape[tape_entry::deps[i]].d)
      cout << i << " " << indep_indx << " " << tape[tape_entry::deps[i]].d << endl;

  cout << endl;
}
/*
 * See active.h
 */
std::ostream& operator<<( std::ostream& out, const active& x) {
  if (TAPE_MODE) {
      out << x.va << "(" << x.v << "): [ " << get_oc(tape[x.va].oc) << ", "
          << tape[x.va].arg1 << ", "
          << tape[x.va].arg2 << ", "
          << tape[x.va].v << ", "
          << tape[x.va].d << " ]" << endl;
  }
  else  out  <<  x.va  <<  " ( "  <<  x.v  <<  " , "  <<  x.d  <<  " )";
  return out;
}
std::ostream& operator<<( std::ostream& out, const tape_entry& x) {
  out << " [" << x.oc << ", " << x.arg1 << ", " << x.arg2
      << ", " << x.v << ", " << x.d << " ]" << endl;
  return out;
}
std::istream& operator>>( std::istream& in, active& x) {
  if (x.va >= 0) in >> x.v;
  return in;
}
/*
 * See active.h
 */
int print_tape_entry_to_dot(ofstream& out, const int& i) {
  out << i << "[label=\""
      << i
      << ":" << tape[i].oc
      << "," << tape[i].v
      << "," << tape[i].d
      << "," << tape[i].arg1
      << "," << tape[i].arg2
      << "\"]"<< endl;

  if (tape[i].arg1 != UNDEF) {
    out  <<  tape[i].arg1  <<  " -> " << i  <<  endl;
    print_tape_entry_to_dot(out, tape[i].arg1);
  }
  if (tape[i].arg2 != UNDEF) {
    out  <<  tape[i].arg2  <<  " -> " << i  <<  endl;
    print_tape_entry_to_dot(out, tape[i].arg2);
  } 
  return 0;
}
/*
 * See active.h
 */
void print_spanning_tree_to_dot(const char* filename, const active& y) {
  cout << "Entering print_spanning_tree_to_dot" << endl;
  if (y.va == UNDEF) {
    cout << "WARN > y.va = " << y.va << endl;
  }
  else {
    ofstream out (filename);
    out<<"digraph {"<<endl<<endl;
    print_tape_entry_to_dot(out, y.va);
/*      if (tape[y.va].arg1 != UNDEF) {
        out  <<  tape[y.va].arg1  <<  " -> " << y.va  <<  endl;
        print_tape_entry_to_dot(out, tape[y.va].arg1);
      }
      if (tape[y.va].arg2 != UNDEF) {
        out  <<  tape[y.va].arg2  <<  " -> " << y.va  <<  endl;
        print_tape_entry_to_dot(out, tape[y.va].arg2);
      } */
    out << "}"  <<  endl;
    out.close();
  }
  cout << "Leaving print_spanning_tree_to_dot" << endl;  
}
/*
 * See active.h
 */
void print_tape_to_dot(const char* filename, int start, int end) {
  ofstream out (filename);
  out << "digraph {" << endl << endl;
  for (int i=start; i<=end; i++) {
    if ( tape[i].arg1 != -1)
    {
      out <<i<<"[label=\""
         <<get_oc(tape[i].oc)
         <<"("<<tape[i].v<<","
         <<tape[i].d<<")\"]"<<endl;

      out<<tape[i].arg1<<"[label=\""
         <<get_oc(tape[tape[i].arg1].oc)
         <<"("<<tape[tape[i].arg1].v<<","
         <<tape[tape[i].arg1].d<<")\"]"  <<  endl;

      out  <<  tape[i].arg1  <<  " -> "<<i  <<  endl;
    }
    if ( tape[i].arg2 != -1)
    {
      out<<tape[i].arg2<<"[label=\""
         <<get_oc(tape[tape[i].arg2].oc)
         <<"("<<tape[tape[i].arg2].v<<","
         <<tape[tape[i].arg2].d<<")\"]"<<endl;
      out<<tape[i].arg2<<" -> "<<i<<endl;
    }
  }
  out << "}"  <<  endl;
  out.close();
}

void tape_to_dot(const char* filename, const int& i, const int& depth) {
  cout<<"print_tape_to_dot:"<<endl;
  ofstream out (filename);
  int ctr = 0;

  out<<"digraph {"<<endl<<endl;

  to_dot(out, i, depth, ctr);

  out << "}"  <<  endl;
  out.close();
}
int to_dot(ofstream& out, const int& i, const int& depth,  int& ctr) {
  out<<i<<"[label=\""
     <<tape[i].oc
     <<":"<<i
     <<","<<tape[i].v
     <<","<<tape[i].d
     <<"\"]"<<endl;

//   print_tape(i);

  if (ctr == depth) return --ctr;
  if (tape[i].arg1 != -1)
  {
    out<<tape[i].arg1<<"[label=\""
//        <<get_oc(tape[tape[i].arg1].oc)
       <<tape[tape[i].arg1].oc
       <<":"<<tape[i].arg1<<","
       <<tape[tape[i].arg1].v<<","
       <<tape[tape[i].arg1].d<< "\"]"<<endl;

    out<<tape[i].arg1<<" -> "<<i<<endl;

    return to_dot(out, tape[i].arg1, depth, ++ctr);
  }
  if ( tape[i].arg2 != -1)
  {
    out<<tape[i].arg2<<"[label=\""
//        <<get_oc(tape[tape[i].arg2].oc)
       <<tape[tape[i].arg2].oc
       <<":"<<tape[i].arg2<<","
       <<tape[tape[i].arg2].v
       <<","<<tape[tape[i].arg2].d
       <<"\"]"<<endl;

    out<<tape[i].arg2<<" -> "<<i<<endl;

    return to_dot(out, tape[i].arg2, depth, ++ctr);
  }
  return --ctr;
}

string get_oc(const int& i) {
  switch (i) {
    case -1 : {
      return "-1";
      break;
    }
    case 0 : {
      return "INDEP";
      break;
    }
    case 1 : {
      return "DEP";
      break;
    }
    case 2 : {
      return "C";
      break;
    }
    case 3 : {
      return "CP";
      break;
    }
    case 4 : {
      return "=";
      break;
    }
    case 5 : {
      return "+=";
      break;
    }
    case 6 : {
      return "-=";
      break;
    }
    case 7 : {
      return "*=";
      break;
    }
    case 8 : {
      return "/=";
      break;
    }
    case 9 : {
      return "+";
      break;
    }
    case 10 : {
      return "-";
      break;
    }
    case 11 : {
      return "=-";
      break;
    }
    case 12 : {
      return "*";
      break;
    }
    case 13 : {
      return "/";
      break;
    }
    case 14 : {
      return "sin";
      break;
    }
    case 15 : {
      return "asin";
      break;
    }
    case 16 : {
      return "sinh";
      break;
    }
    case 17 : {
      return "asinh";
      break;
    }
    case 18 : {
      return "cos";
      break;
    }
    case 19 : {
      return "acos";
      break;
    }
    case 20 : {
      return "cosh";
      break;
    }
    case 21 : {
      return "acosh";
      break;
    }
    case 22 : {
      return "tan";
      break;
    }
    case 23 : {
      return "atan";
      break;
    }
    case 24 : {
      return "tanh";
      break;
    }
    case 25 : {
      return "atanh";
      break;
    }
    case 26 : {
      return "exp";
      break;
    }
    case 27 : {
      return "log";
      break;
    }
    case 28 : {
      return "lop1p";
      break;
    }
    case 29 : {
      return "log10";
      break;
    }
    case 30 : {
      return "pow";
      break;
    }
    case 31 : {
      return "sqrt";
      break;
    }
    case 32 : {
      return "fabs";
      break;
    }
    case 33 : {
      return "max";
      break;
    }
    case 34 : {
      return "min";
      break;
    }
    case 35 : {
      return "hypot";
      break;
    }
    case 36 : {
      return "floor";
      break;
    }
    case 37 : {
      return "ceil";
      break;
    }
    case 38 : {
      return "ldexp";
      break;
    }
    case 39 : {
      return "RECV";
      break;
    }
    case 40 : {
      return "ISEND";
      break;
    }
    case 41 : {
      return "IRECV";
      break;
    }
    case 42 : {
      return "WAIT";
      break;
    }
    case 43 : {
      return "WAITALL";
      break;
    }
    case 44 : {
      return "AWAITALL";
      break;
    }
    case 45 : {
      return "BROADCAST";
      break;
    }
    case 46 : {
      return "REDUCE";
      break;
    }
    case 47 : {
      return "MPI_DUMMY";
      break;
    }
    case 48 : {
      return "AMPI_DUMMY";
      break;
    }
    case 49 : {
      return "SEND";
      break;
    }
    default : {
      return " ";
      break;
    }
  }
  return 0;
}

#if (DEBUG_ACTIVE == 1)
void print_debug_info() {
  cout << "a_ctr         > " << a_ctr << endl
       << "ai_ctr        > " << ai_ctr << endl
       << "al_ctr        > " << al_ctr << endl
       << "ad_ctr        > " << ad_ctr << endl
       << "ald_ctr       > " << ald_ctr << endl       
       << "aa_ctr        > " << aa_ctr << endl       
       << "aeqa_ctr      > " << aeqa_ctr << endl
       << "aeqp_ctr      > " << aeqp_ctr << endl
       << "nsign_ctr     > " << nsign_ctr << endl
       << "ainc_ctr      > " << ainc_ctr << endl
       << "adec_ctr      > " << adec_ctr << endl      
       << "atimesa_ctr   > " << atimesa_ctr << endl       
       << "ptimesa_ctr   > " << ptimesa_ctr << endl
       << "atimesp_ctr   > " << atimesp_ctr << endl
       << "adiva_ctr     > " << adiva_ctr << endl
      
 << "pdiva_ctr     > " << pdiva_ctr << endl
//        << "adivp_ctr     > " << adivp_ctr << endl
       << "aplusa_ctr    > " << aplusa_ctr << endl
       << "pplusa_ctr    > " << pplusa_ctr << endl
       << "aplusp_ctr    > " << aplusp_ctr << endl
       << "aminusa_ctr   > " << aminusa_ctr << endl
       << "pminusa_ctr   > " << pminusa_ctr << endl
       << "aminusp_ctr   > " << aminusp_ctr << endl
       << "atimeseqa_ctr > " << atimeseqa_ctr << endl
       << "atimeseqp_ctr > " << atimeseqp_ctr << endl
       << "adiveqa_ctr   > " << adiveqa_ctr << endl
       << "adiveqp_ctr   > " << adiveqp_ctr << endl
       << "apluseqa_ctr  > " << apluseqa_ctr << endl
       << "apluseqp_ctr  > " << apluseqp_ctr << endl
       << "aminuseqa_ctr > " << aminuseqa_ctr << endl
       << "aminuseqp_ctr > " << aminuseqp_ctr << endl
       << "sin_ctr       > " << sin_ctr << endl
       << "asin_ctr      > " << asin_ctr << endl
       << "sinh_ctr      > " << sinh_ctr << endl
       << "asinh_ctr     > " << asinh_ctr << endl
       << "cos_ctr       > " << cos_ctr << endl
       << "acos_ctr      > " << acos_ctr << endl
       << "cosh_ctr      > " << cosh_ctr << endl
       << "acosh_ctr     > " << acosh_ctr << endl
       << "tan_ctr       > " << tan_ctr << endl
       << "atan_ctr      > " << atan_ctr << endl
       << "tanh_ctr      > " << tanh_ctr  << endl
       << "atanh_ctr     > " << atanh_ctr << endl
       << "aatan2a_ctr   > " << aatan2a_ctr << endl
       << "aatan2p_ctr   > " << aatan2p_ctr << endl
       << "patan2a_ctr   > " << patan2a_ctr << endl
       << "afmaxa_ctr    > " << afmaxa_ctr << endl
       << "afmaxp_ctr    > " << afmaxp_ctr << endl
       << "pfmaxa_ctr    > " << pfmaxa_ctr << endl
       << "afmina_ctr    > " << afmina_ctr << endl
       << "afminp_ctr    > " << afminp_ctr << endl
       << "pfmina_ctr    > " << pfmina_ctr << endl
       << "amaxa_ctr     > " << amaxa_ctr << endl
       << "amaxp_ctr     > " << amaxp_ctr << endl
       << "pmaxa_ctr     > " << pmaxa_ctr << endl
       << "amina_ctr     > " << amina_ctr << endl
       << "aminp_ctr     > " << aminp_ctr << endl
       << "pmina_ctr     > " << pmina_ctr << endl
       << "exp_ctr       > " << exp_ctr << endl
       << "apowa_ctr     > " << apowa_ctr << endl
       << "apowp_ctr     > " << apowp_ctr << endl
       << "ppowa_ctr     > " << ppowa_ctr << endl
       << "sqrt_ctr      > " << sqrt_ctr << endl
       << "ahypota_ctr   > " << ahypota_ctr << endl
       << "ahypotp_ctr   > " << ahypotp_ctr << endl
       << "phypota_ctr   > " << phypota_ctr << endl
       << "log_ctr       > " << log_ctr << endl
       << "log10_ctr     > " << log10_ctr << endl
       << "log1p_ctr     > " << log1p_ctr << endl
       << "fabs_ctr      > " << fabs_ctr << endl
       << "abs_ctr       > " << abs_ctr << endl
       << "aldexpa_ctr   > " << aldexpa_ctr << endl
       << "aldexpp_ctr   > " << aldexpp_ctr << endl
       << "pldexpa_ctr   > " << pldexpa_ctr << endl
       << "afrexpi_ctr   > " << afrexpi_ctr << endl
       << "floor_ctr     > " << floor_ctr << endl
       << "ceil_ctr      > " << ceil_ctr << endl
       << "setv_ctr      > " << setv_ctr << endl << endl
       << "tape_entry::vac > " << tape_entry::vac << endl;
}
#endif
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// #if ENABLE_ACTIVE_OPT
// 
// #include "active.h"
// 
// TEntry *adjointtape=new TEntry[200000000];
// TEntry *ins=adjointtape;
// 
// double active::d=0;
// 
// vector<int> TEntry::deps;
// vector<int> TEntry::indeps;
