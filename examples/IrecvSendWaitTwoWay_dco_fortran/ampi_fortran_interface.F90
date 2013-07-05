!AMPI Fortran to C interface
!---------------------------
!This interface links the AMPI routines directly to
!their AMPI C library counterpart. 

!All routines are prefixed with AMPI_. 


MODULE ampi_fortran_interface

  IMPLICIT NONE

  PRIVATE
  
  PUBLIC   ::    AMPI_SUM, AMPI_PROD, AMPI_MIN, AMPI_MAX

  PUBLIC  :: AMPI_INIT_C, AMPI_FINALIZE_C,           &
       &     AMPI_SEND_C, AMPI_RECV_C, AMPI_ISEND_C, AMPI_IRECV_C,      &
       &     AMPI_REDUCE_C,                                       &
       &     AMPI_ALLREDUCE_C,                                    &
       &     AMPI_WAIT_C,                                         &
       &     AMPI_WAITALL_C,                                      &
       &     AMPI_INTERPRET_TAPE, AMPI_PRINT_TAPE_ENTRY,        &
       &     AMPI_RESET_TAPE


  !AMPI constants

  INTEGER, PARAMETER                          :: AMPI_SUM  = 1 
  INTEGER, PARAMETER                          :: AMPI_PROD = 2 
  INTEGER, PARAMETER                          :: AMPI_MIN  = 3 
  INTEGER, PARAMETER                          :: AMPI_MAX  = 4 


  !AMPI routines which have an MPI counterpart that is adjoined.

  INTERFACE
     SUBROUTINE AMPI_RESET_TAPE (ierr) BIND(c,name='ampi_reset_tape_fort')
       USE iso_c_binding
       INTEGER(c_int)              :: ierr
     END SUBROUTINE AMPI_RESET_TAPE

     SUBROUTINE AMPI_INIT_C (ierr) BIND(c,name='ampi_init_fort')
       USE iso_c_binding
       INTEGER(c_int)              :: ierr
     END SUBROUTINE AMPI_INIT_C

     SUBROUTINE AMPI_FINALIZE_C (ierr) BIND(c,name='ampi_finalize_fort')
       USE iso_c_binding
       INTEGER(c_int)              :: ierr
     END SUBROUTINE AMPI_FINALIZE_C

     !SUBROUTINE AMPI_COMM_RANK_C (comm,myid,ierr) BIND(c,name='ampi_comm_rank_fort')
       !USE iso_c_binding
       !INTEGER(c_int)              :: comm
       !INTEGER(c_int)              :: myid
       !INTEGER(c_int)              :: ierr
     !END SUBROUTINE AMPI_COMM_RANK_C

     !SUBROUTINE AMPI_COMM_SIZE_C (comm,size,ierr) BIND(c,name='ampi_comm_size_fort')
       !USE iso_c_binding
       !INTEGER(c_int)              :: comm
       !INTEGER(c_int)              :: size
       !INTEGER(c_int)              :: ierr
     !END SUBROUTINE AMPI_COMM_SIZE_C

     SUBROUTINE AMPI_SEND_C (buf, count, datatype, dest, tag, & 
          comm, ierr) BIND(c,name='ampi_send_fort')
       USE iso_c_binding
       USE compad_module_adj_tape_common             
       TYPE(compad_type),DIMENSION(*)   :: buf
       INTEGER(c_int)                   :: count
       INTEGER(c_int)                   :: datatype
       INTEGER(c_int)                   :: dest
       INTEGER(c_int)                   :: tag
       INTEGER(c_int)                   :: comm
       INTEGER(c_int)                   :: ierr
     END SUBROUTINE AMPI_SEND_C

     SUBROUTINE AMPI_RECV_C (buf, count, datatype, dest, tag, & 
          comm, status, ierr) BIND(c,name='ampi_recv_fort')
       USE, INTRINSIC              :: iso_c_binding
       USE compad_module_adj_tape_common             
       TYPE(COMPAD_TYPE),DIMENSION(*)   :: buf
       INTEGER(c_int)                   :: count
       INTEGER(c_int)                   :: datatype
       INTEGER(c_int)                   :: dest
       INTEGER(c_int)                   :: tag
       INTEGER(c_int)                   :: comm
       INTEGER(c_int),DIMENSION(*)      :: status
       INTEGER(c_int)                   :: ierr
     END SUBROUTINE AMPI_RECV_C

     SUBROUTINE AMPI_ISEND_C (buf, count, datatype, dest, tag, & 
          comm, request, ierr) BIND(c,name='ampi_isend_fort')
       USE iso_c_binding
       USE compad_module_adj_tape_common             
       TYPE(COMPAD_TYPE),DIMENSION(*)   :: buf
       INTEGER(c_int)                   :: count
       INTEGER(c_int)                   :: datatype
       INTEGER(c_int)                   :: dest
       INTEGER(c_int)                   :: tag
       INTEGER(c_int)                   :: comm
       INTEGER(c_int)              :: request
       INTEGER(c_int)              :: ierr
     END SUBROUTINE AMPI_ISEND_C

     SUBROUTINE AMPI_IRECV_C (buf, count, datatype, dest, tag, & 
          comm, request, ierr) BIND(c,name='ampi_irecv_fort')
       USE iso_c_binding
       USE compad_module_adj_tape_common             
       TYPE(COMPAD_TYPE),DIMENSION(*)   :: buf
       INTEGER(c_int)                   :: count
       INTEGER(c_int)                   :: datatype
       INTEGER(c_int)                   :: dest
       INTEGER(c_int)                   :: tag
       INTEGER(c_int)                   :: comm
       INTEGER(c_int)              :: request
       INTEGER(c_int)              :: ierr
     END SUBROUTINE AMPI_IRECV_C

     SUBROUTINE AMPI_WAIT_C (request, status, ierr) & 
          BIND(c,name='ampi_wait_fort')
       USE, INTRINSIC              :: iso_c_binding
       INTEGER(c_int)              :: request
       INTEGER(c_int), DIMENSION(*):: status
       INTEGER(c_int)              :: ierr
     END SUBROUTINE AMPI_WAIT_C

     SUBROUTINE AMPI_WAITALL_C (count, request, status, ierr) & 
          BIND(c,name='ampi_waitall_fort')
       USE, INTRINSIC                  :: iso_c_binding
       INTEGER(c_int), DIMENSION(*)    :: request
       INTEGER(c_int), DIMENSION(*)  :: status
       INTEGER(c_int)                  :: ierr
       INTEGER(c_int)                  :: count
     END SUBROUTINE AMPI_WAITALL_C

     SUBROUTINE AMPI_REDUCE_C (sendbuf, recvbuf, count, datatype, & 
             op, root, comm, ierr) BIND(c,name='ampi_reduce_fort')
       USE, INTRINSIC              :: iso_c_binding
       USE compad_module_adj_tape_common             
       TYPE(COMPAD_TYPE),DIMENSION(*)   :: sendbuf
       TYPE(COMPAD_TYPE),DIMENSION(*)   :: recvbuf
       INTEGER(c_int)                   :: count
       INTEGER(c_int)                   :: datatype
       INTEGER(c_int)                   :: op
       INTEGER(c_int)                   :: root
       INTEGER(c_int)                   :: tag
       INTEGER(c_int)                   :: comm
       INTEGER(c_int)              :: ierr
     END SUBROUTINE AMPI_REDUCE_C

     SUBROUTINE AMPI_ALLREDUCE_C (sendbuf, recvbuf, count, datatype, & 
             op, comm, ierr) BIND(c,name='ampi_allreduce_fort')
       USE, INTRINSIC              :: iso_c_binding
       USE compad_module_adj_tape_common             
       TYPE(COMPAD_TYPE),DIMENSION(*)   :: sendbuf
       TYPE(COMPAD_TYPE),DIMENSION(*)   :: recvbuf
       INTEGER(c_int)                   :: count
       INTEGER(c_int)                   :: datatype
       INTEGER(c_int)                   :: op
       INTEGER(c_int)                   :: tag
       INTEGER(c_int)                   :: comm
       INTEGER(c_int)              :: ierr
     END SUBROUTINE AMPI_ALLREDUCE_C
     !External AMPI tape interpretation routine called by the COMPAD tape
     !interpreter if an AMPI tape entry is interpreted.

     SUBROUTINE AMPI_INTERPRET_TAPE () BIND(C, name='ampi_interpret_tape')
       USE iso_c_binding
     END SUBROUTINE AMPI_INTERPRET_TAPE

     SUBROUTINE AMPI_PRINT_TAPE_ENTRY (i) BIND(C, name='ampi_print_tape_entry')
       USE iso_c_binding
       INTEGER(c_int)              :: i
     END SUBROUTINE AMPI_PRINT_TAPE_ENTRY
  END INTERFACE

END MODULE ampi_fortran_interface
