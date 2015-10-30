MODULE adj_tape_mpi

  use dco
  use ampi_fortran_interface
  use mpi

  implicit none

  private


  PUBLIC :: AMPI_INIT, AMPI_FINALIZE

  PUBLIC :: AMPI_TYPE
  PUBLIC :: AMPI_OP_SUM
  PUBLIC :: AMPI_OP_MIN
  PUBLIC :: AMPI_OP_MAX

  PUBLIC :: AMPI_REDUCE, AMPI_ALLREDUCE, AMPI_IRECV, AMPI_ISEND
  PUBLIC :: AMPI_WAIT, AMPI_WAITALL, AMPI_SEND, AMPI_RECV


  INTERFACE AMPI_REDUCE
     MODULE PROCEDURE AMPI_REDUCE_0d, AMPI_REDUCE_1d
  END INTERFACE

  INTERFACE AMPI_ALLREDUCE
     MODULE PROCEDURE AMPI_ALLREDUCE_0d, AMPI_ALLREDUCE_1d
  END INTERFACE

  INTEGER   AMPI_TYPE
  INTEGER   AMPI_OP_SUM, AMPI_OP_MIN, AMPI_OP_MAX

CONTAINS
  
  SUBROUTINE AMPI_INIT (ierr)
       INTEGER              :: ierr

       CALL AMPI_INIT_C(ierr)
       AMPI_TYPE = MPI_DOUBLE_PRECISION
       AMPI_OP_SUM = AMPI_SUM
       AMPI_OP_MIN = AMPI_MIN
       AMPI_OP_MAX = AMPI_MAX
  END SUBROUTINE AMPI_INIT

  SUBROUTINE AMPI_FINALIZE (ierr)
      INTEGER              :: ierr

      CALL AMPI_FINALIZE_C(ierr)
  END SUBROUTINE AMPI_FINALIZE

  SUBROUTINE AMPI_SEND (buf, count, datatype, dest, tag, & 
          comm, ierr) 
      TYPE(dco_type),DIMENSION(*)   :: buf
      INTEGER                          :: count
      INTEGER                          :: datatype
      INTEGER                          :: dest
      INTEGER                          :: tag
      INTEGER                          :: comm
      INTEGER                          :: ierr
  END SUBROUTINE AMPI_SEND

  SUBROUTINE AMPI_RECV (buf, count, datatype, dest, tag, & 
          comm, status, ierr)
      TYPE(DCO_TYPE),DIMENSION(*)   :: buf
      INTEGER                          :: count
      INTEGER                          :: datatype
      INTEGER                          :: dest
      INTEGER                          :: tag
      INTEGER                          :: comm
      INTEGER,DIMENSION(*)             :: status
      INTEGER                          :: ierr
  END SUBROUTINE AMPI_RECV

  SUBROUTINE AMPI_IRECV (buf, count, datatype, dest, tag, & 
          comm, request, ierr) 
      TYPE(DCO_TYPE),DIMENSION(*)   :: buf
      INTEGER                          :: count
      INTEGER                          :: datatype
      INTEGER                          :: dest
      INTEGER                          :: tag
      INTEGER                          :: comm
      INTEGER                          :: request
      INTEGER                          :: ierr
      CALL AMPI_IRECV_C(buf,count,datatype,dest,tag,comm,request,ierr)
  END SUBROUTINE AMPI_IRECV

  SUBROUTINE AMPI_ISEND (buf, count, datatype, dest, tag, & 
          comm, request, ierr) 
      TYPE(DCO_TYPE),DIMENSION(*)   :: buf
      INTEGER                          :: count
      INTEGER                          :: datatype
      INTEGER                          :: dest
      INTEGER                          :: tag
      INTEGER                          :: comm
      INTEGER                          :: request
      INTEGER                          :: ierr
      CALL AMPI_ISEND_C(buf,count,datatype,dest,tag,comm,request,ierr)
  END SUBROUTINE AMPI_ISEND

  SUBROUTINE AMPI_WAIT (request, status, ierr)  
      INTEGER                          :: request
      INTEGER, DIMENSION(*)            :: status
      INTEGER                          :: ierr
     CALL AMPI_WAIT_C(request, status, ierr) 
  END SUBROUTINE AMPI_WAIT

  SUBROUTINE AMPI_WAITALL (count, request, status, ierr)  
      INTEGER, DIMENSION(*)            :: request
      INTEGER, DIMENSION(*)            :: status
      INTEGER                          :: ierr
      INTEGER                          :: count
     CALL AMPI_WAITALL_C (count, request, status, ierr) 
  END SUBROUTINE AMPI_WAITALL



  SUBROUTINE AMPI_REDUCE_1d( sendbuf, recvbuf, count, datatype, &
          &                 op,  root, comm , ierror)

      TYPE(DCO_TYPE), DIMENSION(:), INTENT(IN)    :: sendbuf
      TYPE(DCO_TYPE), DIMENSION(:), INTENT(INOUT) :: recvbuf
      INTEGER, INTENT(IN)   :: count, datatype, op,  root, comm
      INTEGER, INTENT(OUT)  :: ierror    


      CALL  AMPI_REDUCE_C( sendbuf, recvbuf, count, datatype, &
          &                 op,  root, comm , ierror)

  END SUBROUTINE AMPI_REDUCE_1d

  SUBROUTINE AMPI_REDUCE_0d( sendbuf, recvbuf, count, datatype, &
       &                    op,  root, comm , ierror)

    TYPE(DCO_TYPE), INTENT(IN)    :: sendbuf
    TYPE(DCO_TYPE), INTENT(INOUT) :: recvbuf
    INTEGER, INTENT(IN)   :: count, datatype, op,  root, comm
    INTEGER, INTENT(OUT)  :: ierror    

    TYPE(DCO_TYPE), DIMENSION(1)  :: sendbuf2
    TYPE(DCO_TYPE), DIMENSION(1)  :: recvbuf2

    sendbuf2(1) = sendbuf
    recvbuf2(1) = recvbuf

    CALL  AMPI_REDUCE_C( sendbuf2, recvbuf2, count, datatype, &
       &                 op,  root, comm , ierror)

    recvbuf = recvbuf2(1)
    
  END SUBROUTINE AMPI_REDUCE_0d




  SUBROUTINE AMPI_ALLREDUCE_1d( sendbuf, recvbuf, count, datatype, &
       &                 op, comm , ierror)

    TYPE(DCO_TYPE), DIMENSION(:), INTENT(IN)    :: sendbuf
    TYPE(DCO_TYPE), DIMENSION(:), INTENT(INOUT) :: recvbuf
    INTEGER, INTENT(IN)   :: count, datatype, op, comm
    INTEGER, INTENT(OUT)  :: ierror    

    REAL, DIMENSION(count)    :: sendbuf_val
    REAL, DIMENSION(count)    :: recvbuf_val

    PRINT *, 'ALLREDUCE 1D IN in mpi.f90'
   
    IF(op .EQ. AMPI_SUM) THEN
    sendbuf_val = sendbuf%v
    CALL  MPI_ALLREDUCE( sendbuf_val, recvbuf_val, count, MPI_DOUBLE_PRECISION, &
       &                 MPI_SUM, comm , ierror)
    recvbuf%v = recvbuf_val
    ELSE
    CALL  AMPI_ALLREDUCE_C( sendbuf, recvbuf, count, datatype, &
       &                 op, comm , ierror)
    ENDIF
    
    PRINT *, 'ALLREDUCE 1D OUT in mpi.F90'
  END SUBROUTINE AMPI_ALLREDUCE_1d

  SUBROUTINE AMPI_ALLREDUCE_0d( sendbuf, recvbuf, count, datatype, &
       &                    op,  comm , ierror)

    REAL, DIMENSION(1)    :: sendbuf_val
    REAL, DIMENSION(1)   :: recvbuf_val
    TYPE(DCO_TYPE), INTENT(IN)    :: sendbuf
    TYPE(DCO_TYPE), INTENT(INOUT) :: recvbuf
    INTEGER, INTENT(IN)   :: count, datatype, op, comm
    INTEGER, INTENT(OUT)  :: ierror    
    INTEGER               :: myid   
    INTEGER               :: ierr  

    TYPE(DCO_TYPE), DIMENSION(1)  :: sendbuf2
    TYPE(DCO_TYPE), DIMENSION(1)  :: recvbuf2

    sendbuf_val(1) = sendbuf%v
    recvbuf_val(1) = recvbuf%v
    sendbuf2(1) = sendbuf
    recvbuf2(1) = recvbuf

    CALL  AMPI_ALLREDUCE_C( sendbuf2, recvbuf2, count, datatype, &
       &                 op, comm , ierror)
    recvbuf = recvbuf2(1)
  END SUBROUTINE AMPI_ALLREDUCE_0d



END MODULE dco_module_adj_tape_mpi
