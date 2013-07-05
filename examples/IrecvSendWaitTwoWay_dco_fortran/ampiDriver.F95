SUBROUTINE HEAD( x, y )

    USE DCO_A1S
    USE AMPI_FORTRAN_INTERFACE

    IMPLICIT NONE

    include 'mpif.h'

    TYPE(DCO_A1S_TYPE), DIMENSION(1)  :: x,y,local
    INTEGER                        :: rank, ier, request

    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ier)

    IF(rank .EQ. 0) THEN
      x(1)=x(1)*2
      CALL AMPI_RECV(y, 1, MPI_DOUBLE_PRECISION, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      !CALL AMPI_IRECV(y, 1, MPI_DOUBLE_PRECISION, 1, 0, MPI_COMM_WORLD, request, ier)
      CALL AMPI_SEND(x, 1, MPI_DOUBLE_PRECISION, 1, 0, MPI_COMM_WORLD,ier)
      !CALL AMPI_WAIT(request, MPI_STATUS_IGNORE, ier)
      y(1)=y(1)*3
    ENDIF
    IF(rank .EQ. 1) THEN
      CALL AMPI_RECV(local, 1, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      !CALL AMPI_IRECV(local, 1, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, request, ier)
      !CALL AMPI_WAIT(request, MPI_STATUS_IGNORE, ier)
      local(1)=SIN(local(1))
      CALL AMPI_SEND(local, 1, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ier)
    ENDIF

END SUBROUTINE HEAD


PROGRAM AMPIDRIVER

  USE DCO_A1S
  USE AMPI_FORTRAN_INTERFACE

  IMPLICIT NONE

  include 'mpif.h'

  TYPE(DCO_A1S_TYPE), DIMENSION(1)             :: xv,x,y
  INTEGER                        :: rank, ier
  DOUBLE PRECISION               :: g

  CALL DCO_A1S_TAPE_CREATE( )
  CALL AMPI_INIT(ier)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ier)

  IF(rank .EQ. 0) THEN
    x(1)=3.5D0
    xv=x
    CALL DCO_A1S_TAPE_REGISTER_VARIABLE( x(1) )
    CALL HEAD(x,y)
    PRINT *, 'Process 0 got number: ', y(1)
    CALL DCO_A1S_SET(y(1), 1.0D0, -1)
    CALL DCO_A1S_TAPE_INTERPRET_ADJOINT
    CALL DCO_A1S_GET( x(1), g, -1)
    PRINT *, 'Process 0 got gradient: ', g
  ENDIF
  IF(rank .EQ. 1) THEN
    CALL HEAD(x,y)
    CALL DCO_A1S_TAPE_INTERPRET_ADJOINT
  ENDIF

  CALL AMPI_FINALIZE(ier)

END PROGRAM AMPIDRIVER
