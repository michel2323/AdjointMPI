SUBROUTINE HEAD( x, y )

    IMPLICIT NONE

    include 'mpif.h'

    DOUBLE PRECISION  :: x,y,local
    INTEGER                        :: rank, ier, request

    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ier)

    IF(rank .EQ. 0) THEN
      x=x*2
      CALL MPI_IRECV(y, 1, MPI_DOUBLE_PRECISION, 1, 0, MPI_COMM_WORLD, request, ier)
      CALL MPI_SEND(x, 1, MPI_DOUBLE_PRECISION, 1, 0, MPI_COMM_WORLD,ier)
      CALL MPI_WAIT(request, MPI_STATUS_IGNORE, ier)
      y=y*3
    ENDIF
    IF(rank .EQ. 1) THEN
      CALL MPI_IRECV(local, 1, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, request, ier)
      CALL MPI_WAIT(request, MPI_STATUS_IGNORE, ier)
      local=SIN(local)
      CALL MPI_SEND(local, 1, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ier)
    ENDIF

END SUBROUTINE HEAD


PROGRAM FDDRIVER

  IMPLICIT NONE

  include 'mpif.h'

  DOUBLE PRECISION               :: xv,x,y,xp,yp,h
  INTEGER                        :: rank, ier

  CALL MPI_INIT(ier)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ier)

  h=1.0D-8

  IF(rank .EQ. 0) THEN
    x=3.5D0
    xv=x
    CALL HEAD(x,y)
    PRINT *, 'Process 0 got number: ', y
  ENDIF
  IF(rank .EQ. 1) THEN
    CALL HEAD(x,y)
  ENDIF

  IF(rank .EQ. 0) THEN
    xp=xv+h
    CALL HEAD(xp,yp)
    PRINT *, 'Process 0 got fd value: ', (yp-y)/h
  ENDIF
  IF(rank .EQ. 1) THEN
    CALL HEAD(xp,yp)
  ENDIF
  CALL MPI_FINALIZE(ier)

END PROGRAM FDDRIVER
