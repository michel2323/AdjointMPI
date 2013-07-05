ifndef MPIF90
  MPIF90=mpif90
endif

FFLAGS=-g -O0
DCO_BUILD=$(DCO_BASE_DIR)/build/dco_fortran

default: driver.out #fdDriver.out ampiDriver.out
	cat $^

driver.out: driver
	mpiexec -n 2 ./$< > $@

driver: driver.F95
	${MPIF90} $(FFLAGS) -o $@ $^

clean:
	rm -f *.out driver fdDriver ampiDriver *.o
