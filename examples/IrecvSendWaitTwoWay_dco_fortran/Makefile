ifndef MPIF90
  MPIF90=mpif90
endif

ifndef MPICXX
  MPICXX=mpicxx
endif

FC=gfortran
FFLAGS=-g -O0 

#DCO_BUILD=$(DCO_BASE_DIR)/build/dco_fortran

DCOF_BASE_DIR=$(DCO_BASE_DIR)/build/dco_fortran/lib/${FC}
DCOCPP_INCLUDE_DIR=$(DCO_BASE_DIR)/build/dco_cpp/include

default: driver.out fdDriver.out ampiDriver.out
	cat $^

driver.out: driver
	mpiexec -n 2 ./$< > $@

fdDriver.out: fdDriver
	mpiexec -n 2 ./$< > $@

ampiDriver.out: ampiDriver
	mpiexec -n 2 ./$< > $@

driver: driver.F95
	${MPIF90} $(FFLAGS) -o $@ $^

fdDriver: fdDriver.F95
	${MPIF90} $(FFLAGS) -o $@ $^

ampiDriver.o: ampiDriver.F95
	${MPIF90} $(FFLAGS) -I$(DCOF_BASE_DIR) -c -o $@ $^

ampi_fortran_interface.o: ampi_fortran_interface.F90
	${MPIF90} $(FFLAGS) -I$(DCOF_BASE_DIR) -c -o $@ $^

ampisupport.o: ampisupport.cpp
	$(MPICXX) $(CXXFLAGS) -I$(DCOCPP_INCLUDE_DIR) -I$(AMPI_DIR)/include -c -o $@ $^

ampi_fortran_stubs.o: ampi_fortran_stubs.c
	$(MPICXX) $(CXXFLAGS) -I$(DCOCPP_INCLUDE_DIR) -I$(AMPI_DIR)/include -c -o $@ $^

ampiDriver: ampiDriver.o ampisupport.o ampi_fortran_stubs.o ampi_fortran_interface.o
	#$(MPICXX) $(CXXFLAGS) -o $@ $^  $(DCOF_BASE_DIR)/libdcof.a -lstdc++ -lrt $(AMPI_DIR)/src/libAMPI.a
	$(MPIF90) $(FFLAGS) -o $@ $^  $(DCOF_BASE_DIR)/libdcof.a -lmpi_cxx -lstdc++ -lrt $(AMPI_DIR)/src/libAMPI.a


clean:
	rm -f *.out driver fdDriver ampiDriver *.o
