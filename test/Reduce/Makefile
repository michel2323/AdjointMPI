include ../Makefile.inc

NPROC=4
PSIZE=10

DCO_BUILD=$(DCO_BASE_DIR)/build/dco_cpp

CXXFLAGS=-g  -O0 -I$(DCO_BUILD)/include -L$(DCO_BUILD)/lib -fopenmp

default: driver.out fdDriver.out ampiDriver.out

driver.out : driver
	mpiexec --output-filename $@ -n $(NPROC)  ./$< $(PSIZE)

fdDriver.out : fdDriver
	mpiexec --output-filename $@ -n $(NPROC)  ./$< $(PSIZE)

ampiDriver.out : ampiDriver
	mpiexec --output-filename $@ -n $(NPROC) ./$< $(PSIZE)

driver: driver.cpp
	${MPICXX} $(CXXFLAGS) -o $@ $^ -lm

fdDriver: fdDriver.cpp
	${MPICXX} $(CXXFLAGS) -o $@ $^ -lm

ampiDriver.o: ampiDriver.cpp
	$(MPICXX) $(CXXFLAGS) -I$(AMPI_DIR)/include -c -o $@ $^

ampi_interface.o: ampi_interface.cpp
	$(MPICXX) $(CXXFLAGS) -I$(AMPI_DIR)/include -c -o $@ $^

ampiDriver: ampiDriver.o ampi_interface.o
	$(MPICXX) $(CXXFLAGS) -o $@ $^ -ldco $(AMPI_DIR)/src/libAMPI.a

clean: 
	rm -f *.out* driver fdDriver ampiDriver *.o 
check:  default
	@echo "SendRecv test..."
	cat driver.out* fdDriver.out* ampiDriver.out* > test.out
	diff test.out all_ok.check	
	$(INDENTATION_OK)
