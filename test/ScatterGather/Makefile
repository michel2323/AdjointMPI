include ../Makefile.inc

NPROC=4
PSIZE=10

DCO_INCLUDE=$(DCO_BASE_DIR)/build



all: driver.out fdDriver.out ampiDriver.out

driver.out: driver
	${MPIEXEC} --output-filename driver.out -n $(NPROC)  ./driver $(PSIZE)

fdDriver.out: fdDriver
	${MPIEXEC} --output-filename fdDriver.out -n $(NPROC)  ./fdDriver $(PSIZE)

ampiDriver.out: ampiDriver
	${MPIEXEC} --output-filename ampiDriver.out -n $(NPROC) ./ampiDriver $(PSIZE)

driver: driver.o
	${MPICXX} $(CXXFLAGS) -o $@ driver.o -lm

driver.o: driver.cpp 
	${MPICXX} $(CXXFLAGS) -c -o $@ driver.cpp 

fdDriver: fdDriver.o
	${MPICXX} $(CXXFLAGS) -o $@ fdDriver.o -lm

fdDriver.o: fdDriver.cpp 
	${MPICXX} $(CXXFLAGS) -c -o $@ fdDriver.cpp 

ampiDriver: ampiDriver.o ampi_interface.o
	$(MPICXX) $(CXXFLAGS) -o $@ ampiDriver.o ampi_interface.o $(AMPI_DIR)/libAMPI.a

ampiDriver.o: ampiDriver.cpp
	$(MPICXX) $(CXXFLAGS) -I$(AMPI_DIR)/include -c -o $@ ampiDriver.cpp

ampi_interface.o: ampi_interface.cpp
	$(MPICXX) $(CXXFLAGS) -I$(AMPI_DIR)/include -c -o $@ ampi_interface.cpp

clean: 
	rm -f *.out* driver fdDriver ampiDriver *.o 
check:  all
	@echo "ScatterGather test..."
	cat driver.out* fdDriver.out* ampiDriver.out* > test.out
	diff test.out all_ok.check	
	$(INDENTATION_OK)
