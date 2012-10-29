include ../Makefile.inc

all: reduce


reduce: reduce.o dco_tape.o
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -L$(AMPI_LIBDIR) -o $@ $^ -lm -lAMPI 
dco_tape.o: dco_tape.cpp
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
reduce.o: reduce.cpp
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
clean:
	rm -f ./*.o reduce
