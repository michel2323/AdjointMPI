include ../Makefile.inc

all: pi


pi: pi.o dco_tape.o
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -L$(AMPI_LIBDIR) -o $@ $^ -lm -lAMPI 
pi.o: pi.cpp
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
dco_tape.o: dco_tape.cpp
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
clean:
	rm -f ./*.o pi
