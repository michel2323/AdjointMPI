include ../Makefile.inc

all: sd 


sd: sd.o f.o
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -I$(DCO_INCLUDE) -L$(AMPI_LIBDIR) -L$(DCO_LIBDIR) -o $@ $^ -lm -lAMPI 
sd.o: sd.cpp
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -I$(DCO_INCLUDE) -c $^  
f.o: f.cpp
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -I$(DCO_INCLUDE) -c $^  
clean:
	rm -f ./*.o sd
