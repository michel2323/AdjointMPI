include ../Makefile.inc

all: sd 


sd: sd.o f.o ampi_interface.o
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -I$(DCO_INCLUDE) -L$(AMPI_LIBDIR) -o $@ $^ $(DCO_LIBDIR)/libdco.a -lm -lAMPI 
sd.o: sd.cpp
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -I$(DCO_INCLUDE) -c $^  
f.o: f.cpp
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -I$(DCO_INCLUDE) -c $^  
ampi_interface.o: ampi_interface.cc
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -I$(DCO_INCLUDE) -c $^  

clean:
	rm -f ./*.o sd
