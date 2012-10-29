include ../Makefile.inc


all: speelpenning


speelpenning: speelpenning.o ampi_interface.o
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -L$(AMPI_LIBDIR) -o $@ $^ -lm -lAMPI 
	mpirun -np 2 ./speelpenning

speelpenning.o: speelpenning.cpp
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
ampi_interface.o: ampi_interface.cpp
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
clean:
	rm -f ./*.o speelpenning
