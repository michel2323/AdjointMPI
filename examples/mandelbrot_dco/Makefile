include ../Makefile.inc

all: mandelbrot 


mandelbrot: mandelbrot.o draw.o complex.o dco_tape.o
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -L$(AMPI_LIBDIR) -o $@ $^ -lm -lAMPI 
mandelbrot.o: mandelbrot.cpp
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
dco_tape.o: dco_tape.cpp
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
draw.o: draw.cpp
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
complex.o: complex.cpp
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
clean:
	rm -f ./*.o mandelbrot
