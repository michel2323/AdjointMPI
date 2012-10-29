include ../Makefile.inc

all: awaitall 


awaitall: awaitall.o
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -L$(AMPI_LIBDIR) -o $@ $^ -lm -lAMPI 
awaitall.o: awaitall.cpp
	$(CXX) $(CFLAGS) -I$(AMPI_INCLUDE) -c $^  
clean:
	rm -f ./*.o awaitall
