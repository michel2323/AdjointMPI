CXX       = mpicxx
CC       = mpicc
CFLAGS   =  -g -O0
#EXAMPLES = heat pi reduce awaitall mandelbrot speelpenning speelpenning_dcc
EXAMPLES = heat_dco speelpenning_dco pi_dco mandelbrot_dco awaitall_dco reduce_dco speelpenning_dcc 
AMPI_LIB_PATH = ../src
 

all: $(EXAMPLES)
heat_dco: force_look 
	cd $@ && $(MAKE) $(MFLAGS)
pi_dco: force_look 
	cd $@ && $(MAKE) $(MFLAGS)
reduce_dco: force_look 
	cd $@ && $(MAKE) $(MFLAGS)
awaitall_dco: force_look 
	cd $@ && $(MAKE) $(MFLAGS)
mandelbrot_dco: force_look 
	cd $@ && $(MAKE) $(MFLAGS)
speelpenning_dco: force_look 
	cd $@ && $(MAKE) $(MFLAGS)
speelpenning_dcc: force_look 
	cd $@ && $(MAKE) $(MFLAGS)
clean:
	-for d in $(EXAMPLES); do (cd $$d; $(MAKE) clean ); done

force_look:
	true
