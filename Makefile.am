CC=@MPICC@
AM_CXXFLAGS = -I$(srcdir)/include  
AM_CFLAGS   = -I$(srcdir)/include $(OPENMP_FLAGS) @MPI_FLAGS@
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

lib_LIBRARIES = libAMPI.a

libAMPI_adir = $(includedir)

libAMPI_a_HEADERS=\
                  include/ampi.h include/ampi_stack.h \
		   \
		  include/ampi_tape.h include/uthash.h \
		   \
		  include/ampi_interface.h include/ampi_tape.hpp 

libAMPI_a_SOURCES = \
                    $(libAMPI_a_HEADERS) \
                    src/ampi.c src/ampi_stack.c \
		    src/ampi_tape.c \
                    include/ampi.h include/ampi_stack.h \
		    \
		    include/ampi_tape.h include/uthash.h \
		    \
		    include/ampi_interface.h include/ampi_tape.hpp 

AUTOMAKE_OPTIONS = subdir-objects

