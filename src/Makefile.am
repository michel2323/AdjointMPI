CC=@MPICC@
AM_CXXFLAGS = -Wall -I$(top_srcdir)/include $(OPENMP_CXXFLAGS)
AM_CFLAGS = -Wall -Wno-unknown-pragmas -std=c89 -I$(top_srcdir)/include $(OPENMP_CFLAGS)
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

lib_LIBRARIES = libAMPI.a

libAMPI_adir = $(includedir)

libAMPI_a_HEADERS=\
                  ../include/ampi.h ../include/ampi_stack.h \
		   \
		  ../include/ampi_tape.h ../include/uthash.h \
		   \
		   ../include/ampi_interface.h ../include/ampi_tape.hpp 

libAMPI_a_SOURCES = \
                    $(libAMPI_a_HEADERS) \
                    ampi.c ampi_stack.c \
		    ampi_tape.c \
                    ../include/ampi.h ../include/ampi_stack.h \
		    \
		    ../include/ampi_tape.h ../include/uthash.h \
		    \
		   ../include/ampi_interface.h ../include/ampi_tape.hpp 
