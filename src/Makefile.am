CC=mpiicc
AM_CXXFLAGS = -Wall -xsse4.1 -O3 -no-prec-dev -I$(top_srcdir)/include
AM_CFLAGS = -Wall -xsse4.1 -ansi-alias  -O3 -Wno-unknown-pragmas -std=c89 -I$(top_srcdir)/include 
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

lib_LIBRARIES = libAMPI.a

libAMPI_adir = $(includedir)

libAMPI_a_HEADERS=\
                  ../include/ampi.h ../include/ampi_stack.h \
		   \
		  ../include/ampi_tape.h

libAMPI_a_SOURCES = \
                    $(libAMPI_a_HEADERS) \
                    ampi.c ampi_stack.c \
		    ampi_tape.c
