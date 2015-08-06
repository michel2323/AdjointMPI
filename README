The AMPI library is used to reverse MPI communication in the context of 
algorithmic differentiation (AD). In principal, the reversal may be 
used for other purposes. The library was closely developed with dco_cpp, 
an overloading AD tool developed at the STCE, RWTH Aachen university. 
The interface is specifically crafted in such a way that coupling 
with any AD tool should be straightforward. Please refer to the more
in depth documentation created by 'make doc'.

Source Files
----------

/src 		        All AMPI library source files.

- ampi_stack.c, 	Dynamic C stack.
- ampi.c, 		Core adjoint MPI routines. All source transformation
  			interfaces access these forward and backward routines.
- ampi_tape.c, 		Generic AMPI tape to integrated into AD overloading
			tools. It works for all tools that use a hash (pointer,
			tape index,...) to link value and adjoint.

/include 		All AMPI library header files.

Interface Files
----------
/interfaces 		AD Tool interfaces. They are not part of the AMPI library
			and should be considered as integration examples for
			other AD Tools.

- ampi_fortran_interface.F90    COMPAD Fortran interface. It integrates 
  ampi_fortran_stubs.c		the generic AMPI tape into the dco_fortran 
  ampi_fortran_stubs.h          overloading tape.

- dcc_alloc.cpp		Adjoint dynamic memory management for source 
  dcc_alloc.hpp  	transformation based second order AD in dcc. 

- dcc_mpi.cpp		Source transformation based second order forward over
  dcc_mpi.hpp		reverse AMPI interface.
- ampi_interface_master.cpp     dco master branch interface.
  
Test
-----------------
/test                   The unit tests currently only work with dco_cpp
