#ifndef AMPI_TAPE_HPP
#define AMPI_TAPE_HPP

/*! \mainpage 
 * The adjoint MPI library served as a prototyping library for developing and applying adjoint MPI
 * patterns in the context of my thesis "Semantics Driven Adjoints of the Message Parsing Interface
 */

#define NO_INCLUDE_MPI
#include <mpi.h>
extern "C" {
#include "ampi_stubs.h"
#include "ampi_tape.h"
}
#endif
