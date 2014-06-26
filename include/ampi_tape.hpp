#ifndef AMPI_TAPE_HPP
#define AMPI_TAPE_HPP

/*! \mainpage 
 * The adjoint MPI library served as a prototyping library for developing and applying adjoint MPI
 * patterns in the context of my disseration "Semantics Driven Adjoints of the Message Parsing Interface".
 * It is provided as is under the MIT licence. In the disseration you will find
 * the theoretical backround for the reversing or adjoining MPI communication.
 * The unique feature of this library is its generic interface in
 * ampi_interface.h that allows an easy integration into other AD tools 
 * not yet aware of MPI calls.
 */

/** \file
 * \brief A C++ wrapper for the C header ampi_tape.h.
 */


#define NO_INCLUDE_MPI
#include <mpi.h>
extern "C" {
#include "ampi_tape.h"
}
#endif
