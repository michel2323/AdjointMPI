#ifndef AMPI_INTERFACE_HPP
#define AMPI_INTERFACE_HPP

/** \file
 * \brief A C++ wrapper for the C header ampi_interface.h.
 */


#define NO_INCLUDE_MPI
#include <mpi.h>
extern "C" {
#include "ampi_interface.h"
}
#endif
