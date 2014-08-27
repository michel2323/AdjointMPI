#ifndef AMPI_FORTRAN_STUBS_H
#define AMPI_FORTRAN_STUBS_H
#include "ampi.h"
#include "ampi_tape.h"
#include <stdio.h>
#include <assert.h>
#define REQUEST_CHUNK_SIZE 1000
static int * AMPI_myid;
static int crequest;
AMPI_Request *request_idx;
typedef struct compad_type {
    double v;
    int j;
} compad_type;
void ampi_check_request_size();
#endif
