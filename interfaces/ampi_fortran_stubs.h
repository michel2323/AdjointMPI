#ifndef AMPI_FORTRAN_STUBS_H
#define AMPI_FORTRAN_STUBS_H
#include "ampi.h"
#include "ampi_tape.h"
#include <stdio.h>
#include <assert.h>
#define REQUEST_IDX_SIZE 10000
#define ASSERT
static int * AMPI_myid;
static int crequest;
AMPI_Request request_idx[REQUEST_IDX_SIZE];
typedef struct compad_type {
    double v;
    int j;
} compad_type;
#endif
