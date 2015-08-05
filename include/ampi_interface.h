#ifndef AMPI_INTERFACE_H
#define AMPI_INTERFACE_H
/** \file 
 * \brief AMPI interface routines which are defined as external. These routines need to be implemented by
 * the external AD tool library. They define the data flow between the external
 * tape and the AMPI tape.
 */

/**
 * \def INT64
 * Defines the type of the tape index as defined by the AD tool. Due to circular
 * dependence, this has to be defined here in addition to the ampi interface.
 */
#define INT64 long int

/** Call AMPI tape interpreter from external tape */
void ampi_interpret_tape(void* handle);
void ampi_reset_entry(void* handle);

/** Get a value *v from a specific tape variable buf[i] */
extern void ampi_get_val(void* buf, int* i, double* v);

/** Set a value *v from a specific tape variable buf[i] */
extern void ampi_set_val(void* buf, int* i, double* v);

/** Get an adjoint a from a specific tape entry with index idx */
extern void ampi_get_adj(INT64* idx, double* a);

/** Set an adjoint a from a specific tape entry with index idx */
extern void ampi_set_adj(INT64*, double*);

/** Get a tape index from a variable buf[i] */
extern void ampi_get_idx(void* buf, int* i, INT64* idx);

/** Create a tape entry in the external tape, indicating an external AMPI call */
extern void ampi_create_tape_entry(void* handle);

/** Create size tape entries to store the values of buf. Refer to receive buffer without
 * initialization */
extern void ampi_create_dummies(void* buf, int *size);

/* Returns 1 if the tape is active, otherwise 0. */
extern int ampi_is_tape_active();
#endif
