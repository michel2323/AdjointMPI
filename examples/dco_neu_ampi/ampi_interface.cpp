#include "active_type.hpp"
//#include "dco.hpp"

#define INT64 double*

using namespace dcoV2::a1s;

extern "C" {
//forward declare von AMPI

//#ifdef DCO_AMPI
    void ampi_interpret_tape();
//#else
    //void ampi_interpret_tape() {}
//#endif



static static_tape *mpi_global_tape=0;

void ampi_set_globale_tape(static_tape *global_tape) {
    mpi_global_tape=global_tape;
}

void ampi_get_val(void *buf, int *i, double *x) {
    *x=static_cast<type*>(buf)[*i]._value;
}
void ampi_set_val(void* buf, int *i, double *v) {
    type &dummy= static_cast<type*>(buf)[*i];
    dummy._value = *v;
}

void ampi_get_idx(void *buf, int *i, INT64 *idx) {
    type &var = static_cast<type*>(buf)[*i];
    static_tape *tape=var._data.tape;
    if(tape==0) {
        *idx=0;
    }
    else {
        *idx=&tape->adjoints[var._data.tape_index];
    }
}


void ampi_get_adj(INT64 *idx, double *x) {
    if(*idx) *x = *(*idx);
    printf("adj: %f\n", *(*idx));
}
void ampi_set_adj(INT64 *idx, double *x) {
    if(*idx) *(*idx)=*x;
}

void ampi_tape_wrapper(static_tape &caller, int mode, void *userdata) {
    ampi_interpret_tape();
}

void ampi_create_tape_entry(int *i) {
    //todo: insert an external function handler!!!
    mpi_global_tape->register_external_function( &ampi_tape_wrapper, 0);
    //this will call ampi_interpret_tape
}

void ampi_create_dummies(void *buf, int *size) {
    type *values=static_cast<type*>(buf);

    for(int i=0;i<*size;++i) {
        type &dummy=values[i];
	//std::cout << "dummy.value=" << dummy._value << " .tapeindex=" << dummy._data.tape_index << std::endl;
	//std::cout << "mpi_global_tape=" << mpi_global_tape << std::endl;
	if(dummy._edgecount==0) {
        	dummy=0;
        	mpi_global_tape->register_variable( dummy);
	}
    }
}


}
