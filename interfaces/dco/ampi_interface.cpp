#include "dco.hpp"
//#include "ampi_tape.hpp"

//#define INT64 int

typedef dco::ga1s<double> AD_MODE;
typedef AD_MODE::type type;

#ifndef DCO_AMPI
void ampi_interpret_tape(long int idx) {}
#endif
//long int ampi_counter=0;

void ampi_get_val(void *buf, int *i, double *x) {
    *x=static_cast<type*>(buf)[*i]._value();
}
void ampi_set_val(void* buf, int *i, double *v) {
    type &dummy= static_cast<type*>(buf)[*i];
    *const_cast<double*>(&(dummy._value())) = *v;
}

void ampi_get_idx(void *buf, int *i, INT64 *idx) {
    type &var = static_cast<type*>(buf)[*i];
    if(!var._data()._is_registered()) {
        *idx=0;
    }
    else {
        *idx = var._data().tape_index();
    }
}

void ampi_get_adj(INT64 *idx, double *x) {
    if(*idx!=0) *x = AD_MODE::global_tape->_adjoint(*idx);
}
void ampi_set_adj(INT64 *idx, double *x) {
    if(*idx!=0) AD_MODE::global_tape->_adjoint(*idx) += *x;
}


/*extern "C" */void ampi_reset_entry(long int idx);

struct AMPI_data : AD_MODE::tape::external_function_base_data  {
    int idx;
    AMPI_data(const int nidx) : idx(nidx) {}
    virtual ~AMPI_data() {
        //std::cout << "ampi_reset_entry with idx=" << idx << std::endl;
        ampi_reset_entry(idx);
    }
};

void ampi_tape_wrapper(AMPI_data *data) {
    ampi_interpret_tape(data->idx);
}


void ampi_create_tape_entry(long int *i) {
    if(AD_MODE::global_tape == NULL || !AD_MODE::global_tape->is_active()) {
        return;
    }
    AD_MODE::global_tape->create_ext_fcn_data<AMPI_data>(&ampi_tape_wrapper, *i);
}

void ampi_create_dummies(void *buf, int *size) {
    type *values=static_cast<type*>(buf);

    for(int i=0; i<*size; ++i) {
        type &dummy=values[i];
        dummy=0;
        AD_MODE::global_tape->register_variable(dummy);
    }
}

int ampi_is_tape_active () {
    if (NULL != AD_MODE::global_tape) {
#ifdef DCO_ALLOW_TAPE_SWITCH_OFF
        return AD_MODE::global_tape->is_active();
#else
        return 1;
#endif
    } else {
        return 0;
    }
}

