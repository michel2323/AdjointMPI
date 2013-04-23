#include "dco.hpp"
#include "ampi_tape.hpp"

//#define INT64 int

using namespace dco::a1s;

extern "C" {
  //forward declare von AMPI

#ifndef DCO_AMPI
//void ampi_interpret_tape() {}
#endif


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
      //      *idx=&(global_tape->_adjoint(var._data().tape_index()));
      *idx = var._data().tape_index();
    }

  }


  void ampi_get_adj(INT64 *idx, double *x) {
    //if(*idx) *x = *(*idx);
    if(*idx!=0) *x = dco::a1s::global_tape->_adjoint(*idx);
    //std::cout << "get adj: " << *(*idx) << " idx=" << *idx << std::endl;
    //std::cout << "AMPI_GET_ADJ: " << *x << std::endl;
  }
  void ampi_set_adj(INT64 *idx, double *x) {
    //if(*idx!=0) const_cast<double&>(dco::a1s::global_tape->_adjoint(*idx)) = *x;
    if(*idx!=0) dco::a1s::global_tape->_adjoint(*idx) += *x;
    //    if(*idx) *(*idx)=*x;
    //std::cout << "set adj: " << *x << " idx=" << *idx << std::endl;
  }

  void ampi_tape_wrapper(tape &caller, const tape::interpretation_settings &settings, dco::a1s::tape::external_function_base_data *userdata) {
    ampi_interpret_tape();
  }

  void ampi_create_tape_entry(int *i) {
    //todo: insert an external function handler!!!
    global_tape->register_external_function(&ampi_tape_wrapper, 0);
    //this will call ampi_interpret_tape
  }

  void ampi_create_dummies(void *buf, int *size) {
    type *values=static_cast<type*>(buf);
    
    for(int i=0;i<*size;++i) {
      type &dummy=values[i];
      //std::cout << "dummy.value=" << dummy._value << " .tapeindex=" << dummy._data.tape_index << std::endl;
      //std::cout << "mpi_global_tape=" << mpi_global_tape << std::endl;
      //if(dummy._edgecount==0) {
      dummy=0;
      global_tape->register_variable(dummy);
    } 
  }
  
  //     const int n=*size;
//     type *actives=static_cast<type*>(buf);

//     int startindex=0;
//     dco::a1s::tape::TAPE_ENTRY *ins=global_tape->_get_insert_ptr_range(n, startindex );

// #ifdef DCO_OPEN_MP
// #pragma omp parallel for
// #endif
//     for(int i=0;i<n;++i) {
//       ins[i].arg=0;   //edgecount =0
//       actives[i] = 0;

//       dco::a1s::data &data = const_cast<dco::a1s::data&>(actives[i]._data());
//       data.register_variable( startindex+i ,global_tape);
//     }

    
  //}


}
