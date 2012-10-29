#ifndef DEBUG_ACTIVE
#define DEBUG_ACTIVE -1
#endif

#ifndef ACTIVE_INCLUDE
#define ACTIVE_INCLUDE

/* AMPI */

#define AMPI_COMM_WORLD MPI_COMM_WORLD
#define AMPI_MAX_PROCESSOR_NAME MPI_MAX_PROCESSOR_NAME
#define AMPI_Status MPI_Status
#define AMPI_DOUBLE MPI_DOUBLE

#include <mpi.h>
#include <set>
#include <map>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <cstring>
#include <limits>
#include <iomanip>
#include <list>
#include <set>                                                                 
extern "C" {
#include "ampi_tape.h"
}

#define TAPE_SIZE 100000000
#define UNDEF -1
#define INDEP 0
#define DEP 1
#define CONST 2
#define COPY 3
#define ASG 4
#define PLUSEQ 5
#define MINUSEQ 6
#define MULTEQ 7
#define DIVEQ 8
#define ADD 9
#define SUB 10
#define NEGSIGN 11
#define MUL 12
#define DIV 13
#define SIN 14
#define ASIN 15
#define SINH 16
#define ASINH 17
#define COS 18
#define ACOS 19
#define COSH 20
#define ACOSH 21
#define TAN 22
#define ATAN 23
#define TANH 24
#define ATANH 25
#define EXP 26
#define LOG 27
#define LOG1P 28
#define LOG10 29
#define POW 30
#define SQRT 31
#define FABS 32
#define MAX 33
#define MIN 34
#define HYPOT 35
#define FLOOR 36  
#define CEIL 37
#define AMPI 38
//#define RECV 39
//#define ISEND 40
//#define IRECV 41
//#define WAIT 42
//#define WAITALL 43
//#define AWAITALL 44
//#define BCAST 45
//#define REDUCE 46
//#define MPI_DUMMY 47
//#define MPI_ADUMMY 48
//#define SEND 49

#define PASSIVE_ACTIVE(OP)                                                                   \
  tmp.va=tape_entry::vac++;                                                                  \
  tape[tmp.va].v=tmp.v;                                                                      \
  tape[tmp.va].d=0;                                                                          \
  tape[tmp.va].oc=OP;                                                                        \
  tape[tmp.va].arg1=UNDEF;                                                                   \
  tape[tmp.va].arg2=x2.va;                                                                   \

#define ACTIVE_PASSIVE(OP)                                                                   \
  tmp.va=tape_entry::vac++;                                                                  \
  tape[tmp.va].v=tmp.v;                                                                      \
  tape[tmp.va].d=0;                                                                          \
  tape[tmp.va].oc=OP;                                                                        \
  tape[tmp.va].arg1=x1.va;                                                                   \
  tape[tmp.va].arg2=UNDEF;                                                                   \

#define ACTIVE_ACTIVE(OP)                                                                    \
  tmp.va=tape_entry::vac++;                                                                  \
  tape[tmp.va].v=tmp.v;                                                                      \
  tape[tmp.va].d=0;                                                                          \
  tape[tmp.va].oc=OP;                                                                        \
  tape[tmp.va].arg1=x1.va;                                                                   \
  tape[tmp.va].arg2=x2.va;                                                                   \

#define ACTIVE(OP)                                                                           \
  tmp.va=tape_entry::vac++;                                                                  \
  tape[tmp.va].v=tmp.v;                                                                      \
  tape[tmp.va].d=0;                                                                          \
  tape[tmp.va].oc=OP;                                                                        \
  tape[tmp.va].arg1=x.va;                                                                    \
  tape[tmp.va].arg2=UNDEF;                                                                   \

#define OPEQ(OP)                                                                             \
  int lhsva=va;                                                                              \
  va=tape_entry::vac++;                                                                      \
  tape[va].v=v;                                                                              \
  tape[va].d=0;                                                                              \
  tape[va].oc=OP;                                                                            \
  tape[va].arg1=lhsva;                                                                       \
  tape[va].arg2=x.va;                                                                        \

#if (DEBUG_ACTIVE == 1)
extern int a_ctr;
extern int ai_ctr;
extern int al_ctr;
extern int ad_ctr;
extern int ald_ctr;
extern int aa_ctr;
extern int aeqa_ctr;
extern int aeqp_ctr;
extern int nsign_ctr;
// extern int aeqp_ctr;
extern int ainc_ctr;
extern int adec_ctr;
extern int atimesa_ctr;
extern int ptimesa_ctr;
extern int atimesp_ctr;
extern int adiva_ctr;
extern int pdiva_ctr;
// extern int adivp_ctr;
extern int aplusa_ctr;
extern int pplusa_ctr;
extern int aplusp_ctr;
extern int aminusa_ctr;
extern int pminusa_ctr;
extern int aminusp_ctr;
extern int atimeseqa_ctr;
extern int atimeseqp_ctr;
extern int adiveqa_ctr;
extern int adiveqp_ctr;
extern int apluseqa_ctr;
extern int apluseqp_ctr;
extern int aminuseqa_ctr;
extern int aminuseqp_ctr;
extern int sin_ctr;
extern int asin_ctr;
extern int sinh_ctr;
extern int asinh_ctr;
extern int cos_ctr;
extern int acos_ctr;
extern int cosh_ctr;
extern int acosh_ctr;
extern int tan_ctr;
extern int atan_ctr;
extern int tanh_ctr;
extern int atanh_ctr;
extern int aatan2a_ctr;
extern int aatan2p_ctr;
extern int patan2a_ctr;
extern int afmaxa_ctr;
extern int afmaxp_ctr;
extern int pfmaxa_ctr;
extern int amaxa_ctr;
extern int amaxp_ctr;
extern int pmaxa_ctr;
extern int afmina_ctr;
extern int afminp_ctr;
extern int pfmina_ctr;
extern int amaxa_ctr;
extern int amaxp_ctr;
extern int pmaxa_ctr;
extern int amina_ctr;
extern int aminp_ctr;
extern int pmina_ctr;
extern int exp_ctr;
extern int apowa_ctr;
extern int apowp_ctr;
extern int ppowa_ctr;
extern int sqrt_ctr;
extern int ahypota_ctr;
extern int ahypotp_ctr;
extern int phypota_ctr;
extern int log_ctr;
extern int log10_ctr;
extern int log1p_ctr;
extern int fabs_ctr;
extern int abs_ctr;
extern int aldexpa_ctr;
extern int aldexpp_ctr;
extern int pldexpa_ctr;
extern int afrexpi_ctr;
extern int floor_ctr;
extern int ceil_ctr;
extern int setv_ctr;
#endif

using namespace std;
extern bool TAPE_MODE;
extern bool DERIV_MODE;
//extern int TAPE_SIZE;
// class tape_entry;
// extern tape_entry* tape;

/* AMPI */



/**
  * augmented double type for computation of directional
  * derivatives by forward mode automatic differentiation
  */
class active {
  public :
  /**
    * function value
    */
  double v;
  /**
    * value of directional derivative
    */
  double d;
  /**
    * tape index of active
    */
  int va;
  active(const int &);
  active(const double &);
  active(const long double &);
  active(const long int &);
  active(const active &);
  active();
  active & operator=(const double &);
  active & operator=(const active &);
  active & operator--();
  active & operator++();
  active & operator+=(const double &);
  active & operator+=(const active &);
  active & operator-=(const double &);
  active & operator-=(const active &);
  active & operator*=(const double &);
  active & operator*=(const active &);
  active & operator/=(const double &);
  active & operator/=(const active &);
  const double& setValue(const double&);
};

/* AMPI */
//typedef struct AMPI_dco_Request{
    //AMPI_Request r;
    //active *buf;
    //int va;
//}AMPI_dco_Request;

class tape_entry {
  public:
    static int vac;
    static vector<int> indeps;
    static vector<int> deps;
    int oc;
    int arg1;
    int arg2;
    double v;
    double d;
    /**
     * AMPI
     */
    //AMPI_dco_Request *request;
    #if (DEBUG_ACTIVE == 3)
    bool ds;
    #endif
    tape_entry() : oc(UNDEF) 
                   , arg1(UNDEF)
                   , arg2(UNDEF)
                   , v(0) 
                   , d(0)
                   #if (DEBUG_ACTIVE == 3)
                   , ds(0)
                   #endif 
                   {};
    inline void reset(double x=0) {
      oc=UNDEF;
      arg1=UNDEF;
      arg2=UNDEF;
      v=x;
      d=0;
      #if (DEBUG_ACTIVE == 3)
      ds=0; 
      #endif
    }
    inline int get_size() {
      int size = 3*sizeof(int) + 2*sizeof(double);
      #if (DEBUG_ACTIVE == 3)
      size += sizeof(bool);
      #endif
      return size;
    }
};

extern tape_entry* tape;
  /**
    * active()
    */

    inline active::active() : v(0), d(0), va(UNDEF) {
#if (DEBUG_ACTIVE == 1)
	a_ctr++;
#endif
    };
  /**
    * active(int)
    */
  inline active::active(const int& x) : v(x), d(0), va(UNDEF) {
    if (DERIV_MODE) {
      if (TAPE_MODE) {
        va=tape_entry::vac++;
        tape[va].v=v;
        tape[va].d=0;
        tape[va].oc=CONST;
        tape[va].arg1=UNDEF;
        tape[va].arg2=UNDEF;
      }
    }
    #if (DEBUG_ACTIVE == 1)
    assert(tape_entry::vac < TAPE_SIZE);
    ai_ctr++;
    #endif
  }
  /**
    * active(long)
    */
  inline active::active(const long& x) : v(double(x)), d(0), va(UNDEF) {
    if (DERIV_MODE) {
      if (TAPE_MODE) {
        va=tape_entry::vac++;
        tape[va].v=v;
        tape[va].d=0;
        tape[va].oc=CONST;
        tape[va].arg1=UNDEF;
        tape[va].arg2=UNDEF;
      }
    }
    #if (DEBUG_ACTIVE == 1)
    assert(tape_entry::vac < TAPE_SIZE);
    al_ctr++;
    #endif
  }
  /**
    * active(double)
    */
  inline active::active(const double& x) : v(x), d(0), va(UNDEF) {
    if (DERIV_MODE) {
      if (TAPE_MODE) {
        va=tape_entry::vac++;
        tape[va].v=v;
        tape[va].d=0;
        tape[va].oc=CONST;
        tape[va].arg1=UNDEF;
        tape[va].arg2=UNDEF;
      }
    }
    #if (DEBUG_ACTIVE == 1)
    assert(++tape_entry::vac<TAPE_SIZE);
    ad_ctr++;
    #endif
  }   
  /**
    * active(double)
    */
  inline active::active(const long double& x) : v (double(x)), d(0), va(UNDEF) {
    if (DERIV_MODE) {
      if (TAPE_MODE) {
        va=tape_entry::vac++;
        tape[va].v=v;
        tape[va].d=0;
        tape[va].oc=CONST;
        tape[va].arg1=UNDEF;
        tape[va].arg2=UNDEF;
      }
    }
    #if (DEBUG_ACTIVE == 1)
    assert(++tape_entry::vac<TAPE_SIZE);
    ald_ctr++;
    #endif
  }
  /**
    * active(active)
    */
  inline active::active(const active& x) : v(x.v), d(x.d), va(UNDEF) {
    if (DERIV_MODE) {
      if (TAPE_MODE && x.va != UNDEF) {
        va=tape_entry::vac++;
        tape[va].v=v;
        tape[va].d=0;
//         if (x.va == UNDEF)  tape[va].oc=CONST;
//         else tape[va].oc=COPY;
        tape[va].oc=COPY;
        tape[va].arg1=x.va;
        tape[va].arg2=UNDEF;
      }
    }
    #if (DEBUG_ACTIVE == 1)
    assert(++tape_entry::vac<TAPE_SIZE);
    aa_ctr++;
    #endif
  }
  /**
    * active = passive  
    */
  inline active& active::operator=(const double& x) {
    if (isnan(x)) v=0;
    else v=x;
    if (TAPE_MODE) {
      va=tape_entry::vac++;
      tape[va].v=v;
      tape[va].d=0;
      tape[va].oc=CONST;
      tape[va].arg1=UNDEF;
      tape[va].arg2=UNDEF;
    }
    else d = 0;
    #if (DEBUG_ACTIVE == 1)
    assert(++tape_entry::vac<TAPE_SIZE);
    aeqp_ctr++;
    #endif
    return *this;
  }
  /**
    * active = active
    */
  inline active& active::operator=(const active& x) {
    v=x.v;
    if (DERIV_MODE) {
      if (TAPE_MODE) {
        va=tape_entry::vac++;
        tape[va].v=v;
        tape[va].d=0;
        if (x.va == UNDEF)  tape[va].oc=CONST;
        else tape[va].oc=ASG;
        tape[va].arg1=x.va;
        tape[va].arg2=UNDEF;
      }
      else d = x.d;
    }
    #if (DEBUG_ACTIVE == 1)
    assert(++tape_entry::vac<TAPE_SIZE);
    if (this==&x) return *this;
    aeqa_ctr++;
    #endif
    return *this;
  }
  /**
    * active ++ 
    */
  inline active& active::operator++() {
    ++v;
    #if (DEBUG_ACTIVE == 1)
    ainc_ctr++;
    #endif        
    return *this;
  }
  /**
    * active --
    */
  inline active& active::operator--() {
    --v;
    #if (DEBUG_ACTIVE == 1)
    adec_ctr++;
    #endif 
    return *this;
  }
  /**
    * active += double
    */
  inline active& active::operator+=(const double& x) {
    v+=x;
    if (DERIV_MODE) {
      if (TAPE_MODE) {
        tape[va].v=v;
      }
    }
    #if (DEBUG_ACTIVE == 1)
    apluseqp_ctr++;
    #endif
    return *this;
  }
  /**
    * active += active
    */
  inline active& active::operator+=(const active& x) {
    v+=x.v;
    if (DERIV_MODE) {
      if (TAPE_MODE) {
        if (va != UNDEF && x.va == UNDEF)  tape[va].v=v;
        else if (va != UNDEF && x.va != UNDEF) {OPEQ(PLUSEQ);}
      }
      else d += x.d;
    }

    #if (DEBUG_ACTIVE == 1)
    assert(++tape_entry::vac<TAPE_SIZE);
    if (this==&x) assert(false);
    apluseqa_ctr++;
    #endif
    return *this;
  }
  /**
    * active -= double
    */
  inline active& active::operator-=(const double& x) {
    v-=x;
    if (DERIV_MODE) {    
      if (TAPE_MODE) {
        tape[va].v=v;
      }
    }
    #if (DEBUG_ACTIVE == 1)
    aminuseqp_ctr++;
    #endif
    return *this;
  }
  /**
    * active -= active
    */
  inline active& active::operator-=(const active& x) {
    v-=x.v;
    if (DERIV_MODE) {
      if (TAPE_MODE) {
        if (va != UNDEF && x.va == UNDEF)  tape[va].v=v;
        else if (va != UNDEF && x.va != UNDEF) {OPEQ(MINUSEQ);}
      }
      else d -= x.d;
    }
    #if (DEBUG_ACTIVE == 1)
    assert(++tape_entry::vac<TAPE_SIZE);

    aminuseqa_ctr++;
    #endif
    return *this;
  }
  /**
    * active *= double
    */
  inline active& active::operator*=(const double& x) {
    v*=x;
    if (DERIV_MODE) {
      if (TAPE_MODE) {
        tape[va].v=v;
      }
    }
    #if (DEBUG_ACTIVE == 1)
    atimeseqp_ctr++;
    #endif
    return *this;    
  }
  /**
    * active *= active
    */
  inline active& active::operator*=(const active& x) {
    double c_v=v;
    v*=x.v;
    if (DERIV_MODE) {
      if (TAPE_MODE) {
        if (va != UNDEF && x.va == UNDEF)  tape[va].v=v;
        else if (va != UNDEF && x.va != UNDEF) {OPEQ(MULTEQ);}
      }
      else d = x.d*c_v+d*x.v;
    }
    #if (DEBUG_ACTIVE == 1)
    assert(++tape_entry::vac<TAPE_SIZE);

    atimeseqa_ctr++;
    #endif
    return *this;    
  }
  /**
    * active /= double
    */
  inline active& active::operator/=(const double& x) {
    v/=x;
    if (DERIV_MODE) {    
      if (TAPE_MODE) {
        tape[va].v=v;
      }
    }
    #if (DEBUG_ACTIVE == 1)
    adiveqp_ctr++;
    #endif
    return *this;
  }
  /**
    * active /= active
    */
  inline active& active::operator/=(const active& x) {
    double c_v = v;
    v/=x.v;
    if (DERIV_MODE) {
      if (TAPE_MODE) {
        if (va != UNDEF && x.va == UNDEF)  tape[va].v=v;
        else if (va != UNDEF && x.va != UNDEF) {OPEQ(DIVEQ);}
      }
      else d = d/x.v-(c_v*x.d)/(x.v*x.v);
    }
    #if (DEBUG_ACTIVE == 1)
    assert(++tape_entry::vac<TAPE_SIZE);
    
    adiveqa_ctr++;
    #endif
    return *this;
  }
  /**
    * sets the value of this active
    */
  inline const double& active::setValue(const double& val) {
    v=val;
    if (DERIV_MODE) {
      if (TAPE_MODE) {
        if (va == UNDEF)  va=tape_entry::vac++;
        tape[va].oc=CONST;
        tape[va].v=v;
        tape[va].d=0;
        tape[va].arg1=UNDEF;
        tape[va].arg2=UNDEF;        
      }
    }
    #if (DEBUG_ACTIVE == 1)
    assert(++tape_entry::vac<TAPE_SIZE);
    setv_ctr++;
    #endif
    return v;    
  }


/**
  * active == != < <= > >= active
  */
inline bool operator==(const active& x1, const active& x2) {return x1.v == x2.v;}
inline bool operator!=(const active& x1, const active& x2) {return x1.v != x2.v;}
inline bool operator<(const active& x1, const active& x2) {return x1.v < x2.v;}
inline bool operator<=(const active& x1, const active& x2) {return x1.v <= x2.v;}
inline bool operator > (const active& x1, const active& x2) {return x1.v > x2.v;}
inline bool operator >= (const active& x1, const active& x2) {return x1.v >= x2.v;}
/**
 * active = double * active
 */
inline active operator*(const double& x1, const active& x2) {
  active tmp;
  tmp.v=x1*x2.v;
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x2.va != UNDEF) {PASSIVE_ACTIVE(MUL);}
    }
    else  tmp.d = x1*x2.d;
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  ptimesa_ctr++;
  #endif
  return tmp;
}
/**
 * active = active * double
 */
inline active operator*(const active& x1, const double& x2) {
  active tmp;
  tmp.v=x1.v*x2;
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x1.va != UNDEF) {ACTIVE_PASSIVE(MUL);}
    }
    else  tmp.d = x1.d*x2;
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  atimesp_ctr++;
  #endif
  return tmp;
}
/**
 * active = active * active
 */
inline active operator*(const active& x1, const active& x2) {
  active tmp;
  tmp.v=x1.v*x2.v;
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if ( x1.va != UNDEF && x2.va != UNDEF) {ACTIVE_ACTIVE(MUL);}
      else if ( x1.va != UNDEF && x2.va == UNDEF) {ACTIVE_PASSIVE(MUL);}
      else if ( x1.va == UNDEF && x2.va != UNDEF) {PASSIVE_ACTIVE(MUL);}
    }
    else  tmp.d = x1.d*x2.v + x1.v*x2.d;
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  atimesa_ctr++;
  #endif
  return tmp;
}
/**
 * active = double / active
 */
inline active operator/(const double x1, const active& x2) {
  active tmp;
  tmp.v=x1/x2.v;
  if (DERIV_MODE) {
    if(TAPE_MODE) {
      if (x2.va != UNDEF) {PASSIVE_ACTIVE(DIV);}
    }
    else tmp.d = -x1/(x2.v*x2.v)*x2.d;
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  pdiva_ctr++;
  #endif
  return tmp;
}
/**
 * active = active / double
 * Numerical instability by gaining 1/x2 from the result using x1,
 * hence commendted out here
 */
// inline active operator/(const active& x1, const double x2) {
//   assert(++tape_entry::vac<TAPE_SIZE);  assert(x2 != 0);
//   active tmp;
//   tmp.v=x1.v/x2;
//   if(TAPE_MODE) {
//     ACTIVE_PASSIVE(DIV);
//   }
//   else tmp.d = (1/x2)*x1.d;
//   #if (DEBUG_ACTIVE == 1)
//   
//   adivp_ctr++;
//   #endif
//   return tmp;
// }
/**
 * active = active / active
 */
inline active operator/(const active& x1, const active& x2) {
  active tmp;
  tmp.v=x1.v/x2.v;
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x1.va != UNDEF && x2.va != UNDEF) {ACTIVE_ACTIVE(DIV);}
      else if (x2.va != UNDEF) {PASSIVE_ACTIVE(DIV);}
    }
    else tmp.d = x1.d/x2.v - x1.v/(x2.v*x2.v)*x2.d;
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);  
  assert( x2 != 0);
  
  adiva_ctr++;
  #endif
  return tmp;
}
/**
 * active = double + active
 */
inline active operator+(const double x1, const active& x2) {
  active tmp;
  tmp.v=x1+x2.v;
  if (DERIV_MODE) {
    if(TAPE_MODE) {
      if (x2.va != UNDEF) {PASSIVE_ACTIVE(ADD);}
    }
    else tmp.d = x2.d;
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  pplusa_ctr++;
  #endif
  return tmp;
}
/**
 * active = active + double
 */
inline active operator+(const active& x1, const double x2) {
  active tmp;
  tmp.v=x1.v+x2;
  if (DERIV_MODE) {
    if(TAPE_MODE) {
      if (x1.va != UNDEF) {ACTIVE_PASSIVE(ADD);}
    }
    else tmp.d = x1.d;
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  aplusp_ctr++;
  #endif
  return tmp;
}
/**
 * active = active + active
 */
inline active operator+(const active& x1, const active& x2) {
  active tmp;
  tmp.v=x1.v+x2.v;
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if ( x1.va != UNDEF && x2.va != UNDEF) {ACTIVE_ACTIVE(ADD);}
      else if ( x1.va != UNDEF && x2.va == UNDEF) {ACTIVE_PASSIVE(ADD);}
      else if ( x1.va == UNDEF && x2.va != UNDEF) {PASSIVE_ACTIVE(ADD);}
    }
    else tmp.d = x1.d+x2.d;
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  aplusa_ctr++;
  #endif
  return tmp;
}
/**
 * active = double - active
 */
inline active operator-(const double x1, const active& x2) {
  active tmp;
  tmp.v=x1-x2.v;
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x2.va != UNDEF) {PASSIVE_ACTIVE(SUB);}
    }
    else tmp.d = -x2.d;
  }
  #if (DEBUG_ACTIVE == 1)
  
  pminusa_ctr++;
  #endif
  return tmp;
}
/**
 * active = active - double
 */
inline active operator-(const active& x1, const double x2) {
  active tmp;
  tmp.v=x1.v-x2;
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x1.va != UNDEF) {ACTIVE_PASSIVE(SUB);}
    }
    else tmp.d = x1.d;
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  aminusp_ctr++;
  #endif
  return tmp;
}
/**
 * active = active - active
 */
inline active operator-(const active& x1, const active& x2) {
  active tmp;
  tmp.v=x1.v-x2.v;
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if ( x1.va != UNDEF && x2.va != UNDEF) {ACTIVE_ACTIVE(SUB);}
      else if (x1.va != UNDEF) {ACTIVE_PASSIVE(SUB);}
      else if (x2.va != UNDEF) {PASSIVE_ACTIVE(SUB);}
    }
    else tmp.d = x1.d-x2.d;
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  aminusa_ctr++;
  #endif
  return tmp;
}
/**
 * active = -active
 */
inline active operator-(const active& x) {
  active tmp;
  tmp.v=-x.v;
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x.va != UNDEF) {ACTIVE(NEGSIGN);}
    }
    else tmp.d = -x.d;
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  nsign_ctr++;
  #endif
  return tmp;
}
/**
 * active = sin(active)
 */
inline active sin(const active& x) {
  active tmp;
  tmp.v=sin(x.v);
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x.va != UNDEF) {ACTIVE(SIN);}
    }
    else tmp.d = cos(x.v)*x.d;
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  sin_ctr++;
  #endif
  return tmp;
}
/**
 * active = asin(active)
 */
inline active asin(const active& x) {
  active tmp;
  tmp.v = asin(x.v);
  if (DERIV_MODE) {
    if (TAPE_MODE){
      if (x.va != UNDEF) {ACTIVE(ASIN);}
    }
    else  {
      double tmp2 = 1.0-x.v*x.v;
      #if (DEBUG_ACTIVE == 1)
      assert(tmp2 > 0);
      #endif
      tmp.d = x.d/sqrt(tmp2);
    }     
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  asin_ctr++;
  #endif        
  return tmp;
}
/**
 * active = sinh(active)
 */
inline active sinh (const active &x) {
  active tmp;
  tmp.v = sinh(x.v);
  if (DERIV_MODE) {
    if (TAPE_MODE){
      if (x.va != UNDEF) {ACTIVE(SINH);}
    }
    else tmp.d = x.d*cosh(x.v);
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  sinh_ctr++;
  #endif        
  return tmp;
}
/**
 * active = asinh(active)
 */
inline active asinh(const active& x) {
  active tmp;
  tmp.v = asinh(x.v);
  if (DERIV_MODE) {
     if (TAPE_MODE){
      if (x.va != UNDEF) {ACTIVE(ASINH);}
    }
    else {
      double tmp2 = x.v*x.v+1.0;
      #if (DEBUG_ACTIVE == 1)
      assert(tmp2 > 0);
      #endif
      tmp.d=x.d/sqrt(tmp2);
    }
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE); 
  
  asinh_ctr++;
  #endif        
  return tmp;
}
/**
 * active = cos(active)
 */
inline active cos(const active& x) {
  active tmp;
  tmp.v=cos(x.v);
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x.va != UNDEF) {ACTIVE(COS);}
    }
    else tmp.d = -sin(x.v)*x.d;
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  cos_ctr++;
  #endif        
  return tmp;
}
/**
 * active = cosh(active)
 */
inline active acos(const active& x) {
  active tmp;
  tmp.v=acos(x.v);
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x.va != UNDEF) {ACTIVE(ACOS);}    
    }
    else {
        double tmp2 = 1.0-x.v*x.v;
        #if (DEBUG_ACTIVE == 1)
        assert(tmp2 > 0);
        #endif
        tmp.d = -x.d/sqrt(tmp2);
    }
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  acos_ctr++;
  #endif        
  return tmp;
}
/**
 * active = acos(active)
 */
inline active cosh(const active& x) {
  active tmp;
  tmp.v=cosh(x.v);
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x.va != UNDEF) {ACTIVE(COSH);}
    }
    else tmp.d = x.d*sinh(x.v);
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  cosh_ctr++;
  #endif
  return tmp;
}
/**
 * active = acosh(active) 
 */
inline active acosh (const active &x) {
  active tmp;
  tmp.v=acosh(x.v);
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x.va != UNDEF) {ACTIVE(ACOSH);}
    }
    else {
      double tmp2 = x.v*x.v-1.0;
      #if (DEBUG_ACTIVE == 1)
      assert(tmp2 > 0);
      #endif
      tmp.d = x.d/sqrt(tmp2);
    }
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  acosh_ctr++;
  
  #endif
  return tmp;
}
/**
 * active = tan(active) 
 */
inline active tan(const active& x) {
  active tmp;
  tmp.v=tan(x.v);
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x.va != UNDEF) {ACTIVE(TAN);}
    }
    else {
      double tmp2 = cos(x.v)*cos(x.v);
      #if (DEBUG_ACTIVE == 1)
      assert(tmp2 != 0);
      #endif
      tmp.d = x.d/tmp2;
    }
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  tan_ctr++;
  
  #endif
  return tmp;
}
/**
 * active = atan(active)
 */
inline active atan(const active& x) {
  active tmp;
  tmp.v=atan(x.v);
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x.va != UNDEF) {ACTIVE(ATAN);}
    }
    else {
      double tmp2 = 1.0+x.v*x.v; 
      #if (DEBUG_ACTIVE == 1)
      assert(tmp2 != 0);
      #endif
      tmp.d=x.d/tmp2;
    }
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  atan_ctr++;
  #endif
  return tmp;
}
/**
 * active = tanh(active)
 */
inline active tanh (const active &x) {
  active tmp;
  tmp.v=tanh(x.v);
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x.va != UNDEF) {ACTIVE(TANH);}
    }
    else {
        double tmp2 = cosh(x.v)*cosh(x.v);
        #if (DEBUG_ACTIVE == 1)
        assert(tmp2 != 0);
        #endif
        tmp.d = x.d/tmp2;
    }
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  tanh_ctr++;
  #endif
  return tmp;
}
/**
 * active = atanh(active)
 */
inline active atanh (const active &x) {
  active tmp;
  tmp.v=atanh(x.v);
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x.va != UNDEF) {ACTIVE(ATANH);}
    }
    else {
      double tmp2 = (1.0-x.v*x.v);
      #if (DEBUG_ACTIVE == 1)
      assert(tmp2 != 0);
      #endif
      tmp.d = x.d/tmp2;
    }
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  atanh_ctr++;
  #endif
  return tmp;
}
/**
 * active = atan2(passive, active)
 */
inline active atan2(const double& x1, const active& x2) {
  #if (DEBUG_ACTIVE == 1)
  patan2a_ctr++;
  atan_ctr--;
  #endif
  return atan(x1/x2);
}
/**
 * active = atan2(active, passive)
 */
inline active atan2(const active& x1, const double& x2) {
  #if (DEBUG_ACTIVE == 1)
  aatan2p_ctr++;
  atan_ctr--;
  #endif
  return atan(x1/x2);
}
/**
 * active = atan2(active,active)
 */
inline active atan2(const active& x1, const active& x2) {
  #if (DEBUG_ACTIVE == 1)
  aatan2a_ctr++;
  atan_ctr--;
  #endif
  return atan(x1/x2);
}
/**
 * active = exp(active)
 */
inline active exp(const active& x) {
  active tmp;
  tmp.v=exp(x.v);
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x.va != UNDEF) {ACTIVE(EXP);}
    }
    else  tmp.d = x.d*tmp.v;
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  exp_ctr++;
  #endif
  return tmp;
}
/**
 * active = log(active)
 */
inline active log(const active& x) {
  active tmp;
  tmp.v=log(x.v);
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x.va != UNDEF) {ACTIVE(LOG);}
    }
    else  tmp.d = x.d/x.v;
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  assert(x.v > 0);
  
  log_ctr++;
  #endif
  return tmp;
}
/**
 * active = log1p(active)
 */
inline active log1p(const active& x) {
  active tmp;
  tmp.v=log1p(x.v);
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x.va != UNDEF) {ACTIVE(LOG1P);}
    }
    else  tmp.d = x.d/(x.v+1.0);
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  assert(x.v >= 0);
  
  log1p_ctr++;
  #endif
  return tmp;
}
/**
 * active = log10(active) (NOT USED IN JURASSIC)
 */
inline active log10(const active& x) {
  active tmp;
  tmp.v=log10(x.v);
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x.va != UNDEF) {ACTIVE(LOG10);}
    }
    else {
        double tmp2=log(10.0)*x.v;
        #if (DEBUG_ACTIVE == 1)
        assert(tmp2 != 0);
        #endif
        tmp.d=x.d/tmp2;
    }
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  log10_ctr++;
  #endif
  return tmp;
}
/**
 * active = pow(double, active)
 */
inline active pow(const double& x1, const active& x2) {
  active tmp;
  tmp.v=pow(x1, x2.v);
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x2.va != UNDEF) {PASSIVE_ACTIVE(POW);}
    }
    else tmp.d = x2.d*log(x1)*tmp.v;
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  ppowa_ctr++;
  #endif
  return tmp;
}
/**
 * active = pow(active, double)
 */
inline active pow(const active& x1, const double& x2) {
  active tmp;
  tmp.v=pow(x1.v, x2);
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x1.va != UNDEF) {ACTIVE_PASSIVE(POW);}
    }
    else tmp.d = x1.d*x2*pow(x1.v, x2-1);
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  apowp_ctr++;
  #endif
  return tmp;
}
/**
 * active = pow(active, active)
 */
inline active pow(const active& x1, const active& x2) {
  active tmp;
  tmp.v=pow(x1.v, x2.v);
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x1.va != UNDEF && x2.va != UNDEF) {ACTIVE_ACTIVE(POW);}
      else if (x1.va != UNDEF) {ACTIVE_PASSIVE(POW);}
      else if (x2.va != UNDEF) {PASSIVE_ACTIVE(POW);}
    }
    else tmp.d = x1.d*(x2.v*pow(x1.v, x2.v-1)) + x2.d*(log(x1.v)*tmp.v);
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  apowa_ctr++;
  #endif
  return tmp;
}
/**
 * active = ldexp(double, active) (NOT USED IN JURASSIC)
 */
inline active ldexp (const double &x1, const active &x2) {
  #if (DEBUG_ACTIVE == 1)
  pldexpa_ctr++;
  ptimesa_ctr--;
  ppowa_ctr--;
  #endif
  return x1*pow(2.0, x2);
}
/**
 * active = ldexp(active, double) (NOT USED IN JURASSIC)
 */
inline active ldexp (const active &x1, const double &x2) {
  #if (DEBUG_ACTIVE == 1)
  aldexpp_ctr++;
  atimesp_ctr--;
  #endif
  return x1*pow(2.0, x2);
}
/**
 * active = ldexp(active, active) (NOT USED IN JURASSIC)
 */
inline active ldexp (const active &x1, const active &x2) {
  #if (DEBUG_ACTIVE == 1)
  aldexpa_ctr++;
  atimesa_ctr--;
  ppowa_ctr--;
  #endif
  return x1*pow(2.0, x2);
}
/**
 * active = frexp(active, active) (NOT USED IN JURASSIC)
 */
inline active frexp (const active &x1, int& x2) {
  #if (DEBUG_ACTIVE == 1)
  afrexpi_ctr++;
  #endif
  return frexp(x1.v, x2);
}
/**
 * active = sqrt(active)
 */
inline active sqrt(const active& x) {
  active tmp;
  tmp.v=sqrt(x.v);
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x.va != UNDEF) {ACTIVE(SQRT);}
    }
    else tmp.d = x.d/(2.0*tmp.v);
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  assert(x.v > 0);
  
  sqrt_ctr++;
  #endif
  return tmp;
}
/**
 * active = fabs(active)
 */
inline active fabs(const active& x) {
  active tmp;
  tmp.v=fabs(x.v);
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x.va != UNDEF) {ACTIVE(FABS);}
    }
    else {
      if (x.v == 0)  tmp.d = 0.;
      else if (x.v < 0)  tmp.d = -x.d;
      else  tmp.d = x.d;
    }
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  fabs_ctr++;
  #endif
  return tmp;
}
/**
 * active = abs(active)
 */
inline active abs(const active& x) {
  #if (DEBUG_ACTIVE == 1)
  abs_ctr++;
  fabs_ctr--;
  #endif
  return fabs(x);
}
/**
 * active = max(double, active)
 */
inline active max(const double& x1, const active& x2) {
  active tmp;
  tmp.v=max(x1,x2.v);
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x2.va != UNDEF) {PASSIVE_ACTIVE(MAX);}
    }
    else {
      if (x1 >= x2.v)  tmp.d = 0;
      else  tmp.d = x2.d;
    }
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  pmaxa_ctr++;
  #endif
  return tmp;
}
/**
 * active = max(active, passive)
 */
inline active max(const active& x1, const double& x2) {
  active tmp;
  tmp.v=max(x1.v,x2);
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x1.va != UNDEF) {ACTIVE_PASSIVE(MAX);}
    }
    else {
      if (x1.v >= x2)  tmp.d = x1.d;
      else  tmp.d = 0;
    }
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  amaxp_ctr++;
  #endif
  return tmp;
}
/**
 * active = max(active, active)
 */
inline active max(const active& x1, const active& x2) {
  active tmp;
  tmp.v=max(x1.v,x2.v);
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x1.va != UNDEF && x2.va != UNDEF) {ACTIVE_ACTIVE(MAX);}
      else if (x1.va != UNDEF) {ACTIVE_PASSIVE(MAX);}
      else if (x2.va != UNDEF) {PASSIVE_ACTIVE(MAX);} 
    }
    else {
      if (x1.v >= x2.v)  tmp.d = x1.d;
      else  tmp.d = x2.d;
    }
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  amaxa_ctr++;
  #endif
  return tmp;
}
/**
 * active = min(passive, active)
 */
inline active min(const double& x1, const active& x2) {
  active tmp;
  tmp.v=min(x1,x2.v);
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x2.va != UNDEF) {PASSIVE_ACTIVE(MIN);}
    }
    else {
      if (x1 < x2.v)  tmp.d = 0;
      else  tmp.d = x2.d;
    }
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  pmina_ctr++;
  #endif
  return tmp;
}
/**
 * active = min(active, passive)
 */
inline active min(const active& x1, const double& x2) {
  active tmp;
  tmp.v=min(x1.v,x2);
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x1.va != UNDEF) {ACTIVE_PASSIVE(MIN);}
    }
    else {
      if (x1.v < x2)  tmp.d = x1.d;
      else  tmp.d = 0;
    }
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  aminp_ctr++;
  #endif
  return tmp;
}
/**
 * active = min(active, active)
 */
inline active min(const active& x1, const active& x2) {
  active tmp;
  tmp.v=min(x1.v,x2.v);
  if (TAPE_MODE) {
    if (x1.va != UNDEF && x2.va != UNDEF) {ACTIVE_ACTIVE(MIN);}
    else if (x1.va != UNDEF) {ACTIVE_PASSIVE(MIN);}
    else if (x2.va != UNDEF) {PASSIVE_ACTIVE(MIN);} 
  }
  else {
    if (x1.v < x2.v)  tmp.d = x1.d;
    else  tmp.d = x2.d;
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  amina_ctr++;
  #endif
  return tmp;
}
/**
 * active = fmax(passive, active)
 */
inline active fmax(const double& x1, const active& x2) {
  #if (DEBUG_ACTIVE == 1)
  pfmaxa_ctr++;
  pmaxa_ctr--;
  #endif
  return max(x1, x2);
}
/**
 * active = fmax(active, passive)
 */
inline active fmax(const active& x1, const double& x2) {
  #if (DEBUG_ACTIVE == 1)
  afmaxp_ctr++;
  amaxp_ctr--;
  #endif
  return max(x1, x2);
}
/**
 * active = fmax(active, active)
 */
inline active fmax(const active& x1, const active& x2) {
  #if (DEBUG_ACTIVE == 1)
  afmaxa_ctr++;
  amaxa_ctr--;
  #endif
  return max(x1, x2);
}
/**
 * active = fmin(passive, active)
 */
inline active fmin(const double& x1, const active& x2) {
  #if (DEBUG_ACTIVE == 1)
  pfmina_ctr++;
  pmina_ctr--;
  #endif
  return min(x1, x2);
}
/**
 * active = fmin(active, passive)
 */
inline active fmin(const active& x1, const double& x2) {
  #if (DEBUG_ACTIVE == 1)
  afminp_ctr++;
  aminp_ctr--;
  #endif
  return min(x1, x2);
}
/**
 * active = fmin(active, active)
 */
inline active fmin(const active& x1, const active& x2) {
  #if (DEBUG_ACTIVE == 1)
  afmina_ctr++;
  amina_ctr--;
  #endif
  return min(x1, x2);
}
/**
 * active = hypot(double, active)
 * sqrt(x1^2 + x2^2)
 */
inline active hypot(const double& x1, const active& x2) {
  active tmp;
  tmp.v=hypot(x1, x2.v);
  assert(tmp.v != 0);
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x2.va != UNDEF) {PASSIVE_ACTIVE(HYPOT);} 
    }
    else  tmp.d=(x2.v/tmp.v)*x2.d;
  }
  #if (DEBUG_ACTIVE == 1)
  assert(tmp.v != 0);
  assert(++tape_entry::vac<TAPE_SIZE);
  
  phypota_ctr++;
  #endif
  return tmp;
}
/**
 * active = hypot(active, double)
 */
inline active hypot(const active& x1, const double& x2) {
  active tmp;
  tmp.v=hypot(x1.v,x2);
  assert(tmp.v != 0);
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x1.va != UNDEF) {ACTIVE_PASSIVE(HYPOT);}
    }
    else  tmp.d=(x1.v/tmp.v)*x1.d;
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  assert(tmp.v != 0);
  
  ahypotp_ctr++;
  #endif
  return tmp;
}
/**
 * active = hypot(active, active)
 */
inline active hypot(const active& x1, const active& x2) {
  active tmp;
  tmp.v=hypot(x1.v,x2.v);
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x1.va != UNDEF && x2.va != UNDEF) {ACTIVE_ACTIVE(HYPOT);}
      else if (x1.va != UNDEF) {ACTIVE_PASSIVE(HYPOT);}
      else if (x2.va != UNDEF) {PASSIVE_ACTIVE(HYPOT);} 
    }
    else  tmp.d=(x1.v/tmp.v)*x1.d+(x2.v/tmp.v)*x2.d;
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  assert(tmp.v != 0);
  
  ahypota_ctr++;
  #endif
  return tmp;
}
/**
 * active = floor(active)
 */
inline active floor (const active& x) {
  active tmp;
  tmp.v=floor(x.v);
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x.va != UNDEF) {
        ACTIVE(FLOOR);
      }
    }
    else  tmp.d = 0.;
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE);
  
  floor_ctr++;
  #endif
  return tmp;  
}
/**
 * active = ceil(active)
 */
inline active ceil (const active &x) {
  active tmp;
  tmp.v=ceil(x.v);
  if (DERIV_MODE) {
    if (TAPE_MODE) {
      if (x.va != UNDEF) {
        ACTIVE(CEIL);
      }
    }
    else  tmp.d = 0.;
  }
  #if (DEBUG_ACTIVE == 1)
  assert(++tape_entry::vac<TAPE_SIZE); 
  
  ceil_ctr++;
  #endif

  return tmp;  
}
/**
 * returns pure value withour any derivative.
 */
inline const double& value(const active& val) {
  return val.v;
}
inline const double& value(const double& val) {
  return val;
}
/**
  * isnan(active)
  */
inline bool isnan(active a) {
  if (isnan(a.v) ||  isnan(a.d))  return true;
  return  false;
}
/**
 * enables activate evaluation per default on tape 
 * by setting both DERIV_MDOE to true, 
 * where TAPE_MODE is initially true.  
 */
inline void trace_on(bool on_tape=true, int tape_size=10000000) {
  #if (DEBUG_ACTIVE == 0)
  cout << "INFO::active.h::trace_on > Openning active session ..." << endl; 
  #endif

  if (!tape) {
      //TAPE_SIZE = tape_size;

    #if (DEBUG_ACTIVE == 0)
    cout << "  Tape allocation of size " <<  TAPE_SIZE  << " runs ..." << endl;
    #endif

    tape = new tape_entry[TAPE_SIZE];
  }

  DERIV_MODE = true;
  TAPE_MODE = on_tape;
}
/**
 * disabels the active evaluation by setting DERIV_MDOE to false 
 */
inline void trace_off() {
  DERIV_MODE = false;
  TAPE_MODE = false;

  #if (DEBUG_ACTIVE == 0)
  cout << "INFO::active.h:: trace_off < Closing active session" << endl; 
  #endif
}
/**
 * disabels the active evaluation by setting DERIV_MDOE to false 
 */
inline void destroy_tape() {
  if (tape) {
    #if (DEBUG_ACTIVE == 0)
    cout << "INFO (active.h::destroy_tape > tape deallocation of size " 
         <<  TAPE_SIZE << " runs ..." << endl;
    #endif
    delete [] tape;
    tape = 0;

    #if (DEBUG_ACTIVE == 0)
    cout << "INFO (active.h::destroy_tape < tape deallocation of size " 
         <<  TAPE_SIZE << " terminates!" << endl;
    #endif
  }
}
/**
* set ad_mode to 0 for tapeless forward and to 1 for tape based adjoint mode.
* set deriv_mode to 0 for computation of function values
* to 1 for propagation of derivative.
*/
inline void init_mode(bool ad_mode, bool deriv_mode) {
  TAPE_MODE = ad_mode;
  DERIV_MODE = deriv_mode;
  if (!TAPE_MODE) cout << "TAPELESS MODE AD" << endl;
  else cout << "TAPED MODE AD" << endl;
}
/**
 * marks input x as independent active variable
 */
inline void independent(active& x) {
  #if (DEBUG_ACTIVE > 0)
  cout << "INFO::active.h > independent(active& x(" << x.va << ", " << x.v << ")" << endl; 
  #endif

  if (TAPE_MODE) {
//     if (x.va == UNDEF) {
      x.va=tape_entry::vac++;
      tape[x.va].v = x.v;
      tape[x.va].d = 0;
      tape[x.va].arg1 = UNDEF;
      tape[x.va].arg2 = UNDEF;
//     }
    tape[x.va].oc=INDEP;   
    tape_entry::indeps.push_back(x.va);
  }

  #if (DEBUG_ACTIVE > 0)
  cout << "INFO::active.h < independent(active& x(" << x.va << ", " << x.v << ")" << endl; 
  #endif
}
inline void independent(active& x, std::vector<int>& indeps) {
  #if (DEBUG_ACTIVE > 0)
  cout << "INFO::active.h > independent(active& x(" << x.va << ", " << x.v << ")" << endl; 
  #endif

  if (TAPE_MODE) {
//     if (x.va == UNDEF) {
      x.va=tape_entry::vac++;
      tape[x.va].v = x.v;
      tape[x.va].d = 0;
      tape[x.va].arg1 = UNDEF;
      tape[x.va].arg2 = UNDEF;
//     }
    tape[x.va].oc=INDEP;
  }
  indeps.push_back(x.va);

  #if (DEBUG_ACTIVE > 0)
  cout << "INFO::active.h < independent(active& x(" << x.va << ", " << x.v << ")" << endl; 
  #endif
}
inline void independent(active& x, std::vector<active*>& indeps) {
  #if (DEBUG_ACTIVE > 0)
  cout << "INFO::active.h > independent(active& x(" << x.va << ", " << x.v << ")" << endl; 
  #endif

  if (TAPE_MODE) {
//     if (x.va == UNDEF) {
      x.va=tape_entry::vac++;
      tape[x.va].v = x.v;
      tape[x.va].d = 0;
      tape[x.va].arg1 = UNDEF;
      tape[x.va].arg2 = UNDEF;
//     }
    tape[x.va].oc=INDEP;
  }
  indeps.push_back(&x);

  #if (DEBUG_ACTIVE > 0)
  cout << "INFO::active.h < independent(active& x(" << x.va << ", " << x.v << ")" << endl; 
  #endif
}
/**
 * marks x as dependent active variable
 */
inline void dependent(active& x) {
  #if (DEBUG_ACTIVE > 0)
  cout << "dependent > va (" << x.va <<"),  v (" << x.v << "),   d (" << x.d << ")" << endl;
  #endif

  if (TAPE_MODE) {
//     tape[x.va].oc=DEP;
    tape_entry::deps.push_back(x.va);
  }

  #if (DEBUG_ACTIVE > 0)
  cout << "dependent < va (" << x.va <<"),  v (" << x.v << "),   d (" << x.d << ")" << endl;
  #endif
}
inline void dependent(active& x, std::vector<int>& deps) {
  #if (DEBUG_ACTIVE > 0)
  cout << "dependent > va (" << x.va <<"),  v (" << x.v << "),   d (" << x.d << ")" << endl;
  #endif

  if (TAPE_MODE) {
//     tape[x.va].oc=DEP;
    deps.push_back(x.va);
  }

  #if (DEBUG_ACTIVE > 0)
  cout << "dependent < va(" << x.va <<"),  d (" << x.d << ")" << endl;
  #endif
}
inline void dependent(const active& x, std::vector<int>& deps) {
  #if (DEBUG_ACTIVE > 0)
    cout << "dependent > va (" << x.va <<"),  v (" << x.v << "),   d (" << x.d << ")" << endl;
  #endif

  if (TAPE_MODE) {
//     tape[x.va].oc=DEP;
    deps.push_back(x.va);
  }

  #if (DEBUG_ACTIVE > 0)
  cout << "dependent < va(" << x.va <<"),  d (" << x.d << ")" << endl;
  #endif
}

inline void dependent(active& x, std::vector<active*>& deps) {
  #if (DEBUG_ACTIVE > 0)
  cout << "dependent > va (" << x.va <<"),  v (" << x.v << "),   d (" << x.d << ")" << endl;
  #endif

  if (TAPE_MODE) {
//     tape[x.va].oc=DEP;
    tape_entry::deps.push_back(x.va);
  }
  deps.push_back(&x);

  #if (DEBUG_ACTIVE > 0)
  cout << "dependent < va(" << x.va <<"),  d (" << x.d << ")" << endl;
  #endif
}
/**
 * set the derivative component of the tape entry corresponding to 
 * ith input to 1 and all others to 0.
 */
inline void seed_independent(const int& i) {
  for (int j=0; j<tape_entry::vac; j++)  tape[j].d=0;
  tape[tape_entry::indeps[i]].d=1;
}
inline void seed_independent(const int& i, std::vector<active*>& indeps) {
  for (int j=0; j<int(indeps.size()); j++)
    if (i==j)  indeps[j]->d=1;
    else  indeps[j]->d=0;
}
inline void seed_independent(const int& i, std::vector<int> indeps) {
  for (int j=0; j<tape_entry::vac; j++)
    if (indeps[i] == j) tape[j].d=1;
    else tape[j].d=0;
}
/**
 * set the derivative component of the tape entry corresponding to 
 * jth output to 1 and all others to 0.
 */
inline void seed_dependent(const int& i) {
  for (int j=0; j<tape_entry::vac; j++) tape[j].d=0;
  tape[tape_entry::deps[i]].d=1;
}
inline void seed_dependent(const int& i, const double& d) {
//   for (int j=0; j<tape_entry::vac; j++) tape[j].d=0;
  tape[tape_entry::deps[i]].d=d;
}
inline void seed_dependent(const active& x) {
  for (int j=0; j<tape_entry::vac; j++) tape[j].d=0;
  tape[x.va].d=1;
}
inline void seed_dependent(const int& i, std::vector<int>& deps) {
  for (int j=0; j<tape_entry::vac; j++) tape[j].d=0;
  tape[deps[i]].d=1;
}
inline void seed_dependent(const int& i, std::vector<active*>& deps) {
  for (int j=0; j<int(deps.size()); j++)
    if (i==j)  deps[j]->d=1;
    else  deps[j]->d=0;
}
/**
 * interprets the tape entry i in forward mode
 */
void forward_tape_interpreter(int i);
/**
 * forward tape interpretation
 */
inline void forward_tape_interpreter() {
  #if (DEBUG_ACTIVE == 2)
  cout << "INFO::(active.h) > Forward Tape Interpretation runs ..."
       << endl;
  #endif

  for (int i=0; i<tape_entry::vac; i++)
    forward_tape_interpreter(i);

  #if (DEBUG_ACTIVE == 2)
  cout << "INFO::(active.h) < Forward Tape Interpretation terminates ..."
       << endl;
  #endif
}
/**
 * interprets the tape entry i in reverse mode
 */
void reverse_tape_interpreter(int i);
/**
 * reverse tape interpretation
 */
inline void reverse_tape_interpreter() {
  #if (DEBUG_ACTIVE == 2)
  cout << "INFO::(active.h) > Reverse Tape Interpretation runs ..."
       << endl;
  #endif

  for (int i=tape_entry::vac-1; i>=0; i--)
    reverse_tape_interpreter(i);

  #if (DEBUG_ACTIVE == 2)
  cout << "INFO::(active.h) < Reverse Tape Interpretation terminates ..."
       << endl;
  #endif
}
/**
 * reverse tape interpretation with d of tape entries of actives y initialized to by
 */
inline void reverse_tape_interpreter(active* y, const std::vector<double> by, int len) {
  for (int i=0; i<tape_entry::vac; i++)  tape[i].d=0;

  for (int i=0; i<len; i++)  tape[y[i].va].d=by[i];

  for (int i=tape_entry::vac-1; i>=0; i--)   reverse_tape_interpreter(i);
}
/**
 * reverse tape interpretation with d of tape entries of actives y initialized to by
 */
inline void reverse_tape_interpreter(active* y, const double* by, int len) {
  for (int i=0; i<tape_entry::vac; i++)
    tape[i].d=0;

  for (int i=0; i<len; i++)
    tape[y[i].va].d=by[i];

  for (int i=tape_entry::vac; i>=0; i--)
    reverse_tape_interpreter(i);
}
/**
 * sets tape counter back to start index
 */
inline void cleanup_tape(int start) {
  #if (DEBUG_ACTIVE == 2)
  int vac_before=tape_entry::vac;
  #endif

  tape_entry::deps.clear();
  tape_entry::vac=start;

  #if (DEBUG_ACTIVE == 2)
  cerr << "   > start            : " << start << endl
       << "   > tape size before : " << vac_before << endl
       << "   > tape size after  : " << tape_entry::vac << endl;
  #endif
}
inline void reset_tape() {

  #if (DEBUG_ACTIVE == 2)
  cout << "INFO active.h >  reset_tape runs ..." << endl;
  #endif

  tape_entry::indeps.clear();
  tape_entry::deps.clear();

  #if (DEBUG_ACTIVE == 2)
  int vac_before=tape_entry::vac;
  #endif

  // optional, since tape entries are overwritten each time
  for (int i=0; i<tape_entry::vac; i++)
    tape[i].reset();

  tape_entry::vac=0;

  #if (DEBUG_ACTIVE == 2)
  cerr << "   > tape size before : " << vac_before << endl
       << "   > tape size after  : " << tape_entry::vac << endl;

  cout << "INFO active.h >  reset_tape terminates" << endl;
  #endif
}
/**
* resets the independent tape entries by
* by preserving their values.
*/
inline void reset_indep_tape_entries() {

  #if (DEBUG_ACTIVE == 2)
  cout << "INFO active.h >  reset_indep_tape_entries runs ..." << endl;
  #endif

  vector<int>::iterator it;
  for (it=tape_entry::indeps.begin(); it !=tape_entry::indeps.end(); it++)
    tape[*it].d=0;

  #if (DEBUG_ACTIVE == 2)
  cout << "INFO active.h >  reset_indep_tape_entries terminates" << endl;
  #endif
}

inline void tape_stat () {
  cout << " Tape Status : " << endl
       << "  tape            : " << tape << endl
       << "  tape_size       : " <<  TAPE_SIZE << endl
       << "  tape_entry::vac : " << tape_entry::vac << endl
       << "  indeps.size()   : " <<  tape_entry::indeps.size() << endl
       << "  deps.size()     : " <<  tape_entry::deps.size() << endl;
}
/**
* return the static counter
*/
inline int get_ctr() {
  return tape_entry::vac;
}
/**
 * sets the static counter to n
 */
inline void set_counter(const int & n) {
  tape_entry::vac=n;
}
/**
 *  return the partial derivative if ith tape element
 */
inline double get_partial(const int& i) {
  return tape[i].d;
}
/**
 *  return a reference to the ith tape element
 */
inline tape_entry& get_tape_entry(const int& i) {
  return tape[i];
}
#if (DEBUG_ACTIVE == 3)
/**
* computes the pattern of the Jacobian on tape
*/
void compute_jacobian_pattern(unsigned int**& P, int& m);
#endif
/**
 * prints tape onto screen
 */
void print_tape(int i, string kind="INFO");
/**
 * prints tape onto screen
 */
void print_tape();
/**
 * prints tape onto screen
 */
void print_tape_indeps();
/**
 * prints tape onto screen
 */
void print_tape_deps();
/**
 * prints tape onto screen
 */
void print_tape_extremes();
/**
 * prints tape onto screen
 */
void print_deriv_of_indeps(int dep_indx=0);
void print_deriv_of_deps(int indep_indx=0);
/**
* return the operation code corresponding to i
*/
std::string get_oc(const int& i);
/**
* standard output / Input
*/
std::ostream& operator<<( std::ostream& out, const active& x);
std::ostream& operator<<( std::ostream& out, const tape_entry& x);
std::istream& operator>>( std::istream& in, active& x);
/**
 * prints the dot version of ith tape entry 
 */
int print_tape_entry_to_dot(ofstream& out, const int& i);
/**
 * prints graphviz file dagj.dot representing
 * the spanning tree of active variable y
 */
void print_spanning_tree_to_dot(const char* filename, const active& y);
/**
 * prints graphviz file dag.dot
 */
void print_tape_to_dot(const char* filename, int start, int end);
/**
 * prints dag.dot starting at vertex corresponding to tape entry i backwards
 * with spanned graph of path length depths.
 */
void tape_to_dot(const char* filename, const int& i, const int& depth);
/**
 * recursive to dot of a vertex corresponding to tape i and its dependent vertices 
 */
int to_dot(std::ofstream& out, const int& i, const int& depth, int& ctr);

#if (DEBUG_ACTIVE == 1)
/**
 * prints the tape status 
 */
void print_debug_info();
#endif
// /*//////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////////////
// #if ENABLE_ACTIVE_OPT   // to differ optimized version from previous version of active tape
// 
// #ifndef ACTIVE_H_
// #define ACTIVE_H_
// 
// #include <stdio.h>
// #include <cmath>
// #include <stdio.h>
// #include <iostream>
// #include <vector>
// 
// using namespace std;
// 
// #ifdef TAPECOUNTER
//     extern unsigned int tapecounter;
//     #define INCTAPE ins++; tapecounter++;
// #else
//     #define INCTAPE   printf("INCTAPE!\n"); ins++;
//     #define INCTAPE ins++;
// #endif
// 
// typedef struct tape_entry
// {
//    double* arg1adj;
//    double* arg2adj;
//    double pval1;
//    double pval2;
//    double d;
// 
//    static vector<int> deps;
//    static vector<int> indeps;
// } TEntry;
// 
// 
// extern TEntry *adjointtape;
// extern TEntry *ins; 
// 
// class active
// {
// public:
//     double v;
//     TEntry* node;
//     bool isactive;
//     static double d;
// 
//     
//     
//     active() __attribute__((always_inline));
//     active(const double& x) __attribute__((always_inline));
//     active(const active& x) __attribute__((always_inline));
//     
//     active& operator=(const active &vneu) __attribute__((always_inline)) {
//         active ret;
//         this->node=vneu.node;
//         this->v=vneu.v;
//         this->isactive=vneu.isactive;   
//         return *this;
//     }
// 
//     /*active& operator=(const active &vneu) {
//         active ret;
//         if(vneu.isactive==false) {
//             this->node=ins; ins++;
//             this->v=vneu.v;
//             this->isactive=true;
//         }
//         else {
//             this->node=vneu.node;
//             this->v=vneu.v;
//             this->isactive=vneu.isactive;   
//         }
//             
//         return *this;
//     }*/
// 
//     inline void setValue(const double &x) __attribute__((always_inline)) {
//         this->v=x;
//     }
// 
//     inline double getderv() __attribute__((always_inline)) {
//         if(node==0) return 0;
//         return node->d; 
//     }   
//     
//     inline active& operator+=(const active& x) __attribute__((always_inline));
//     inline active& operator-=(const active& x) __attribute__((always_inline));
//     inline active& operator*=(const active& x) __attribute__((always_inline));
//     inline active& operator/=(const active& x) __attribute__((always_inline));
//     
// };
// 
// inline active::active(): v(0.), node(0), isactive(false) {
// }
// 
// inline active::active(const double& x) : v(x), node(ins), isactive(false) {
// }
// 
// inline active::active(const active& x) : v(x.v), node(x.node), isactive(x.isactive) {
// }
// 
// 
// inline bool operator == (const active& x1, const active& x2) __attribute__((always_inline));
// inline bool operator != (const active& x1, const active& x2) __attribute__((always_inline));
// inline bool operator <  (const active& x1, const active& x2) __attribute__((always_inline));
// inline bool operator <= (const active& x1, const active& x2) __attribute__((always_inline));
// inline bool operator >  (const active& x1, const active& x2) __attribute__((always_inline));
// inline bool operator >= (const active& x1, const active& x2) __attribute__((always_inline));
// 
// 
// inline bool operator == (const active& x1, const active& x2) {  return x1.v==x2.v; }
// inline bool operator != (const active& x1, const active& x2) { return x1.v!=x2.v; }
// inline bool operator <  (const active& x1, const active& x2) { return x1.v<x2.v; }
// inline bool operator <= (const active& x1, const active& x2) { return x1.v<=x2.v; }
// inline bool operator >  (const active& x1, const active& x2) { return x1.v>x2.v; }
// inline bool operator >= (const active& x1, const active& x2) { return x1.v>=x2.v; }
// 
// 
// 
// 
// #define BINOPERATOR(operatorX, pval1X, pval2X)  inline active operator operatorX(const active& x1, const active& x2) __attribute__((always_inline)); inline active operator operatorX(const active& x1, const active& x2) { active ret(x1.v operatorX x2.v);    if(x1.isactive==false && x2.isactive==false) return ret;    ret.isactive=true; ret.node=ins;    if(x1.isactive) ins->arg1adj=&x1.node->d;           if(x2.isactive) ins->arg2adj=&x2.node->d;           ins->pval1=pval1X;          ins->pval2=pval2X;                  INCTAPE return ret;}
// 
// BINOPERATOR(+,1,1)
// BINOPERATOR(-,1,-1)
// BINOPERATOR(*,x2.v,x1.v)
// BINOPERATOR(/,1/x2.v,-x1.v/x2.v/x2.v)
// 
// 
// #define UNNOPERATOR(operatorX, pval1X, pval2X) inline active& active::operator operatorX (const active& x2) {       if(this->isactive==false && x2.isactive==false) {       this->v operatorX x2.v;     return *this;   }   if(this->isactive) { ins->arg1adj=&this->node->d; ins->pval1=pval1X; }  if(x2.isactive) { ins->arg2adj=&x2.node->d; ins->pval2=pval2X; }    this->v operatorX x2.v;         this->isactive=true; this->node=ins; INCTAPE return *this; }
// 
// UNNOPERATOR(+=,1,1)
// UNNOPERATOR(-=,1,-1)
// UNNOPERATOR(*=,x2.v,this->v)
// UNNOPERATOR(/=,1/x2.v,-this->v/x2.v/x2.v)
// 
// inline active operator-(const active& x1) __attribute__((always_inline));
// inline active operator-(const active& x1)  {
//     active ret(-x1.v);
//     if(x1.isactive==false) return ret;
//     ret.isactive=true;
//     ret.node=ins;
//     ins->arg1adj=&x1.node->d;
//     ins->pval1=-1;
//     ins++;
//     return ret;
// }
// 
// 
// #define UNARYFUNCTION(fuu, dfuu) inline active fuu(const active &x) __attribute__((always_inline)); inline active fuu(const active &x) {    active ret( fuu(x.v) ); if(x.isactive==false) return ret;   if(x.node->arg2adj==0) {        ret.isactive=true; ret.node=x.node;     ret.node->pval1*=dfuu;      return ret; }   else  {     ret.isactive=true; ret.node=ins;        ins->arg1adj=&x.node->d;        ins->pval1=dfuu;            INCTAPE return ret; } }
// 
// 
// #define UNARYFUNCTION(fuu, dfuu) inline active fuu(const active &x) __attribute__((always_inline)); inline active fuu(const active &x) {  active ret( fuu(x.v) ); if(x.isactive==false) return ret;   ret.isactive=true; ret.node=ins;    ins->arg1adj=&x.node->d;    ins->pval1=dfuu;        INCTAPE return ret;}
// 
// UNARYFUNCTION(sin, cos(x.v))
// UNARYFUNCTION(asin, 1/sqrt(1.0-x.v*x.v))
// UNARYFUNCTION(cos, -sin(x.v))
// UNARYFUNCTION(acos, -1/sqrt(1.0-x.v*x.v))
// UNARYFUNCTION(exp, exp(x.v))
// UNARYFUNCTION(atan, 1/(1+x.v*x.v))
// UNARYFUNCTION(tanh, 1/cosh(x.v)/cosh(x.v))
// UNARYFUNCTION(sqrt, 1/(2.0*sqrt(x.v)))
// UNARYFUNCTION(log1p, 1/(x.v+1))
// UNARYFUNCTION(log, 1/x.v)
// 
// 
// 
// #define BINARYFUNCTION(fuu, pval1X, pval2X) active inline fuu(const active &x1, const active &x2) __attribute__((always_inline)); active inline fuu(const active &x1, const active &x2) {   active ret( fuu(x1.v,x2.v) );   if(x1.isactive==false && x2.isactive==false) return ret;    ret.node=ins;   if(x1.isactive) { ins->arg1adj=&x1.node->d; ins->pval1=pval1X; }    if(x2.isactive) { ins->arg2adj=&x2.node->d; ins->pval2=pval2X; }    INCTAPE return ret;}
// 
// BINARYFUNCTION(pow, x2.v*pow(x1.v,x2.v-1.0), log(x1.v)*ret.v)
// BINARYFUNCTION(hypot, x1.v/ret.v, x2.v/ret.v)
// 
// active inline atan2(const active &x1, const active &x2) __attribute__((always_inline));
// active inline atan2(const active &x1, const active &x2) {
//     return atan(x1/x2);
// }
// 
// inline active abs(const active& x)  __attribute__((always_inline));
// inline active abs(const active& x) {
//     if(x.v>=0) return x;
//     return -x;
// }
// 
// inline active max(const active& x1, const active& x2) __attribute__((always_inline));
// inline active max(const active& x1, const active& x2) {
//     if(x1.v>=x2.v) return x1;
//     return x2;
// }
// 
// inline active fmax(const active& x1, const active& x2) __attribute__((always_inline));
// inline active fmax(const active& x1, const active& x2) {
//     if(x1.v>=x2.v) return x1;
//     return x2;
// }
// 
// inline active min(const active& x1, const active& x2) __attribute__((always_inline));
// inline active min(const active& x1, const active& x2) {
//     if(x1.v<=x2.v) return x1;
//     return x2;
// }
// 
// inline active fmin(const active& x1, const active& x2) __attribute__((always_inline));
// inline active fmin(const active& x1, const active& x2) {
//     if(x1.v<=x2.v) return x1;
//     return x2;
// }
// 
// inline bool isnan(const active& a) {
//     return isnan(a.v);  
// } 
// 
// inline ostream& operator<<( ostream& out, const active& x) {
//     
//     out << x.v << endl;
//     return out; 
// }
// 
// inline ostream& operator<<( ostream& out, const tape_entry& x) {
//     out << "adj=" << x.d << endl;
//     return out; 
// }
// 
// inline istream& operator>>( istream& in, active& x) {
//   if (x.va >= 0) in >> x.v;
//   if (DEBUG_ACTIVE && x.va >= 0) cout << "in >> tape[" << x.va << "].v" << x << endl;
//   cout << "in >> tape " << x << endl;
//   return in;
// }
// 
// inline const double& value(const active& val) {
//     return val.v;
// }
// 
// 
// inline void reverse_tape_interpreter() {
//     TEntry *act=ins;
//     
//     printf("tapecounter=%li\n",(act-adjointtape));
//     
//     while(act!=adjointtape) {
//         act--;
//         if(act->d!=0) {
//             if(act->arg1adj!=0) *act->arg1adj+= act->d*act->pval1;
//             if(act->arg2adj!=0) *act->arg2adj+= act->d*act->pval2;
//         }
//     }
// }
// 
// 
// inline TEntry&  get_tape_entry(const int& i) {
//     return adjointtape[i];  
// }
// 
// inline void seed_independent(const int& i){
//     cout << "seed indep " << i << endl;
// }
// 
// inline void seed_dependent(const int& i) 
// {
//     printf("seed with i=%i\n",i);
//     
//     TEntry *actNode2Interpret=ins;
//     while(actNode2Interpret!=adjointtape) {
//         actNode2Interpret--;
//         actNode2Interpret->d=0;
//     }   
//     if(i>0) adjointtape[i].d=1;
// }
// 
// 
// //////////////////////////////////////////
// inline void forward_tape_interpreter() {
// 
// }
// 
// inline void print_tape() {
//     
// }
// 
// inline void print_tape(int i) {
//     cout << " print_tape with i=" << i << endl;
// }
// 
// inline void independent(active& x) {
//     x.isactive=true;
//     x.node=ins; ins++;
// }
// 
// inline void dependent(const active& x, vector<int>& deps) {
//   unsigned long int idx=(x.node-adjointtape);
//   deps.push_back(int(idx));
// }
// 
// 
// inline void independent(active& x, vector<int>& deps) {
//   x.isactive=true;
//   x.node=ins; ins++;
//   unsigned long int idx=(x.node-adjointtape);
//   deps.push_back(int(idx));
// }
// 
// inline void independent(active& x, vector<active*>& deps) {
//     x.isactive=true;
//     x.node=ins; ins++;
//     deps.push_back(&x);
// }
// 
// inline void init_mode(int mode, int b) {
//     cout << " initmode << " << mode << "b=" << b << endl;
// }
// #endif /*ACTIVE_H_*/*/
// #endif

/* AMPI */

//void ampi_get_val(int64_t, double);
void ampi_get_val(void*, int*, double*);
void ampi_set_val(void*, int*, double*);
void ampi_get_adj(int64_t*, double*);
void ampi_set_adj(int64_t*, double*);
void ampi_get_idx(void *, int *, int64_t*);
void ampi_create_tape_entry(int*);

//int AMPI_Send(active*, int, MPI_Datatype, int, int, MPI_Comm);
//int AMPI_Recv(active*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status*);

//int AMPI_Isend(active*, int, MPI_Datatype, int, int, MPI_Comm, AMPI_dco_Request *);
//int AMPI_Irecv(active*, int, MPI_Datatype, int, int, MPI_Comm, AMPI_dco_Request *);

//int AMPI_Wait(AMPI_dco_Request *, MPI_Status *);
//int AMPI_Waitall(int , AMPI_dco_Request *, MPI_Status *);
//int AMPI_Awaitall(int , AMPI_dco_Request *, MPI_Status *);

//int AMPI_Bcast(active *buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm);
//int AMPI_Reduce(active *sendbuf, active *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
#endif // ACTIVE_INCLUDE
