#include <iostream>
#include <stdlib.h>
#include <cmath>

#define DCO_INLINE __attribute__((always_inline))


/////////////////////////
    #define TEMPLATE_BUILDER_CONSTUCTOR(ARG1) template<ARG1 >
    #define TEMPLATE_BUILDER_1(ARG1) template<class DCO_TAPE_REAL, ARG1 >
    #define TEMPLATE_BUILDER_2(ARG1, ARG2) template<class DCO_TAPE_REAL, ARG1, ARG2 >
    /////////////////////////
    #define TEMPLATE_ARG1_ACTIVE_CLASSES  class DCO_HANDLER_DATA_1, class DCO_HANDLER_1
    #define TEMPLATE_ARG2_ACTIVE_CLASSES  class DCO_HANDLER_DATA_2, class DCO_HANDLER_2

    #define TEMPLATE_ARG1_UN_CLASSES class A1_T, class A1_OP
    #define TEMPLATE_ARG2_UN_CLASSES class A2_T, class A2_OP

    #define TEMPLATE_ARG1_BIN_AA_CLASSES class A1_T1, class A1_T2, class A1_OP
    #define TEMPLATE_ARG2_BIN_AA_CLASSES class A2_T1, class A2_T2, class A2_OP

    #define TEMPLATE_ARG1_BIN_AP_CLASSES class A1_T1, class A1_OP
    #define TEMPLATE_ARG2_BIN_AP_CLASSES class A2_T1, class A2_OP

    #define TEMPLATE_ARG1_BIN_PA_CLASSES class A1_T2, class A1_OP
    #define TEMPLATE_ARG2_BIN_PA_CLASSES class A2_T2, class A2_OP

    /////////////////////////
    #define DCO_ACTIVE_TYPE_1 active_type<DCO_TAPE_REAL, DCO_HANDLER_DATA_1, DCO_HANDLER_1>
    #define DCO_ACTIVE_TYPE_2 active_type<DCO_TAPE_REAL, DCO_HANDLER_DATA_2, DCO_HANDLER_2>

    #define DCO_UN_TYPE_1 unary_intermediate<DCO_TAPE_REAL, A1_T, A1_OP>
    #define DCO_UN_TYPE_2 unary_intermediate<DCO_TAPE_REAL, A2_T, A2_OP>

    #define DCO_BIN_AA_TYPE_1 binary_intermediate_aa<DCO_TAPE_REAL, A1_T1, A1_T2, A1_OP>
    #define DCO_BIN_AA_TYPE_2 binary_intermediate_aa<DCO_TAPE_REAL, A2_T1, A2_T2, A2_OP>

    #define DCO_BIN_AP_TYPE_1 binary_intermediate_ap<DCO_TAPE_REAL, A1_T1, A1_OP>
    #define DCO_BIN_AP_TYPE_2 binary_intermediate_ap<DCO_TAPE_REAL, A2_T1, A2_OP>

    #define DCO_BIN_PA_TYPE_1 binary_intermediate_pa<DCO_TAPE_REAL, A1_T2, A1_OP>
    #define DCO_BIN_PA_TYPE_2 binary_intermediate_pa<DCO_TAPE_REAL, A2_T2, A2_OP>

/////////////////////////

namespace dcoV2
{
    typedef     unsigned char   DCO_EDGE_COUNT_TYPE;

    #include "tlm_type.hpp"

#define DCO_ACTIVITY 1

#ifdef DCO_ACTIVITY
#define DCO_ACTIVITY_STATEMENT(a) a
#define DCO_ACTIVITY_STATEMENT_CONSTRUCT(a) a,
#else
#define DCO_ACTIVITY_STATEMENT(a)
#define DCO_ACTIVITY_STATEMENT_CONSTRUCT(a)
#endif

    template<class DCO_TAPE_REAL, class DCO_ARG, class DCO_OPERATION> struct unary_intermediate{
        const   DCO_TAPE_REAL   _value;
        DCO_ACTIVITY_STATEMENT(const   DCO_EDGE_COUNT_TYPE     _edgecount;)
        const   DCO_ARG         &_arg;

        unary_intermediate(const DCO_ARG &arg) :
                _value(DCO_OPERATION::eval(arg)),
                DCO_ACTIVITY_STATEMENT_CONSTRUCT(_edgecount(arg._edgecount))
                _arg(arg) {}

        static inline const DCO_TAPE_REAL pval(const DCO_TAPE_REAL &_value, const DCO_ARG &x1, const DCO_TAPE_REAL &pval) {
            return DCO_OPERATION::calc_partial(_value, x1, pval);
        }
    };

    template<class DCO_TAPE_REAL, class DCO_ARG1, class DCO_ARG2, class DCO_OPERATION> struct binary_intermediate_aa{
        const   DCO_TAPE_REAL   _value;
        DCO_ACTIVITY_STATEMENT(const   DCO_EDGE_COUNT_TYPE     _edgecount;)
        const DCO_ARG1          &_arg1;
        const DCO_ARG2          &_arg2;

        binary_intermediate_aa(const DCO_ARG1 &arg1, const DCO_ARG2 &arg2) :
                _value(DCO_OPERATION::eval(arg1,arg2)),
                DCO_ACTIVITY_STATEMENT_CONSTRUCT(_edgecount(arg1._edgecount + arg2._edgecount))
                _arg1(arg1),
                _arg2(arg2) {}

        static inline const DCO_TAPE_REAL pval1(const DCO_TAPE_REAL &_value, const DCO_ARG1 &x1, const DCO_ARG2 &x2, const DCO_TAPE_REAL &pval) {
            return DCO_OPERATION::calc_partial1(_value, x1,x2, pval);
        }

        static inline const DCO_TAPE_REAL pval2(const DCO_TAPE_REAL &_value, const DCO_ARG1 &x1, const DCO_ARG2 &x2, const DCO_TAPE_REAL &pval) {
            return DCO_OPERATION::calc_partial2(_value, x1,x2, pval);
        }
    };

    template<class DCO_TAPE_REAL, class DCO_ARG1, class DCO_OPERATION> struct binary_intermediate_ap{
        const   DCO_TAPE_REAL   _value;
        DCO_ACTIVITY_STATEMENT(const   DCO_EDGE_COUNT_TYPE     _edgecount;)
        const DCO_ARG1          &_arg1;
        const double            _arg2;

        binary_intermediate_ap(const DCO_ARG1 &arg1, const double &arg2) :
                _value(DCO_OPERATION::eval(arg1,arg2)),
                DCO_ACTIVITY_STATEMENT_CONSTRUCT(_edgecount(arg1._edgecount))
                _arg1(arg1),
                _arg2(arg2) {}

        static inline const DCO_TAPE_REAL pval1(const DCO_TAPE_REAL &_value, const DCO_ARG1 &x1, const double &x2, const DCO_TAPE_REAL &pval) {
            return DCO_OPERATION::calc_partial1(_value,x1,x2, pval);
        }

    };

    template<class DCO_TAPE_REAL, class DCO_ARG2, class DCO_OPERATION> struct binary_intermediate_pa{
        const   DCO_TAPE_REAL   _value;
        DCO_ACTIVITY_STATEMENT(const   DCO_EDGE_COUNT_TYPE     _edgecount;)
        const double            _arg1;
        const DCO_ARG2          &_arg2;

        binary_intermediate_pa(const double &arg1, const DCO_ARG2 &arg2) :
                _value(DCO_OPERATION::eval(arg1,arg2)),
                DCO_ACTIVITY_STATEMENT_CONSTRUCT(_edgecount(arg2._edgecount))
                _arg1(arg1),
                _arg2(arg2) {}

        inline DCO_TAPE_REAL pval1Ex(const DCO_TAPE_REAL &pval) {
            return DCO_OPERATION::calc_partial1(_value,_arg1,_arg2, pval);
        }

        static inline DCO_TAPE_REAL pval1(const DCO_TAPE_REAL &_value, const double &x1, const DCO_ARG2 &x2, const DCO_TAPE_REAL &pval) {
            return DCO_OPERATION::calc_partial1(_value,x1,x2, pval);
        }

        static inline DCO_TAPE_REAL pval2(const DCO_TAPE_REAL &_value, const double &x1, const DCO_ARG2 &x2, const DCO_TAPE_REAL &pval) {
            return DCO_OPERATION::calc_partial2(_value,x1,x2, pval);
        }
    };

#define DCO_UNUSED __attribute__((unused))

    namespace internal {
        #define BIN_OPERATOR_AA_CLASS(OPNAME, OPERATORX, PVAL1, PVAL2) \
        template<class DCO_TAPE_REAL>struct dco_##OPNAME##_aa{ \
            template<class T1, class T2>static inline const DCO_TAPE_REAL eval(const T1 &arg1, const T2 &arg2) { \
                return arg1._value OPERATORX arg2._value; \
            } \
            template<class T1,class T2>static inline const DCO_TAPE_REAL calc_partial1(const DCO_TAPE_REAL &_value DCO_UNUSED, const T1 &arg1 DCO_UNUSED, const T2 &arg2 DCO_UNUSED, const DCO_TAPE_REAL &pval) { \
                return PVAL1; \
            } \
            template<class T1,class T2>static inline const DCO_TAPE_REAL calc_partial2(const DCO_TAPE_REAL &_value DCO_UNUSED, const T1 &arg1 DCO_UNUSED, const T2 &arg2 DCO_UNUSED, const DCO_TAPE_REAL &pval) { \
                return PVAL2; \
            } \
        };

        BIN_OPERATOR_AA_CLASS(add,+,pval,pval)
        BIN_OPERATOR_AA_CLASS(sub,-,pval,-pval)
        BIN_OPERATOR_AA_CLASS(mul,*,arg2._value*pval,arg1._value*pval)
        BIN_OPERATOR_AA_CLASS(div,/,pval/arg2._value,-pval*_value/arg2._value)

        #define BIN_OPERATOR_AP_PA_CLASS(OPNAME, OPERATORX, PVAL1, PVAL2) \
        template<class DCO_TAPE_REAL>struct dco_##OPNAME##_ap{ \
            template<class T1>static inline const DCO_TAPE_REAL eval(const T1 &arg1, const double &arg2) { \
                return arg1._value OPERATORX arg2; \
            } \
            template<class T1>static inline const DCO_TAPE_REAL calc_partial1(const DCO_TAPE_REAL &_value DCO_UNUSED, const T1 &arg1 DCO_UNUSED, const double &arg2 DCO_UNUSED, const DCO_TAPE_REAL &pval) { \
                return PVAL1; \
            } \
        }; \
        template<class DCO_TAPE_REAL>struct dco_##OPNAME##_pa{ \
            template<class T2>static inline const DCO_TAPE_REAL eval(const double &arg1, const T2 &arg2) { \
                return arg1 OPERATORX arg2._value; \
            } \
            template<class T2>static inline const DCO_TAPE_REAL calc_partial2(const DCO_TAPE_REAL &_value DCO_UNUSED, const double &arg1 DCO_UNUSED, const T2 &arg2 DCO_UNUSED, const DCO_TAPE_REAL &pval) { \
                return PVAL2; \
            } \
        };

        BIN_OPERATOR_AP_PA_CLASS(add,+,pval,pval)
        BIN_OPERATOR_AP_PA_CLASS(sub,-,pval,-pval)
        BIN_OPERATOR_AP_PA_CLASS(mul,*,arg2*pval,arg1*pval)
        BIN_OPERATOR_AP_PA_CLASS(div,/,pval/arg2,-pval*_value/arg2._value)

        #define UNARYFUNCTION_CLASS(FUU,DFUU,PREPROCESS, DEFINE) \
        using DEFINE::FUU; \
        template<class DCO_TAPE_REAL>struct dco_##FUU{ \
            template<class T>static inline const DCO_TAPE_REAL eval(const T &arg) { \
                return FUU(arg._value); \
            } \
            template<class T>static inline const DCO_TAPE_REAL calc_partial(const DCO_TAPE_REAL &_value DCO_UNUSED, const T &x, const DCO_TAPE_REAL &pval) { \
                return (DFUU) * pval;\
            }\
        };

        UNARYFUNCTION_CLASS(sin, cos(x._value),,std)
        UNARYFUNCTION_CLASS(cos, -sin(x._value),,std)
        UNARYFUNCTION_CLASS(tan, (1.0+(tan(x._value)*tan(x._value))),,std)
        UNARYFUNCTION_CLASS(cosh, sinh(x._value),,std)
        UNARYFUNCTION_CLASS(sinh, cosh(x._value),,std)
        UNARYFUNCTION_CLASS(asin, 1/sqrt(1.0-x._value*x._value),,std)
        UNARYFUNCTION_CLASS(acos, -1/sqrt(1.0-x._value*x._value),,std)
        UNARYFUNCTION_CLASS(exp, exp(x._value),,std)
        UNARYFUNCTION_CLASS(atan, 1.0/(1.0+x._value*x._value),,std)
        UNARYFUNCTION_CLASS(tanh, 1.0/cosh(x._value)/cosh(x._value),,std)
        UNARYFUNCTION_CLASS(sqrt, 1.0/(2.0*sqrt(x._value)),,std)
        UNARYFUNCTION_CLASS(log, 1.0/x._value,,std)
        #ifndef MSVC
        UNARYFUNCTION_CLASS(log1p, 1.0/(x._value+1),,)
        UNARYFUNCTION_CLASS(log10, 1.0/(x._value*log(10)),,)
        #endif
        #undef UNARYFUNCTION

        template<class DCO_TAPE_REAL>struct dco_minus{
            template<class T>static inline const DCO_TAPE_REAL eval(const T &arg1) {
                return -arg1._value;
            }
            template<class T>static inline const DCO_TAPE_REAL calc_partial(const DCO_TAPE_REAL &_value DCO_UNUSED, const T &arg1 DCO_UNUSED, const DCO_TAPE_REAL &pval) { \
                return -pval;
            }
        };

        using ::fabs;
        template<class DCO_TAPE_REAL>struct dco_fabs{
            template<class T>static inline const DCO_TAPE_REAL eval(const T &arg1) {
                return fabs(arg1._value);
            }
            template<class T>static inline const DCO_TAPE_REAL calc_partial(const DCO_TAPE_REAL &_value DCO_UNUSED, const T &arg1 DCO_UNUSED, const DCO_TAPE_REAL &pval) { \
                if(_value<0)
                    return -pval;
                else
                    return pval;
            }
        };
        
        template<class DCO_TAPE_REAL>struct dco_abs{
            template<class T>static inline const DCO_TAPE_REAL eval(const T &arg1) {
                return fabs(arg1._value);
            }
            template<class T>static inline const DCO_TAPE_REAL calc_partial(const DCO_TAPE_REAL &_value DCO_UNUSED, const T &arg1 DCO_UNUSED, const DCO_TAPE_REAL &pval) { \
                if(_value<0)
                    return -pval;
                else
                    return pval;
            }
        };

        #define BIN_OPERATION_AA_CLASS(OPNAME, FUNCTION, PVAL1, PVAL2) \
        template<class DCO_TAPE_REAL>struct dco_##OPNAME##_aa{ \
            template<class T1, class T2>static inline const DCO_TAPE_REAL eval(const T1 &arg1, const T2 &arg2) { \
                return FUNCTION(arg1._value,arg2._value); \
            } \
            template<class T1,class T2>static inline const DCO_TAPE_REAL calc_partial1(const DCO_TAPE_REAL _value DCO_UNUSED, const T1 &arg1, const T2 &arg2, const DCO_TAPE_REAL &pval) { \
                return PVAL1; \
            } \
            template<class T1,class T2>static inline const DCO_TAPE_REAL calc_partial2(const DCO_TAPE_REAL _value DCO_UNUSED, const T1 &arg1, const T2 &arg2, const DCO_TAPE_REAL &pval) { \
                return PVAL2; \
            } \
        };

        #define BIN_OPERATION_AP_PA_CLASS(OPNAME, FUNCTION, PVAL1, PVAL2) \
        template<class DCO_TAPE_REAL>struct dco_##OPNAME##_ap{ \
            template<class T>static inline const DCO_TAPE_REAL eval(const T &arg1, const double &arg2) { \
                return FUNCTION(arg1._value,arg2); \
            } \
            template<class T>static inline const DCO_TAPE_REAL calc_partial1(const DCO_TAPE_REAL _value DCO_UNUSED, const T &arg1, const double &arg2, const DCO_TAPE_REAL &pval) { \
                return PVAL1; \
            } \
        }; \
        template<class DCO_TAPE_REAL>struct dco_##OPNAME##_pa{ \
            template<class T>static inline const DCO_TAPE_REAL eval(const double &arg1, const T &arg2) { \
                return FUNCTION(arg1,arg2._value); \
            } \
            template<class T>static inline const DCO_TAPE_REAL calc_partial2(const DCO_TAPE_REAL _value DCO_UNUSED, const double &arg1, const T &arg2, const DCO_TAPE_REAL &pval) { \
                return PVAL2; \
            } \
        };

        #define BIN_OPERATION_CLASS(OPNAME, FUNCTION, PVAL1AA, PVAL2AA, PVAL1AP, PVAL2PA) \
        BIN_OPERATION_AA_CLASS(OPNAME, FUNCTION, PVAL1AA, PVAL2AA) \
        BIN_OPERATION_AP_PA_CLASS(OPNAME, FUNCTION, PVAL1AP, PVAL2PA)

        BIN_OPERATION_CLASS(atan2,atan2,
                            pval/(1.0+(arg1._value/arg2._value)*(arg1._value/arg2._value))/arg2._value,
                            -pval/(1.0+(arg1._value/arg2._value)*(arg1._value/arg2._value))/arg2._value *arg1._value/arg2._value,
                            pval/(1.0+(arg1._value/arg2)*(arg1._value/arg2))/arg2,
                            -pval/(1.0+(arg1/arg2._value)*(arg1/arg2._value))/arg2._value *arg1/arg2._value)

        BIN_OPERATION_CLASS(pow,pow,
                            arg2._value * pval*_value / arg1._value,
                            pval*_value / arg1._value,
                            arg2 * pval*_value / arg1._value,
                            pval*_value / arg1)

        BIN_OPERATION_CLASS(hypot,hypot,
                            pval*arg1._value/_value,
                            pval*arg2._value/_value,
                            pval*arg1._value/_value,
                            pval*arg2._value/_value)

    }

    ///// Active type for partial, bit pattern etc.

    template<class DCO_TAPE_REAL, class DCO_HANDLER_DATA, class DCO_HANDLER> struct active_type{
        DCO_TAPE_REAL       _value;
        DCO_EDGE_COUNT_TYPE _edgecount;
        DCO_HANDLER_DATA    _data;

        active_type() : _value(0), _edgecount(0) {}
        active_type(const double &val) : _value(val), _edgecount(0) {}

        active_type(const active_type &x) : _value(x._value), _edgecount(x._edgecount) { if(_edgecount>0) _data = x._data;}

        active_type& operator = (const active_type &x) {
            this->_value = x._value;
            this->_edgecount = x._edgecount;
            if(this->_edgecount>0) {
                this->_data = x._data;
            }
            return *this;
        }

#ifdef DCO_ACTIVITY
        #define DCO_ACTIVE_TYPE_KONSTRUTOR(PREFIX, DT) \
        PREFIX inline void DCO_INLINE build_from(const DT &x) { \
            if(x._edgecount>0) { \
                DCO_HANDLER handler (x, *this); \
                _edgecount=1; \
                /*handler.interpret(x, 1.0);*/ \
                /*this->_edgecount=handler.finalize(this->_data)?1:0;*/ \
            } \
            else { \
                if(_edgecount!=0) _edgecount=0; \
            } \
            this->_value = x._value; \
        } \
        PREFIX active_type(const DT &x) { \
            build_from(x); \
        } \
        PREFIX inline active_type& operator = (const DT &x) { \
            build_from(x); return *this; \
        }
#else
        #define DCO_ACTIVE_TYPE_KONSTRUTOR(PREFIX, DT) \
        PREFIX inline void build_from(const DT &x) { \
            DCO_HANDLER handler(0); \
            handler.interpret(x, 1.0); \
            this->_edgecount=handler.finalize(this->_data)?1:0; \
            this->_value = x._value; \
        } \
        PREFIX active_type(const DT &x) { \
            build_from(x); \
        } \
        PREFIX inline active_type& operator = (const DT &x) { \
            build_from(x); return *this; \
        }

#endif


        DCO_ACTIVE_TYPE_KONSTRUTOR( TEMPLATE_BUILDER_CONSTUCTOR(TEMPLATE_ARG1_BIN_AA_CLASSES), DCO_BIN_AA_TYPE_1)
        DCO_ACTIVE_TYPE_KONSTRUTOR( TEMPLATE_BUILDER_CONSTUCTOR(TEMPLATE_ARG1_BIN_AP_CLASSES), DCO_BIN_AP_TYPE_1)
        DCO_ACTIVE_TYPE_KONSTRUTOR( TEMPLATE_BUILDER_CONSTUCTOR(TEMPLATE_ARG1_BIN_PA_CLASSES), DCO_BIN_PA_TYPE_1)
        DCO_ACTIVE_TYPE_KONSTRUTOR( TEMPLATE_BUILDER_CONSTUCTOR(TEMPLATE_ARG1_UN_CLASSES), DCO_UN_TYPE_1)

        #define UNARY_OPERATOR_TEMPLATE(OPERATORX, PUREOPERATOR) \
        inline active_type& operator OPERATORX (const active_type &x){ *this = *this PUREOPERATOR x; return *this; } \
        TEMPLATE_BUILDER_CONSTUCTOR( TEMPLATE_ARG1_BIN_AA_CLASSES) inline active_type& operator OPERATORX (const DCO_BIN_AA_TYPE_1 &x){ *this = *this PUREOPERATOR x; return *this; } \
        TEMPLATE_BUILDER_CONSTUCTOR( TEMPLATE_ARG1_BIN_AP_CLASSES) inline active_type& operator OPERATORX (const DCO_BIN_AP_TYPE_1 &x){ *this = *this PUREOPERATOR x; return *this; } \
        TEMPLATE_BUILDER_CONSTUCTOR( TEMPLATE_ARG1_BIN_PA_CLASSES) inline active_type& operator OPERATORX (const DCO_BIN_PA_TYPE_1 &x){ *this = *this PUREOPERATOR x; return *this; } \
        TEMPLATE_BUILDER_CONSTUCTOR( TEMPLATE_ARG1_UN_CLASSES) inline active_type& operator OPERATORX (const DCO_UN_TYPE_1 &x){ *this = *this PUREOPERATOR x; return *this; } \
        inline active_type& operator OPERATORX (const double &x) { *this = *this PUREOPERATOR x; return *this; }

        UNARY_OPERATOR_TEMPLATE(+=,+)
        UNARY_OPERATOR_TEMPLATE(-=,-)
        UNARY_OPERATOR_TEMPLATE(*=,*)
        UNARY_OPERATOR_TEMPLATE(/=,/)

        inline active_type& operator ++() { ++this->_value; return *this; }
        inline active_type& operator --() { --this->_value; return *this; }

        inline active_type operator ++(int) {
            active_type ret(*this); ++ret._value; return ret;
        }
        inline active_type operator --(int) {
            active_type ret(*this); --ret._value; return ret;
        }
    };




    #define DCO_UNARY_OPERATION_TEMPLATE(PREFIX, DT, FUNCTION, OPNAME) \
    PREFIX unary_intermediate<DCO_TAPE_REAL, DT, internal::dco_##OPNAME<DCO_TAPE_REAL> > \
    FUNCTION( const DT &x1) { return unary_intermediate<DCO_TAPE_REAL, DT,internal::dco_##OPNAME<DCO_TAPE_REAL> >(x1); }

    #define DCO_UNARY_OPERATION(FUNCTION, OPNAME) \
    DCO_UNARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_ACTIVE_CLASSES), DCO_ACTIVE_TYPE_1,FUNCTION,OPNAME) \
    DCO_UNARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_BIN_AA_CLASSES), DCO_BIN_AA_TYPE_1,FUNCTION,OPNAME) \
    DCO_UNARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_BIN_AP_CLASSES), DCO_BIN_AP_TYPE_1,FUNCTION,OPNAME) \
    DCO_UNARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_BIN_PA_CLASSES), DCO_BIN_PA_TYPE_1,FUNCTION,OPNAME) \
    DCO_UNARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_UN_CLASSES), DCO_UN_TYPE_1,FUNCTION,OPNAME)

    DCO_UNARY_OPERATION( operator -, minus)

    #define DCO_UNARY_FUNCTION(NAME) DCO_UNARY_OPERATION(NAME,NAME)

    DCO_UNARY_FUNCTION(sin)
    DCO_UNARY_FUNCTION(cos)
    DCO_UNARY_FUNCTION(tan)
    DCO_UNARY_FUNCTION(cosh)
    DCO_UNARY_FUNCTION(sinh)
    DCO_UNARY_FUNCTION(asin)
    DCO_UNARY_FUNCTION(acos)
    DCO_UNARY_FUNCTION(exp)
    DCO_UNARY_FUNCTION(atan)
    DCO_UNARY_FUNCTION(tanh)
    DCO_UNARY_FUNCTION(sqrt)
    DCO_UNARY_FUNCTION(log)
    #ifndef MSVC
    DCO_UNARY_FUNCTION(log1p)
    DCO_UNARY_FUNCTION(log10)
    #endif
    #undef UNARYFUNCTION
    DCO_UNARY_FUNCTION(fabs)
    DCO_UNARY_FUNCTION(abs)

    #define DCO_BINARY_OPERATION_TEMPLATE(PREFIX, DT1, DT2, FUNCTION, OPNAME) \
    PREFIX \
    binary_intermediate_aa<DCO_TAPE_REAL, DT1, DT2, internal::dco_##OPNAME##_aa<DCO_TAPE_REAL> > \
    FUNCTION (const DT1 &x1, const DT2 &x2) { \
        return binary_intermediate_aa<DCO_TAPE_REAL, DT1,DT2, internal::dco_##OPNAME##_aa<DCO_TAPE_REAL> >(x1,x2); \
    }

    #define DCO_BINARY_AP_PA_OPERATION_TEMPLATE(PREFIX, DT, FUNCTION, OPNAME) \
    PREFIX \
    binary_intermediate_ap<DCO_TAPE_REAL, DT, internal::dco_##OPNAME##_ap<DCO_TAPE_REAL> > \
    FUNCTION (const DT &x1, const double &x2) { \
        return binary_intermediate_ap<DCO_TAPE_REAL, DT, internal::dco_##OPNAME##_ap<DCO_TAPE_REAL> >(x1,x2); \
    } \
    PREFIX \
    binary_intermediate_pa<DCO_TAPE_REAL, DT, internal::dco_##OPNAME##_pa<DCO_TAPE_REAL> > \
    FUNCTION (const double &x1, const DT &x2) { \
        return binary_intermediate_pa<DCO_TAPE_REAL, DT, internal::dco_##OPNAME##_pa<DCO_TAPE_REAL> >(x1,x2); \
    }

    #define BIN_AP_PA_OPERATION(OPERATION, OPNAME) \
    DCO_BINARY_AP_PA_OPERATION_TEMPLATE( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_ACTIVE_CLASSES), DCO_ACTIVE_TYPE_1, OPERATION, OPNAME) \
    DCO_BINARY_AP_PA_OPERATION_TEMPLATE( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_UN_CLASSES), DCO_UN_TYPE_1, OPERATION, OPNAME) \
    DCO_BINARY_AP_PA_OPERATION_TEMPLATE( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_BIN_AA_CLASSES), DCO_BIN_AA_TYPE_1, OPERATION, OPNAME) \
    DCO_BINARY_AP_PA_OPERATION_TEMPLATE( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_BIN_AP_CLASSES), DCO_BIN_AP_TYPE_1, OPERATION, OPNAME) \
    DCO_BINARY_AP_PA_OPERATION_TEMPLATE( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_BIN_PA_CLASSES), DCO_BIN_PA_TYPE_1, OPERATION, OPNAME)

    #define BIN_OPERATION(OPERATION,OPNAME) \
    DCO_BINARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_ACTIVE_CLASSES), DCO_ACTIVE_TYPE_1, DCO_ACTIVE_TYPE_1, OPERATION, OPNAME) \
    DCO_BINARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_ACTIVE_CLASSES, TEMPLATE_ARG2_UN_CLASSES), DCO_ACTIVE_TYPE_1, DCO_UN_TYPE_2, OPERATION, OPNAME) \
    DCO_BINARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_ACTIVE_CLASSES, TEMPLATE_ARG2_BIN_AA_CLASSES), DCO_ACTIVE_TYPE_1, DCO_BIN_AA_TYPE_2, OPERATION, OPNAME) \
    DCO_BINARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_ACTIVE_CLASSES, TEMPLATE_ARG2_BIN_AP_CLASSES), DCO_ACTIVE_TYPE_1, DCO_BIN_AP_TYPE_2, OPERATION, OPNAME) \
    DCO_BINARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_ACTIVE_CLASSES, TEMPLATE_ARG2_BIN_PA_CLASSES), DCO_ACTIVE_TYPE_1, DCO_BIN_PA_TYPE_2, OPERATION, OPNAME) \
    \
    DCO_BINARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_UN_CLASSES, TEMPLATE_ARG2_ACTIVE_CLASSES), DCO_UN_TYPE_1, DCO_ACTIVE_TYPE_2, OPERATION, OPNAME) \
    DCO_BINARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_UN_CLASSES, TEMPLATE_ARG2_UN_CLASSES), DCO_UN_TYPE_1, DCO_UN_TYPE_2, OPERATION, OPNAME) \
    DCO_BINARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_UN_CLASSES, TEMPLATE_ARG2_BIN_AA_CLASSES), DCO_UN_TYPE_1, DCO_BIN_AA_TYPE_2, OPERATION, OPNAME) \
    DCO_BINARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_UN_CLASSES, TEMPLATE_ARG2_BIN_AP_CLASSES), DCO_UN_TYPE_1, DCO_BIN_AP_TYPE_2, OPERATION, OPNAME) \
    DCO_BINARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_UN_CLASSES, TEMPLATE_ARG2_BIN_PA_CLASSES), DCO_UN_TYPE_1, DCO_BIN_PA_TYPE_2, OPERATION, OPNAME) \
    \
    DCO_BINARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_AA_CLASSES, TEMPLATE_ARG2_ACTIVE_CLASSES), DCO_BIN_AA_TYPE_1, DCO_ACTIVE_TYPE_2, OPERATION, OPNAME) \
    DCO_BINARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_AA_CLASSES, TEMPLATE_ARG2_UN_CLASSES), DCO_BIN_AA_TYPE_1, DCO_UN_TYPE_2, OPERATION, OPNAME) \
    DCO_BINARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_AA_CLASSES, TEMPLATE_ARG2_BIN_AA_CLASSES), DCO_BIN_AA_TYPE_1, DCO_BIN_AA_TYPE_2, OPERATION, OPNAME) \
    DCO_BINARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_AA_CLASSES, TEMPLATE_ARG2_BIN_AP_CLASSES), DCO_BIN_AA_TYPE_1, DCO_BIN_AP_TYPE_2, OPERATION, OPNAME) \
    DCO_BINARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_AA_CLASSES, TEMPLATE_ARG2_BIN_PA_CLASSES), DCO_BIN_AA_TYPE_1, DCO_BIN_PA_TYPE_2, OPERATION, OPNAME) \
    \
    DCO_BINARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_AP_CLASSES, TEMPLATE_ARG2_ACTIVE_CLASSES), DCO_BIN_AP_TYPE_1, DCO_ACTIVE_TYPE_2, OPERATION, OPNAME) \
    DCO_BINARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_AP_CLASSES, TEMPLATE_ARG2_UN_CLASSES), DCO_BIN_AP_TYPE_1, DCO_UN_TYPE_2, OPERATION, OPNAME) \
    DCO_BINARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_AP_CLASSES, TEMPLATE_ARG2_BIN_AA_CLASSES), DCO_BIN_AP_TYPE_1, DCO_BIN_AA_TYPE_2, OPERATION, OPNAME) \
    DCO_BINARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_AP_CLASSES, TEMPLATE_ARG2_BIN_AP_CLASSES), DCO_BIN_AP_TYPE_1, DCO_BIN_AP_TYPE_2, OPERATION, OPNAME) \
    DCO_BINARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_AP_CLASSES, TEMPLATE_ARG2_BIN_PA_CLASSES), DCO_BIN_AP_TYPE_1, DCO_BIN_PA_TYPE_2, OPERATION, OPNAME) \
    \
    DCO_BINARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_PA_CLASSES, TEMPLATE_ARG2_ACTIVE_CLASSES), DCO_BIN_PA_TYPE_1, DCO_ACTIVE_TYPE_2, OPERATION, OPNAME) \
    DCO_BINARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_PA_CLASSES, TEMPLATE_ARG2_UN_CLASSES), DCO_BIN_PA_TYPE_1, DCO_UN_TYPE_2, OPERATION, OPNAME) \
    DCO_BINARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_PA_CLASSES, TEMPLATE_ARG2_BIN_AA_CLASSES), DCO_BIN_PA_TYPE_1, DCO_BIN_AA_TYPE_2, OPERATION, OPNAME) \
    DCO_BINARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_PA_CLASSES, TEMPLATE_ARG2_BIN_AP_CLASSES), DCO_BIN_PA_TYPE_1, DCO_BIN_AP_TYPE_2, OPERATION, OPNAME) \
    DCO_BINARY_OPERATION_TEMPLATE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_PA_CLASSES, TEMPLATE_ARG2_BIN_PA_CLASSES), DCO_BIN_PA_TYPE_1, DCO_BIN_PA_TYPE_2, OPERATION, OPNAME) \
    \
    DCO_BINARY_AP_PA_OPERATION_TEMPLATE( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_ACTIVE_CLASSES), DCO_ACTIVE_TYPE_1, OPERATION, OPNAME) \
    DCO_BINARY_AP_PA_OPERATION_TEMPLATE( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_UN_CLASSES), DCO_UN_TYPE_1, OPERATION, OPNAME) \
    DCO_BINARY_AP_PA_OPERATION_TEMPLATE( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_BIN_AA_CLASSES), DCO_BIN_AA_TYPE_1, OPERATION, OPNAME) \
    DCO_BINARY_AP_PA_OPERATION_TEMPLATE( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_BIN_AP_CLASSES), DCO_BIN_AP_TYPE_1, OPERATION, OPNAME) \
    DCO_BINARY_AP_PA_OPERATION_TEMPLATE( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_BIN_PA_CLASSES), DCO_BIN_PA_TYPE_1, OPERATION, OPNAME)

    BIN_OPERATION(operator +,add)
    BIN_OPERATION(operator -,sub)
    BIN_OPERATION(operator *,mul)
    BIN_OPERATION(operator /,div)

    BIN_OPERATION(atan2,atan2)
    BIN_OPERATION(pow,pow)
    BIN_OPERATION(hypot,hypot)


    #define COMPARE_AA_DATATYPE(PREFIX, DT1, DT2, OPERATOR) \
    PREFIX static inline bool operator OPERATOR (const DT1& x1, const DT2& x2) { return x1._value OPERATOR x2._value; }

    #define COMPARE_AP_PA_DATATYPE(PREFIX, DT, OPERATOR) \
    PREFIX static inline bool operator OPERATOR (const DT& x1, const double& x2) { return x1._value OPERATOR x2; } \
    PREFIX static inline bool operator OPERATOR (const double& x1, const DT& x2) { return x1 OPERATOR x2._value; }

    #define COMPARE_AA(OPERATOR) \
    COMPARE_AA_DATATYPE( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_ACTIVE_CLASSES), DCO_ACTIVE_TYPE_1, DCO_ACTIVE_TYPE_1, OPERATOR) \
    COMPARE_AA_DATATYPE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_ACTIVE_CLASSES, TEMPLATE_ARG2_UN_CLASSES), DCO_ACTIVE_TYPE_1, DCO_UN_TYPE_2, OPERATOR) \
    COMPARE_AA_DATATYPE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_ACTIVE_CLASSES, TEMPLATE_ARG2_BIN_AA_CLASSES), DCO_ACTIVE_TYPE_1, DCO_BIN_AA_TYPE_2, OPERATOR) \
    COMPARE_AA_DATATYPE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_ACTIVE_CLASSES, TEMPLATE_ARG2_BIN_AP_CLASSES), DCO_ACTIVE_TYPE_1, DCO_BIN_AP_TYPE_2, OPERATOR) \
    COMPARE_AA_DATATYPE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_ACTIVE_CLASSES, TEMPLATE_ARG2_BIN_PA_CLASSES), DCO_ACTIVE_TYPE_1, DCO_BIN_PA_TYPE_2, OPERATOR) \
    \
    COMPARE_AA_DATATYPE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_UN_CLASSES, TEMPLATE_ARG2_ACTIVE_CLASSES), DCO_UN_TYPE_1, DCO_ACTIVE_TYPE_2, OPERATOR) \
    COMPARE_AA_DATATYPE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_UN_CLASSES, TEMPLATE_ARG2_UN_CLASSES), DCO_UN_TYPE_1, DCO_UN_TYPE_2, OPERATOR) \
    COMPARE_AA_DATATYPE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_UN_CLASSES, TEMPLATE_ARG2_BIN_AA_CLASSES), DCO_UN_TYPE_1, DCO_BIN_AA_TYPE_2, OPERATOR) \
    COMPARE_AA_DATATYPE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_UN_CLASSES, TEMPLATE_ARG2_BIN_AP_CLASSES), DCO_UN_TYPE_1, DCO_BIN_AP_TYPE_2, OPERATOR) \
    COMPARE_AA_DATATYPE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_UN_CLASSES, TEMPLATE_ARG2_BIN_PA_CLASSES), DCO_UN_TYPE_1, DCO_BIN_PA_TYPE_2, OPERATOR) \
    \
    COMPARE_AA_DATATYPE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_AA_CLASSES, TEMPLATE_ARG2_ACTIVE_CLASSES), DCO_BIN_AA_TYPE_1, DCO_ACTIVE_TYPE_2, OPERATOR) \
    COMPARE_AA_DATATYPE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_AA_CLASSES, TEMPLATE_ARG2_UN_CLASSES), DCO_BIN_AA_TYPE_1, DCO_UN_TYPE_2, OPERATOR) \
    COMPARE_AA_DATATYPE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_AA_CLASSES, TEMPLATE_ARG2_BIN_AA_CLASSES), DCO_BIN_AA_TYPE_1, DCO_BIN_AA_TYPE_2, OPERATOR) \
    COMPARE_AA_DATATYPE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_AA_CLASSES, TEMPLATE_ARG2_BIN_AP_CLASSES), DCO_BIN_AA_TYPE_1, DCO_BIN_AP_TYPE_2, OPERATOR) \
    COMPARE_AA_DATATYPE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_AA_CLASSES, TEMPLATE_ARG2_BIN_PA_CLASSES), DCO_BIN_AA_TYPE_1, DCO_BIN_PA_TYPE_2, OPERATOR) \
    \
    COMPARE_AA_DATATYPE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_AP_CLASSES, TEMPLATE_ARG2_ACTIVE_CLASSES), DCO_BIN_AP_TYPE_1, DCO_ACTIVE_TYPE_2, OPERATOR) \
    COMPARE_AA_DATATYPE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_AP_CLASSES, TEMPLATE_ARG2_UN_CLASSES), DCO_BIN_AP_TYPE_1, DCO_UN_TYPE_2, OPERATOR) \
    COMPARE_AA_DATATYPE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_AP_CLASSES, TEMPLATE_ARG2_BIN_AA_CLASSES), DCO_BIN_AP_TYPE_1, DCO_BIN_AA_TYPE_2, OPERATOR) \
    COMPARE_AA_DATATYPE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_AP_CLASSES, TEMPLATE_ARG2_BIN_AP_CLASSES), DCO_BIN_AP_TYPE_1, DCO_BIN_AP_TYPE_2, OPERATOR) \
    COMPARE_AA_DATATYPE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_AP_CLASSES, TEMPLATE_ARG2_BIN_PA_CLASSES), DCO_BIN_AP_TYPE_1, DCO_BIN_PA_TYPE_2, OPERATOR) \
    \
    COMPARE_AA_DATATYPE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_PA_CLASSES, TEMPLATE_ARG2_ACTIVE_CLASSES), DCO_BIN_PA_TYPE_1, DCO_ACTIVE_TYPE_2, OPERATOR) \
    COMPARE_AA_DATATYPE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_PA_CLASSES, TEMPLATE_ARG2_UN_CLASSES), DCO_BIN_PA_TYPE_1, DCO_UN_TYPE_2, OPERATOR) \
    COMPARE_AA_DATATYPE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_PA_CLASSES, TEMPLATE_ARG2_BIN_AA_CLASSES), DCO_BIN_PA_TYPE_1, DCO_BIN_AA_TYPE_2, OPERATOR) \
    COMPARE_AA_DATATYPE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_PA_CLASSES, TEMPLATE_ARG2_BIN_AP_CLASSES), DCO_BIN_PA_TYPE_1, DCO_BIN_AP_TYPE_2, OPERATOR) \
    COMPARE_AA_DATATYPE( TEMPLATE_BUILDER_2(TEMPLATE_ARG1_BIN_PA_CLASSES, TEMPLATE_ARG2_BIN_PA_CLASSES), DCO_BIN_PA_TYPE_1, DCO_BIN_PA_TYPE_2, OPERATOR) \
    \
    COMPARE_AP_PA_DATATYPE( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_ACTIVE_CLASSES), DCO_ACTIVE_TYPE_1, OPERATOR) \
    COMPARE_AP_PA_DATATYPE( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_UN_CLASSES), DCO_UN_TYPE_1, OPERATOR) \
    COMPARE_AP_PA_DATATYPE( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_BIN_AA_CLASSES), DCO_BIN_AA_TYPE_1, OPERATOR) \
    COMPARE_AP_PA_DATATYPE( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_BIN_AP_CLASSES), DCO_BIN_AP_TYPE_1, OPERATOR) \
    COMPARE_AP_PA_DATATYPE( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_BIN_PA_CLASSES), DCO_BIN_PA_TYPE_1, OPERATOR)


    COMPARE_AA(==) COMPARE_AA(!=) COMPARE_AA(<) COMPARE_AA(<=) COMPARE_AA(>) COMPARE_AA(>=)

    using ::isnan;
    #define DCO_ISNAN_OPERATION(PREFIX, DT) \
    PREFIX static inline bool isnan(const DT& x) { \
        return isnan(x._value); \
    }
    DCO_ISNAN_OPERATION( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_ACTIVE_CLASSES), DCO_ACTIVE_TYPE_1) \
    DCO_ISNAN_OPERATION( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_BIN_AA_CLASSES), DCO_BIN_AA_TYPE_1) \
    DCO_ISNAN_OPERATION( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_BIN_AP_CLASSES), DCO_BIN_AP_TYPE_1) \
    DCO_ISNAN_OPERATION( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_BIN_PA_CLASSES), DCO_BIN_PA_TYPE_1) \
    DCO_ISNAN_OPERATION( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_UN_CLASSES), DCO_UN_TYPE_1)



    #define DCO_OUTSTREAM_OPERATION(PREFIX, DT) \
    PREFIX static inline std::ostream& operator << (std::ostream& out, const DT& x) { \
        out << x._value; \
        return out; \
    }

    DCO_OUTSTREAM_OPERATION( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_ACTIVE_CLASSES), DCO_ACTIVE_TYPE_1) \
    DCO_OUTSTREAM_OPERATION( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_BIN_AA_CLASSES), DCO_BIN_AA_TYPE_1) \
    DCO_OUTSTREAM_OPERATION( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_BIN_AP_CLASSES), DCO_BIN_AP_TYPE_1) \
    DCO_OUTSTREAM_OPERATION( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_BIN_PA_CLASSES), DCO_BIN_PA_TYPE_1) \
    DCO_OUTSTREAM_OPERATION( TEMPLATE_BUILDER_1(TEMPLATE_ARG1_UN_CLASSES), DCO_UN_TYPE_1)

    struct POSITION {
      int stackcounter;
      int progvarcounter;

      POSITION() : stackcounter(0), progvarcounter(0) {}
      POSITION(const int nstackcounter, const int nprogvarcounter) : stackcounter(nstackcounter), progvarcounter(nprogvarcounter) {}
    };

    template<class DCO_TAPE_REAL>
    struct static_tape {
      struct TAPE_ENTRY {
        int arg;
        DCO_TAPE_REAL  pval;
      };
      TAPE_ENTRY    *stack;
      TAPE_ENTRY    *topOfStack;
      DCO_TAPE_REAL *adjoints;
      int           progvarcounter;
      int           local_edgecount;

      typedef     void (*EXTERNALFUNCTION)(static_tape &caller, int mode, void *userdata);

      struct external_functionhandler {
         POSITION            position;
         EXTERNALFUNCTION    handler;
         void                *userdata;

         external_functionhandler() : position(0,0) {}
      };

      external_functionhandler    external_functions[1000000];
      int                         external_function_count;

      static_tape(const unsigned int size, int numprogvars=0) : stack(0), topOfStack(0), adjoints(0), progvarcounter(0), local_edgecount(0) {
        stack=new TAPE_ENTRY[size];
        topOfStack=stack;
        progvarcounter=0;
        if(numprogvars==0)
          numprogvars=size/2;
        adjoints=new DCO_TAPE_REAL[numprogvars];

        for(unsigned int i=0;i<size;++i) {
          stack[i].arg=0;
          stack[i].pval=0;
        }

        for(int i=0;i<numprogvars;++i)
          adjoints[i]=0;

        external_function_count=0;

        //cout << "Size of one stack-entry = " << sizeof(TAPE_ENTRY) << "; Memory usage = " << size*sizeof(TAPE_ENTRY)/1024/1024 << " MB" << endl;
        //cout << "Size of one adjoint-entry = " << sizeof(double) << "; Memory usage = " << numprogvars*sizeof(double)/1024/1024 << " MB" << endl;

      }

      inline ~static_tape() {
        delete [] stack;
        delete [] adjoints;
        stack=topOfStack=0;
        adjoints=0;
      }

      POSITION    get_position() {
        return POSITION( static_cast<int>(topOfStack-stack), progvarcounter);
      }

      inline void interpret_reverse_internal_plain(const POSITION &from, const POSITION &to, const bool shrinkTape=true) {
        TAPE_ENTRY *last = stack + to.stackcounter;
        TAPE_ENTRY *cur = stack + from.stackcounter;
        //register DCO_TAPE_REAL *adjoints = this->adjoints;
        register int progvaridx = from.progvarcounter;

        while (cur-- != last) {
          const DCO_TAPE_REAL &adj = adjoints[progvaridx--];
          const int &edgecount = cur->arg;
          if(false)  {
              switch(edgecount) {
              case 0:
                  break;
              case 1:
                  --cur; adjoints[cur->arg] += adj * cur->pval;
                  break;
              case 2:
                  --cur; adjoints[cur->arg] += adj * cur->pval;
                  --cur; adjoints[cur->arg] += adj * cur->pval;
                  break;
              case 3:
                  adjoints[cur[0].arg] += adj * cur[0].pval;
                  adjoints[cur[1].arg] += adj * cur[1].pval;
                  adjoints[cur[2].arg] += adj * cur[2].pval;
                  /*--cur; adjoints[cur->arg] += adj * cur->pval;
                  --cur; adjoints[cur->arg] += adj * cur->pval;*/
                  cur -= 3;
                  break;
              }
          }
          else {
              for (int j = 0; j < edgecount; j++) {
                --cur;
                adjoints[cur->arg] += adj * cur->pval;
              }
          }
        }
        if(shrinkTape) {
          progvarcounter = progvaridx;
          topOfStack = last;
        }
      }

      inline void interpret_reverse_internal(const POSITION &from, const POSITION &to, const bool shrinkTape=true) {
           int external_first=-1;
           int external_count=0;

         for(int i=external_function_count-1;i>=0;--i) {
             if( (external_functions[i].position.stackcounter <= from.stackcounter) &&
                 (external_functions[i].position.stackcounter >= to.stackcounter) ) {
                 if(external_first==-1) external_first=i;
                 ++external_count;
             }
         }

        POSITION myfrom=from;
         POSITION myto(0,0);

        for(int i=external_first;external_count>0;--external_count) {
             myto = external_functions[i].position;
             interpret_reverse_internal_plain( myfrom, myto, shrinkTape);
             external_functions[i].handler(*this, 0, external_functions[i].userdata);
             myfrom=myto;
        }
        interpret_reverse_internal_plain(myfrom,to,shrinkTape);
      }

      inline void register_external_function(EXTERNALFUNCTION handler, void *userdata) {
         external_functions[external_function_count].position = get_position();
         external_functions[external_function_count].handler = handler;
         external_functions[external_function_count].userdata=userdata;
         ++external_function_count;
	std::cout << "external_function_count=" << external_function_count << std::endl;
      }

      inline void interpret_reverse() {
        POSITION to(0,0);
        interpret_reverse_internal(get_position(),to,false);
      }
      inline void interpret_reverse_to(const POSITION &to) {
        interpret_reverse_internal(get_position(),to,false);
      }
      inline void interpret_reverse_from(const POSITION &from) {
        POSITION to(0,0);
        interpret_reverse_internal(from,to,false);
      }
      inline void interpret_reverse_from_to(const POSITION &from, const POSITION &to) { interpret_reverse_internal(from,to,false); }

      inline void interpret_and_reset_reverse() {
        POSITION save(get_position());
        POSITION to(0,0);
        interpret_reverse_internal(save,to,true);
        zero_adjoints_internal(save, to);
      }
      inline void interpret_and_reset_reverse_to(const POSITION &to) {
        POSITION save(get_position());
      interpret_reverse_internal(save, to, true);
        zero_adjoints_internal(save, to);
      }

      inline void reset_to(const POSITION &to) {
        zero_adjoints_internal(get_position(), to);
        this->progvarcounter=to.progvarcounter;
        this->topOfStack = this->stack + to.stackcounter;
      }

      inline void reset() {
        reset_to(POSITION(0,0));
      }

      inline void zero_adjoints_internal(const POSITION &from, const POSITION &to) {
        for(int i=from.progvarcounter;i>to.progvarcounter;--i) {
          adjoints[i]=0;
        }
      }

      inline void zero_adjoints() { POSITION to(0,0); zero_adjoints_internal(get_position(), to); }
      inline void zero_adjoints_to(const POSITION &to) { zero_adjoints_internal(get_position(),to); }
      inline void zero_adjoints_from(const POSITION &from) { POSITION to(0,0); zero_adjoints_internal(from,to); }
      inline void zero_adjoints_from_to(const POSITION &from, const POSITION &to) { zero_adjoints_internal(from,to); }

        template<class DCO_DATA, class DCO_HANDLER>
        inline void register_variable( active_type<DCO_TAPE_REAL, DCO_DATA, DCO_HANDLER> &x) {
            x._data.tape_index=++progvarcounter;
            x._data.tape=this;
            x._edgecount=1;
            topOfStack->arg=0;
            topOfStack->pval=0;
            ++topOfStack;
        }

        template<class DCO_DATA>
        inline void new_edge(const DCO_DATA &x, const DCO_TAPE_REAL &pval) {
            topOfStack->arg=x.tape_index;
            topOfStack->pval=pval;
            ++topOfStack;
            ++local_edgecount;
        }

        template<class DCO_DATA>
        inline void new_edge_intern(const DCO_DATA &x, const DCO_TAPE_REAL &pval) {
            topOfStack->arg=x.tape_index;
            topOfStack->pval=pval;
            ++topOfStack;
        }


        template<class DCO_DATA>
        inline bool new_progvar(DCO_DATA &x) {
            x.tape_index=++progvarcounter;
            x.tape=this;
            topOfStack->arg=local_edgecount; local_edgecount=0;
            topOfStack->pval=0;
            ++topOfStack;
            return true;
        }

        template<class DCO_DATA>
        inline bool new_progvar(DCO_DATA &x, const DCO_EDGE_COUNT_TYPE &edgecount) {
            x.tape_index=++progvarcounter;
            x.tape=this;
            topOfStack->arg=edgecount;
            topOfStack->pval=0;
            ++topOfStack;
            return true;
        }


    };

    template<class DCO_TAPE_REAL, class DCO_DATA, class DCO_SHARED_DATA>struct statement_handler {
        typedef active_type<DCO_TAPE_REAL, DCO_DATA, statement_handler<DCO_TAPE_REAL, DCO_DATA, DCO_SHARED_DATA> > type;
        //const DCO_EDGE_COUNT_TYPE &edgecount;
        DCO_SHARED_DATA    *shared_object;

        /*statement_handler(const DCO_EDGE_COUNT_TYPE &nedgecount) : edgecount(nedgecount), shared_object(0)  {
        }*/

        template<class T, class T2> inline statement_handler(const T &x, T2 &ret) :  shared_object(0) {
            this->interpret(x,1.0);
            shared_object->new_progvar(ret._data, x._edgecount);
        }

        /*DCO_SHARED_DATA    *shared_object;
        statement_handler(const DCO_EDGE_COUNT_TYPE &nedgecount) : shared_object(0)  {}*/

        inline void interpret(const type &x, const DCO_TAPE_REAL &pval) {
#ifdef DCO_ACTIVITY
            if(x._edgecount>0) shared_object = DCO_DATA::new_edge(x._data, pval);
#else
            shared_object = DCO_DATA::new_edge(x._data, pval);
#endif
        }

        //TEMPLATE_BUILDER_CONSTUCTOR( TEMPLATE_ARG1_BIN_AA_CLASSES)
        //inline void interpret(const DCO_BIN_AA_TYPE_1 &x, const DCO_TAPE_REAL &pval) __attribute__((always_inline));
        TEMPLATE_BUILDER_CONSTUCTOR( TEMPLATE_ARG1_BIN_AA_CLASSES)
        inline void interpret(const DCO_BIN_AA_TYPE_1 &x, const DCO_TAPE_REAL &pval) {
            this->interpret(x._arg1, x.pval1(x._value, x._arg1,x._arg2, pval));
            this->interpret(x._arg2, x.pval2(x._value, x._arg1,x._arg2, pval));
        }

        TEMPLATE_BUILDER_CONSTUCTOR( TEMPLATE_ARG1_UN_CLASSES)
        inline void interpret(const DCO_UN_TYPE_1 &x, const DCO_TAPE_REAL &pval) {
            this->interpret(x._arg, x.pval(x._value,x._arg, pval));
        }

        TEMPLATE_BUILDER_CONSTUCTOR( TEMPLATE_ARG1_BIN_AP_CLASSES)
        inline void interpret(const DCO_BIN_AP_TYPE_1 &x, const DCO_TAPE_REAL &pval) {
            this->interpret(x._arg1, x.pval1(x._value,x._arg1,x._arg2, pval));
        }

        TEMPLATE_BUILDER_CONSTUCTOR( TEMPLATE_ARG1_BIN_PA_CLASSES)
        inline void interpret(const DCO_BIN_PA_TYPE_1 &x, const DCO_TAPE_REAL &pval) {
            this->interpret(x._arg2, x.pval2(x._value,x._arg1,x._arg2, pval));
        }

    };



    template<class DCO_TAPE_REAL>struct data_object_tapeonly {
        int         tape_index;
        static_tape<DCO_TAPE_REAL> *tape;

        static inline static_tape<DCO_TAPE_REAL>* new_edge(const data_object_tapeonly &x, const DCO_TAPE_REAL &pval) {
            x.tape->new_edge_intern(x,pval);
            return x.tape;
        }
    };

    /////////////////////////
    template<class DCO_TAPE_REAL, const int vecsize>    struct tlm_buffer;
    template<class DCO_TAPE_REAL, const int vecsize>    struct data_object_tlmVecktor;

    template<class DCO_TAPE_REAL, const int vecsize>
    struct tlm_buffer {
        typedef data_object_tlmVecktor<DCO_TAPE_REAL, vecsize> data;
        struct edge {
            const data_object_tlmVecktor<DCO_TAPE_REAL, vecsize> *node;
            DCO_TAPE_REAL pval;
        };
        edge    edges[40];
        int     edgecount;

        tlm_buffer() : edgecount(0) {}

        inline void new_edge(const data_object_tlmVecktor<DCO_TAPE_REAL, vecsize> &node, const DCO_TAPE_REAL &pval) {
            edges[edgecount].node = &node;
            edges[edgecount].pval = pval;
            ++edgecount;
        }


        template<class T>
        inline void activate_variable(T &x) {
            x._edgecount=1;
            x._data.buffer=this;
            for(int i=0;i<vecsize;++i)
                x._data.tlms[i]=0;
        }

        //__attribute__((always_inline))
        inline bool new_progvar(data &x, const int nedgecount DCO_UNUSED)  {
            if(edgecount==0) return false;
            x.buffer=this;
            //cout << "edgecount=" << edgecount << endl;

            for(int i=0;i<vecsize;++i) {
                switch(edgecount) {
                    case 1:
                        x.tlms[i]= edges[0].node->tlms[i] * edges[0].pval;
                        break;
                    case 2:
                        x.tlms[i]= edges[0].node->tlms[i] * edges[0].pval +
                                   edges[1].node->tlms[i] * edges[1].pval;
                        break;

                    case 3:
                        x.tlms[i]= edges[0].node->tlms[i] * edges[0].pval +
                                   edges[1].node->tlms[i] * edges[1].pval +
                                   edges[2].node->tlms[i] * edges[2].pval;
                        break;

                    case 4:
                        x.tlms[i]= edges[0].node->tlms[i] * edges[0].pval +
                                   edges[1].node->tlms[i] * edges[1].pval +
                                   edges[2].node->tlms[i] * edges[2].pval +
                                   edges[3].node->tlms[i] * edges[3].pval;
                        break;

                    default:
                        //cout << "edgecount=" << edgecount << endl;
                        DCO_TAPE_REAL d=0;
                        for(int j=0;j<edgecount;++j) {
                            d += edges[j].node->tlms[i] * edges[j].pval;
                        }
                        x.tlms[i]=d;
                }

            }
            edgecount=0;
            return true;
        }

    };

    template<class DCO_TAPE_REAL, const int vecsize>
    struct tlm_buffer_new {
        typedef data_object_tlmVecktor<DCO_TAPE_REAL, vecsize> data;
        DCO_TAPE_REAL   buffer[vecsize];

        tlm_buffer_new()  {
            for(int i=0;i<vecsize;++i) buffer[i]=0;
        }

        inline void new_edge(const data_object_tlmVecktor<DCO_TAPE_REAL, vecsize> &node, const DCO_TAPE_REAL &pval) {
            for(int i=0;i<vecsize;++i)
                buffer[i] += pval * node.tlms[i];
        }

        template<class T>
        inline void activate_variable(T &x) {
            x._edgecount=1;
            x._data.buffer=this;
            for(int i=0;i<vecsize;++i)
                x._data.tlms[i]=0;
        }

        //__attribute__((always_inline))
        inline bool new_progvar(data &x, const int edgecount DCO_UNUSED)  {
            if(edgecount==0) return false;
            x.buffer=this;
            for(int i=0;i<vecsize;++i) {
                x.tlms[i] = buffer[i]; buffer[i]=0;
            }
            return true;
        }
    };


    template<class DCO_TAPE_REAL, const int vecsize>struct data_object_tlmVecktor {
        double      tlms[vecsize];
        tlm_buffer_new<DCO_TAPE_REAL, vecsize> *buffer;

        static inline tlm_buffer_new<DCO_TAPE_REAL, vecsize>* new_edge(const data_object_tlmVecktor<DCO_TAPE_REAL, vecsize> &x, const DCO_TAPE_REAL &pval) {
            x.buffer->new_edge(x,pval);
            return x.buffer;
        }
    };

    template<class T>
    struct access_exception{
        const T    &argument;

        access_exception(const T &arg) : argument(arg) {

        }
    };


  namespace t1s {
    typedef tlm_type<double> type;
  }

  namespace t1v {

    //typedef tlm_buffer<double, 1> buffer;
    typedef tlm_buffer_new<double, 1> buffer;
    typedef data_object_tlmVecktor<double, 1> data;
    typedef statement_handler<double, data, buffer> statement_handler;
    typedef active_type<double, data, statement_handler > type;

    static inline void set(type &x, const double val, const int k0, const int subidx=0) {
      if(k0==0)
	x._value=val;
      else if(k0==1)
	{
	  x._data.tlms[subidx]=val;
	}
    }

    static inline void get(const type &x, double &val, const int k0,const int subidx=0) {
      if(k0==0)
	val=x._value;
      else if(k0==1)
	{
	  val=x._data.tlms[subidx];
	}
    }

  }
  
  namespace a1s {
    typedef static_tape<double> static_tape;
    typedef data_object_tapeonly<double> data;
    typedef statement_handler<double, data, static_tape> statement_handler;
    
    typedef active_type<double, data, statement_handler > type;
    
    static inline void set(type &x, const double &val, const int k1) {
      if(k1==-1) {
	if(x._edgecount>0) {
	  x._data.tape->adjoints[x._data.tape_index]=val;
	}
      }
    }
    static inline void set(type &x, const double &val) {
      x._value=val;
    }
    static inline void get(const type &x, double &val, const int k1) {
      if(k1==-1){
	if(x._edgecount>0) {
	  val=x._data.tape->adjoints[x._data.tape_index];
	}
      }
      else {
	throw new access_exception<type>(x);
      }
    }
    static inline void get(const type &x, double &val) {
      val=x._value;
    }
  }
  
  

    namespace t2s_a1s {
        typedef static_tape<t1s::type> tape;
        typedef data_object_tapeonly<t1s::type> data;
        typedef statement_handler<t1s::type, data, tape> statement_handler;
        typedef active_type<t1s::type, data, statement_handler > type;
    }
}

