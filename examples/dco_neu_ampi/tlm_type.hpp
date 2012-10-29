    // TLM Type für pure tangent mode or tangent over adjoint etc.
    template<class DCO_TLMREAL>
    struct tlm_type{
        DCO_TLMREAL v;
        DCO_TLMREAL d;

        tlm_type() : v(0), d(0) {}
        tlm_type(const double &vneu) : v(vneu), d(0) {}
        //tlm_type(const tlm_type &vneu) : v(vneu.v), d(vneu.d) {}

        tlm_type(const DCO_TLMREAL &vneu, const DCO_TLMREAL &dneu) : v(vneu), d(dneu) {}

        inline const tlm_type& operator = (const DCO_TLMREAL &vneu) {
            v=vneu; d=0;
            return *this;
        }

        tlm_type& operator +=(const double &vneu) { this->v += vneu; return *this; }
        tlm_type& operator -=(const double &vneu) { this->v -= vneu; return *this; }

        tlm_type& operator *=(const double &vneu)   { this->d *= vneu;                                          this->v*=vneu;  return *this; }
        tlm_type& operator /=(const double &vneu)   { this->d /= vneu;                                          this->v/=vneu;  return *this; }

        tlm_type& operator +=(const tlm_type &vneu) { this->v += vneu.v; this->d += vneu.d; return *this; }
        tlm_type& operator -=(const tlm_type &vneu) { this->v -= vneu.v; this->d -= vneu.d; return *this; }

        tlm_type& operator *=(const tlm_type &vneu) { this->d = this->d * vneu.v + vneu.d * this->v;            this->v*=vneu.v; return *this; }
        tlm_type& operator /=(const tlm_type &vneu) { this->d = this->d/vneu.v-vneu.d*this->v/vneu.v/vneu.v;    this->v/=vneu.v; return *this; }

        tlm_type& operator --() { --this->v; return *this; }
        tlm_type operator --(int) { tlm_type ret(this->v--,this->d); return ret; }
        tlm_type& operator ++() { ++this->v; return *this; }
        tlm_type operator ++(int) { tlm_type ret(this->v++,this->d); return ret; }
    };

    template<class DCO_TLMREAL>
    static inline tlm_type<DCO_TLMREAL> operator -(const tlm_type<DCO_TLMREAL> &x) {
        return tlm_type<DCO_TLMREAL>(-x.v,-x.d);
    }

    //binary pure active operations
    #define BINOPERATOR_ACTIVE(OPERATORX, PVAL) \
    template<class DCO_TLMREAL> \
    static inline tlm_type<DCO_TLMREAL> operator OPERATORX (const tlm_type<DCO_TLMREAL> &x1, const tlm_type<DCO_TLMREAL> &x2) { \
        return tlm_type<DCO_TLMREAL>(x1.v OPERATORX x2.v, PVAL); \
    }
    BINOPERATOR_ACTIVE(+,x1.d + x2.d)
    BINOPERATOR_ACTIVE(-,x1.d - x2.d)
    BINOPERATOR_ACTIVE(*,x1.d * x2.v + x2.d * x1.v)
    BINOPERATOR_ACTIVE(/,x1.d/x2.v-x2.d*x1.v/x2.v/x2.v)
    #undef BINOPERATOR_ACTIVE

    #define BINOPERATOR_DATATYPE(OPERATORX, PVAL1, PVAL2) \
    template<class DCO_TLMREAL> \
    static inline tlm_type<DCO_TLMREAL> operator OPERATORX(const tlm_type<DCO_TLMREAL> &x1, const double x2) { return tlm_type<DCO_TLMREAL>(x1.v OPERATORX x2, PVAL1); } \
    template<class DCO_TLMREAL> \
    static inline tlm_type<DCO_TLMREAL> operator OPERATORX(const double &x1, const tlm_type<DCO_TLMREAL> x2) { return tlm_type<DCO_TLMREAL>(x1 OPERATORX x2.v, PVAL2); }

    BINOPERATOR_DATATYPE( +, x1.d, x2.d)
    BINOPERATOR_DATATYPE( -, x1.d, -x2.d)
    BINOPERATOR_DATATYPE( *, x1.d*x2,x2.d*x1)
    BINOPERATOR_DATATYPE( /, x1.d/x2, -x2.d*x1/x2.v/x2.v)
    #undef BINOPERATOR_DATATYPE

    //compare operator <,<= == != etc.
    #define COMPRARE_AA(OPERATOR) \
    template<class DCO_TLMREAL> static inline bool operator OPERATOR (const tlm_type<DCO_TLMREAL>& x1, const tlm_type<DCO_TLMREAL>& x2) { return x1.v OPERATOR x2.v; } \
    template<class DCO_TLMREAL> static inline bool operator OPERATOR (const tlm_type<DCO_TLMREAL>& x1, const double& x2) { return x1.v OPERATOR x2; } \
    template<class DCO_TLMREAL> static inline bool operator OPERATOR (const double& x1, const tlm_type<DCO_TLMREAL>& x2) { return x1 OPERATOR x2.v; }
    COMPRARE_AA(==) COMPRARE_AA(!=) COMPRARE_AA(<=) COMPRARE_AA(<) COMPRARE_AA(>) COMPRARE_AA(>=)
    #undef COMPRARE_AA

    //now sin,cos etc.
    #define UNARYFUNCTION(FUU,DFUU,PREPROCESS, DEFINE)				\
      using DEFINE::FUU;							\
      template<class DCO_TLMREAL>static inline tlm_type<DCO_TLMREAL> FUU(const tlm_type<DCO_TLMREAL> &x) { PREPROCESS return tlm_type<DCO_TLMREAL>(FUU(x.v),DFUU*x.d); }

    UNARYFUNCTION(sin, cos(x.v),,std)
    UNARYFUNCTION(cos, -sin(x.v),,std)
    UNARYFUNCTION(tan, (1.0+(tan(x.v)*tan(x.v))),,std)
    UNARYFUNCTION(cosh, sinh(x.v),,std)
    UNARYFUNCTION(sinh, cosh(x.v),,std)
    UNARYFUNCTION(asin, 1/sqrt(1.0-x.v*x.v),,std)
    UNARYFUNCTION(acos, -1/sqrt(1.0-x.v*x.v),,std)
    UNARYFUNCTION(exp, exp(x.v),,std)
    UNARYFUNCTION(atan, 1.0/(1.0+x.v*x.v),,std)
    UNARYFUNCTION(tanh, 1.0/cosh(x.v)/cosh(x.v),,std)
    UNARYFUNCTION(sqrt, 1.0/(2.0*sqrt(x.v)),,std)
    UNARYFUNCTION(log, 1.0/x.v,,std)
    #ifndef MSVC
    UNARYFUNCTION(log1p, 1.0/(x.v+1),,)
    UNARYFUNCTION(log10, 1.0/(x.v*log(10)),,)
    #endif
    #undef UNARYFUNCTION

    #define DUMMYFUNCTION(name) \
    template<class DCO_TLMREAL>static inline double name(const tlm_type<DCO_TLMREAL> &x) { \
        return name(x.v); \
    }
    DUMMYFUNCTION(ceil)
    DUMMYFUNCTION(floor)
    #undef DUMMYFUNCTION

    using std::isnan;
    template<class DCO_TLMREAL> static inline bool isnan(const tlm_type<DCO_TLMREAL> &x) {
      return isnan(x.v);
    }

    template<class DCO_TLMREAL> static inline const tlm_type<DCO_TLMREAL> fabs(const tlm_type<DCO_TLMREAL> &x) {
        if(x.v<0) return -x;
        return x;
    }

    template<class DCO_TLMREAL> static inline const tlm_type<DCO_TLMREAL> abs(const tlm_type<DCO_TLMREAL> &x) {
      return fabs(x);
    }

    template<class DCO_TLMREAL>
    static inline std::ostream& operator<<(std::ostream& out, const tlm_type<DCO_TLMREAL>& x) {
        out << x.v;
        return out;
    }

    using std::atan2;

    template<class DCO_TLMREAL> static inline tlm_type<DCO_TLMREAL> atan2(const tlm_type<DCO_TLMREAL> &y, const tlm_type<DCO_TLMREAL> &x) {
        if(x.v==0) return tlm_type<DCO_TLMREAL>( atan2(y.v,x.v),0);
        DCO_TLMREAL tmp=1/(1+(y.v/x.v)*(y.v/x.v));
        return tlm_type<DCO_TLMREAL>(atan2(y.v,x.v), y.d * tmp/x.v - x.d*tmp*y.v/x.v/x.v);
    }
    template<class DCO_TLMREAL> static inline tlm_type<DCO_TLMREAL> atan2(const tlm_type<DCO_TLMREAL> &y, const double &x) {
        if(x==0) return tlm_type<DCO_TLMREAL>( atan2(y.v,x),0);
        DCO_TLMREAL tmp=1/(1+(y.v/x)*(y.v/x));
        return tlm_type<DCO_TLMREAL>(atan2(y.v,x), y.d * tmp/x);
    }
    template<class DCO_TLMREAL> static inline tlm_type<DCO_TLMREAL> atan2(const double &y, const tlm_type<DCO_TLMREAL> &x) {
        if(x.v==0) return tlm_type<DCO_TLMREAL>( atan2(y,x.v),0);
        DCO_TLMREAL tmp=1/(1+(y/x.v)*(y/x.v));
        return tlm_type<DCO_TLMREAL>(atan2(y,x.v), - x.d*tmp*y/x.v/x.v);
    }

    using std::pow;
    template<class DCO_TLMREAL> static inline tlm_type<DCO_TLMREAL> pow(const tlm_type<DCO_TLMREAL> &x, const tlm_type<DCO_TLMREAL> &y) {
        DCO_TLMREAL tmp=pow(x.v,y.v);
        return tlm_type<DCO_TLMREAL>( tmp, x.d * y.v * tmp/x.v + y.d * tmp * log(x.v));
    }
    template<class DCO_TLMREAL> static inline tlm_type<DCO_TLMREAL> pow(const tlm_type<DCO_TLMREAL> &x, const double &y) {
        DCO_TLMREAL tmp=pow(x.v,y);
        return tlm_type<DCO_TLMREAL>( tmp, x.d * y * tmp/x.v );
    }
    template<class DCO_TLMREAL> static inline tlm_type<DCO_TLMREAL> pow(const double &x, const tlm_type<DCO_TLMREAL> &y) {
        DCO_TLMREAL tmp=pow(x,y.v);
        return tlm_type<DCO_TLMREAL>( tmp, y.d * tmp * log(x));
    }

    using ::hypot;
    template<class DCO_TLMREAL> static inline tlm_type<DCO_TLMREAL> hypot(const tlm_type<DCO_TLMREAL> &x1, const tlm_type<DCO_TLMREAL> &x2) {
        DCO_TLMREAL tmp=hypot(x1.v,x2.v);
        return tlm_type<DCO_TLMREAL>( tmp, x1.d * x1.v/tmp + x2.d * x2.v/tmp);
    }
    template<class DCO_TLMREAL> static inline tlm_type<DCO_TLMREAL> hypot(const tlm_type<DCO_TLMREAL> &x1, const double &x2) {
        DCO_TLMREAL tmp=hypot(x1.v,x2);
        return tlm_type<DCO_TLMREAL>( tmp, x1.d * x1.v/tmp);
    }
    template<class DCO_TLMREAL> static inline tlm_type<DCO_TLMREAL> hypot(const double &x1, const tlm_type<DCO_TLMREAL> &x2) {
        DCO_TLMREAL tmp=hypot(x1,x2.v);
        return tlm_type<DCO_TLMREAL>( tmp, x2.d * x2.v/tmp);
    }

