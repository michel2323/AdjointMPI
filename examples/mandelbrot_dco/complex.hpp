#ifndef COMPLEX_H
#define COMPLEX_H
#include <math.h>
#include "dco_tape.hpp"

class complex {
    public:
	active r;
	active i;
	complex(const active&, const active&);
	complex();
	complex& operator=(const complex);
};


complex complex_pow(complex&);
active complex_norm(complex&);

complex operator+(const complex&, const complex&);
complex operator-(const complex&, const complex&);
complex operator*(const complex&, const complex&);
#endif
