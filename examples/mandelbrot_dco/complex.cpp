#include "complex.hpp"
#include "dco_tape.hpp"

complex complex_pow(complex& x1) {
    complex res;
    active two;
    two = 2.0;
    res.r = (x1.r * x1.r) - (x1.i * x1.i);
    res.i = two * (x1.r * x1.i);
    return res;
}

active complex_norm(complex& x1) {
    return (x1.r * x1.r) + (x1.i * x1.i);
}

complex operator+(const complex& x1, const complex& x2) {
    complex  res;
    res.r = (x1.r + x2.r);
    res.i = (x1.i + x2.i);
    return res;
}

complex operator-(const complex& x1, const complex& x2) {
    complex res;
    res.r = (x1.r - x2.r);
    res.i = (x1.i - x2.i);
    return res;
}

complex operator*(const complex& x1, const complex& x2) {
    complex res;
    res.r = (x1.r * x2.r) - (x1.i * x2.i);
    //	res.r = (x1.r * x2.r) - (x1.i * x2.i);
    //	res.i = (x1.r * x2.i) + (x2.r * x1.i);
    res.i = (x1.r * x2.i) + (x2.r * x1.i);
    return res;
}

complex::complex(const active& r, const active& i) {
    this->r = r;
    this->i = i;
}

complex& complex::operator=(const complex y) {
    this->r = y.r;
    this->i = y.i;
    return *this;
}

complex::complex() {
    this->r = 0;
    this->i = 0;
}
