#ifndef F_HPP
#define F_HPP
#include <dco.hpp>
#include <ampi_tape.hpp>
typedef dco::a1s::type active;



void f(int& nx, 
	int& nt, 
	int& stride, 
	active& delta_t, 
	active& a,
	active& b,
	active& c,
	active* temp,
	active* temp_obs,
	active& cost); 
void sf(int& nx, 
	int& nt, 
	int& stride, 
	active& delta_t, 
	active& a,
	active& b,
	active& c,
	active* temp,
	active* temp_obs,
	active& cost); 
//$ad indep temp
//$ad dep cost
#endif
