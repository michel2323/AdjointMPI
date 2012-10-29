#include "dcc_mpi.hpp"
#include <iostream>


void alloc(double *myres, int &size, double *x) 
//$ad indep x 
//$ad dep myres 
{
    int i=0;
    double *buf;
    double *buf2;
    myres[0]=1;
    dcc_new(buf,size); 
    for(i=0; i<size; i=i+1){ buf[i]=x[i]; }
    for(i=0; i<size; i=i+1) {
	buf[i]=i+2.0;
	print_num(buf[i]);
	myres[0]=myres[0]*buf[i];
    }
    dcc_delete(buf);
    size=size+5;
    dcc_new(buf2,size); 
    for(i=0; i<size; i=i+1){ buf2[i]=x[i]; }
    for(i=0; i<size; i=i+1) {
	buf2[i]=i+4.0;
	print_num(buf2[i]);
	myres[0]=myres[0]*buf2[i];
    }
    dcc_delete(buf2);
}
