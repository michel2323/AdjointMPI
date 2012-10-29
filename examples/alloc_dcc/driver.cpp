#include <mpi.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ampi.h>
using namespace std;

extern void alloc(double *, int &, double *); 
extern void b1_alloc(int &, double *, double *, int &, double *, double *);
extern void d2_b1_alloc(int &, double *, double *, double *, double *, int &, double *, double *, double *, double*);

int main(int argc, char** argv){
    int size=2;
    int bmode=1;
    double myres[1];
    double b1_myres[1];
    double *x;
    double *b1_x;
    x=new double[size];
    b1_x=new double[size];
    b1_myres[0]=1;
    b1_alloc(bmode,myres,b1_myres,size,x,b1_x); 

} 

