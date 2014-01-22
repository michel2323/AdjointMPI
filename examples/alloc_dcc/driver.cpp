#include <mpi.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ampi.h>
using namespace std;

extern void alloc(double *, int &, double *); 
extern void a1_alloc(int &, double *, double *, int &, double *, double *);
extern void t2_a1_alloc(int &, double *, double *, double *, double *, int &, double *, double *, double *, double*);

int main(int argc, char** argv){
    int size=2;
    int bmode=1;
    double myres[1];
    double a1_myres[1];
    double *x;
    double *a1_x;
    x=new double[size];
    a1_x=new double[size];
    a1_myres[0]=1;
    a1_alloc(bmode,myres,b1_myres,size,x,a1_x); 

} 

