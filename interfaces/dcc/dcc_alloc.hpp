#ifndef DCC_ALLOC_H
#define DCC_ALLOC_H
#include <math.h>
#include <cassert>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <stack>
//#include <unordered_map>

using namespace std;

// Memory allocator

// hash table entry

//class mentry{
    //public:
    ////stac< data;
    //int size;
    //vector<double> data;
    ////vector<double> a1_data;
    //double *hash;
    //double *depth;
//};

//typedef std::unordered_map<double*, mentry> dcc_mmap_id;
// Buffer allocation

void dcc_new(double *&buf, int size);
void dcc_delete(double *&buf);

void print_num(double &num);

// First order

// Buffer allocation

void a1_dcc_new(int bmode, double *&buf, double *&a1_buf, int &size);
void a1_dcc_delete(int bmode, double *&buf, double *&a1_buf);

void a1_print_num(int bmode, double &num, double &a1_num);

// Forward over Reverse

// Buffer allocation

void t2_a1_dcc_new(int bmode, double *&buf, double *&t2_buf, double *&a1_buf, double *&t2_a1_buf, int &size);
void t2_a1_dcc_delete(int bmode, double *&buf, double *&t2_buf, double *&a1_buf, double *&t2_a1_buf);

void t2_a1_print_num(int bmode, double &num, double &t2_num, double &a1_num, double &t2_a1_num);
#endif 
