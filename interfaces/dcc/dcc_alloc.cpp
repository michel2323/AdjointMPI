#include "dcc_alloc.hpp"

//dcc_mmap_id dcc_mmap;

stack<double> dcc_mstack;


//double *buf_ampi_rst[200];
//double **ptr_buf_ampi_rst[200];

void dcc_new(double *&buf, int size) {
    buf = new double[size];
}

void dcc_delete(double *&buf) {
    delete [] buf;
}

// FIRST ORDER ADJOINTS

// Buffer allocation

void a1_dcc_new(int bmode, double *&buf, double *&a1_buf, int &size) {
    if(bmode == 2) {
	buf = new double[size];
	a1_buf = new double[1];
	for(int i = size ; i > 0 ; i--) buf[i-1] = 0;
	a1_buf[0] = size;

    }
    if(bmode == 1) {
	delete [] buf;
	delete [] a1_buf;
    }
}

void a1_dcc_delete(int bmode, double *&buf, double *&a1_buf) {
    int size=0;
    if(bmode == 2) {
	size=a1_buf[0];
	for(int i = 0 ; i < size ; i++) dcc_mstack.push(buf[i]);
	dcc_mstack.push((double) size);
    	delete [] buf;
    	delete [] a1_buf;
    }
    if(bmode == 1) {
	size=(int) dcc_mstack.top();
	dcc_mstack.pop();
	buf=new double[size];
	for(int i = size-1 ; i >=0 ; i--) {
	    buf[i] = dcc_mstack.top();
	    dcc_mstack.pop();
	    printf("buf[%d]: %f\n",i,buf[i]);
	}
	a1_buf=new double[size];
	for(int i = 0 ; i < size ; i++) a1_buf[i] = 0;
    }
}

// FORWARD OVER REVERSE

void t2_a1_dcc_new(int bmode, double *&buf, double *&t2_buf, double *&a1_buf, double *&t2_a1_buf, int &size) {
    if(bmode == 2) {
	buf = new double[size];
	t2_buf = new double[size];
	a1_buf = new double[1];
	for(int i = size ; i > 0 ; i--) buf[i-1] = 0;
	for(int i = size ; i > 0 ; i--) t2_buf[i-1] = 0;
	a1_buf[0] = size;

    }
    if(bmode == 1) {
	delete [] buf;
	delete [] t2_buf;
	delete [] a1_buf;
	delete [] t2_a1_buf;
    }
}

void t2_a1_dcc_delete(int bmode, double *&buf, double *&t2_buf, double *&a1_buf, double *&t2_a1_buf) {
    int size=0;
    if(bmode == 2) {
	size=a1_buf[0];
	for(int i = 0 ; i < size ; i++) dcc_mstack.push(t2_buf[i]);
	for(int i = 0 ; i < size ; i++) dcc_mstack.push(buf[i]);
	dcc_mstack.push((double) size);
    	delete [] buf;
    	delete [] t2_buf;
    	delete [] a1_buf;
    }
    if(bmode == 1) {
	size = (int) dcc_mstack.top();
	dcc_mstack.pop();
	buf=new double[size];
	t2_buf=new double[size];
	for(int i = size-1 ; i >=0 ; i--) {
	    buf[i] = dcc_mstack.top();
	    dcc_mstack.pop();
	}
	for(int i = size-1 ; i >=0 ; i--) {
	    t2_buf[i] = dcc_mstack.top();
	    dcc_mstack.pop();
	}
	a1_buf=new double[size];
	t2_a1_buf=new double[size];
	for(int i = 0 ; i < size ; i++) a1_buf[i] = 0;
	for(int i = 0 ; i < size ; i++) t2_a1_buf[i] = 0;
    }
}

