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

void a1_dcc_new(int bmode, double *&buf, double *&b1_buf, int &size) {
    if(bmode == 2) {
	buf = new double[size];
	b1_buf = new double[1];
	for(int i = size ; i > 0 ; i--) buf[i-1] = 0;
	b1_buf[0] = size;

    }
    if(bmode == 1) {
	delete [] buf;
	delete [] b1_buf;
    }
}

void a1_dcc_delete(int bmode, double *&buf, double *&b1_buf) {
    int size=0;
    if(bmode == 2) {
	size=b1_buf[0];
	for(int i = 0 ; i < size ; i++) dcc_mstack.push(buf[i]);
	dcc_mstack.push((double) size);
    	delete [] buf;
    	delete [] b1_buf;
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
	b1_buf=new double[size];
	for(int i = 0 ; i < size ; i++) b1_buf[i] = 0;
    }
}

// FORWARD OVER REVERSE

void t2_a1_dcc_new(int bmode, double *&buf, double *&d2_buf, double *&b1_buf, double *&d2_b1_buf, int &size) {
    if(bmode == 2) {
	buf = new double[size];
	d2_buf = new double[size];
	b1_buf = new double[1];
	for(int i = size ; i > 0 ; i--) buf[i-1] = 0;
	for(int i = size ; i > 0 ; i--) d2_buf[i-1] = 0;
	b1_buf[0] = size;

    }
    if(bmode == 1) {
	delete [] buf;
	delete [] d2_buf;
	delete [] b1_buf;
	delete [] d2_b1_buf;
    }
}

void t2_a1_dcc_delete(int bmode, double *&buf, double *&d2_buf, double *&b1_buf, double *&d2_b1_buf) {
    int size=0;
    if(bmode == 2) {
	size=b1_buf[0];
	for(int i = 0 ; i < size ; i++) dcc_mstack.push(d2_buf[i]);
	for(int i = 0 ; i < size ; i++) dcc_mstack.push(buf[i]);
	dcc_mstack.push((double) size);
    	delete [] buf;
    	delete [] d2_buf;
    	delete [] b1_buf;
    }
    if(bmode == 1) {
	size = (int) dcc_mstack.top();
	dcc_mstack.pop();
	buf=new double[size];
	d2_buf=new double[size];
	for(int i = size-1 ; i >=0 ; i--) {
	    buf[i] = dcc_mstack.top();
	    dcc_mstack.pop();
	}
	for(int i = size-1 ; i >=0 ; i--) {
	    d2_buf[i] = dcc_mstack.top();
	    dcc_mstack.pop();
	}
	b1_buf=new double[size];
	d2_b1_buf=new double[size];
	for(int i = 0 ; i < size ; i++) b1_buf[i] = 0;
	for(int i = 0 ; i < size ; i++) d2_b1_buf[i] = 0;
    }
}

