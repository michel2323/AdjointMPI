#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "dco_tape.hpp"
#include "draw.hpp"
#include "complex.hpp"

#define TAPE_MODE 1
#define DERIV_MODE 1

using std::vector;


int main(int argc, char** argv){

    // FORWARD
    int  myid, numprocs;
    int namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    init_mode(TAPE_MODE, DERIV_MODE);
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Get_processor_name(processor_name,&namelen);
    MPI_Status status;
    //	print_tape();
    //	interpret_tape();
    //	print_tape();

    // fractal
    //
    int maxiterations = 1000;
    int iterations = 0;
    complex z;
    complex c;
    int j = 0;
    int height = 1000;
    int width = 1000;
    int i = 0;
    // image
    int *** image;
    double ** mpi_image;
    double * mpi_row;
    if(myid == 0) {
	image = new int**[height];
	for(i = 0 ; i < height ; i++) {
	    image[i] = new int*[width];
	    for(j = 0 ; j < width ; j++) {
		image[i][j] = new int[3];
	    }
	}
    }

    mpi_image = new double*[height];
    for(j = 0 ; j < height ; j++) {
	mpi_image[j] = new double[width];
    }

    mpi_row = new double[width];

    int r_va;
    int i_va;
    complex d;
    active resd;
    // master/slave
    int buffer = 0;
    if(myid == 0) {
	reset_tape();
	for(i = 0 ; i < numprocs - 1 ; i++) {
	    printf("i: %d\n", (int) ((((double) height / 100.0))* (double) i));
	    MPI_Send(&i, 1 , MPI_INT, i+1 , 0,MPI_COMM_WORLD);
	}
	for(i = numprocs-1 ; i < height ; i++) {
	    printf("%d\n", (int) (((100.0 / (double) height)) * (double) i));
	    MPI_Recv(&buffer, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,&status);
	    MPI_Recv(mpi_row, width, MPI_DOUBLE, status.MPI_SOURCE, 0, MPI_COMM_WORLD,&status);
	    for(j = 0 ; j < width ; j++) {
		mpi_image[buffer][j] = mpi_row[j];
	    }
	    MPI_Send(&i, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
	}
	for(i = 1 ; i < numprocs ; i++){
	    MPI_Recv(&buffer, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,&status);
	    MPI_Recv(mpi_row, width, MPI_DOUBLE, status.MPI_SOURCE, 0, MPI_COMM_WORLD,&status);
	    for(j = 0 ; j < width ; j++) {
		mpi_image[buffer][j] = mpi_row[j];
	    }
	    buffer = -1;
	    MPI_Send(&buffer, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
	}
	reverse_tape_interpreter();
	print_tape();
    }
    else {
	while(buffer != -1) {
	    MPI_Recv(&buffer, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,&status);
	    if(buffer != -1) {
		for(i = 0 ; i < width ; i++) {
	    	reset_tape();
//		    			c.r = 2.0*((double) i/((double) width * ((double) height/ (double) width))) - 2.0;
//		    			c.i = 2.0*((double) buffer /(double) height) - 1.0;
		    c.r = 0.8*0.01*((double) i/((double) width * ((double) height/ (double) width))) - 0.84;
		    c.i = 0.8*0.01* ((double) buffer /(double) height) - 0.21;
		    independent(c.r);
		    independent(c.i);
		    z = c;
		    iterations = 0;
		    while(iterations < maxiterations){
			z = complex_pow(z);
			z = z + c;
			if(complex_norm(z) > 2.0)
			    break;
			iterations++;
		    }
		   dependent(z.r); 
		   dependent(z.i); 
		   seed_dependent(z.r); 
		   seed_dependent(z.i); 
		   reverse_tape_interpreter();
		   d.r = tape[tape_entry::indeps[0]].d;
		   d.i = tape[tape_entry::indeps[1]].d;
//		   mpi_row[i] = (double) iterations;
		   resd = complex_norm(d);
		   mpi_row[i] = resd.v;
		}
		MPI_Send(&buffer, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		MPI_Send(mpi_row, width, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
//		print_tape();
	    }
	}
    }
    double maxd = 0;
    if(myid == 0) {
/*	for(j = 0 ; j < height ; j++) {
	    for(i = 0 ; i < width ; i++) {
		iterations = (int) mpi_image[j][i];			

		//			printf("iterations: %d\n", iterations);
		if(iterations == maxiterations){
		    image[j][i][0] = 0;
		    image[j][i][1] = 0;
		    image[j][i][2] = 0;
		}	
		else {	
		    image[j][i][0] = 255 - (int) (((double) iterations / (double) maxiterations)*255.0);	
		    image[j][i][1] = (int) ((((double) iterations / (double) maxiterations)*255.0) + 256.0);
		    //		printf("old: %d ",image[j][i][1]);
		    image[j][i][1] = image[j][i][1] % 256;
		    //		printf("new: %d\n",image[j][i][1]);
		    image[j][i][2] = (int) (((double) iterations / (double) maxiterations)*255.0);	
		}
	    }
	}	
*/	for(j = 0 ; j < height ; j++) {
	    for(i = 0 ; i < width ; i++) {
	         if(mpi_image[j][i] > 1.0){
		     mpi_image[j][i] = log(mpi_image[j][i]);
		 }
	    }
	}
	for(j = 0 ; j < height ; j++) {
	    for(i = 0 ; i < width ; i++) {
	         if(mpi_image[j][i] > 1.0){
		     mpi_image[j][i] = log(mpi_image[j][i]);
		 }
	    }
	}

	for(j = 0 ; j < height ; j++) {
	    for(i = 0 ; i < width ; i++) {
	         if(mpi_image[j][i] > maxd){
		     maxd = mpi_image[j][i];
		 }
	    }
	}
		printf(" max: %f \n", maxd);
	for(j = 0 ; j < height ; j++) {
	    for(i = 0 ; i < width ; i++) {
		    image[j][i][0] = 255 - (int) ((mpi_image[j][i] / maxd)*255.0);	
//		    printf(" %d ", image[j][i][0]);
		    image[j][i][1] = (int) (((mpi_image[j][i] / maxd)*255.0) + 256.0);
		    //		printf("old: %d ",image[j][i][1]);
		    image[j][i][1] = image[j][i][1] % 256;
		    //		printf("new: %d\n",image[j][i][1]);
		    image[j][i][2] = (int) ((mpi_image[j][i] / maxd)*255.0);	
	    }
	}		
	draw(image, height, width); 
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
} 



