#include <omp.h>
#include <iostream>




int main (int argc, char *argv[]) {
    int n=1000000;
    int myid, num_threads;
    double x[n];
    double res=0;
    double temp=0;
    for(int i=0 ; i<n ; i++) {
	x[i] = i+2;
    }
#pragma omp parallel private(myid)
    {
	myid = omp_get_thread_num();
	std::cout << "Hello World from thread" << myid << std::endl;
#pragma omp barrier
	if ( myid == 0 ) {
	    num_threads = omp_get_num_threads();
	    std::cout << "There are " << num_threads << " threads" << std::endl;
	}
    }
#pragma omp parallel 
    {
#pragma omp for private(temp)
	for(int i=0 ; i<n ; i++) {
	    temp = temp+x[i]*x[i];
	}
#pragma omp atomic
	    res+=temp;
    }
    std::cout << "Res: " << res << std::endl;
    return 0;
}
