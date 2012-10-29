#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include "dco_tape.hpp"

#define TAPE_MODE 1
#define DERIV_MODE 1

using namespace std;

extern active CalcPi (int n, int iRank, int iNumProc);

int main(int argc, char **argv)
{
    int n = 1500;
    int iMyRank, iNumProcs;
    const double fPi25DT = 3.141592653589793238462643;

    init_mode(TAPE_MODE, DERIV_MODE);
    active fLocalPi;
    /* MPI Initialization */
    AMPI_Init(&argc, &argv);
    AMPI_Comm_size(MPI_COMM_WORLD, &iNumProcs);
    AMPI_Comm_rank(MPI_COMM_WORLD, &iMyRank);
    MPI_Status status;
#ifdef READ_INPUT  
    if (iMyRank == 0) 
    {
        cout << "Enter the number of intervals: ";
        cin >> n;
    }
#endif

    double fTimeStart = MPI_Wtime();
    
    /* broadcast n */
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (n <= 0 || n > 2147483647)
    {
       cout << "\n given value has to be between 0 and 2147483647 ..\n";
       return 1;
    }

    /* the calculation is done here*/
    active fPi = 0;
    fPi = CalcPi(n, iMyRank, iNumProcs);
    independent(fPi);
    fLocalPi = fPi;
    AMPI_Reduce(&fLocalPi, &fPi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    double fTimeEnd = MPI_Wtime();
if(iMyRank == 0)
    dependent(fPi);
    seed_dependent(fPi);
    if (iMyRank == 0)
    {
        cout << "\npi is approximately = " << setprecision(20) << fPi.v<< "\n"
             << "Error is            = " << fabs(fPi.v - fPi25DT) << "\n"
             << "wall clock time     = " << fTimeEnd - fTimeStart << "\n";
    }
    reverse_tape_interpreter();

    // PRINT TAPE AFTER REVERSE MODE

    for(int i = 0 ; i < iNumProcs ; i++) {
	if(iMyRank == i)
	    print_tape_indeps();
	MPI_Barrier(MPI_COMM_WORLD);
    }
    AMPI_Finalize();
    return 0;
}

active f(active a)
{
    return (4.0 / (1.0 + a*a));
}

active CalcPi (int n, int iRank, int iNumProcs)
{
    active fH   = 1.0 / (double) n;
    active fSum = 0.0;

    for (int i = iRank; i < n; i += iNumProcs)
    {
        active fX = fH * ((double)i + 0.5);
        fSum = fSum + f(fX);
    }  
    return (fH * fSum);
}
