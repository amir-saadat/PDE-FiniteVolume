#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h" 
#include "Routines.h"

int main(int argc, char **argv)
{
    /*>>>mpi>>>*/
    //...startup:
    int nPROC, myID, ierr;
    // arrays to store integer and double parameters for master and all workers
    int Niparms = 2; int Nparms = 9;
	int *iparms = malloc((Niparms)*sizeof(int));
	double *parms = malloc((Nparms)*sizeof(double));
	
	ierr = MPI_Init(&argc,&argv);
	//....... nPROC is specified at mpirun or mpiexec, see Makefile....
    ierr = MPI_Comm_size(MPI_COMM_WORLD,&nPROC);
    int mster = 0;		// master has rank=0
    int nWRs  = nPROC - 1;	// =number of workers
    /*----------------- start 0, ... ,nWRs tasks ---------------*/
    ierr = MPI_Comm_rank(MPI_COMM_WORLD,&myID);    //..assigns myID
    //cout<<">>>main>>> running on " << nWRs <<"WRs"<<endl;
    if (myID == mster){
        double tt0 = MPI_Wtime();       //...start CPU timer on MR
        MASTER( nWRs, mster, iparms, parms);
        double tt1 = MPI_Wtime();       //...end timer
		printf(">>main>> MR timing= %lf sec on %i workers\n",tt1-tt0,nWRs);
    }
    else {
        WORKER( nWRs, myID, iparms, parms  );  /*... now MPI is running ...*/
		printf("Bye from WR: %i with ierr: %i\n",myID,ierr);
		if( ierr != 0 ){
			printf(">>>> worker: %i ended with ierr= %i\n",myID,ierr);
        }
    }
    ierr = MPI_Finalize();
 
    free(iparms); free(parms);
	return 0;
}



