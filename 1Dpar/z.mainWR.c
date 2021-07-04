#include  <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Routines.h"
#include "mpi.h"

int WORKER(int nWRs, int Me, int *iparms, double *parms)
{
	int MM, M;
	double tend, dtout, dtfactor, D, a, b;
	double dx, dtexpl, dt;
	
    int Niparms = 1; int Nparms = 7;

	/*Unpack initial iparms and parms arrays, local M = M / nWRs*/
    MPI_Bcast( iparms, Niparms, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast( parms, Nparms, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	//calculating M, a, and b of each process
	M = iparms[0]/ nWRs;
	double alfa = (parms[3]-parms[2])/nWRs;
	a = parms[2]+(Me-1)*alfa; b = a + alfa;
	dt = parms[0]; D = parms[1]; dx = parms[4]; 
	tend = parms[5]; dtout = parms[6]; 

	double *x = malloc((M+2)*sizeof(double));
	double *U = malloc((M+2)*sizeof(double));
	double *F = malloc((M+2)*sizeof(double));

	//Calling MESH and INIT by each worker for its own chunk
    MESH(x, M, a, b, dx);

	// .......... initialize ........... //
    int nsteps = 0; double time = 0.0; double tout = dtout; int maxsteps = (double)(tend/dt) + 1;

	INIT( Me, U, x, M);
	
	// ........ begin timestepping ......... //

	int NodeUP = Me+1; int NodeDN = Me-1;

	SEND_output_MPI( Me, nWRs, NodeUP, NodeDN, M, U );

	for (nsteps=1; nsteps<=maxsteps; nsteps++){
	    
		//............synchronize everyone.............//
	    MPI_Barrier( MPI_COMM_WORLD );

		//........... Exchange "boundary" values ...........//
        EXCHANGE_bry_MPI( nWRs, Me, NodeUP, NodeDN, M, U );
		
		FLUX(nWRs, Me, M, D, parms[3], dx, U, F, time);
        PDE( M, dt, dx, U, F);
   	    
		time = time + dt;
   	    
		if( time >= tout ){
			//............ send  output to the  MASTER  ...........//
			SEND_output_MPI( Me, nWRs, NodeUP, NodeDN, M, U );
   		    tout = tout + dtout;
	    }
	}
	// ........... end of timestepping ..............//
	
	free(x); free(U); free(F);
	return 0;
}
