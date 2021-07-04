#include  <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Routines.h"
#include "mpi.h"

#define dim(i,j) ( (int)(i)*(Mr+2) + (int)(j) )

int WORKER(int nWRs, int Me, int *iparms, double *parms)
{
    int i, j, Mr, Mz;
    double tend, dtout, dtfactor, D, Rin, Rout, Z, a, b;
    double dr, dz, dtexpl, dt;
	
    int Niparms = 2; int Nparms = 9;
    
    /*Unpack initial iparms and parms arrays, local Mz = Mz / nWRs*/
    MPI_Bcast( iparms, Niparms, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast( parms, Nparms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    //calculating Mz, a, and b of each process
    Mr = iparms[0]; Mz = iparms[1]/ nWRs;
	Rin=parms[2]; Rout=parms[3];
    double alfa = parms[4]/nWRs;
    a = (Me-1)*alfa; b = a + alfa;
    dt = parms[0]; D = parms[1]; dr = parms[5]; dz = parms[6];  
    tend = parms[7]; dtout = parms[8]; 
   
	double *r, *z, *U;
	double **Fr, **Fz;
    r = (double *)malloc((Mr+2)*sizeof(double));
    z = (double *)malloc((Mz+2)*sizeof(double));
    U = (double *)malloc(((Mz+2)*(Mr+2))*sizeof(double));

    Fr = (double **)malloc((Mz+2)*sizeof(double *));
	for (i=1; i<(Mz+2); i++)
		Fr[i] = (double *)malloc((Mr+2)*sizeof(double));

    Fz = (double **)malloc((Mz+2)*sizeof(double *));
	for (i=1; i<(Mz+2); i++)
		Fz[i] = (double *)malloc((Mr+2)*sizeof(double));

	MESH(a,b,Mz,z);
    MESH(Rin,Rout,Mr,r);

    /// --------- initialize -----------------
    int nsteps = 0;
    double time = 0.0;
    double tout = dtout;
    int maxsteps = tend/dt + 1;

    INIT(Me,nWRs,Mr,Mz,U,r,z);

    /// --------- begin timestepping -------------

    int NodeUP = Me+1; int NodeDN = Me-1;

    EXCHANGE_bry_MPI( nWRs, Me, NodeUP, NodeDN, Mr, Mz, U );
    
	SEND_output_MPI( Me, nWRs, NodeUP, NodeDN, Mr, Mz, U );

    for (nsteps=1; nsteps<=maxsteps; nsteps++){
        
	    //............synchronize everyone.............//
        MPI_Barrier( MPI_COMM_WORLD );

        time = time + dt;
	    //--------------- Exchange "boundary" values ----------------//
        EXCHANGE_bry_MPI( nWRs, Me, NodeUP, NodeDN, Mr, Mz, U );
	        
	    FLUXr(Mr, Mz, D, dr, Rout, z, U, Fr, time);
	    FLUXz(Me, nWRs, Mr, Mz, D, dz, U, Fz, time);
        
		PDE( Mr, Mz, dt, dr, dz, r, U, Fr, Fz);
   	
        if( time >= tout ){
            //-------------- send  output to the  MASTER ----------//
            SEND_output_MPI( Me, nWRs, NodeUP, NodeDN, Mr, Mz, U );
            tout = tout + dtout;
        }
    }
    /// --------- end of timestepping ----------------

	// free the memory
	free(r);free(z);
	free(U);
	for (j=1; j<(Mz+2); j++)
		free(Fr[j]);
	free(Fr);
	for (j=1; j<(Mz+2); j++)
		free(Fz[j]);
	free(Fz);

    return 0;
}


