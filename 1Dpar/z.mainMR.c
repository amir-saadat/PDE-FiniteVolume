#include  <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Routines.h"
#include "mpi.h"

int MASTER(int nWRs, int mster, int *iparms, double *parms)
{
	int MM, M, ierr;
	double tend, dtout, dtfactor, D, a, b;
	double dx, dtexpl, dt;
	int Niparms = 1; int Nparms = 7;

	/*Read run-time parameters from data file, print them in o.out*/
	FILE *data;

	data = fopen("data.dat","r");
	// reading MM, t_end, dtfactor and dtout
	fscanf(data,"%i\t %lf\t %lf\t %lf\n", &MM, &tend, &dtfactor, &dtout);
	// reading D, a and b
	fscanf(data,"%lf\t %lf\t %lf", &D, &a, &b);
	fclose(data);

	printf("MM=%i\ttend=%lf\tdtfactor=%lf\tdtout=%lf\tD=%lf\ta=%lf\tb=%lf\n",
			MM,tend,dtfactor,dtout,D,a,b);
	
	// character to specify if we want to compare with an exact solution
	const char *compare="Yes";

	dx = 1.0/MM;
	M = (b-a)*MM;
	dtexpl = dx*dx/(2*D);
	dt = dtfactor*dtexpl;

	/*Pack integers in iparms array, reals in parms array, and send to all*/
	iparms[0]=M;
	parms[0]=dt;parms[1]=D;parms[2]=a;parms[3]=b;
	parms[4]=dx;parms[5]=tend;parms[6]=dtout;

	MPI_Bcast( iparms, Niparms, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast( parms, Nparms, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	double *x = malloc((M+2)*sizeof(double));
	double *U = malloc((M+2)*sizeof(double));

	// ............... initialize .............. //
	int nsteps = 0; double time = 0.0; double tout = dtout; int maxsteps = (double)(tend/dt) + 1;
	
	MESH(x, M, a, b, dx);

	RECV_output_MPI( nWRs, M/nWRs , U );

	OUTPUT(x, U, M, time);

	// .......... begin timestepping ........... //

	for (nsteps=1; nsteps<=maxsteps; nsteps++){
	
		//............synchronize everyone.............//
		MPI_Barrier( MPI_COMM_WORLD );

   	    time = time + dt;
   	    
		if( time >= tout ){
			
			RECV_output_MPI( nWRs, M/nWRs , U );
            OUTPUT(x, U, M, time);
			tout = tout + dtout;
   	    }
		if ( (compare=="Yes") && (time >= tend) ){
			COMPARISON(x, U, M, D, time);}
	}

	// ....... end of timestepping ....... //
	
	printf("DONE, exiting at t=%lf after %i steps\n", time, nsteps-1);
		
	free(x); free(U);
	return 0;
}
