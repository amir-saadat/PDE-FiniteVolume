#include  <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "Routines.h"

#define dim(i,j) ( (int)(i)*(Mr+2) + (int)(j) ) 

int MASTER(int nWRs, int mster, int *iparms, double *parms)
{
    int i, j, Mr, Mz, MMr, MMz, ierr, err;
    double tend, dtout, dtfactor, D, Rin, Rout, Z;
    double dr, dz, dtexpl, dt;
    int Niparms = 2; int Nparms = 9;
    const double PI = asin(1.0)*2.0;

	/*Read run-time parameters from data file, print them in o.out*/
	FILE *data;
	
	data = fopen("data.dat","r");
	// reading MM, t_end, dtfactor and dtout
	fscanf(data,"%i\t %i\t %lf\t %lf\t %lf\n", &Mr, &Mz, &tend, &dtfactor, &dtout);
	// reading D, a and b
	fscanf(data,"%lf\t %lf\t %lf", &D, &Rin, &Rout);
	fclose(data);
	

	Z=PI;
	//dr = 1.0/MMr; Mr = (int) ((Rout-Rin)*MMr); 
	//dz = 1.0/MMz; Mz = (int) (Z*MMz);
	dr = (Rout-Rin)/Mr;dz=Z/Mz;
    //cout<<"dr:"<<dr<<"; dz:"<<dz<<endl;

	printf("Mr=%i\tMz=%i\ttend=%lf\tdtfactor=%lf\tdtout=%lf\tD=%lf\tRin=%lf\tRout=%lf\n",
			Mr,Mz,tend,dtfactor,dtout,D,Rin,Rout);

	if ((Mz%nWRs) != 0){
		printf("please select Mz to be divisible by nWRs\n");
		printf("Mz:%i\tMr:%i\tnWRs:%i\n",Mz,Mr,nWRs);
		MPI_Abort(MPI_COMM_WORLD,err);
	}

    dtexpl = 1.0/(2.0*D*(1.0/(dr*dr)+1.0/(dz*dz))) ; dt = dtfactor*dtexpl;


    /*Pack integers in iparms array, reals in parms array, and send to all*/
    iparms[0]=Mr;iparms[1]=Mz;
    parms[0]=dt;parms[1]=D;parms[2]=Rin;parms[3]=Rout;parms[4]=Z;
    parms[5]=dr;parms[6]=dz;parms[7]=tend;parms[8]=dtout;
    
    MPI_Bcast( iparms, Niparms, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast( parms, Nparms, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	double *r = malloc((Mr+2)*sizeof(double));
	double *z = malloc((Mz+2)*sizeof(double));
	//   Here in lab9, U indecis were selected as U[j][i] j=0~Mz+1, i=0~Mr+1. 
	//   In order to transferr the U variable in contiguous fashion, it is necessary to define it this
	//      way. 
	//   Also, in order for U to be contiguous,  Macro (dim(j,i)) to facilitate the understanding of 
	//      dimensions of it as a matrix. So, U[DN(j,i)]=U[j][i].

    double *U = malloc(((Mz+2)*(Mr+2))*sizeof(double));
        
    // --------- initialize -----------------
    int nsteps = 0;	
    double time = 0.0;
    double tout = dtout;
    int maxsteps = tend/dt + 1;

    MESH(Rin,Rout,Mr,r);
    MESH(0,Z,Mz,z);
     
    RECV_output_MPI( nWRs, Mz/nWRs , Mr , U );

    OUTPUT(Mr,Mz,U,r,z,time);

    /// --------- begin timestepping -------------
    
    for (nsteps=1; nsteps<=maxsteps; nsteps++){
	
        //--------synchronize everyone--------//
        MPI_Barrier( MPI_COMM_WORLD );
        
        time = time + dt;
        
        if( time >= tout ){
            RECV_output_MPI( nWRs, Mz/nWRs , Mr , U );
            OUTPUT(Mr,Mz,U,r,z,time);
            tout = tout + dtout;
        }
    }    
    /// --------- end of timestepping ----------------
    COMPARISON(Mr,Mz,r,z,U,D,time);
   
    printf("DONE, exiting at t=%lf after %i steps\n", time, nsteps-1);

    free(r); free(z);
    free(U);
    return 0;
}


