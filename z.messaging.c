#include <math.h>
#include "mpi.h"

#define dim(i,j) ( (int)(i)*(Mr+2) + (int)(j) )

void EXCHANGE_bry_MPI( int nWRs, int Me, int NodeUP, int NodeDN, int Mr, int Mz, double *U){

    /*..........Exchange "boundary" values btn neighbors.........*/
    /*.................... every WR does this ...................*/
    int Jup = Mz;
    int Jup1 = Jup + 1;
    int msgUP = 10;
    int msgDN = 20;
    int msgtag;
    MPI_Status status;
	MPI_Datatype Mytype;

	int I2 = Mr + 2; //size of U(j,:) array

	// U has to be one contiguous, So I defined U as U[dim(j,i)] = U[i][j] in previous case (lab 5).

	MPI_Type_contiguous(I2,MPI_DOUBLE,&Mytype);
	MPI_Type_commit(&Mytype);
    
	//.................send bottom row to neighbor down:
    if ( Me != 1 ) {
        msgtag = msgDN;
        MPI_Send(&U[dim(1,0)],1,Mytype,NodeDN,msgtag,MPI_COMM_WORLD);
    }

    //.....receive bottom row from neighbor up and save as upper bry:
    if ( Me != nWRs ) {
        msgtag = msgDN;
        MPI_Recv(&U[dim(Jup1,0)],I2,MPI_DOUBLE,NodeUP,msgtag,MPI_COMM_WORLD,&status);
    }

    //...................send the top row to neighbor up:
    if ( Me != nWRs ) {
        msgtag = msgUP;
        MPI_Send(&U[dim(Jup,0)],1,Mytype,NodeUP,msgtag,MPI_COMM_WORLD);
    }

    //......receive top row from neighbor down and save as lower bry:
    if ( Me != 1 ) {
        msgtag = msgUP;
        MPI_Recv( &U[dim(0,0)],I2,MPI_DOUBLE,NodeDN,msgtag,MPI_COMM_WORLD,&status);
    }

    MPI_Type_free(&Mytype);
}

void RECV_output_MPI( int nWRs, int Mz, int Mr, double *U ){
	
    /*............ only MR does this ..........*/

    //.........receive values from everybody for output...........

    int i,j,J ;	//Size of each process's chunk
    int Jme; //Offset from the starting point in the global U array
    int msgtag;
    MPI_Status status;
    MPI_Datatype Mytype;
    int number_amount;
    
	int I2 = Mr + 2; //size of U(j,:) array
	MPI_Type_contiguous(I2,MPI_DOUBLE,&Mytype);
	MPI_Type_commit(&Mytype);

    for (i = 1; i <= nWRs; i++){
        //Jme and J are selected based on the Worker number as:

		if (nWRs == 1){
			Jme=0; J=Mz+2;
		}
		else{
            if (i == 1 ) { Jme = 0; J = Mz+1;}
            else if ( i == nWRs ) { Jme = (i-1)*Mz + 1; J = Mz+1;}
            else { Jme = (i-1)*Mz + 1; J = Mz;}
		}

        msgtag = 1000 + i;
        
		MPI_Recv(&U[dim(Jme,0)],J*I2,MPI_DOUBLE,MPI_ANY_SOURCE,msgtag,MPI_COMM_WORLD,&status);
	}

    MPI_Type_free(&Mytype);
}

void SEND_output_MPI( int Me, int nWRs, int NodeUP, int NodeDN, int Mr, int Mz, double *U ){
    /*.................. every WR does this .................*/
    int I2 = Mr + 2;
    MPI_Datatype Mytype;
    MPI_Status status;
    MPI_Type_contiguous(I2,MPI_DOUBLE,&Mytype);
    MPI_Type_commit (&Mytype);
    //.....everybody  sends values to the Master for output.....
    int mster = 0;
    
    int J;
    int msgtag = 1000 + Me;
    
	//Start and End points of U array of each one of the workers is based on their position as:
    int start, end;
	if (nWRs == 1) {
		start=0; end=Mz+1; J=Mz+2;
	}
	else {
        if (Me == 1 ) { start = 0; end = Mz; J = Mz+1;}
        else if ( Me == nWRs ) { start = 1; end = Mz+1; J = Mz+1;}
        else { start = 1; end = Mz; J = Mz;}
	}

	int i,j;

    MPI_Send(&U[dim(start,0)], J, Mytype, mster, msgtag, MPI_COMM_WORLD);
    
    MPI_Type_free(&Mytype);
}


