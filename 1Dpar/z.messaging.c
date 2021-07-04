#include <math.h>
#include "mpi.h"

void EXCHANGE_bry_MPI( int nWRs, int Me, int NodeUP, int NodeDN, int M, double *U){

	/*..........Exchange "boundary" values btn neighbors.........*/
	/*.................... every WR does this ...................*/
	int Jup = M;
	int Jup1 = Jup + 1;
	int msgUP = 10;
	int msgDN = 20;
	int msgtag;
	MPI_Status status;

	//.................send bottom row to neighbor down:
	if ( Me != 1 ) {
		msgtag = msgDN;
		MPI_Send(&U[1],1,MPI_DOUBLE,NodeDN,msgtag,MPI_COMM_WORLD);
	}

	//.....receive bottom row from neighbor up and save as upper bry:
	if ( Me != nWRs ) {
		msgtag = msgDN;
		MPI_Recv(&U[Jup1],1,MPI_DOUBLE,NodeUP,msgtag,MPI_COMM_WORLD,&status);
	}

	//...................send the top row to neighbor up:
	if ( Me != nWRs ) {
		msgtag = msgUP;
		MPI_Send(&U[Jup],1,MPI_DOUBLE,NodeUP,msgtag,MPI_COMM_WORLD);
    }

	//......receive top row from neighbor down and save as lower bry:
	if ( Me != 1 ) {
		msgtag = msgUP;
		MPI_Recv( &U[0],1,MPI_DOUBLE,NodeDN,msgtag,MPI_COMM_WORLD,&status);
	}

}

void RECV_output_MPI( int nWRs, int M, double *U ){
	
	/*............ only MR does this ..........*/
		
	//.........receive values from everybody for output...........

	int i;
	int J ;	//size of U(:) array
	int Jme; //Offset from the first element of the global U array
	int msgtag;
	MPI_Status status;

	for (i = 1; i <= nWRs; i++){

		//Jme and J are calculated according to Worker number:

		if (i == 1 ) { 
			Jme = 0; 
			J = M+1;
		}
	    else if ( i == nWRs ) { 
			Jme = (i-1)*M + 1; 
			J = M+1;
		}
	    else { 
			Jme = (i-1)*M + 1; 
			J = M;
		}

		msgtag = 1000 + i;
		MPI_Recv(&U[Jme],J,MPI_DOUBLE,MPI_ANY_SOURCE,msgtag,MPI_COMM_WORLD,&status);
	}
}

void SEND_output_MPI( int Me, int nWRs, int NodeUP, int NodeDN, int M, double *U ){
	/*.................. every WR does this .................*/
	
	//.....everybody  sends values to the Master for output.....
	int mster = 0;
	int J;
	int msgtag = 1000 + Me;
    
	//First and last elements of U array of each one of the workers is as follows:
	int first, last;
	if (Me == 1 ) { 
		first = 0; 
		last = M;//Not needed 
		J = M+1;
	}
	else if ( Me == nWRs ) { 
		first = 1; 
		last = M+1;//Not needed 
		J = M+1;
	}
	else { 
		first = 1; 
		last = M;//Not needed 
		J = M;
	}

	MPI_Send(&U[first],J,MPI_DOUBLE, mster,msgtag,MPI_COMM_WORLD);
}
