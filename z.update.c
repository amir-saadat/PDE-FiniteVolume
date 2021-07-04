#include  <stdio.h>
#include <stdlib.h>
#include <math.h>

#define dim(i,j) ( (int)(i)*(Mr+2) + (int)(j) )

void FLUXr(int Mr, int Mz, double D, double dr, double Rout, double *z, double *U, double **Fr, double time){
	int i, j;
		
	for(j=0;j<(Mz+2);j++){
		U[dim(j,0)]=0.0;
	}
	
	for(j=0;j<(Mz+2);j++){
		U[dim(j,Mr+1)]=exp(-time)*log(Rout)*sin(z[j]);
	}
	
	for(j=1;j<(Mz+1);j++){
		Fr[j][1]=-D*(U[dim(j,1)]-U[dim(j,0)])/(dr/2);
		for(i=2;i<(Mr+1);i++){
			Fr[j][i]=-D*(U[dim(j,i)]-U[dim(j,i-1)])/dr;
		}
		Fr[j][Mr+1]=-D*(U[dim(j,Mr+1)]-U[dim(j,Mr)])/(dr/2);
	}
}

void FLUXz(int Me, int nWRs, int Mr, int Mz, double D, double dz, double *U, double **Fz, double time){
    int i, j;

	if (nWRs == 1) {
		// As the problem is boundary driven, we write as:
		for(i=1;i<(Mr+1);i++){ 
			U[dim(0,i)]=0.0;
		}
    	for(i=1;i<(Mr+1);i++){
    	    U[dim(Mz+1,i)]=0.0;
    	}
    	    
		for(i=1;i<(Mr+1);i++){
			Fz[1][i]=-D*(U[dim(1,i)]-U[dim(0,i)])/(dz/2);
			for(j=2;j<(Mz+1);j++){
				Fz[j][i]= -D*(U[dim(j,i)]-U[dim(j-1,i)])/dz;
			}
    		Fz[Mz+1][i]=-D*(U[dim(Mz+1,i)]-U[dim(Mz,i)])/(dz/2);
		}
	}

	else {
		
		if (Me == 1) {
			
			// As the problem is boundary driven, we write as:
    		for(i=1;i<(Mr+1);i++){ 
    			U[dim(0,i)]=0.0;
    		}
    
    	    for(i=1;i<(Mr+1);i++){
    		   Fz[1][i]=-D*(U[dim(1,i)]-U[dim(0,i)])/(dz/2);
    		   for(j=2;j<(Mz+2);j++){
    			   Fz[j][i]= -D*(U[dim(j,i)]-U[dim(j-1,i)])/dz;
    		   }
    		}
    	}
    		
    	else if (Me == nWRs) {
            // As the problem is boundary driven, we write as:
    		for(i=1;i<(Mr+1);i++){
    			U[dim(Mz+1,i)]=0.0;
    		}
    	
    		for(i=1;i<(Mr+1);i++){
    			for(j=1;j<(Mz+1);j++){
    				Fz[j][i]= -D*(U[dim(j,i)]-U[dim(j-1,i)])/dz;
    			}
    			Fz[Mz+1][i]=-D*(U[dim(Mz+1,i)]-U[dim(Mz,i)])/(dz/2);
    		}
    	}
    	
    	else { //Not Me=1 or Me=nWRs
    		for(i=1;i<(Mr+1);i++){
    			for(j=1;j<(Mz+2);j++){
    				Fz[j][i]= -D*(U[dim(j,i)]-U[dim(j-1,i)])/dz;
    			}
    		}
    	}
	}

}

void PDE(int Mr, int Mz, double dt, double dr, double dz, double *r, double *U, double **Fr, double **Fz){
	int i,j;
	
	for(j=1;j<(Mz+1);j++){
		for(i=1;i<(Mr+1);i++){
			U[dim(j,i)] = U[dim(j,i)]+dt/(r[i]*dr*dz)*((r[i]-dr/2)*dz*Fr[j][i]-(r[i]+dr/2)*dz*Fr[j][i+1]+
					r[i]*dr*Fz[j][i]-r[i]*dr*Fz[j+1][i]);
			
		}
	}
}

void COMPARISON(int Mr, int Mz, double *r, double *z, double *U, double D, double time){

	int i, j;

	// allocate u_exact and error
	double **u_exact = malloc((Mz+2)*sizeof(double *));
	for (j=0; j<(Mz+2); j++)
		u_exact[j] = (double *)malloc((Mr+2)*sizeof(double));
	
	double **u_error = malloc((Mz+2)*sizeof(double *));
	for (j=0; j<(Mz+2); j++)
		u_error[j] = (double *)malloc((Mr+2)*sizeof(double));
	
    for (i = 0; i<(Mr+2); i++){
		for (j = 0; j<(Mz+2); j++){
			u_exact[j][i] = exp(-time)*log(r[i])*sin(z[j]);
			u_error[j][i] = u_exact[j][i]-U[dim(j,i)];
		}
	}

    FILE * result;
	result = fopen("result.dat","wb");	    
	for (i = 0; i<(Mr+2); i++){
		for (j = 0; j<(Mz+2); j++){
			fprintf(result,"%lf\t%lf\t%lf\t%lf\t%lf\n",r[i],z[j],U[dim(j,i)],u_exact[j][i],u_error[j][i]);}}
	
	double MaxAbsErr = fabs(u_error[0][0]);
	for (i = 0; i<(Mr+2); i++) {
		for (j = 0; j<(Mz+2); j++) {
			if ( fabs(u_error[j][i]) > MaxAbsErr ) {
				MaxAbsErr = fabs(u_error[j][i]);
			}
		}
	}
	printf("Maximum Absolute Error is:%2.10f\n", MaxAbsErr);
	
	//deallocate arrays
	for (j=0; j<(Mr+2); j++)
		free(u_exact[j]);
	free(u_exact);
	for (j=0; j<(Mr+2); j++)
		free(u_error[j]);
	free(u_error);
}

