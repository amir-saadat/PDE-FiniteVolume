#include  <stdio.h>
#include <stdlib.h>
#include <math.h>

void FLUX(int nWRs, int Me, int M, double D, double b, double dx, double *U, 
		double *F, double time){
	int i;
	if (Me == 1) {
		// As the problem is boundary driven, we write as:
		U[0] = 1.0;

		F[1] = -D*(U[1]-U[0])/(dx/2.0);
		for(i=2;i<=M+1;i++){
			F[i]=-D*(U[i]-U[i-1])/dx;
		}
	}

	else if (Me == nWRs) {
		// As the problem is boundary driven, we write as:
		U[M+1] = 1.0 - erf(b/(2.0*sqrt(D*(time))));

		for(i=1;i<M+1;i++){
			F[i]=-D*(U[i]-U[i-1])/dx;
		}
		F[M+1]=-D*(U[M+1]-U[M])/(dx/2.0);
	}

	else {
		for(i=1;i<=M+1;i++){
			F[i]=-D*(U[i]-U[i-1])/dx;
		}
	}
}

void PDE(int M, double dt, double dx, double *U, double *F){
    int i;
	
	for(i=1;i<M+1;i++){
		U[i]=U[i]+dt/dx*(F[i] - F[i+1]);
	}
}

void COMPARISON(double *x, double *U, int M, double D, double time){
	int i;
	FILE * result;
	double *u_exact,*u_error;
	
	// allocate u_exact and error
	u_exact = (double *)malloc((M+2)*sizeof(double));
	u_error = (double *)malloc((M+2)*sizeof(double));
	
	for (i = 0; i<=(M+1); i++){
		u_exact[i] = erfc(x[i]/(2.0*sqrt(D*time)));
		u_error[i] = u_exact[i]-U[i];
	}
	
	result = fopen("result.dat","wb");
	for (i = 0; i<(M+2); i++)
		fprintf(result,"%lf\t%lf\t%lf\t%lf\n",x[i],U[i],u_exact[i],u_error[i]);
	
	double MaxAbsErr = fabs(u_error[0]);
	for (i = 1; i<=(M+1); i++) {
		if ( fabs(u_error[i]) > MaxAbsErr ) {
			MaxAbsErr = fabs(u_error[i]);
		}
	}
	printf("Maximum Absolute Error is:%lf\n", MaxAbsErr);
	
	//deallocate arrays
	free(u_exact);free(u_error);
}

double TRAPZRULE(int M, double *U,double dx){
	int i;
	double area;
	area = (U[0]+3*(U[1]+U[M])+U[M+1])/4.0;
	for(i=2;i<M;i++){
		area = area + U[i];
	}
	area = area * dx;
	return area;
}
