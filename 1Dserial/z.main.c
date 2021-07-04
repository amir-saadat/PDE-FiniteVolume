#include  <stdio.h>
#include <stdlib.h>
#include <math.h>

// prototype of the Subroutines and Function

void MESH(double *x, int M, double a, double b, double dx);
void INIT(double *U, double *x, int M);
void OUTPUT(double *x, double *U, int M, double time);
void FLUX(double *U, double *F, double dx, double D, double b, int M, double time);
void PDE(double *U, double *F, double dx, double dt, int M);
void COMPARISON(double *x, double *U, int M, double D, double time);
double TrapzRule (int M, double dx, double *U);

/////// MAIN PROGRAM ///////

int main(int argc, char *argv[]){
	// declaration of parameters and files
	int i, MM, M;
	double tend, dtfactor, dtout, D, a, b;
	double dx, dtEXPL, dt;
	FILE * data;

	// declaration of arrays
	double *x, *U, *F;
    data = fopen("data.dat","r");
	// reading MM, t_end, dtfactor and dtout
	fscanf(data,"%i\t %lf\t %lf\t %lf\n", &MM, &tend, &dtfactor, &dtout);
	// reading D, a and b
	fscanf(data,"%lf\t %lf\t %lf", &D, &a, &b);
	fclose(data);

	// character to specify if we want to compare with an exact solution
	const char *compare="Yes";

	//calculating other parameters
    dx = 1.0/MM ; M = (b-a)*MM;
    dtEXPL = dx*dx/(2*D) ; dt = dtfactor*dtEXPL;

    // allocate x
    x = (double *)malloc((M+2)*sizeof(double));
 
	MESH(x, M, a, b, dx);

    // initialize //
	int nsteps = 0; double time = 0.0; double tout = dtout; int maxSteps = (double)(tend/dt) + 1;

	// allocate U
	U = (double *)malloc((M+2)*sizeof(double));

	INIT(U, x, M);

	OUTPUT(x, U, M, time);

	//calculating area under the curve
	double area = TrapzRule(M,dx,U);
	printf("Area = %lf, at t = %lf\n", area, time);
	printf("dt = %lf\n", dt);


	// Begin Time Stepping //
    // allocate F
    // In C index starts from 0. To keep the proposed indexing, I generate an array with M+2 elements.
    //   But I don't use the first element.

    F = (double *)malloc((M+2)*sizeof(double));

    for (nsteps=1; nsteps <= maxSteps; nsteps++){

		FLUX(U,F,dx,D,b,M,time);
		PDE(U,F,dx,dt,M);
		time = time + dt;
		if ( time >= tout ){
			OUTPUT(x,U,M,time);
			tout = tout + dtout;
	    }
	    if ( time >= 40.0){
			//calculating area under the curve
		    double area = TrapzRule(M,dx,U);
		    printf("Area = %lf, at t = %lf\n", area, time);
		}
		if ( (compare=="Yes") && (time >= tend) ){
			COMPARISON(x, U, M, D, time);}

	}

	printf("DONE, exiting at t=%lf after %i steps\n", time, nsteps-1);

	// free the memory
	free(x);free(U);free(F);
	return 0;
}

