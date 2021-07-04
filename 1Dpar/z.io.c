#include  <stdio.h>
#include <stdlib.h>
#include <math.h>

void OUTPUT(double *x, double *U, int M, double time){
	int i;
	FILE * output;
	char filename[32];
	sprintf( filename, "o.prof%i",(int)time );
	output = fopen(filename,"wb");
	for (i = 0; i<(M+2); i++)
		fprintf(output,"%lf\t%lf\n",x[i],U[i]);
}

