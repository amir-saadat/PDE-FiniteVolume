#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define dim(i,j) ( (int)(i)*(Mr+2) + (int)(j) )

void OUTPUT(int Mr, int Mz, double *U, double *r, double *z, double time){
	int i,j;
	FILE * output;
	char filename[32];
	sprintf( filename, "o.prof%i",(int)time);
	output = fopen(filename,"wb");
	for (i = 0; i<(Mr+2); i++){
		for (j = 0; j<(Mz+2); j++){
			fprintf(output,"%lf\t%lf\t%2.14f\n",z[j],r[i],U[dim(j,i)]);
		}
	}
	fclose(output);
}



