#include  <stdio.h>
#include <stdlib.h>
#include <math.h>

#define dim(i,j) ( (int)(i)*(Mr+2) + (int)(j) )

void MESH(double a, double b, int M, double *x){
    int i;
	double dx = (b-a)/M;
    x[0]=a;
    x[1]=a + dx/2.0;
    for (i=2; i<(M+1);i++ )
    {
            x[i] = x[i-1]+dx;
    }
    x[M+1] = b;
}

void INIT(int Me, int nWRs, int Mr, int Mz, double *U, double *r, double *z){
    int i, j, start, end;
	if (nWRs == 1) {
		start=0;end=Mz+2;
	}
	else {
       if (Me == 1){start=0;end=Mz+1;}
       else if (Me == nWRs) {start=1;end=Mz+2;}
       else {start=1;end=Mz+1;}
	}
    for(j=start; j<end; j++){
       for(i=0; i<(Mr+2); i++){
          U[dim(j,i)] = sin(z[j])*log(r[i]);
       }
    }
}

