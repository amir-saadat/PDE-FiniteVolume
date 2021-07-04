void MESH(double *x, int M, double a, double b, double dx){
	int i;
	x[0] = a; x[1] = a + dx/2.0;
	for (i = 2; i<(M+1); i++)
		x[i] = x[1] + (i-1)*dx;
	x[M+1] = b;
}

void INIT(int Me, double *U, double *x, int M){
	int i;
	if (Me == 1) {U[0]=1.0;}
	else {U[0]=0.0;}

	for (i = 1; i<=(M+1); i++){
		U[i] = 0.0;}
}

