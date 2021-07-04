void MESH(double *x, int M, double a, double b, double dx);
void INIT(int Me, double *U, double *x, int M);
void OUTPUT(double *x, double *U, int M, double time);
void FLUX(int nWRs, int Me, int M, double D, double b, double dx, double *U, 
		double *F, double time);
void PDE(int M, double dt, double dx, double *U, double *F);
double TRAPZRULE(int M, double *U, double dx);
int MASTER(int nWRS, int mster, int *iparms, double *parms);
int WORKER(int nWRS, int Me, int *iparms, double *parms);
void EXCHANGE_bry_MPI(int nWRS, int Me, int NodeUP, int NodeDN, int M, double *U);
void SEND_output_MPI(int Me, int nWRs, int NodeUP, int NodeDN, int M, double *U);
void RECV_output_MPI(int nWRs, int M, double *U);
