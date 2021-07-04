
void MESH(double a, double b, int M, double *x);
void INIT(int Me, int nWRs, int Mr, int Mz, double *U, double *r, double *z);
void OUTPUT(int Mr, int Mz, double *U, double *r, double *z,double time);
void FLUXr(int Mr, int Mz, double D, double dr, double Rout, double *z, double *U, double **Fr, double time);
void FLUXz(int Me, int nWRs, int Mr, int Mz, double D, double dz, double *U, double **Fz, double time);
void PDE(int Mr, int Mz, double dt, double dr, double dz, double *r, double *U, double **Fr, double **Fz);
void COMPARISON(int Mr, int Mz, double *r, double *z, double *U, double D, double time);
int MASTER(int nWRS, int mster, int *iparms, double *parms);
int WORKER(int nWRS, int Me, int *iparms, double *parms);
void EXCHANGE_bry_MPI(int nWRS, int Me, int NodeUP, int NodeDN, int Mr, int Mz, double *U);
void SEND_output_MPI(int Me, int nWRs, int NodeUP, int NodeDN, int Mr, int Mz, double *U);
void RECV_output_MPI(int nWRs, int Mz, int Mr, double *U);



