#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

double vec_inner_prod1(double* XX, double* theta, int t, int j, int M)
{
	double ans = 0.0;
	for(int k=0; k<M; k++){
		ans += XX[j*M+k]*theta[t*M+k];
	}
	ans -= XX[j*M+j]*theta[t*M+j];
	return(ans);
}

double beta_j_norm1(double* vec, int P)
{
	double ans = 0;
	for(int t=0; t<P; t++){
		ans += vec[t]*vec[t];
	}
	return(sqrt(ans));
}

void inner_iter2(SEXP R_XX, double* XY, double* theta, int M, int P, double* beta_j_lasso, double* lambda1, double lambda2, double* Xnorm, double* bannot)
{
	double* tmp_XX;
	double tmp_XYj = 0.0;
	double bnorm = 0.0;
	
	for(int j=0; j<M; j++){
        for (int t = 0; t < P; t++) {
			tmp_XX = REAL(VECTOR_ELT(R_XX, t));
		    tmp_XYj = XY[t*M+j] - vec_inner_prod1(tmp_XX, theta, t, j, M);
		    beta_j_lasso[t] = fmax(fabs(tmp_XYj) - lambda1[t], 0)*copysign(1.0, tmp_XYj);
        }
        
		bnorm = beta_j_norm1(beta_j_lasso, P);
		
		if (bnorm > 0.0) {
            for (int t = 0; t < P; t++) {
                theta[t*M+j] = fmax(1 - (lambda2*bannot[j]) / bnorm, 0)*beta_j_lasso[t]/Xnorm[j];
                // NB: the 'correct' formula should be (1 - lambda2/bnorm)*beta_j_lasso[t]/Xnorm[j];
			    // but we use ()+ instead to gain stability in practice
            }
        }else{
            for (int t=0; t<P; t++) {
                theta[t*M+j] = 0.0;
            }
        }
	}
}

extern SEXP wrapper2(SEXP R_XX, SEXP R_XY, SEXP R_theta, SEXP R_M, SEXP R_P, SEXP R_beta_j_lasso, SEXP R_lambda1, SEXP R_lambda2, SEXP R_Xnorm, SEXP R_bannot)
{
	SEXP answer;
	
	double* XY = REAL(R_XY);
    double* theta = REAL(R_theta);
    double* beta_j_lasso = REAL(R_beta_j_lasso);
    double* Xnorm = REAL(R_Xnorm);
    double* lambda1 = REAL(R_lambda1);
    double  lambda2 = REAL(R_lambda2)[0];
    int M = INTEGER(R_M)[0];
    int P = INTEGER(R_P)[0];
    double* bannot = REAL(R_bannot);
    
	inner_iter2(R_XX, XY, theta, M, P, beta_j_lasso, lambda1, lambda2, Xnorm, bannot);
	
	PROTECT(answer = allocVector(REALSXP,1));
	REAL(answer)[0] = 1;
  	UNPROTECT(1);
  	return(answer);
}


static const R_CallMethodDef CallEntries[] = {
    {"wrapper2", (DL_FUNC) &wrapper2, 10},
    {NULL, NULL, 0}
};

// Init function called when package loads
void R_init_yourpackagename(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);  // Optional: enforce strict symbol use
}

