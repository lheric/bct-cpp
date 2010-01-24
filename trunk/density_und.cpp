#include "bct.h"
#include <gsl/gsl_matrix.h>

/*
 * Computes the connection density of an undirected matrix.  Connection weights
 * are ignored.
 */
double bct::density_und(const gsl_matrix* m) {
	if (safe_mode) check_status(m, UNDIRECTED, "density_und");
	
	// N = size(CIJ,1);
	int N = m->size1;
	
	// K = nnz(triu(CIJ));
	gsl_matrix* triu_m = triu(m);
	int K = nnz(triu_m);
	gsl_matrix_free(triu_m);
	
	// kden = K/((N^2-N)/2);
	return (double)K / ((double)(N * N - N) / 2.0);
}
