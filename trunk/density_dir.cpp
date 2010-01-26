#include "bct.h"
#include <gsl/gsl_matrix.h>

/*
 * Computes density for a directed graph.  Connection weights are ignored.
 */
double bct::density_dir(const gsl_matrix* m) {
	if (safe_mode) check_status(m, DIRECTED, "density_dir");
	
	// N = size(CIJ,1);
	int N = m->size1;
	
	// K = nnz(CIJ);
	int K = nnz(m);
	
	// kden = K/(N^2-N);
	return (double)K / (double)(N * N - N);
}
