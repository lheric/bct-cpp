#include "bct.h"
#include <gsl/gsl_matrix.h>

/*
 * Computes the connection density of the directed graph whose
 * adjaceny matrix is m. The total number of possible connections 
 * is given by the form N(N-1). N is the number of nodes.
 */
double bct::density_dir(const gsl_matrix* m) {
	int N,K;
	double kden;
	N = m->size1;
	K = nnz(m);
	kden = (double)K/(double)(N*N-N);
	return kden;
}

