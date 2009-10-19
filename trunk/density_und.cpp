#include "bct.h"
#include <gsl/gsl_matrix.h>

/*
 * Computes the connection density of the undirected graph whose
 * adjaceny matrix is m. The total number of possible connections 
 * is given by the form N(N-1)/2. N is the number of nodes.
 */
double bct::density_und(const gsl_matrix* m) {
	int N,K;
	double kden;
	N = m->size1;
	gsl_matrix* triu_m = triu(m,0);
	K = nnz(triu_m);
	kden = (double)K/((double)(N*N-N)/2.0);
	return kden;
}

