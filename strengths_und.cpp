#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
 * Computes the strength for each node in an undirected matrix.
 */
gsl_vector* bct::strengths_und(const gsl_matrix* m) {
	if (safe_mode) check_status(m, UNDIRECTED, "strengths_und");
	
	// str = sum(CIJ);
	return sum(m);
}
