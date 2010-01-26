#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
 * Computes degree for an undirected graph.  Connection weights are ignored.
 */
gsl_vector* bct::degrees_und(const gsl_matrix* CIJ) {
	if (safe_mode) check_status(CIJ, UNDIRECTED, "degrees_und");
	
	// CIJ = double(CIJ~=0);
	// deg = sum(CIJ);
	gsl_vector* deg = gsl_vector_alloc(CIJ->size2);
	for (int i = 0; i < CIJ->size2; i++) {
		gsl_vector_const_view column = gsl_matrix_const_column(CIJ, i);
		gsl_vector_set(deg, i, nnz(&column.vector));
	}
	return deg;
}
