#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
 * Computes the degree for each node in an undirected matrix.
 */
gsl_vector* bct::degrees_und(const gsl_matrix* m) {
	if (safe_mode) check_status(m, UNDIRECTED, "degrees_und");
	
	// CIJ = double(CIJ~=0);
	// deg = sum(CIJ);
	gsl_vector* deg = gsl_vector_alloc(m->size2);
	for (int i = 0; i < m->size2; i++) {
		gsl_vector_const_view column = gsl_matrix_const_column(m, i);
		gsl_vector_set(deg, i, nnz(&column.vector));
	}
	return deg;
}
