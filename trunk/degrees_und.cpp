#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
 * Computes the degree for each node in an undirected binary matrix.  Weights
 * are discarded.  Results are returned in a vector where each element is the
 * degree of the corresponding node.
 */
gsl_vector* bct::degrees_und(const gsl_matrix* m) {
	gsl_vector* degrees = gsl_vector_calloc(m->size2);
	for (int i = 0; i < m->size2; i++) {
		gsl_vector_const_view column = gsl_matrix_const_column(m, i);
		gsl_vector_set(degrees, i, bct::nnz(&column.vector));
	}
	return degrees;
}
