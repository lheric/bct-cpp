#include "bct.h"
#include <cassert>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
 * Computes the clustering coefficient for a binary undirected graph.  Results
 * are returned in a vector where each element is the clustering coefficient of
 * the corresponding node.
 */
gsl_vector* bct::clustering_coef_bu(const gsl_matrix* m) {
	assert(m->size1 == m->size2);
	gsl_vector* clustering_coef = gsl_vector_calloc(m->size1);
	for (int i = 0; i < m->size1; i++) {
		gsl_vector_const_view row = gsl_matrix_const_row(m, i);
		int k = bct::nnz(&row.vector);
		if (k >= 2) {
			gsl_vector* neighbors = bct::find(&row.vector);
			gsl_matrix* s = bct::submatrix(m, neighbors, neighbors);
			int actual_connections = bct::nnz(s);
			int possible_connections = k * (k - 1);
			gsl_vector_set(clustering_coef, i, (double)actual_connections / (double)possible_connections);
			gsl_matrix_free(s);
			gsl_vector_free(neighbors);
		}
	}
	return clustering_coef;
}
