#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
 * Computes the clustering coefficient for a binary undirected graph.  Results
 * are returned in a vector where each element is the clustering coefficient of
 * the corresponding node.
 */
gsl_vector* bct::clustering_coef_bu(const gsl_matrix* m) {
	if (safe_mode) check_status(m, BINARY | UNDIRECTED, "clustering_coef_bu");
	if (m->size1 != m->size2) {
		throw size_exception();
	}
	gsl_vector* clustering_coef = gsl_vector_calloc(m->size1);
	for (int i = 0; i < m->size1; i++) {
		gsl_vector_const_view row = gsl_matrix_const_row(m, i);
		int k = nnz(&row.vector);
		if (k >= 2) {
			gsl_vector* neighbors = find(&row.vector);
			gsl_matrix* neighbors_m = ordinal_index(m, neighbors, neighbors);
			gsl_vector* connections_v = sum(neighbors_m);
			int connections = sum(connections_v);
			int possible_connections = k * (k - 1);
			gsl_vector_set(clustering_coef, i, (double)connections / (double)possible_connections);
			gsl_vector_free(connections_v);
			gsl_matrix_free(neighbors_m);
			gsl_vector_free(neighbors);
		}
	}
	return clustering_coef;
}
