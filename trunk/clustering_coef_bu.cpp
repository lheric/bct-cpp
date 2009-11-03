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
	gsl_matrix* m_loopless = remove_loops(m);
	for (int i = 0; i < m_loopless->size1; i++) {
		gsl_vector_const_view row = gsl_matrix_const_row(m_loopless, i);
		int k = nnz(&row.vector);
		if (k >= 2) {
			gsl_vector* neighbors = find(&row.vector);
			gsl_matrix* s = submatrix(m_loopless, neighbors, neighbors);
			int actual_connections = nnz(s);
			int possible_connections = k * (k - 1);
			gsl_vector_set(clustering_coef, i, (double)actual_connections / (double)possible_connections);
			gsl_matrix_free(s);
			gsl_vector_free(neighbors);
		}
	}
	gsl_matrix_free(m_loopless);
	return clustering_coef;
}
