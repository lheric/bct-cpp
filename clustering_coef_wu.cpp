#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
 * Computes the clustering coefficient for a weighted undirected graph.  Results
 * are returned in a vector where each element is the clustering coefficient of
 * the corresponding node.
 */
gsl_vector* bct::clustering_coef_wu(const gsl_matrix* m) {
	if (safe_mode) check_status(m, WEIGHTED | UNDIRECTED, "clustering_coef_wu");
	if (m->size1 != m->size2) {
		throw size_exception();
	}
	
	// K=sum(W~=0,2);
	gsl_matrix* binary_m = to_binary(m);
	gsl_vector* k = sum(binary_m, 2);
	
	// cyc3=diag((W.^(1/3))^3);
	gsl_matrix* root3_m = pow_elements(m, 1.0 / 3.0);
	gsl_matrix* cyc3_m = pow(root3_m, 3);
	gsl_vector_const_view cyc3 = gsl_matrix_const_diagonal(cyc3_m);
	
	// C=cyc3./(K.*(K-1));
	gsl_vector* k_less_1 = copy(k);
	gsl_vector_add_constant(k_less_1, -1.0);
	gsl_vector_mul(k, k_less_1);
	gsl_vector* clustering_coef = copy(&cyc3.vector);
	gsl_vector_div(clustering_coef, k);
	
	// If no 3-cycles exist, set clustering coefficient to 0
	for (int i = 0; i < m->size1; i++) {
		if (fp_zero(gsl_vector_get(&cyc3.vector, i))) {
			gsl_vector_set(clustering_coef, i, 0.0);
		}
	}

	gsl_vector_free(k_less_1);
	gsl_matrix_free(cyc3_m);
	gsl_matrix_free(root3_m);
	gsl_vector_free(k);
	gsl_matrix_free(binary_m);
	
	return clustering_coef;
}
