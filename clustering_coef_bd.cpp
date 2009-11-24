#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
 * Computes the clustering coefficient for a binary directed graph.  Results are
 * returned in a vector where each element is the clustering coefficient of the
 * corresponding node.
 */
gsl_vector* bct::clustering_coef_bd(const gsl_matrix* m) {
	if (safe_mode) check_status(m, BINARY | DIRECTED, "clustering_coef_bd");
	if (m->size1 != m->size2) {
		throw size_exception();
	}
	
	// S=A+A.';
	gsl_matrix* s = gsl_matrix_alloc(m->size1, m->size2);
	gsl_matrix_transpose_memcpy(s, m);
	gsl_matrix_add(s, m);
	
	// K=sum(S,2);
	gsl_vector* k = sum(s, 2);
	
	// cyc3=diag(S^3)/2;
	gsl_matrix* cyc3_m = pow(s, 3);
	gsl_matrix_scale(cyc3_m, 0.5);
	gsl_vector_const_view cyc3 = gsl_matrix_const_diagonal(cyc3_m);
	
	// CYC3=K.*(K-1)-2*diag(A^2);
	gsl_vector* k_less_1 = copy(k);
	gsl_vector_add_constant(k_less_1, -1.0);
	gsl_matrix* pow2_m = pow(m, 2);
	gsl_matrix_scale(pow2_m, 2.0);
	gsl_vector_const_view pow2 = gsl_matrix_const_diagonal(pow2_m);
	gsl_vector* possible_cyc3 = copy(k);
	gsl_vector_mul(possible_cyc3, k_less_1);
	gsl_vector_sub(possible_cyc3, &pow2.vector);
	
	// C=cyc3./CYC3
	gsl_vector* clustering_coef = copy(&cyc3.vector);
	gsl_vector_div(clustering_coef, possible_cyc3);
	
	// If no 3-cycles exist, set clustering coefficient to 0
	for (int i = 0; i < m->size1; i++) {
		if (fp_zero(gsl_vector_get(&cyc3.vector, i))) {
			gsl_vector_set(clustering_coef, i, 0.0);
		}
	}
	
	gsl_vector_free(possible_cyc3);
	gsl_matrix_free(pow2_m);
	gsl_vector_free(k_less_1);
	gsl_matrix_free(cyc3_m);
	gsl_vector_free(k);
	gsl_matrix_free(s);
	
	return clustering_coef;
}
