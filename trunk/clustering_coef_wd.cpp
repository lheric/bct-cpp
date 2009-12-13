#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
 * Computes the clustering coefficient for a weighted directed graph.  Results
 * are returned in a vector where each element is the clustering coefficient of
 * the corresponding node.
 */
gsl_vector* bct::clustering_coef_wd(const gsl_matrix* m) {
	if (safe_mode) check_status(m, WEIGHTED | DIRECTED, "clustering_coef_wd");
	if (m->size1 != m->size2) {
		throw size_exception();
	}
	
	// A=W~=0;
	gsl_matrix* a = to_binary(m);
	
	// S=W.^(1/3)+(W.').^(1/3);
	gsl_matrix* s = pow_elements(m, 1.0 / 3.0);
	gsl_matrix* transpose_m = gsl_matrix_alloc(m->size2, m->size1);
	gsl_matrix_transpose_memcpy(transpose_m, m);
	gsl_matrix* root3_transpose_m = pow_elements(transpose_m, 1.0 / 3.0);
	gsl_matrix_add(s, root3_transpose_m);
	
	// K=sum(A+A.',2);
	gsl_matrix* k_m = gsl_matrix_alloc(a->size1, a->size2);
	gsl_matrix_transpose_memcpy(k_m, a);
	gsl_matrix_add(k_m, a);
	gsl_vector* k = sum(k_m, 2);
	
	// cyc3=diag(S^3)/2;
	gsl_matrix* cyc3_m = pow(s, 3);
	gsl_matrix_scale(cyc3_m, 0.5);
	gsl_vector_const_view cyc3 = gsl_matrix_const_diagonal(cyc3_m);
	
	// CYC3=K.*(K-1)-2*diag(A^2);
	gsl_vector* k_less_1 = copy(k);
	gsl_vector_add_constant(k_less_1, -1.0);
	gsl_matrix* pow2_a = pow(a, 2);
	gsl_matrix_scale(pow2_a, 2.0);
	gsl_vector_const_view pow2 = gsl_matrix_const_diagonal(pow2_a);
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
	gsl_matrix_free(pow2_a);
	gsl_vector_free(k_less_1);
	gsl_matrix_free(cyc3_m);
	gsl_vector_free(k);
	gsl_matrix_free(k_m);
	gsl_matrix_free(root3_transpose_m);
	gsl_matrix_free(transpose_m);
	gsl_matrix_free(s);
	gsl_matrix_free(a);
	
	return clustering_coef;
}
