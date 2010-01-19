#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>

/*
 * Computes the distance matrix for binary undirected graph m
 * Mean distance (excluding the main diagonal) equals the characteristic path length
 * Algebraic shortest path algorithm.
 */

gsl_matrix* bct::distance_bin(const gsl_matrix* m) {
	gsl_matrix* D = eye(m->size1);
	int n = 1;
	gsl_matrix* npath = copy(m);
	gsl_matrix* L = compare_elements(npath, fp_not_equal, 0.0);
	gsl_vector* L_non_zero;
	while((L_non_zero = find(L, 1)) != NULL) {
		gsl_vector_free(L_non_zero);
		//D=D+n.*L;
		gsl_matrix_scale(L, n);
		gsl_matrix_add(D, L);
		n++;
		gsl_matrix* npath_temp = mul(npath, m);
		gsl_matrix_free(npath);
		npath = npath_temp;
		gsl_matrix* npath_ne_zero = compare_elements(npath, fp_not_equal, 0.0);
		gsl_matrix* D_eq_zero = compare_elements(D, fp_equal, 0.0);
		gsl_matrix_mul_elements(npath_ne_zero, D_eq_zero);
		gsl_matrix_free(L);
		L = npath_ne_zero;
		gsl_matrix_free(D_eq_zero);
	}
	if(L_non_zero != NULL) {
		gsl_vector_free(L_non_zero);
	}
	gsl_matrix* not_D = logical_not(D);
	logical_index_assign(D, not_D, GSL_POSINF);
	gsl_matrix* m_eye = eye(m->size1);
	gsl_matrix_sub(D, m_eye);
	
	gsl_matrix_free(npath);
	gsl_matrix_free(L);
	gsl_matrix_free(m_eye);
	
	return D;
}
	
