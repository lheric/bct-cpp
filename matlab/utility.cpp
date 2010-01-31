#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector.h>
#include "matlab.h"

/*
 * Permutes the columns of a matrix.
 */
gsl_matrix* matlab::permute_columns(const gsl_permutation* p, const gsl_matrix* m) {
	if (p->size != m->size2) return NULL;
	gsl_matrix* permuted_m = gsl_matrix_alloc(m->size1, m->size2);
	for (int i = 0; i < (int)p->size; i++) {
		int i_col = gsl_permutation_get(p, i);
		gsl_vector_const_view m_col_i_col = gsl_matrix_const_column(m, i_col);
		gsl_matrix_set_col(permuted_m, i, &m_col_i_col.vector);
	}
	return permuted_m;
}

/*
 * Permutes the rows of a matrix.
 */
gsl_matrix* matlab::permute_rows(const gsl_permutation* p, const gsl_matrix* m) {
	if (p->size != m->size1) return NULL;
	gsl_matrix* permuted_m = gsl_matrix_alloc(m->size1, m->size2);
	for (int i = 0; i < (int)p->size; i++) {
		int i_row = gsl_permutation_get(p, i);
		gsl_vector_const_view m_row_i_row = gsl_matrix_const_row(m, i_row);
		gsl_matrix_set_row(permuted_m, i, &m_row_i_row.vector);
	}
	return permuted_m;
}
