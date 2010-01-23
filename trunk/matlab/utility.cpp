#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector.h>
#include "matlab.h"

/*
 * Permutes the columns of a matrix.
 */
gsl_matrix* matlab::permute_columns(const gsl_permutation* p, const gsl_matrix* m) {
	if (p->size != m->size2) {
		return NULL;
	}
	gsl_matrix* permute_m = gsl_matrix_alloc(m->size1, m->size2);
	for (int i = 0; i < p->size; i++) {
		int i_column = gsl_permutation_get(p, i);
		gsl_vector_const_view column = gsl_matrix_const_column(m, i_column);
		gsl_matrix_set_col(permute_m, i, &column.vector);
	}
	return permute_m;
}

/*
 * Permutes the rows of a matrix.
 */
gsl_matrix* matlab::permute_rows(const gsl_permutation* p, const gsl_matrix* m) {
	if (p->size != m->size1) {
		return NULL;
	}
	gsl_matrix* permute_m = gsl_matrix_alloc(m->size1, m->size2);
	for (int i = 0; i < p->size; i++) {
		int i_row = gsl_permutation_get(p, i);
		gsl_vector_const_view row = gsl_matrix_const_row(m, i_row);
		gsl_matrix_set_row(permute_m, i, &row.vector);
	}
	return permute_m;
}
