#include "bct.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector.h>

/*
 * Catches GSL errors and throws BCT exceptions.
 */
void bct::gsl_error_handler(const char* reason, const char* file, int line, int gsl_errno) {
	if (gsl_errno == GSL_ENOMEM) {
		throw out_of_memory_exception();
	} else {
		throw gsl_exception();
	}
}

/*
 * Overloaded convenience function for freeing GSL vectors and matrices.
 */
void bct::gsl_free(gsl_vector* v) { gsl_vector_free(v); }
void bct::gsl_free(gsl_matrix* m) { gsl_matrix_free(m); }

/*
 * Initializes the BCT library for external use.
 */
void bct::init() {
	gsl_set_error_handler(gsl_error_handler);
}

/*
 * Permutes the columns of a matrix.
 */
gsl_matrix* bct::permute_columns(const gsl_permutation* p, const gsl_matrix* m) {
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
gsl_matrix* bct::permute_rows(const gsl_permutation* p, const gsl_matrix* m) {
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

/*
 * Allocates a square matrix and initializes all elements to the given value.
 */
gsl_matrix* bct::yens(int size, double value) {
	return yens(size, size, value);
}

/*
 * Allocates a matrix and initializes all elements to the given value.
 */
gsl_matrix* bct::yens(int size1, int size2, double value) {
	gsl_matrix* m = gsl_matrix_alloc(size1, size2);
	gsl_matrix_set_all(m, value);
	return m;
}
