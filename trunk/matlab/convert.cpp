#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "matlab.h"

/*
 * Emulates matrix-to-vector conversion.  The vector is constructed by
 * consecutively appending columns.
 */
gsl_vector* matlab::to_vector(const gsl_matrix* m) {
	gsl_vector* v = gsl_vector_alloc(m->size1 * m->size2);
	for (int j = 0; j < m->size2; j++) {
		for (int i = 0; i < m->size1; i++) {
			double value = gsl_matrix_get(m, i, j);
			gsl_vector_set(v, j * m->size1 + i, value);
		}
	}
	return v;
}

/*
 * Converts a vector to a column matrix
 */
gsl_matrix* matlab::to_column_matrix(const gsl_vector* v) {
	gsl_matrix* m = gsl_matrix_alloc(v->size, 1);
	for(int i=0;i < v->size;i++) {
		gsl_matrix_set(m, i, 0, gsl_vector_get(v, i));
	}
	return m;
}

/*
 * Converts a vector to a row matrix
 */
gsl_matrix* matlab::to_row_matrix(const gsl_vector* v) {
	gsl_matrix* m = gsl_matrix_alloc(1, v->size);
	for(int i=0;i < v->size;i++) {
		gsl_matrix_set(m, 0, i, gsl_vector_get(v, i));
	}
	return m;
}
