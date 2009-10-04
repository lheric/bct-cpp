#include "bct_test.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <octave/oct.h>

/*
 * Converts an Octave matrix to a GSL matrix.
 */
gsl_matrix* bct_test::to_gsl(const Matrix m) {
	int rows = m.dims()(0);
	int columns = m.dims()(1);
	gsl_matrix* gslm = gsl_matrix_alloc(rows, columns);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++) {
			gsl_matrix_set(gslm, i, j, m(i, j));
		}
	}
	return gslm;
}

/*
 * Converts a GSL matrix to an Octave matrix.
 */
Matrix bct_test::from_gsl(const gsl_matrix* gslm) {
	int rows = gslm->size1;
	int columns = gslm->size2;
	Matrix m = Matrix(rows, columns);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++) {
			m(i, j) = gsl_matrix_get(gslm, i, j);
		}
	}
	return m;
}

/*
 * Converts a GSL vector to an Octave matrix.
 */
Matrix bct_test::from_gsl(const gsl_vector* gslv) {
	gsl_matrix_const_view gslmv = gsl_matrix_const_view_vector(gslv, 1, gslv->size);
	return from_gsl(&gslmv.matrix);
}
