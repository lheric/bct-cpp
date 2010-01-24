#ifndef BCT_TEST_CPP
#define BCT_TEST_CPP

#include "bct_test.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <octave/oct.h>

/*
 * Converts an Octave array to a GSL vector.
 */
gsl_vector* bct_test::to_gsl(const Array<int> a) {
	int length = a.length();
	gsl_vector* gslv = gsl_vector_alloc(length);
	for (int i = 0; i < length; i++) {
		gsl_vector_set(gslv, i, a.elem(i));
	}
	return gslv;
}

/*
 * Converts an Octave matrix to a GSL matrix.
 */
gsl_matrix* bct_test::to_gsl(const Matrix m) {
	int n_rows = m.dims()(0);
	int n_columns = m.dims()(1);
	gsl_matrix* gslm = gsl_matrix_alloc(n_rows, n_columns);
	for (int i = 0; i < n_rows; i++) {
		for (int j = 0; j < n_columns; j++) {
			gsl_matrix_set(gslm, i, j, m(i, j));
		}
	}
	return gslm;
}

/*
 * Converts an Octave 3D matrix to a GSL 3D matrix.
 */
gsl_matrix** bct_test::to_gsl(const NDArray a) {
	dim_vector dims = a.dims();
	int n_rows = dims(0);
	int n_columns = dims(1);
	int depth = dims(2);
	gsl_matrix** gslm = new gsl_matrix*[depth];
	for (int k = 0; k < depth; k++) {
		gslm[k] = gsl_matrix_alloc(n_rows, n_columns);
		for (int i = 0; i < n_rows; i++) {
			for (int j = 0; j < n_columns; j++) {
				gsl_matrix_set(gslm[k], i, j, a(i, j, k));
			}
		}
	}
	return gslm;
}

/*
 * Converts a GSL vector to an Octave matrix.
 */
Matrix bct_test::from_gsl(const gsl_vector* gslv) {
	gsl_matrix_const_view gslmv = gsl_matrix_const_view_vector(gslv, 1, gslv->size);
	return from_gsl(&gslmv.matrix);
}

/*
 * Converts a GSL matrix to an Octave matrix.
 */
Matrix bct_test::from_gsl(const gsl_matrix* gslm) {
	int n_rows = gslm->size1;
	int n_columns = gslm->size2;
	Matrix m = Matrix(n_rows, n_columns);
	for (int i = 0; i < n_rows; i++) {
		for (int j = 0; j < n_columns; j++) {
			m(i, j) = gsl_matrix_get(gslm, i, j);
		}
	}
	return m;
}

/*
 * Converts a GSL 3D matrix to an Octave 3D matrix.
 */
NDArray bct_test::from_gsl(const gsl_matrix** gslm, int depth) {
	dim_vector dims(3);
	dims(0) = gslm[0]->size1;
	dims(1) = gslm[0]->size2;
	dims(2) = depth;
	NDArray m(dims);
	for (int k = 0; k < depth; k++) {
		for (int i = 0; i < gslm[0]->size1; i++) {
			for (int j = 0; j < gslm[0]->size2; j++) {
				m(i, j, k) = gsl_matrix_get(gslm[k], i, j);
			}
		}
	}
	return m;
}

#endif
