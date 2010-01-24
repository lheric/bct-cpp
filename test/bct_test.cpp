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
	gsl_vector* gslv = gsl_vector_alloc(a.length());
	for (int i = 0; i < a.length(); i++) {
		gsl_vector_set(gslv, i, a.elem(i));
	}
	return gslv;
}

/*
 * Converts an Octave matrix to a GSL matrix.
 */
gsl_matrix* bct_test::to_gsl(const Matrix m) {
	gsl_matrix* gslm = gsl_matrix_alloc(m.dims()(0), m.dims()(1));
	for (int i = 0; i < m.dims()(0); i++) {
		for (int j = 0; j < m.dims()(1); j++) {
			gsl_matrix_set(gslm, i, j, m(i, j));
		}
	}
	return gslm;
}

/*
 * Converts an Octave 3D matrix to a GSL 3D matrix.
 */
gsl_matrix** bct_test::to_gsl(const NDArray a) {
	gsl_matrix** gslm = new gsl_matrix*[a.dims()(2)];
	for (int k = 0; k < a.dims()(2); k++) {
		gslm[k] = gsl_matrix_alloc(a.dims()(0), a.dims()(1));
		for (int i = 0; i < a.dims()(0); i++) {
			for (int j = 0; j < a.dims()(1); j++) {
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
	Matrix m = Matrix(1, gslv->size);
	for (int i = 0; i < gslv->size; i++) {
		m(0, i) = gsl_vector_get(gslv, i);
	}
	return m;
}

/*
 * Converts a GSL matrix to an Octave matrix.
 */
Matrix bct_test::from_gsl(const gsl_matrix* gslm) {
	Matrix m = Matrix(gslm->size1, gslm->size2);
	for (int i = 0; i < gslm->size1; i++) {
		for (int j = 0; j < gslm->size2; j++) {
			m(i, j) = gsl_matrix_get(gslm, i, j);
		}
	}
	return m;
}

/*
 * Converts a GSL 3D matrix to an Octave 3D matrix.
 */
NDArray bct_test::from_gsl(const gsl_matrix* const* gslm, int depth) {
	dim_vector dim_v(3);
	dim_v.resize(3);
	dim_v(0) = gslm[0]->size1;
	dim_v(1) = gslm[0]->size2;
	dim_v(2) = depth;
	NDArray m(dim_v);
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
