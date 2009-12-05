#include "bct_test.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <octave/oct.h>
#include <cstdio>

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
 * Converts an Octave array to a GSL vector.
 */
gsl_vector* bct_test::to_gsl(const Array<int> a) {
	int len = a.length();
	gsl_vector* gslv = gsl_vector_alloc(len);
	for (int i = 0; i < len; i++) {
		double elem = a.elem(i);
		gsl_vector_set(gslv, i, elem);
	}
	return gslv;
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
 * Converts a GSL 3D matrix to an Octave 3D matrix.
 */
NDArray bct_test::from_gsl(gsl_matrix** gslm, int gslm_len) {
	//First identify the dimensions of the Nd-array to be filled
	dim_vector dimv(3);
	dimv.resize(3);
	dimv(0) = (*gslm)->size1;
	dimv(1) = (*gslm)->size2;
	dimv(2) = gslm_len;
	NDArray nd_oct(dimv);
	int rows = dimv(0);
	int columns = dimv(1);
	for(int i = 0;i < gslm_len;i++) {
		for (int j = 0; j < rows; j++) {
			for (int k = 0; k < columns; k++) {
				nd_oct(j, k, i) = gsl_matrix_get(gslm[i], j, k);
			}
		}
	}
	return nd_oct;
}

/*
 * Converts a Octave 3D matrix to GSL 3D matrix.
 */
gsl_matrix** bct_test::to_gsl(const NDArray arr) {
	dim_vector dimv = arr.dims();
	int depth = dimv(2);
	int rows = dimv(0);
	int columns = dimv(1);
	gsl_matrix** gslm = new gsl_matrix* [depth];
	for(int i = 0;i < depth;i++) {
		gslm[i] = gsl_matrix_alloc(rows, columns);
		for (int j = 0; j < rows; j++) {
			for (int k = 0; k < columns; k++) {
				gsl_matrix_set(gslm[i], j, k, arr(j, k, i));
			}
		}
	}
	return gslm;
}

/*
 * Converts a GSL vector to an Octave matrix.
 */
Matrix bct_test::from_gsl(const gsl_vector* gslv) {
	/*
	gsl_matrix_const_view gslmv = gsl_matrix_const_view_vector(gslv, 1, gslv->size);
	return from_gsl(&gslmv.matrix);
	*/
	int columns = gslv->size;
	Matrix m = Matrix(1, columns);
	for (int i = 0;i < columns;i++)
		m(0, i) = gsl_vector_get(gslv, i);
	return m;
}




