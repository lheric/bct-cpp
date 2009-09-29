#include "bct.h"
#include <cmath>
#include <cstdio>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
 * Returns a binary copy of the given vector.
 */
gsl_vector* bct::binary(const gsl_vector* v) {
	gsl_vector* bv = gsl_vector_calloc(v->size);
	for (int i = 0; i < v->size; i++) {
		if (std::abs(gsl_vector_get(v, i)) > EPSILON) {
			gsl_vector_set(bv, i, 1.0);
		}
	}
	return bv;
}

/*
 * Returns a binary copy of the given matrix.
 */
gsl_matrix* bct::binary(const gsl_matrix* m) {
	gsl_matrix* bm = gsl_matrix_calloc(m->size1, m->size2);
	for (int i = 0; i < m->size1; i++) {
		for (int j = 0; j < m->size2; j++) {
			if (std::abs(gsl_matrix_get(m, i, j)) > EPSILON) {
				gsl_matrix_set(bm, i, j, 1.0);
			}
		}
	}
	return bm;
}

/*
 * Returns the number of nonzero elements in the given vector.
 */
int bct::nnz(const gsl_vector* v) {
	int count = 0;
	for (int i = 0; i < v->size; i++) {
		if (std::abs(gsl_vector_get(v, i)) > EPSILON) {
			count++;
		}
	}
	return count;
}

/*
 * Returns the number of nonzero elements in the given matrix.
 */
int bct::nnz(const gsl_matrix* m) {
	int count = 0;
	for (int i = 0; i < m->size1; i++) {
		for (int j = 0; j < m->size2; j++) {
			if (std::abs(gsl_matrix_get(m, i, j)) > EPSILON) {
				count++;
			}
		}
	}
	return count;
}

/*
 * Returns a vector of indices of nonzero elements in the given vector.
 */
gsl_vector* bct::find(const gsl_vector* v) {
	gsl_vector* nz = gsl_vector_alloc(nnz(v));
	int pos = 0;
	for (int i = 0; i < v->size; i++) {
		if (std::abs(gsl_vector_get(v, i)) > EPSILON) {
			gsl_vector_set(nz, pos, i);
			pos++;
		}
	}
	return nz;
}

/*
 * Returns a matrix of size (rows->size, columns->size).  Elements are copied
 * from the specified rows and columns of m.
 */
gsl_matrix* bct::submatrix(const gsl_matrix* m, const gsl_vector* rows, const gsl_vector* columns) {
	gsl_matrix* s = gsl_matrix_alloc(rows->size, columns->size);
	for (int i = 0; i < rows->size; i++) {
		for (int j = 0; j < columns->size; j++) {
			gsl_matrix_set(s, i, j, gsl_matrix_get(m, gsl_vector_get(rows, i), gsl_vector_get(columns, j)));
		}
	}
	return s;
}

/*
 * Prints a vector using the given format for each element.  This is only
 * provided for debugging purposes.  In other cases, use gsl_vector_fprintf.
 */
void bct::printf(const gsl_vector* v, const char* format) {
	for (int i = 0; i < v->size; i++) {
		std::printf(format, gsl_vector_get(v, i));
		std::printf(" ");
	}
	std::printf("\n");
}

/*
 * Prints a matrix using the given format for each element.  This is only
 * provided for debugging purposes.  In other cases, use gsl_matrix_fprintf.
 */
void bct::printf(const gsl_matrix* m, const char* format) {
	for (int i = 0; i < m->size1; i++) {
		for (int j = 0; j < m->size2; j++) {
			std::printf(format, gsl_matrix_get(m, i, j));
			std::printf(" ");
		}
		std::printf("\n");
	}
}
