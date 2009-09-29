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
