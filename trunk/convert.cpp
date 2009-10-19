#include "bct.h"
#include <cmath>
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
 * Returns a copy of the given matrix with all elements along the diagonal set
 * to zero.
 */
gsl_matrix* bct::zero_diagonal(const gsl_matrix* m) {
	gsl_matrix* zdm = gsl_matrix_alloc(m->size1, m->size2);
	gsl_matrix_memcpy(zdm, m);
	gsl_vector_view diagonal = gsl_matrix_diagonal(zdm);
	gsl_vector_set_zero(&diagonal.vector);
	return zdm;
}
