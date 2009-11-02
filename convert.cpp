#include "bct.h"
#include <cmath>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

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
 * Returns a positive copy of the given matrix.
 */
gsl_matrix* bct::positive(const gsl_matrix* m) {
	gsl_matrix* pm = gsl_matrix_alloc(m->size1, m->size2);
	for (int i = 0; i < m->size1; i++) {
		for (int j = 0; j < m->size2; j++) {
			gsl_matrix_set(pm, i, j, std::abs(gsl_matrix_get(m, i, j)));
		}
	}
	return pm;
}

/*
 * Returns a copy of the given matrix with no loops (all elements along the
 * diagonal set to zero).
 */
gsl_matrix* bct::remove_loops(const gsl_matrix* m) {
	gsl_matrix* m_loopless = gsl_matrix_alloc(m->size1, m->size2);
	gsl_matrix_memcpy(m_loopless, m);
	gsl_vector_view diagonal = gsl_matrix_diagonal(m_loopless);
	gsl_vector_set_zero(&diagonal.vector);
	return m_loopless;
}

/*
 * Returns an undirected copy of the given matrix by mirroring either the upper
 * or lower triangle across the diagonal.
 */
gsl_matrix* bct::undirected(const gsl_matrix* m, bool upper) {
	if (m->size1 != m->size2) {
		throw matrix_size_exception();
	}
	gsl_matrix* um = gsl_matrix_alloc(m->size1, m->size2);
	for (int i = 0; i < m->size1; i++) {
		for (int j = 0; j < m->size2; j++) {
			if ((upper && i <= j) || (!upper && i >= j)) {
				gsl_matrix_set(um, i, j, gsl_matrix_get(m, i, j));
			} else {
				gsl_matrix_set(um, i, j, gsl_matrix_get(m, j, i));
			}
		}
	}
	return um;
}
