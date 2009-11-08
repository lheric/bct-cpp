#include "bct.h"
#include <cmath>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
 * Returns a binary copy of the given matrix.
 */
gsl_matrix* bct::binary(const gsl_matrix* m) {
	gsl_matrix* binary_m = gsl_matrix_calloc(m->size1, m->size2);
	for (int i = 0; i < m->size1; i++) {
		for (int j = 0; j < m->size2; j++) {
			if (is_nonzero(gsl_matrix_get(m, i, j))) {
				gsl_matrix_set(binary_m, i, j, 1.0);
			}
		}
	}
	return binary_m;
}

/*
 * Returns a copy of the given matrix with no loops.
 */
gsl_matrix* bct::no_loops(const gsl_matrix* m) {
	gsl_matrix* no_loops_m = gsl_matrix_alloc(m->size1, m->size2);
	gsl_matrix_memcpy(no_loops_m, m);
	gsl_vector_view diagonal = gsl_matrix_diagonal(no_loops_m);
	gsl_vector_set_zero(&diagonal.vector);
	return no_loops_m;
}

/*
 * Returns a positive copy of the given matrix.
 */
gsl_matrix* bct::positive(const gsl_matrix* m) {
	gsl_matrix* positive_m = gsl_matrix_alloc(m->size1, m->size2);
	for (int i = 0; i < m->size1; i++) {
		for (int j = 0; j < m->size2; j++) {
			gsl_matrix_set(positive_m, i, j, std::abs(gsl_matrix_get(m, i, j)));
		}
	}
	return positive_m;
}

/*
 * Returns an undirected copy of the given matrix.  For every pair of nodes:
 * - If m(i, j) and m(j, i) are the same, no change is made.
 * - If m(i, j) is zero and m(j, i) is nonzero, m(i, j) is set to m(j, i).
 * - If m(i, j) is nonzero and m(j, i) is zero, m(j, i) is set to m(i, j).
 * - If m(i, j) and m(j, i) are different nonzero values, they are both set to
 *   the average of these values.
 */
gsl_matrix* bct::undirected(const gsl_matrix* m, bool upper) {
	if (m->size1 != m->size2) {
		throw size_exception();
	}
	gsl_matrix* undirected_m = gsl_matrix_alloc(m->size1, m->size2);
	gsl_matrix_memcpy(undirected_m, m);
	for (int i = 0; i < m->size1; i++) {
		for (int j = i; j < m->size2; j++) {
			double value_ij = gsl_matrix_get(m, i, j);
			double value_ji = gsl_matrix_get(m, j, i);
			if (is_zero(value_ij) && is_nonzero(value_ji)) {
				gsl_matrix_set(undirected_m, i, j, value_ji);
			} else if (is_nonzero(value_ij) && is_zero(value_ji)) {
				gsl_matrix_set(undirected_m, j, i, value_ij);
			} else if (is_not_equal(value_ij, value_ji)) {
				double average = (value_ij + value_ji) / 2.0;
				gsl_matrix_set(undirected_m, i, j, average);
				gsl_matrix_set(undirected_m, j, i, average);
			}
		}
	}
	return undirected_m;
}
