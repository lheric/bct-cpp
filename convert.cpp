#include "bct.h"
#include <cmath>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
 * Returns a copy of the given matrix with no loops.
 */
gsl_matrix* bct::remove_loops(const gsl_matrix* m) {
	gsl_matrix* nl_m = copy(m);
	gsl_vector_view diag_nl_m = gsl_matrix_diagonal(nl_m);
	gsl_vector_set_zero(&diag_nl_m.vector);
	return nl_m;
}

/*
 * Returns a binary copy of the given matrix.
 */
gsl_matrix* bct::to_binary(const gsl_matrix* m) {
	return compare_elements(m, fp_not_equal, 0.0);
}

/*
 * Returns a positive copy of the given matrix.
 */
gsl_matrix* bct::to_positive(const gsl_matrix* m) {
	gsl_matrix* pos_m = gsl_matrix_alloc(m->size1, m->size2);
	for (int i = 0; i < (int)m->size1; i++) {
		for (int j = 0; j < (int)m->size2; j++) {
			gsl_matrix_set(pos_m, i, j, std::abs(gsl_matrix_get(m, i, j)));
		}
	}
	return pos_m;
}

/*
 * Returns an undirected copy of the given matrix.  For every pair of nodes:
 * - If m(i, j) and m(j, i) are the same, no change is made.
 * - If m(i, j) is zero and m(j, i) is nonzero, m(i, j) is set to m(j, i).
 * - If m(i, j) is nonzero and m(j, i) is zero, m(j, i) is set to m(i, j).
 * - If m(i, j) and m(j, i) are different nonzero values, they are both set to
 *   the average of these values.
 */
gsl_matrix* bct::to_undirected(const gsl_matrix* m) {
	if (m->size1 != m->size2) throw size_exception();
	gsl_matrix* und_m = copy(m);
	for (int i = 0; i < (int)m->size1; i++) {
		for (int j = i; j < (int)m->size2; j++) {
			double value_ij = gsl_matrix_get(m, i, j);
			double value_ji = gsl_matrix_get(m, j, i);
			if (fp_zero(value_ij) && fp_nonzero(value_ji)) {
				gsl_matrix_set(und_m, i, j, value_ji);
			} else if (fp_nonzero(value_ij) && fp_zero(value_ji)) {
				gsl_matrix_set(und_m, j, i, value_ij);
			} else if (fp_not_equal(value_ij, value_ji)) {
				double average = (value_ij + value_ji) / 2.0;
				gsl_matrix_set(und_m, i, j, average);
				gsl_matrix_set(und_m, j, i, average);
			}
		}
	}
	return und_m;
}
