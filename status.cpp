#include "bct.h"
#include <cmath>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
 * Returns whether a matrix matches the given status flags.
 */
bool bct::check_status(const gsl_matrix* m, int flags) {
	bool matches = true;
	if (flags & UNDIRECTED) matches = matches && is_undirected(m);
	if (flags & DIRECTED) matches = matches && is_directed(m);
	if (flags & BINARY) matches = matches && is_binary(m);
	if (flags & WEIGHTED) matches = matches && is_weighted(m);
	if (flags & POSITIVE) matches = matches && is_positive(m);
	if (flags & SIGNED) matches = matches && is_signed(m);
	if (flags & NO_LOOPS) matches = matches && has_no_loops(m);
	if (flags & LOOPS) matches = matches && has_loops(m);
	return matches;
}

/*
 * Returns whether the given matrix has loops (nonzero elements along the
 * diagonal).
 */
bool bct::has_loops(const gsl_matrix* m) {
	if (m->size1 != m->size2) {
		throw size_exception();
	}
	for (int i = 0; i < m->size1; i++) {
		if (std::abs(gsl_matrix_get(m, i, i)) > EPSILON) {
			return true;
		}
	}
	return false;
}

/*
 * Returns whether the given matrix has no loops (nonzero elements along the
 * diagonal).
 */
bool bct::has_no_loops(const gsl_matrix* m) {
	return !has_loops(m);
}

/*
 * Returns whether the given matrix is binary.
 */
bool bct::is_binary(const gsl_matrix* m) {
	return !is_weighted(m);
}

/*
 * Returns whether the given matrix is directed.
 */
bool bct::is_directed(const gsl_matrix* m) {
	if (m->size1 != m->size2) {
		throw size_exception();
	}
	for (int i = 0; i < m->size1; i++) {
		for (int j = 0; j < m->size2; j++) {
			if (std::abs(gsl_matrix_get(m, i, j) - gsl_matrix_get(m, j, i)) > EPSILON) {
				return true;
			}
		}
	}
	return false;
}

/*
 * Returns whether the given matrix is positive.
 */
bool bct::is_positive(const gsl_matrix* m) {
	return !is_signed(m);
}

/*
 * Returns whether the given matrix is signed.
 */
bool bct::is_signed(const gsl_matrix* m) {
	for (int i = 0; i < m->size1; i++) {
		for (int j = 0; j < m->size2; j++) {
			if (gsl_matrix_get(m, i, j) < -EPSILON) {
				return true;
			}
		}
	}
	return false;
}

/*
 * Returns whether the given matrix is undirected.
 */
bool bct::is_undirected(const gsl_matrix* m) {
	return !is_directed(m);
}

/*
 * Returns whether the given matrix is weighted.
 */
bool bct::is_weighted(const gsl_matrix* m) {
	for (int i = 0; i < m->size1; i++) {
		for (int j = 0; j < m->size2; j++) {
			if (std::abs(gsl_matrix_get(m, i, j)) > EPSILON &&
				std::abs(gsl_matrix_get(m, i, j) - 1) > EPSILON) {
				return true;
			}
		}
	}
	return false;
}
