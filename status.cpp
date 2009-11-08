#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <iostream>

bool bct::safe_mode = false;

/*
 * Returns whether a matrix matches the given status flags.  If the check fails,
 * a message is printed to stderr starting with the given text.
 */
bool bct::check_status(const gsl_matrix* m, int flags, const char* text) {
	bool ret = true;
	if (flags & UNDIRECTED) ret = ret && is_undirected(m);
	if (flags & DIRECTED) ret = ret && is_directed(m);
	if (flags & BINARY) ret = ret && is_binary(m);
	if (flags & WEIGHTED) ret = ret && is_weighted(m);
	if (flags & POSITIVE) ret = ret && is_positive(m);
	if (flags & SIGNED) ret = ret && is_signed(m);
	if (flags & NO_LOOPS) ret = ret && has_no_loops(m);
	if (flags & LOOPS) ret = ret && has_loops(m);
	if (!ret) {
		if (text != NULL) std::cerr << text << ": ";
		std::cerr << "Matrix status check failed." << std::endl;
	}
	return ret;
}

/*
 * Returns whether the given matrix is undirected.
 */
bool bct::is_undirected(const gsl_matrix* m) {
	return !is_directed(m);
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
			double value_ij = gsl_matrix_get(m, i, j);
			double value_ji = gsl_matrix_get(m, j, i);
			if (is_not_equal(value_ij, value_ji)) {
				return true;
			}
		}
	}
	return false;
}

/*
 * Returns whether the given matrix is binary.
 */
bool bct::is_binary(const gsl_matrix* m) {
	return !is_weighted(m);
}

/*
 * Returns whether the given matrix is weighted.
 */
bool bct::is_weighted(const gsl_matrix* m) {
	for (int i = 0; i < m->size1; i++) {
		for (int j = 0; j < m->size2; j++) {
			double value = gsl_matrix_get(m, i, j);
			if (is_nonzero(value) && is_not_equal(value, 1.0)) {
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
			double value = gsl_matrix_get(m, i, j);
			if (is_negative(value)) {
				return true;
			}
		}
	}
	return false;
}

/*
 * Returns whether the given matrix has loops.
 */
bool bct::has_loops(const gsl_matrix* m) {
	return !has_no_loops(m);
}

/*
 * Returns whether the given matrix has no loops.
 */
bool bct::has_no_loops(const gsl_matrix* m) {
	if (m->size1 != m->size2) {
		throw size_exception();
	}
	gsl_vector_const_view diagonal = gsl_matrix_const_diagonal(m);
	return gsl_vector_isnull(&diagonal.vector) == 1;
}
