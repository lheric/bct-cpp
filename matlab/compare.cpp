#include <cmath>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "matlab.h"

/*
 * Compares two floating-point numbers.
 */
int matlab::fp_compare(double x, double y) {
	
	// gsl_fcmp is not suitable for testing whether a value is approximately zero
	if (fp_zero(x) && fp_zero(y)) {
		return 0;
	} else {
		return gsl_fcmp(x, y, EPSILON);
	}
}

bool matlab::fp_equal(double x, double y) { return fp_compare(x, y) == 0; }
bool matlab::fp_not_equal(double x, double y) { return fp_compare(x, y) != 0; }
bool matlab::fp_zero(double x) { return std::abs(x) < EPSILON; }
bool matlab::fp_nonzero(double x) { return std::abs(x) > EPSILON; }
bool matlab::fp_positive(double x) { return x > EPSILON; }
bool matlab::fp_negative(double x) { return x < -EPSILON; }

/*
 * Compares two vectors lexicographically, returning -1, 0, or 1 if the first
 * vector is less than, equal to, or greater than the second.
 */
int matlab::compare_vectors(const gsl_vector* v1, const gsl_vector* v2) {
	for (int i = 0; i < v1->size; i++) {
		if (i >= v2->size) {
			return 1;
		}
		int result = fp_compare(gsl_vector_get(v1, i), gsl_vector_get(v2, i));
		if (result != 0) {
			return result;
		}
	}
	if (v1->size < v2->size) {
		return -1;
	} else {
		return 0;
	}
}

/*
 * For use with gsl_heapsort.
 */
int matlab::compare_vectorps(const gsl_vector** v1, const gsl_vector** v2) {
	return compare_vectors(*v1, *v2);
}

/*
 * Compares two matrices lexicographically, returning -1, 0, or 1 if the first
 * matrix is less than, equal to, or greater than the second.
 */
int matlab::compare_matrices(const gsl_matrix* m1, const gsl_matrix* m2) {
	int size1 = (int)m1->size1 * (int)m1->size2;
	int size2 = (int)m2->size1 * (int)m2->size2;
	for (int i = 0; i < size1; i++) {
		if (i >= size2) {
			return 1;
		}
		int result = fp_compare(ordinal_index(m1, i), ordinal_index(m2, i));
		if (result != 0) {
			return result;
		}
	}
	if (size1 < size2) {
		return -1;
	} else {
		return 0;
	}
}

/*
 * For use with gsl_heapsort.
 */
int matlab::compare_matrixps(const gsl_matrix** m1, const gsl_matrix** m2) {
	return compare_matrices(*m1, *m2);
}

/*
 * Emulates (v op x), where op is a binary comparison operator.
 */
gsl_vector* matlab::compare_elements(const gsl_vector* v, fp_cmp_fn compare, double x) {
	gsl_vector* cmp_v = gsl_vector_alloc(v->size);
	for (int i = 0; i < v->size; i++) {
		double value = gsl_vector_get(v, i);
		gsl_vector_set(cmp_v, i, (double)compare(value, x));
	}
	return cmp_v;
}

/*
 * Emulates (v1 op v2), where op is a binary comparison operator.
 */
gsl_vector* matlab::compare_elements(const gsl_vector* v1, fp_cmp_fn compare, const gsl_vector* v2) {
	if (v1->size != v2->size) {
		return NULL;
	}
	gsl_vector* cmp_v = gsl_vector_alloc(v1->size);
	for (int i = 0; i < v1->size; i++) {
		double value1 = gsl_vector_get(v1, i);
		double value2 = gsl_vector_get(v2, i);
		gsl_vector_set(cmp_v, i, (double)compare(value1, value2));
	}
	return cmp_v;
}

/*
 * Emulates (m op x), where op is a binary comparison operator.
 */
gsl_matrix* matlab::compare_elements(const gsl_matrix* m, fp_cmp_fn compare, double x) {
	gsl_matrix* cmp_m = gsl_matrix_alloc(m->size1, m->size2);
	for (int i = 0; i < m->size1; i++) {
		for (int j = 0; j < m->size2; j++) {
			double value = gsl_matrix_get(m, i, j);
			gsl_matrix_set(cmp_m, i, j, (double)compare(value, x));
		}
	}
	return cmp_m;
}

/*
 * Emulates (m1 op m2), where op is a binary comparison operator.
 */
gsl_matrix* matlab::compare_elements(const gsl_matrix* m1, fp_cmp_fn compare, const gsl_matrix* m2) {
	if (m1->size1 != m2->size1 || m1->size2 != m2->size2) {
		return NULL;
	}
	gsl_matrix* cmp_m = gsl_matrix_alloc(m1->size1, m1->size2);
	for (int i = 0; i < m1->size1; i++) {
		for (int j = 0; j < m1->size2; j++) {
			double value1 = gsl_matrix_get(m1, i, j);
			double value2 = gsl_matrix_get(m2, i, j);
			gsl_matrix_set(cmp_m, i, j, (double)compare(value1, value2));
		}
	}
	return cmp_m;
}

bool matlab::cmp_equal(double x, double y) { return fp_equal(x, y); }
bool matlab::cmp_not_equal(double x, double y) { return fp_not_equal(x, y); }
bool matlab::cmp_greater(double x, double y) { return fp_compare(x, y) == 1; }
bool matlab::cmp_greater_or_equal(double x, double y) { return fp_compare(x, y) == 1 || fp_equal(x, y); }
bool matlab::cmp_less(double x, double y) { return fp_compare(x, y) == -1; }
bool matlab::cmp_less_or_equal(double x, double y) { return fp_compare(x, y) == -1 || fp_equal(x, y); }
