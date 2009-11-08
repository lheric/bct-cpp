#include <cmath>
#include <cstring>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "matlab.h"

/*
 * See MATLAB documentation for descriptions of these functions.  We will
 * document instances where our version differs from the MATLAB version.
 */

gsl_vector* matlab::find(const gsl_vector* v, int n, const char* direction) {
	int size = nnz(v);
	if (size == 0) {
		return NULL;
	}
	gsl_vector* found = gsl_vector_alloc((n < size) ? n : size);
	if (direction == NULL || std::strcmp(direction, "first") == 0) {
		int position = 0;
		for (int i = 0; i < v->size && position < found->size; i++) {
			if (fp_nonzero(gsl_vector_get(v, i))) {
				gsl_vector_set(found, position, i);
				position++;
			}
		}
		return found;
	} else if (std::strcmp(direction, "last") == 0) {
		int position = found->size - 1;
		for (int i = v->size - 1; i >= 0 && position >= 0; i--) {
			if (fp_nonzero(gsl_vector_get(v, i))) {
				gsl_vector_set(found, position, i);
				position--;
			}
		}
		return found;
	} else {
		gsl_vector_free(found);
		return NULL;
	}
}

gsl_vector* matlab::find(const gsl_matrix* m, int n, const char* direction) {
	gsl_vector* v = to_vector(m);
	gsl_vector* ret = find(v, n, direction);
	gsl_vector_free(v);
	return ret;
}

int matlab::nnz(const gsl_vector* v) {
	int nnz = 0;
	for (int i = 0; i < v->size; i++) {
		if (fp_nonzero(gsl_vector_get(v, i))) {
			nnz++;
		}
	}
	return nnz;
}

int matlab::nnz(const gsl_matrix* m) {
	gsl_vector* v = to_vector(m);
	int ret = nnz(v);
	gsl_vector_free(v);
	return ret;
}

double matlab::sum(const gsl_vector* v) {
	double sum = 0.0;
	for (int i = 0; i < v->size; i++) {
		sum += gsl_vector_get(v, i);
	}
	return sum;
}

gsl_vector* matlab::sum(const gsl_matrix* m, int dim) {
	if (dim == 1) {
		gsl_vector* sum = gsl_vector_calloc(m->size2);
		for (int i = 0; i < m->size1; i++) {
			gsl_vector_const_view row = gsl_matrix_const_row(m, i);
			gsl_vector_add(sum, &row.vector);
		}
		return sum;
	} else if (dim == 2) {
		gsl_vector* sum = gsl_vector_calloc(m->size1);
		for (int i = 0; i < m->size2; i++) {
			gsl_vector_const_view column = gsl_matrix_const_column(m, i);
			gsl_vector_add(sum, &column.vector);
		}
		return sum;
	} else {
		return NULL;
	}
}

gsl_matrix* matlab::tril(const gsl_matrix* m, int k) {
	if (k <= -(int)m->size1 || k >= (int)m->size2) {
		return NULL;
	}
	gsl_matrix* tril = gsl_matrix_alloc(m->size1, m->size2);
	gsl_matrix_memcpy(tril, m);
	for (int i = 0; i < m->size1; i++) {
		for (int j = i + k + 1; j < m->size2; j++) {
			if (j >= 0) {
				gsl_matrix_set(tril, i, j, 0.0);
			}
		}
	}
	return tril;
}

gsl_matrix* matlab::triu(const gsl_matrix* m, int k) {
	if (k <= -(int)m->size1 || k >= (int)m->size2) {
		return NULL;
	}
	gsl_matrix* triu = gsl_matrix_alloc(m->size1, m->size2);
	gsl_matrix_memcpy(triu, m);
	for (int i = 0; i < m->size1; i++) {
		for (int j = i + k - 1; j >= 0; j--) {
			if (j < m->size2) {
				gsl_matrix_set(triu, i, j, 0.0);
			}
		}
	}
	return triu;
}

/*
 * Emulates (v1 & v2).
 */
gsl_vector* matlab::logical_and(const gsl_vector* v1, const gsl_vector* v2) {
	if (v1->size != v2->size) {
		return NULL;
	}
	gsl_vector* and_v = gsl_vector_alloc(v1->size);
	for (int i = 0; i < v1->size; i++) {
		bool nz1 = fp_nonzero(gsl_vector_get(v1, i));
		bool nz2 = fp_nonzero(gsl_vector_get(v2, i));
		gsl_vector_set(and_v, i, (double)(nz1 && nz2));
	}
	return and_v;
}

/*
 * Emulates (~v).
 */
gsl_vector* matlab::logical_not(const gsl_vector* v) {
	gsl_vector* not_v = gsl_vector_alloc(v->size);
	for (int i = 0; i < v->size; i++) {
		bool z = fp_zero(gsl_vector_get(v, i));
		gsl_vector_set(not_v, i, (double)z);
	}
	return not_v;
}

/*
 * Emulates (v1 | v2).
 */
gsl_vector* matlab::logical_or(const gsl_vector* v1, const gsl_vector* v2) {
	if (v1->size != v2->size) {
		return NULL;
	}
	gsl_vector* or_v = gsl_vector_alloc(v1->size);
	for (int i = 0; i < v1->size; i++) {
		bool nz1 = fp_nonzero(gsl_vector_get(v1, i));
		bool nz2 = fp_nonzero(gsl_vector_get(v2, i));
		gsl_vector_set(or_v, i, (double)(nz1 || nz2));
	}
	return or_v;
}

/*
 * Emulates (v .^ power).
 */
gsl_vector* matlab::pow_elements(const gsl_vector* v, double power) {
	gsl_vector* pow_v = gsl_vector_alloc(v->size);
	for (int i = 0; i < v->size; i++) {
		double value = std::pow(gsl_vector_get(v, i), power);
		gsl_vector_set(pow_v, i, value);
	}
	return pow_v;
}

/*
 * Emulates (v .^ powers).
 */
gsl_vector* matlab::pow_elements(const gsl_vector* v, const gsl_vector* powers) {
	if (v->size != powers->size) {
		return NULL;
	}
	gsl_vector* pow_v = gsl_vector_alloc(v->size);
	for (int i = 0; i < v->size; i++) {
		double value = std::pow(gsl_vector_get(v, i), gsl_vector_get(powers, i));
		gsl_vector_set(pow_v, i, value);
	}
	return pow_v;
}

/*
 * Emulates (m .^ power).
 */
gsl_matrix* matlab::pow_elements(const gsl_matrix* m, double power) {
	gsl_matrix* pow_m = gsl_matrix_alloc(m->size1, m->size2);
	for (int i = 0; i < m->size1; i++) {
		for (int j = 0; j < m->size2; j++) {
			double value = std::pow(gsl_matrix_get(m, i, j), power);
			gsl_matrix_set(pow_m, i, j, value);
		}
	}
	return pow_m;
}

/*
 * Emulates (m .^ powers).
 */
gsl_matrix* matlab::pow_elements(const gsl_matrix* m, const gsl_matrix* powers) {
	if (m->size1 != powers->size1 || m->size2 != powers->size2) {
		return NULL;
	}
	gsl_matrix* pow_m = gsl_matrix_alloc(m->size1, m->size2);
	for (int i = 0; i < m->size1; i++) {
		for (int j = 0; j < m->size2; j++) {
			double value = std::pow(gsl_matrix_get(m, i, j), gsl_matrix_get(powers, i, j));
			gsl_matrix_set(pow_m, i, j, value);
		}
	}
	return pow_m;
}

/*
 * Emulates (m ^ power).
 */
gsl_matrix* matlab::pow(const gsl_matrix* m, int power) {
	if (m->size1 != m->size2 || power < 1) {
		return NULL;
	}
	gsl_matrix* pow_m = gsl_matrix_alloc(m->size1, m->size2);
	gsl_matrix_memcpy(pow_m, m);
	for (int i = 2; i <= power; i++) {
		gsl_matrix* temp = gsl_matrix_alloc(m->size1, m->size2);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, pow_m, m, 0.0, temp);
		gsl_matrix_free(pow_m);
		pow_m = temp;
	}
	return pow_m;
}

gsl_vector* matlab::compare_elements(const gsl_vector* v, compare_fn cmp_fn, double x) {
	gsl_vector* cmp_v = gsl_vector_alloc(v->size);
	for (int i = 0; i < v->size; i++) {
		double value = gsl_vector_get(v, i);
		gsl_vector_set(cmp_v, i, (double)cmp_fn(value, x));
	}
	return cmp_v;
}

gsl_vector* matlab::compare_elements(const gsl_vector* v1, compare_fn cmp_fn, const gsl_vector* v2) {
	if (v1->size != v2->size) {
		return NULL;
	}
	gsl_vector* cmp_v = gsl_vector_alloc(v1->size);
	for (int i = 0; i < v1->size; i++) {
		double value1 = gsl_vector_get(v1, i);
		double value2 = gsl_vector_get(v2, i);
		gsl_vector_set(cmp_v, i, (double)cmp_fn(value1, value2));
	}
	return cmp_v;
}

gsl_matrix* matlab::compare_elements(const gsl_matrix* m, compare_fn cmp_fn, double x) {
	gsl_matrix* cmp_m = gsl_matrix_alloc(m->size1, m->size2);
	for (int i = 0; i < m->size1; i++) {
		for (int j = 0; j < m->size2; j++) {
			double value = gsl_matrix_get(m, i, j);
			gsl_matrix_set(cmp_m, i, j, (double)cmp_fn(value, x));
		}
	}
	return cmp_m;
}

gsl_matrix* matlab::compare_elements(const gsl_matrix* m1, compare_fn cmp_fn, const gsl_matrix* m2) {
	if (m1->size1 != m2->size1 || m1->size2 != m2->size2) {
		return NULL;
	}
	gsl_matrix* cmp_m = gsl_matrix_alloc(m1->size1, m1->size2);
	for (int i = 0; i < m1->size1; i++) {
		for (int j = 0; j < m1->size2; j++) {
			double value1 = gsl_matrix_get(m1, i, j);
			double value2 = gsl_matrix_get(m2, i, j);
			gsl_matrix_set(cmp_m, i, j, (double)cmp_fn(value1, value2));
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

int matlab::fp_compare(double x, double y) {
	if (fp_zero(x) || fp_zero(y)) {
		
		// gsl_fcmp is not suitable for testing whether a value is approximately zero
		if (fp_zero(x) && fp_zero(y)) {
			return 0;
		} else {
			return (x < y) ? -1 : 1;
		}
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
 * Emulates vector indexing by another vector.
 */
gsl_vector* matlab::index(const gsl_vector* v, const gsl_vector* indices) {
	gsl_vector* indexed = gsl_vector_alloc(indices->size);
	for (int i = 0; i < indices->size; i++) {
		int position = (int)gsl_vector_get(indices, i);
		double value = gsl_vector_get(v, position);
		gsl_vector_set(indexed, i, value);
	}
	return indexed;
}

/*
 * Emulates matrix indexing by a scalar.
 */
double matlab::index(const gsl_matrix* m, int index) {
	int row = index % (int)m->size1;
	int column = index / (int)m->size1;
	return gsl_matrix_get(m, row, column);
}

/*
 * Emulates matrix indexing by two vectors.
 */
gsl_matrix* matlab::index(const gsl_matrix* m, const gsl_vector* rows, const gsl_vector* columns) {
	gsl_matrix* indexed = gsl_matrix_alloc(rows->size, columns->size);
	for (int i = 0; i < rows->size; i++) {
		for (int j = 0; j < columns->size; j++) {
			int row = (int)gsl_vector_get(rows, i);
			int column = (int)gsl_vector_get(columns, j);
			double value = gsl_matrix_get(m, row, column);
			gsl_matrix_set(indexed, i, j, value);
		}
	}
	return indexed;
}

/*
 * Emulates matrix indexing by another matrix.
 */
gsl_matrix* matlab::index(const gsl_matrix* m, const gsl_matrix* indices) {
	gsl_matrix* indexed = gsl_matrix_alloc(indices->size1, indices->size2);
	for (int i = 0; i < indices->size1; i++) {
		for (int j = 0; j < indices->size2; j++) {
			double value = index(m, gsl_matrix_get(indices, i, j));
			gsl_matrix_set(indexed, i, j, value);
		}
	}
	return indexed;
}

/*
 * Emulates vector logical indexing by another vector.
 */
gsl_vector* matlab::logical_index(const gsl_vector* v, const gsl_vector* lv) {
	int size = nnz(lv);
	if (size == 0 || v->size < lv->size) {
		return NULL;
	}
	gsl_vector* indexed = gsl_vector_alloc(size);
	int position = 0;
	for (int i = 0; i < lv->size; i++) {
		if (fp_nonzero(gsl_vector_get(lv, i))) {
			double value = gsl_vector_get(v, i);
			gsl_vector_set(indexed, position, value);
			position++;
		}
	}
	return indexed;
}

/*
 * Emulates matrix logical indexing by a vector.
 */
gsl_vector* matlab::logical_index(const gsl_matrix* m, const gsl_vector* lv) {
	int size = nnz(lv);
	int m_size = (int)m->size1 * (int)m->size2;
	if (size == 0 || m_size < lv->size) {
		return NULL;
	}
	gsl_vector* indexed = gsl_vector_alloc(size);
	int position = 0;
	for (int i = 0; i < lv->size; i++) {
		if (fp_nonzero(gsl_vector_get(lv, i))) {
			double value = index(m, i);
			gsl_vector_set(indexed, position, value);
			position++;
		}
	}
	return indexed;
}

/*
 * Emulates matrix logical indexing by another matrix.
 */
gsl_vector* matlab::logical_index(const gsl_matrix* m, const gsl_matrix* lm) {
	int size = nnz(lm);
	int m_size = (int)m->size1 * (int)m->size2;
	int lm_size = (int)lm->size1 * (int)lm->size2;
	if (size == 0 || m_size < lm_size) {
		return NULL;
	}
	gsl_vector* indexed = gsl_vector_alloc(size);
	int position = 0;
	for (int j = 0; j < lm->size2; j++) {
		for (int i = 0; i < lm->size1; i++) {
			if (fp_nonzero(gsl_matrix_get(lm, i, j))) {
				double value = index(m, j * lm->size1 + i);
				gsl_vector_set(indexed, position, value);
				position++;
			}
		}
	}
	return indexed;
}

/*
 * Emulates matrix-to-vector conversion.  The vector is constructed by
 * consecutively appending columns.
 */
gsl_vector* matlab::to_vector(const gsl_matrix* m) {
	gsl_vector* v = gsl_vector_alloc(m->size1 * m->size2);
	for (int j = 0; j < m->size2; j++) {
		for (int i = 0; i < m->size1; i++) {
			double value = gsl_matrix_get(m, i, j);
			gsl_vector_set(v, j * m->size2 + i, value);
		}
	}
	return v;
}
