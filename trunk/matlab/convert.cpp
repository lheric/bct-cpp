#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
 * Converts a vector to an array.
 */
void matlab::to_array(const VECTOR_TYPE* v, FP_TYPE* array) {
	for (int i = 0; i < (int)v->size; i++) {
		array[i] = VECTOR_ID(get)(v, i);
	}
}

/*
 * Converts a vector to a boolean: true if all elements are nonzero, false
 * otherwise.
 */
bool matlab::to_bool(const VECTOR_TYPE* v) {
	return all(v) == 1;
}

/*
 * Converts a matrix to a boolean: true if all elements are nonzero, false
 * otherwise.
 */
bool matlab::to_bool(const MATRIX_TYPE* m) {
	VECTOR_TYPE* all_v = all(m);
	bool ret = all(all_v);
	VECTOR_ID(free)(all_v);
	return ret;
}

/*
 * Converts a vector to float precision.
 */
gsl_vector_float* matlab::to_vector_float(const VECTOR_TYPE* v) {
	gsl_vector_float* float_v = gsl_vector_float_alloc(v->size);
	for (int i = 0; i < (int)v->size; i++) {
		gsl_vector_float_set(float_v, i, (float)VECTOR_ID(get)(v, i));
	}
	return float_v;
}

/*
 * Converts a vector to double precision.
 */
gsl_vector* matlab::to_vector_double(const VECTOR_TYPE* v) {
	gsl_vector* double_v = gsl_vector_alloc(v->size);
	for (int i = 0; i < (int)v->size; i++) {
		gsl_vector_set(double_v, i, (double)VECTOR_ID(get)(v, i));
	}
	return double_v;
}

/*
 * Converts a vector to long double precision.
 */
gsl_vector_long_double* matlab::to_vector_long_double(const VECTOR_TYPE* v) {
	gsl_vector_long_double* long_double_v = gsl_vector_long_double_alloc(v->size);
	for (int i = 0; i < (int)v->size; i++) {
		gsl_vector_long_double_set(long_double_v, i, (long double)VECTOR_ID(get)(v, i));
	}
	return long_double_v;
}

/*
 * Converts a matrix to a vector.  The vector is constructed by consecutively
 * appending columns.
 */
VECTOR_TYPE* matlab::to_vector(const MATRIX_TYPE* m) {
	VECTOR_TYPE* v = VECTOR_ID(alloc)(m->size1 * m->size2);
	for (int j = 0; j < (int)m->size2; j++) {
		for (int i = 0; i < (int)m->size1; i++) {
			FP_TYPE value = MATRIX_ID(get)(m, i, j);
			VECTOR_ID(set)(v, j * m->size1 + i, value);
		}
	}
	return v;
}

/*
 * Converts a matrix to float precision.
 */
gsl_matrix_float* matlab::to_matrix_float(const MATRIX_TYPE* m) {
	gsl_matrix_float* float_m = gsl_matrix_float_alloc(m->size1, m->size2);
	for (int i = 0; i < (int)m->size1; i++) {
		for (int j = 0; j < (int)m->size2; j++) {
			gsl_matrix_float_set(float_m, i, j, (float)MATRIX_ID(get)(m, i, j));
		}
	}
	return float_m;
}

/*
 * Converts a matrix to float precision.
 */
gsl_matrix* matlab::to_matrix_double(const MATRIX_TYPE* m) {
	gsl_matrix* double_m = gsl_matrix_alloc(m->size1, m->size2);
	for (int i = 0; i < (int)m->size1; i++) {
		for (int j = 0; j < (int)m->size2; j++) {
			gsl_matrix_set(double_m, i, j, (double)MATRIX_ID(get)(m, i, j));
		}
	}
	return double_m;
}

/*
 * Converts a matrix to float precision.
 */
gsl_matrix_long_double* matlab::to_matrix_long_double(const MATRIX_TYPE* m) {
	gsl_matrix_long_double* long_double_m = gsl_matrix_long_double_alloc(m->size1, m->size2);
	for (int i = 0; i < (int)m->size1; i++) {
		for (int j = 0; j < (int)m->size2; j++) {
			gsl_matrix_long_double_set(long_double_m, i, j, (long double)MATRIX_ID(get)(m, i, j));
		}
	}
	return long_double_m;
}

/*
 * Converts a vector to a single-column matrix.
 */
MATRIX_TYPE* matlab::to_column_matrix(const VECTOR_TYPE* v) {
	MATRIX_TYPE* m = MATRIX_ID(alloc)(v->size, 1);
	for (int i = 0; i < (int)v->size; i++) {
		MATRIX_ID(set)(m, i, 0, VECTOR_ID(get)(v, i));
	}
	return m;
}

/*
 * Converts a vector to a single-row matrix.
 */
MATRIX_TYPE* matlab::to_row_matrix(const VECTOR_TYPE* v) {
	MATRIX_TYPE* m = MATRIX_ID(alloc)(1, v->size);
	for (int i = 0; i < (int)v->size; i++) {
		MATRIX_ID(set)(m, 0, i, VECTOR_ID(get)(v, i));
	}
	return m;
}

/*
 * Converts a vector to a permutation.
 */
gsl_permutation* matlab::to_permutation(const VECTOR_TYPE* v) {
	gsl_permutation* p = gsl_permutation_alloc(v->size);
	for (int i = 0; i < (int)v->size; i++) {
		p->data[i] = (int)VECTOR_ID(get)(v, i);
	}
	if (gsl_permutation_valid(p) == 1) {
		gsl_permutation_free(p);
		return NULL;
	} else {
		return p;
	}
}
