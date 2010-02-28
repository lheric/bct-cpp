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
