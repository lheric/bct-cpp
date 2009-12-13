#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "matlab.h"

/*
 * Our indexing-with-assignment functions do not automatically resize vectors
 * and matrices as MATLAB's do.  Thus, while the MATLAB code
 * (A = [1 2 3]; A([4 5]) = 5) yields (A = [1 2 3 5 5]), the equivalent C++ code
 * would result in an error.  In general, it is the caller's responsibility to
 * ensure that a call to one of these functions does not attempt to index past
 * the end of a vector or matrix.
 */

/*
 * Emulates (v[indices]) using ordinal indexing.
 */
gsl_vector* matlab::ordinal_index(const gsl_vector* v, const gsl_vector* indices) {
	gsl_vector* indexed = gsl_vector_alloc(indices->size);
	for (int i = 0; i < indices->size; i++) {
		int index = (int)gsl_vector_get(indices, i);
		double value = gsl_vector_get(v, index);
		gsl_vector_set(indexed, i, value);
	}
	return indexed;
}

/*
 * Emulates (v[indices] = value) using ordinal indexing.
 */
void matlab::ordinal_index_assign(gsl_vector* v, const gsl_vector* indices, double value) {
	for (int i = 0; i < indices->size; i++) {
		int index = (int)gsl_vector_get(indices, i);
		gsl_vector_set(v, index, value);
	}
}

/*
 * Emulates (v[indices] = values) using ordinal indexing.
 */
void matlab::ordinal_index_assign(gsl_vector* v, const gsl_vector* indices, const gsl_vector* values) {
	for (int i = 0; i < indices->size; i++) {
		int index = (int)gsl_vector_get(indices, i);
		double value = gsl_vector_get(values, i);
		gsl_vector_set(v, index, value);
	}
}

/*
 * Emulates (v[log_v]) using logical indexing.
 */
gsl_vector* matlab::logical_index(const gsl_vector* v, const gsl_vector* log_v) {
	gsl_vector* indexed = gsl_vector_alloc(nnz(log_v));
	for (int i = 0, index = 0; i < log_v->size; i++) {
		if (fp_nonzero(gsl_vector_get(log_v, i))) {
			double value = gsl_vector_get(v, i);
			gsl_vector_set(indexed, index++, value);
		}
	}
	return indexed;
}

/*
 * Emulates (v[log_v] = value) using logical indexing.
 */
void matlab::logical_index_assign(gsl_vector* v, const gsl_vector* log_v, double value) {
	for (int i = 0; i < log_v->size; i++) {
		if (fp_nonzero(gsl_vector_get(log_v, i))) {
			gsl_vector_set(v, i, value);
		}
	}
}

/*
 * Emulates (v[log_v] = values) using logical indexing.
 */
void matlab::logical_index_assign(gsl_vector* v, const gsl_vector* log_v, const gsl_vector* values) {
	for (int i = 0, index = 0; i < log_v->size; i++) {
		if (fp_nonzero(gsl_vector_get(log_v, i))) {
			double value = gsl_vector_get(values, index++);
			gsl_vector_set(v, i, value);
		}
	}
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
		int row = (int)gsl_vector_get(rows, i);
		for (int j = 0; j < columns->size; j++) {
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
 * Emulates matrix indexing by vectors with row and column indices, followed by assignment.
 * m([r], [c])=x
 */
void matlab::index_assign(gsl_matrix* m, const gsl_vector* row_indices, const gsl_vector* col_indices, double x) {
	for (int i = 0;i < row_indices->size; i++) {
		int row = gsl_vector_get(row_indices, i);
		for(int j = 0;j < col_indices->size;j++) {
			int col = gsl_vector_get(col_indices, j);
			int index = col*m->size1 + row;
			index_assign(m, index, x);
		}
	}
}

/*
 * Emulates matrix indexing by vectors with row and column indices, followed by assignment.
 * m([r], [c])=m2
 */
void matlab::index_assign(gsl_matrix* m, const gsl_vector* row_indices, const gsl_vector* col_indices, const gsl_matrix* source) {
	for (int i = 0, source_i = 0;i < row_indices->size && source_i < source->size1; i++, source_i++) {
		int row = gsl_vector_get(row_indices, i);
		for(int j = 0, source_j = 0;j < col_indices->size && source_j < source->size2; j++, source_j++) {
			int col = gsl_vector_get(col_indices, j);
			int index = col*m->size1 + row;
			double source_val = gsl_matrix_get(source, source_i, source_j);
			index_assign(m, index, source_val);
		}
	}
}

/*
 * Emulates matrix indexing and assignment (m(index) = x).
 */
void matlab::index_assign(gsl_matrix* m, int index, double x) {
	int row = index % (int)m->size1;
	int column = index / (int)m->size1;
	if (row < m->size1 && column < m->size2) {
		gsl_matrix_set(m, row, column, x);
	}
}

/*
 * Emulates matrix indexing and assignment (m(indices) = x).
 */
void matlab::index_assign(gsl_matrix* m, const gsl_vector* indices, double x) {
	for (int i = 0; i < indices->size; i++) {
		int index = (int)gsl_vector_get(indices, i);
		index_assign(m, index, x);
	}
}

/*
 * Emulates matrix indexing and assignment (m(indices) = x).
 */
void matlab::index_assign(gsl_matrix* m, const gsl_matrix* indices, double x) {
	for (int i = 0; i < indices->size1; i++) {
		for (int j = 0; j < indices->size2; j++) {
			int index = (int)gsl_matrix_get(indices, i, j);
			index_assign(m, index, x);
		}
	}
}

/*
 * Emulates m(r,c) where r is a vector of actual row indices and c is a vector of 1's and 0's
 */
gsl_matrix* matlab::mixed_logical_index(const gsl_matrix* m, const gsl_vector* rows, const gsl_vector* logical_cols) {
	int columns = nnz(logical_cols);
	if (columns == 0 || rows == NULL || m == NULL) {
		return NULL;
	}
	gsl_matrix* indexed = gsl_matrix_alloc(rows->size, columns);
	int column = 0;
	for (int j = 0; j < logical_cols->size; j++) {
		if (fp_nonzero(gsl_vector_get(logical_cols, j))) {
			for (int i = 0; i < rows->size; i++) {
				int row = (int)gsl_vector_get(rows, i);
				double value = gsl_matrix_get(m, row, j);
				gsl_matrix_set(indexed, i, column, value);
			}
			column += 1;
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
 * Emulates matrix logical indexing and assignment (m(lv) = x).
 */
void matlab::logical_index_assign(gsl_matrix* m, const gsl_vector* lv, double x) {
	for (int i = 0; i < lv->size; i++) {
		if (fp_nonzero(gsl_vector_get(lv, i))) {
			index_assign(m, i, x);
		}
	}
}

/*
 * Emulates matrix logical indexing and assignment (m(lm) = x).
 */
void matlab::logical_index_assign(gsl_matrix* m, const gsl_matrix* lm, double x) {
	for (int j = 0; j < lm->size2; j++) {
		for (int i = 0; i < lm->size1; i++) {
			if (fp_nonzero(gsl_matrix_get(lm, i, j))) {
				int index = j * lm->size1 + i;
				index_assign(m, index, x);
			}
		}
	}
}
