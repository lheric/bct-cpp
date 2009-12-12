#include <cmath>
#include <cstring>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "matlab.h"

/*
 * See MATLAB documentation for descriptions of these functions.  We will
 * document instances where our version differs from the MATLAB version.
 */
 
int matlab::all(const gsl_vector* v) {
	for (int i = 0; i < v->size; i++) {
		if (fp_zero(gsl_vector_get(v, i))) {
			return 0;
		}
	}
	return 1;
}

gsl_vector* matlab::all(const gsl_matrix* m, int dim) {
	if (dim == 1) {
		gsl_vector* all_v = gsl_vector_alloc(m->size2);
		for (int i = 0; i < m->size2; i++) {
			gsl_vector_const_view column = gsl_matrix_const_column(m, i);
			gsl_vector_set(all_v, i, all(&column.vector));
		}
		return all_v;
	} else if (dim == 2) {
		gsl_vector* all_v = gsl_vector_alloc(m->size1);
		for (int i = 0; i < m->size1; i++) {
			gsl_vector_const_view row = gsl_matrix_const_row(m, i);
			gsl_vector_set(all_v, i, all(&row.vector));
		}
		return all_v;
	} else {
		return NULL;
	}
}

int matlab::any(const gsl_vector* v) {
	for (int i = 0; i < v->size; i++) {
		if (fp_nonzero(gsl_vector_get(v, i))) {
			return 1;
		}
	}
	return 0;
}

gsl_vector* matlab::any(const gsl_matrix* m, int dim) {
	if (dim == 1) {
		gsl_vector* any_v = gsl_vector_alloc(m->size2);
		for (int i = 0; i < m->size2; i++) {
			gsl_vector_const_view column = gsl_matrix_const_column(m, i);
			gsl_vector_set(any_v, i, any(&column.vector));
		}
		return any_v;
	} else if (dim == 2) {
		gsl_vector* any_v = gsl_vector_alloc(m->size1);
		for (int i = 0; i < m->size1; i++) {
			gsl_vector_const_view row = gsl_matrix_const_row(m, i);
			gsl_vector_set(any_v, i, any(&row.vector));
		}
		return any_v;
	} else {
		return NULL;
	}
}

void matlab::dec2bin(int n, char* bin) {
	if (n <= 0) {
		bin[0] = '0';
		bin[1] = '\0';
		return;
	}
	int len = (int)(std::floor(1.0 + std::log(n) / std::log(2)));
	bin[len] = '\0';
	for (int i = len - 1; i >= 0; i--) {
		int remainder = n % 2;
		if (remainder) {
			bin[i] = '1';
		} else {
			bin[i] = '0';
		}
		n >>= 1;
	}
}

void matlab::dec2bin(int n, int len, char* bin) {
	dec2bin(n, bin);
	int strlen = std::strlen(bin);
	if (len > strlen) {
		for (int i = strlen, j = len; j >= 0; i--, j--) {
			if (i >= 0) {
				bin[j] = bin[i];
			} else {
				bin[j] = '0';
			}
		}
	}
}

gsl_matrix* matlab::eye(int size) {
	return eye(size, size);
}

gsl_matrix* matlab::eye(int size1, int size2) {
	gsl_matrix* m = gsl_matrix_calloc(size1, size2);
	gsl_vector_view diagonal = gsl_matrix_diagonal(m);
	gsl_vector_set_all(&diagonal.vector, 1.0);
	return m;
}

gsl_vector* matlab::find(const gsl_vector* v, int n, const char* direction) {
	int size = nnz(v);
	if (size == 0 || n < 1) {
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

/*
 * Emulates the two-return version of find.
 */
gsl_matrix* matlab::find_ij(const gsl_matrix* m, int n, const char* direction) {
	gsl_vector* found_v = find(m, n, direction);
	gsl_matrix* found_m;
	if (found_v == NULL) {
		found_m = NULL;
	} else {
		found_m = gsl_matrix_alloc(found_v->size, 2);
		for (int i = 0; i < found_v->size; i++) {
			int index = (int)gsl_vector_get(found_v, i);
			int row = index % (int)m->size1;
			int column = index / (int)m->size1;
			gsl_matrix_set(found_m, i, 0, (double)row);
			gsl_matrix_set(found_m, i, 1, (double)column);
		}
		gsl_vector_free(found_v);
	}
	return found_m;
}

gsl_vector* matlab::min(const gsl_matrix* m) {
	if(m == NULL) {
		return NULL;
	}
	gsl_vector* minv = gsl_vector_alloc(m->size2);
	for(int col = 0;col < m->size2;col++) {
		gsl_vector_const_view column = gsl_matrix_const_column(m, col);
		double minval = gsl_vector_min(&column.vector);
		gsl_vector_set(minv, col, minval);
	}
	return minv;
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

gsl_vector* matlab::nonzeros(const gsl_matrix* m) {
	gsl_vector* found = find(m);
	if (found != NULL) {
		for (int i = 0; i < found->size; i++) {
			int m_index = (int)gsl_vector_get(found, i);
			double value = index(m, m_index);
			gsl_vector_set(found, i, value);
		}
	}
	return found;
}

gsl_matrix* matlab::ones(int size) {
	return ones(size, size);
}

gsl_matrix* matlab::ones(int size1, int size2) {
	gsl_matrix* m = gsl_matrix_alloc(size1, size2);
	gsl_matrix_set_all(m, 1);
	return m;
}

/*
 * Emulates ones with a scale
 */
gsl_matrix* matlab::yens(int size, double n) {
	return yens(size, size, n);
}

gsl_matrix* matlab::yens(int size1, int size2, double n) {
	gsl_matrix* m = gsl_matrix_alloc(size1, size2);
	gsl_matrix_set_all(m, n);
	return m;
}

// TODO: Add the two-argument version
// TODO: Replace bubble sort with something from <algorithm>?
gsl_matrix* matlab::sortrows(const gsl_matrix* m) {
	gsl_matrix* sorted = copy(m);
	int swaps;
	do {
		swaps = 0;
		for (int i = 0; i < sorted->size1 - 1; i++) {
			gsl_vector_view row1 = gsl_matrix_row(sorted, i);
			gsl_vector_view row2 = gsl_matrix_row(sorted, i + 1);
			if (compare(&row1.vector, &row2.vector) == 1) {
				gsl_matrix_swap_rows(sorted, i, i + 1);
				swaps++;
			}
		}
	} while (swaps > 0);
	return sorted;
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
	gsl_matrix* tril = copy(m);
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
	gsl_matrix* triu = copy(m);
	for (int i = 0; i < m->size1; i++) {
		for (int j = i + k - 1; j >= 0; j--) {
			if (j < m->size2) {
				gsl_matrix_set(triu, i, j, 0.0);
			}
		}
	}
	return triu;
}

// TODO: Implement other "unique" variants?
/*
 * Emulates unique(m, "rows").
 */
gsl_matrix* matlab::unique_rows(const gsl_matrix* m, gsl_vector* i, gsl_vector* j) {
	gsl_matrix* unique_m = gsl_matrix_alloc(m->size1, m->size2);
	gsl_vector* temp_i = gsl_vector_alloc(m->size1);
	gsl_vector* temp_j = gsl_vector_alloc(m->size1);
	int uniqueSize = 0;
	for (int rowIndex = 0; rowIndex < m->size1; rowIndex++) {
		gsl_vector_const_view row = gsl_matrix_const_row(m, rowIndex);
		bool found = false;
		for (int uniqueIndex = 0; uniqueIndex < uniqueSize; uniqueIndex++) {
			gsl_vector_view uniqueRow = gsl_matrix_row(unique_m, uniqueIndex);
			if (compare(&row.vector, &uniqueRow.vector) == 0) {
				found = true;
				gsl_vector_set(temp_j, rowIndex, uniqueIndex);
				break;
			}
		}
		if (!found) {
			gsl_matrix_set_row(unique_m, uniqueSize, &row.vector);
			gsl_vector_set(temp_i, uniqueSize, rowIndex);
			gsl_vector_set(temp_j, rowIndex, uniqueSize);
			uniqueSize++;
		}
	}
	
	// TODO: This should be sorted before returning
	gsl_matrix* unique_m_sized = gsl_matrix_alloc(uniqueSize, m->size2);
	gsl_matrix_view unique_mv = gsl_matrix_submatrix(unique_m, 0, 0, uniqueSize, m->size2);
	gsl_matrix_memcpy(unique_m_sized, &unique_mv.matrix);
	gsl_matrix_free(unique_m);
	
	if (i != NULL) {
		gsl_vector_free(i);
		i = gsl_vector_alloc(uniqueSize);
		gsl_vector_view temp_iv = gsl_vector_subvector(temp_i, 0, uniqueSize);
		gsl_vector_memcpy(i, &temp_iv.vector);
	}
	gsl_vector_free(temp_i);
	
	if (j != NULL) {
		gsl_vector_free(j);
		j = temp_j;
	} else {
		gsl_vector_free(temp_j);
	}
	
	return unique_m_sized;
}

gsl_matrix* matlab::zeros(int size) {
	return gsl_matrix_calloc(size, size);
}

gsl_matrix* matlab::zeros(int size1, int size2) {
	return gsl_matrix_calloc(size1, size2);
}
