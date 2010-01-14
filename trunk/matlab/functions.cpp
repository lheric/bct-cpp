#include <cmath>
#include <cstring>
#include <gsl/gsl_heapsort.h>
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
	int strlen = (int)(std::floor(1.0 + std::log(n) / std::log(2)));
	bin[strlen] = '\0';
	for (int i = strlen - 1; i >= 0; i--) {
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
	gsl_matrix* eye_m = gsl_matrix_calloc(size1, size2);
	gsl_vector_view diagonal = gsl_matrix_diagonal(eye_m);
	gsl_vector_set_all(&diagonal.vector, 1.0);
	return eye_m;
}

gsl_vector* matlab::find(const gsl_vector* v, int n, const char* direction) {
	int n_find = nnz(v);
	if (n_find == 0 || n < 1) {
		return NULL;
	}
	gsl_vector* find_v = gsl_vector_alloc((n < n_find) ? n : n_find);
	if (direction == NULL || std::strcmp(direction, "first") == 0) {
		int position = 0;
		for (int i = 0; i < v->size && position < find_v->size; i++) {
			if (fp_nonzero(gsl_vector_get(v, i))) {
				gsl_vector_set(find_v, position, i);
				position++;
			}
		}
		return find_v;
	} else if (std::strcmp(direction, "last") == 0) {
		int position = find_v->size - 1;
		for (int i = v->size - 1; i >= 0 && position >= 0; i--) {
			if (fp_nonzero(gsl_vector_get(v, i))) {
				gsl_vector_set(find_v, position, i);
				position--;
			}
		}
		return find_v;
	} else {
		gsl_vector_free(find_v);
		return NULL;
	}
}

gsl_vector* matlab::find(const gsl_matrix* m, int n, const char* direction) {
	gsl_vector* v = to_vector(m);
	gsl_vector* find_v = find(v, n, direction);
	gsl_vector_free(v);
	return find_v;
}

/*
 * Emulates the two-return version of "find".
 */
gsl_matrix* matlab::find_ij(const gsl_matrix* m, int n, const char* direction) {
	gsl_vector* find_v = find(m, n, direction);
	if (find_v == NULL) {
		return NULL;
	} else {
		gsl_matrix* found_m = gsl_matrix_alloc(find_v->size, 2);
		for (int i = 0; i < find_v->size; i++) {
			int index = (int)gsl_vector_get(find_v, i);
			int row = index % (int)m->size1;
			int column = index / (int)m->size1;
			gsl_matrix_set(found_m, i, 0, (double)row);
			gsl_matrix_set(found_m, i, 1, (double)column);
		}
		gsl_vector_free(find_v);
		return found_m;
	}
}

// TODO: Implement other "max" and "min" variants?

double matlab::max(const gsl_vector* v) {
	return gsl_vector_max(v);
}

gsl_vector* matlab::max(const gsl_matrix* m) {
	gsl_vector* max_v = gsl_vector_alloc(m->size2);
	for (int i = 0; i < m->size2; i++) {
		gsl_vector_const_view column = gsl_matrix_const_column(m, i);
		double value = gsl_vector_max(&column.vector);
		gsl_vector_set(max_v, i, value);
	}
	return max_v;
}

gsl_vector* matlab::max(const gsl_matrix* m, int dim) { //no default value for 'dim'
	if(dim == 1) {
		gsl_vector* v = gsl_vector_alloc(m->size1);
		for(int i = 0;i < m->size2;i++) {
			gsl_vector_const_view column = gsl_matrix_const_column(m, i);
			double col_max = gsl_vector_max(&column.vector);
			gsl_vector_set(v, i, col_max);
		}
		return v;
	}
	if(dim == 2) {
		gsl_vector* v = gsl_vector_alloc(m->size1);
		for(int i = 0;i < m->size1;i++) {
			gsl_vector_const_view row = gsl_matrix_const_row(m, i);
			double row_max = gsl_vector_max(&row.vector);
			gsl_vector_set(v, i, row_max);
		}
		return v;
	}
}

double matlab::min(const gsl_vector* v) {
	return gsl_vector_min(v);
}

gsl_vector* matlab::min(const gsl_matrix* m) {
	gsl_vector* min_v = gsl_vector_alloc(m->size2);
	for (int i = 0; i < m->size2; i++) {
		gsl_vector_const_view column = gsl_matrix_const_column(m, i);
		double value = gsl_vector_min(&column.vector);
		gsl_vector_set(min_v, i, value);
	}
	return min_v;
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
	int nnz_v = nnz(v);
	gsl_vector_free(v);
	return nnz_v;
}

gsl_vector* matlab::nonzeros(const gsl_matrix* m) {
	gsl_vector* nz_v = find(m);
	if (nz_v != NULL) {
		for (int i = 0; i < nz_v->size; i++) {
			int i_m = (int)gsl_vector_get(nz_v, i);
			double value = ordinal_index(m, i_m);
			gsl_vector_set(nz_v, i, value);
		}
	}
	return nz_v;
}

gsl_matrix* matlab::ones(int size) {
	return ones(size, size);
}

gsl_matrix* matlab::ones(int size1, int size2) {
	gsl_matrix* ones_m = gsl_matrix_alloc(size1, size2);
	gsl_matrix_set_all(ones_m, 1);
	return ones_m;
}

gsl_vector* matlab::reverse(gsl_vector* v) {
	gsl_vector* rev_v = gsl_vector_alloc(v->size);
	for(int i = (v->size-1),j = 0;i >= 0;i--,j++) {
		gsl_vector_set(rev_v, j, gsl_vector_get(v, i));
	}
	return rev_v;
}

gsl_matrix* matlab::sortrows(const gsl_matrix* m, gsl_vector* ind) {
	gsl_vector* rows[m->size1];
	for (int i = 0; i < m->size1; i++) {
		rows[i] = gsl_vector_alloc(m->size2);
		gsl_matrix_get_row(rows[i], m, i);
	}
	size_t indices[m->size1];
	quicksort_index(indices, rows, m->size1, (cmp_fn)compare_vectorps);
	for (int i = 0; i < m->size1; i++) {
		gsl_vector_free(rows[i]);
	}
	gsl_matrix* sorted_m = gsl_matrix_alloc(m->size1, m->size2);
	for (int i = 0; i < m->size1; i++) {
		int index = indices[i];
		gsl_vector_const_view row = gsl_matrix_const_row(m, index);
		gsl_matrix_set_row(sorted_m, i, &row.vector);
		if (ind != NULL) {
			gsl_vector_set(ind, i, index);
		}
	}
	return sorted_m;
}

double matlab::sum(const gsl_vector* v) {
	double total = 0.0;
	for (int i = 0; i < v->size; i++) {
		total += gsl_vector_get(v, i);
	}
	return total;
}

gsl_vector* matlab::sum(const gsl_matrix* m, int dim) {
	if (dim == 1) {
		gsl_vector* sum_v = gsl_vector_calloc(m->size2);
		for (int i = 0; i < m->size1; i++) {
			gsl_vector_const_view row = gsl_matrix_const_row(m, i);
			gsl_vector_add(sum_v, &row.vector);
		}
		return sum_v;
	} else if (dim == 2) {
		gsl_vector* sum_v = gsl_vector_calloc(m->size1);
		for (int i = 0; i < m->size2; i++) {
			gsl_vector_const_view column = gsl_matrix_const_column(m, i);
			gsl_vector_add(sum_v, &column.vector);
		}
		return sum_v;
	} else {
		return NULL;
	}
}

gsl_matrix* matlab::tril(const gsl_matrix* m, int k) {
	if (k <= -(int)m->size1 || k >= (int)m->size2) {
		return NULL;
	}
	gsl_matrix* tril_m = copy(m);
	for (int i = 0; i < m->size1; i++) {
		for (int j = i + k + 1; j < m->size2; j++) {
			if (j >= 0) {
				gsl_matrix_set(tril_m, i, j, 0.0);
			}
		}
	}
	return tril_m;
}

gsl_matrix* matlab::triu(const gsl_matrix* m, int k) {
	if (k <= -(int)m->size1 || k >= (int)m->size2) {
		return NULL;
	}
	gsl_matrix* triu_m = copy(m);
	for (int i = 0; i < m->size1; i++) {
		for (int j = i + k - 1; j >= 0; j--) {
			if (j < m->size2) {
				gsl_matrix_set(triu_m, i, j, 0.0);
			}
		}
	}
	return triu_m;
}

/*
 * Emulates (unique(m, "rows")).
 */
gsl_matrix* matlab::unique_rows(const gsl_matrix* m, gsl_vector* i, gsl_vector* j) {
	return unique_rows(m, "last", i, j);
}

/*
 * Emulates (unique(m, "rows", "first")) or (unique(m, "rows", "last")).
 */
gsl_matrix* matlab::unique_rows(const gsl_matrix* m, const char* first_or_last, gsl_vector* i, gsl_vector* j) {
	if (std::strcmp(first_or_last, "first") != 0 && std::strcmp(first_or_last, "last") != 0) {
		return NULL;
	}
	
	// Find unique rows
	gsl_matrix* unique_m = gsl_matrix_alloc(m->size1, m->size2);
	int n_unique = 0;
	for (int i_m = 0; i_m < m->size1; i_m++) {
		gsl_vector_const_view row_m = gsl_matrix_const_row(m, i_m);
		bool found = false;
		for (int i_unique = 0; i_unique < n_unique; i_unique++) {
			gsl_vector_view row_unique = gsl_matrix_row(unique_m, i_unique);
			if (compare_vectors(&row_m.vector, &row_unique.vector) == 0) {
				found = true;
				break;
			}
		}
		if (!found) {
			gsl_matrix_set_row(unique_m, n_unique, &row_m.vector);
			n_unique++;
		}
	}
	
	// Resize and sort unique matrix
	gsl_matrix* unique_m_sized = gsl_matrix_alloc(n_unique, unique_m->size2);
	gsl_matrix_view unique_subm = gsl_matrix_submatrix(unique_m, 0, 0, n_unique, unique_m->size2);
	gsl_matrix_memcpy(unique_m_sized, &unique_subm.matrix);
	gsl_matrix_free(unique_m);
	gsl_matrix* unique_m_sorted = sortrows(unique_m_sized);
	gsl_matrix_free(unique_m_sized);
	
	// Calculate additional returns if necessary
	if (i != NULL) {
		
		// TODO: Can't do this, but how does the caller know the right size?
		if (i->size != n_unique) {
			gsl_vector_free(i);
			i = gsl_vector_alloc(n_unique);
		}
		for (int i_unique = 0; i_unique < n_unique; i_unique++) {
			gsl_vector_view row_unique = gsl_matrix_row(unique_m_sorted, i_unique);
			for (int i_m = 0; i_m < m->size1; i_m++) {
				gsl_vector_const_view row_m = gsl_matrix_const_row(m, i_m);
				if (compare_vectors(&row_unique.vector, &row_m.vector) == 0) {
					gsl_vector_set(i, i_unique, i_m);
					if (std::strcmp(first_or_last, "first") == 0) {
						break;
					}
				}
			}
		}
	}
	if (j != NULL) {
		for (int i_m = 0; i_m < m->size1; i_m++) {
			gsl_vector_const_view row_m = gsl_matrix_const_row(m, i_m);
			for (int i_unique = 0; i_unique < n_unique; i_unique++) {
				gsl_vector_view row_unique = gsl_matrix_row(unique_m_sorted, i_unique);
				if (compare_vectors(&row_m.vector, &row_unique.vector) == 0) {
					gsl_vector_set(j, i_m, i_unique);
					break;
				}
			}
		}
	}
	
	return unique_m_sorted;
}

gsl_matrix* matlab::zeros(int size) {
	return gsl_matrix_calloc(size, size);
}

gsl_matrix* matlab::zeros(int size1, int size2) {
	return gsl_matrix_calloc(size1, size2);
}
