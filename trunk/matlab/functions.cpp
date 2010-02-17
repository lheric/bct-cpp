#include <cmath>
#include <cstring>
#include <ctime>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include "matlab.h"

/*
 * See MATLAB documentation for descriptions of these functions.  We will
 * document instances where our version differs from the MATLAB version.
 */
 
int matlab::all(const gsl_vector* v) {
	for (int i = 0; i < (int)v->size; i++) {
		if (fp_zero(gsl_vector_get(v, i))) {
			return 0;
		}
	}
	return 1;
}

gsl_vector* matlab::all(const gsl_matrix* m, int dim) {
	if (dim == 1) {
		gsl_vector* all_v = gsl_vector_alloc(m->size2);
		for (int i = 0; i < (int)m->size2; i++) {
			gsl_vector_const_view m_col_i = gsl_matrix_const_column(m, i);
			gsl_vector_set(all_v, i, all(&m_col_i.vector));
		}
		return all_v;
	} else if (dim == 2) {
		gsl_vector* all_v = gsl_vector_alloc(m->size1);
		for (int i = 0; i < (int)m->size1; i++) {
			gsl_vector_const_view m_row_i = gsl_matrix_const_row(m, i);
			gsl_vector_set(all_v, i, all(&m_row_i.vector));
		}
		return all_v;
	} else {
		return NULL;
	}
}

int matlab::any(const gsl_vector* v) {
	for (int i = 0; i < (int)v->size; i++) {
		if (fp_nonzero(gsl_vector_get(v, i))) {
			return 1;
		}
	}
	return 0;
}

gsl_vector* matlab::any(const gsl_matrix* m, int dim) {
	if (dim == 1) {
		gsl_vector* any_v = gsl_vector_alloc(m->size2);
		for (int i = 0; i < (int)m->size2; i++) {
			gsl_vector_const_view m_col_i = gsl_matrix_const_column(m, i);
			gsl_vector_set(any_v, i, any(&m_col_i.vector));
		}
		return any_v;
	} else if (dim == 2) {
		gsl_vector* any_v = gsl_vector_alloc(m->size1);
		for (int i = 0; i < (int)m->size1; i++) {
			gsl_vector_const_view m_row_i = gsl_matrix_const_row(m, i);
			gsl_vector_set(any_v, i, any(&m_row_i.vector));
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
	gsl_vector_view diag_eye_m = gsl_matrix_diagonal(eye_m);
	gsl_vector_set_all(&diag_eye_m.vector, 1.0);
	return eye_m;
}

gsl_vector* matlab::find(const gsl_vector* v, int n, const char* direction) {
	int n_find = nnz(v);
	if (n_find == 0 || n < 1) {
		return NULL;
	}
	gsl_vector* find_v = gsl_vector_alloc(n < n_find ? n : n_find);
	if (std::strcmp(direction, "first") == 0) {
		int position = 0;
		for (int i = 0; i < (int)v->size && position < (int)find_v->size; i++) {
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
		gsl_matrix* find_m = gsl_matrix_alloc(find_v->size, 2);
		for (int i = 0; i < (int)find_v->size; i++) {
			int index = (int)gsl_vector_get(find_v, i);
			int row = index % (int)m->size1;
			int column = index / (int)m->size1;
			gsl_matrix_set(find_m, i, 0, (double)row);
			gsl_matrix_set(find_m, i, 1, (double)column);
		}
		gsl_vector_free(find_v);
		return find_m;
	}
}

gsl_vector* matlab::hist(const gsl_vector* v, int n) {
	gsl_vector* centers = gsl_vector_alloc(n);
	double min = gsl_vector_min(v);
	double max = gsl_vector_max(v);
	double width = (max - min) / (double)n;
	for (int i = 0; i < n; i++) {
		gsl_vector_set(centers, i, min + (i + 0.5) * width);
	}
	gsl_vector* hist_v = hist(v, centers);
	gsl_vector_free(centers);
	return hist_v;
}

gsl_vector* matlab::hist(const gsl_vector* v, const gsl_vector* centers) {
	int n = centers->size;
	gsl_vector* hist_v = gsl_vector_calloc(n);
	for (int i = 0; i < (int)v->size; i++) {
		double value = gsl_vector_get(v, i);
		int index = n - 1;
		for (int j = 0; j < n - 1; j++) {
			double left = gsl_vector_get(centers, j);
			double right = gsl_vector_get(centers, j + 1);
			if (value < left) {
				index = j;
				break;
			} else if (value < right) {
				double middle = (left + right) / 2.0;
				if (fp_compare(value, middle) <= 0) {
					index = j;
				} else {
					index = j + 1;
				}
				break;
			}
		}
		gsl_vector_set(hist_v, index, gsl_vector_get(hist_v, index) + 1.0);
	}
	return hist_v;
}

int matlab::length(const gsl_vector* v) {
	return v->size;
}

int matlab::length(const gsl_matrix* m) {
	return m->size1 > m->size2 ? m->size1 : m->size2;
}

double matlab::max(double x, double y) {
	return x > y ? x : y;
}

double matlab::max(const gsl_vector* v) {
	return gsl_vector_max(v);
}

/*
 * Emulates (max(m)) or (max(m')).
 */
gsl_vector* matlab::max(const gsl_matrix* m, int dim) {
	if (dim == 1) {
		gsl_vector* max_v = gsl_vector_alloc(m->size2);
		for (int i = 0; i < (int)m->size2; i++) {
			gsl_vector_const_view m_col_i = gsl_matrix_const_column(m, i);
			double value = gsl_vector_max(&m_col_i.vector);
			gsl_vector_set(max_v, i, value);
		}
		return max_v;
	} else if (dim == 2) {
		gsl_vector* max_v = gsl_vector_alloc(m->size1);
		for (int i = 0; i < (int)m->size1; i++) {
			gsl_vector_const_view m_row_i = gsl_matrix_const_row(m, i);
			double value = gsl_vector_max(&m_row_i.vector);
			gsl_vector_set(max_v, i, value);
		}
		return max_v;
	} else {
		return NULL;
	}
}

double matlab::min(double x, double y) {
	return x < y ? x : y;
}

double matlab::min(const gsl_vector* v) {
	return gsl_vector_min(v);
}

/*
 * Emulates (min(m)) or (min(m')).
 */
gsl_vector* matlab::min(const gsl_matrix* m, int dim) {
	if (dim == 1) {
		gsl_vector* min_v = gsl_vector_alloc(m->size2);
		for (int i = 0; i < (int)m->size2; i++) {
			gsl_vector_const_view m_col_i = gsl_matrix_const_column(m, i);
			double value = gsl_vector_min(&m_col_i.vector);
			gsl_vector_set(min_v, i, value);
		}
		return min_v;
	} else if (dim == 2) {
		gsl_vector* min_v = gsl_vector_alloc(m->size1);
		for (int i = 0; i < (int)m->size1; i++) {
			gsl_vector_const_view m_row_i = gsl_matrix_const_row(m, i);
			double value = gsl_vector_min(&m_row_i.vector);
			gsl_vector_set(min_v, i, value);
		}
		return min_v;
	} else {
		return NULL;
	}
}

int matlab::nnz(const gsl_vector* v) {
	int nnz = 0;
	for (int i = 0; i < (int)v->size; i++) {
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
		for (int i = 0; i < (int)nz_v->size; i++) {
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

double matlab::prod(const gsl_vector* v) {
	double prod = 1.0;
	for (int i = 0; i < (int)v->size; i++) {
		prod *= gsl_vector_get(v, i);
	}
	return prod;
}

/*
 * Generates a permutation of the integers 0 to (size - 1), whereas the MATLAB
 * version uses the integers 1 to size.
 */
gsl_permutation* matlab::randperm(int size) {
	gsl_rng_default_seed = std::time(NULL);
	gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);
	double values[size];
	for (int i = 0; i < size; i++) {
		values[i] = (double)i;
	}
	gsl_ran_shuffle(rng, values, size, sizeof(double));
	gsl_vector_view values_vv = gsl_vector_view_array(values, size);
	gsl_permutation* values_p = to_permutation(&values_vv.vector);
	gsl_rng_free(rng);
	return values_p;
}

gsl_vector* matlab::prod(const gsl_matrix* m, int dim) {
	if (dim == 1) {
		gsl_vector* prod_v = gsl_vector_alloc(m->size2);
		gsl_vector_set_all(prod_v, 1.0);
		for (int i = 0; i < (int)m->size1; i++) {
			gsl_vector_const_view m_row_i = gsl_matrix_const_row(m, i);
			gsl_vector_mul(prod_v, &m_row_i.vector);
		}
		return prod_v;
	} else if (dim == 2) {
		gsl_vector* prod_v = gsl_vector_alloc(m->size1);
		gsl_vector_set_all(prod_v, 1.0);
		for (int i = 0; i < (int)m->size2; i++) {
			gsl_vector_const_view m_col_i = gsl_matrix_const_column(m, i);
			gsl_vector_mul(prod_v, &m_col_i.vector);
		}
		return prod_v;
	} else {
		return NULL;
	}
}

gsl_vector* matlab::reverse(const gsl_vector* v) {
	gsl_vector* rev_v = gsl_vector_alloc(v->size);
	for(int i = (v->size-1),j = 0;i >= 0;i--,j++) {
		gsl_vector_set(rev_v, j, gsl_vector_get(v, i));
	}
	return rev_v;
}

gsl_vector* matlab::setxor(const gsl_vector* v1, const gsl_vector* v2) {
	if (v1 == NULL && v2 == NULL) {
		return NULL;
	}
	if (v1 == NULL) {
		return unique(v2);
	}
	if (v2 == NULL) {
		return unique(v1);
	}
	gsl_vector* unique_v1 = unique(v1);
	gsl_vector* unique_v2 = unique(v2);
	gsl_vector* unsized_v = gsl_vector_alloc(v1->size + v2->size);
	int n = 0;
	for (int i = 0; i < (int)unique_v1->size; i++) {
		bool found = false;
		double v1_value = gsl_vector_get(unique_v1, i);
		for (int j = 0; j < (int)unique_v2->size; j++) {
			double v2_value = gsl_vector_get(unique_v2, j);
			if (fp_equal(v1_value, v2_value)) {
				found = true;
				break;
			}
		}
		if (!found) {
			gsl_vector_set(unsized_v, n++, v1_value);
		}
	}
	for (int i = 0; i < (int)unique_v2->size; i++) {
		bool found = false;
		double v2_value = gsl_vector_get(unique_v2, i);
		for (int j = 0; j < (int)unique_v1->size; j++) {
			double v1_value = gsl_vector_get(unique_v1, j);
			if (fp_equal(v2_value, v1_value)) {
				found = true;
				break;
			}
		}
		if (!found) {
			gsl_vector_set(unsized_v, n++, v2_value);
		}
	}
	gsl_vector_free(unique_v1);
	gsl_vector_free(unique_v2);
	if (n > 0) {
		gsl_vector* unsorted_v = gsl_vector_alloc(n);
		gsl_vector_view unsized_subv = gsl_vector_subvector(unsized_v, 0, n);
		gsl_vector_memcpy(unsorted_v, &unsized_subv.vector);
		gsl_vector_free(unsized_v);
		gsl_vector* setxor_v = sort(unsorted_v);
		gsl_vector_free(unsorted_v);
		return setxor_v;
	} else {
		gsl_vector_free(unsized_v);
		return NULL;
	}
}

gsl_vector* matlab::sort(const gsl_vector* v, const char* mode, gsl_vector** ind) {
	if (std::strcmp(mode, "ascend") != 0 && std::strcmp(mode, "descend") != 0) {
		return NULL;
	}
	double elements[v->size];
	to_array(v, elements);
	std::size_t indices[v->size];
	if (std::strcmp(mode, "ascend") == 0) {
		stable_sort_index(indices, elements, v->size);
	} else {
		stable_sort_index(indices, elements, v->size, fp_greater);
	}
	gsl_vector* sort_v = gsl_vector_alloc(v->size);
	if (ind != NULL) {
		*ind = gsl_vector_alloc(v->size);
	}
	for (int i = 0; i < (int)v->size; i++) {
		int index = indices[i];
		gsl_vector_set(sort_v, i, elements[index]);
		if (ind != NULL) {
			gsl_vector_set(*ind, i, (double)index);
		}
	}
	return sort_v;
}

gsl_matrix* matlab::sort(const gsl_matrix* m, int dim, const char* mode, gsl_matrix** ind) {
	if (std::strcmp(mode, "ascend") != 0 && std::strcmp(mode, "descend") != 0) {
		return NULL;
	}
	if (dim == 1) {
		gsl_matrix* sort_m = gsl_matrix_alloc(m->size1, m->size2);
		if (ind != NULL) {
			*ind = gsl_matrix_alloc(m->size1, m->size2);
		}
		for (int i = 0; i < (int)m->size2; i++) {
			gsl_vector_const_view m_col_i = gsl_matrix_const_column(m, i);
			gsl_vector* sort_m_col_i;
			if (ind == NULL) {
				sort_m_col_i = sort(&m_col_i.vector, mode);
			} else {
				gsl_vector* ind_col_i;
				sort_m_col_i = sort(&m_col_i.vector, mode, &ind_col_i);
				gsl_matrix_set_col(*ind, i, ind_col_i);
				gsl_vector_free(ind_col_i);
			}
			gsl_matrix_set_col(sort_m, i, sort_m_col_i);
			gsl_vector_free(sort_m_col_i);
		}
		return sort_m;
	} else if (dim == 2) {
		gsl_matrix* sort_m = gsl_matrix_alloc(m->size1, m->size2);
		if (ind != NULL) {
			*ind = gsl_matrix_alloc(m->size1, m->size2);
		}
		for (int i = 0; i < (int)m->size1; i++) {
			gsl_vector_const_view m_row_i = gsl_matrix_const_row(m, i);
			gsl_vector* sort_m_row_i;
			if (ind == NULL) {
				sort_m_row_i = sort(&m_row_i.vector, mode);
			} else {
				gsl_vector* ind_row_i;
				sort_m_row_i = sort(&m_row_i.vector, mode, &ind_row_i);
				gsl_matrix_set_row(*ind, i, ind_row_i);
				gsl_vector_free(ind_row_i);
			}
			gsl_matrix_set_row(sort_m, i, sort_m_row_i);
			gsl_vector_free(sort_m_row_i);
		}
		return sort_m;
	} else {
		return NULL;
	}
}

/*
 * Emulates (sortrows(v)) for a column vector.
 */
gsl_vector* matlab::sortrows(const gsl_vector* v, gsl_vector** ind) {
	return sort(v, "ascend", ind);
}

gsl_matrix* matlab::sortrows(const gsl_matrix* m, gsl_vector** ind) {
	gsl_vector* rows[m->size1];
	for (int i = 0; i < (int)m->size1; i++) {
		rows[i] = gsl_vector_alloc(m->size2);
		gsl_matrix_get_row(rows[i], m, i);
	}
	std::size_t indices[m->size1];
	stable_sort_index(indices, rows, m->size1, vector_less);
	for (int i = 0; i < (int)m->size1; i++) {
		gsl_vector_free(rows[i]);
	}
	gsl_matrix* sort_m = gsl_matrix_alloc(m->size1, m->size2);
	if (ind != NULL) {
		*ind = gsl_vector_alloc(m->size1);
	}
	for (int i = 0; i < (int)m->size1; i++) {
		int index = indices[i];
		gsl_vector_const_view m_row_index = gsl_matrix_const_row(m, index);
		gsl_matrix_set_row(sort_m, i, &m_row_index.vector);
		if (ind != NULL) {
			gsl_vector_set(*ind, i, (double)index);
		}
	}
	return sort_m;
}

double matlab::sum(const gsl_vector* v) {
	double sum = 0.0;
	for (int i = 0; i < (int)v->size; i++) {
		sum += gsl_vector_get(v, i);
	}
	return sum;
}

gsl_vector* matlab::sum(const gsl_matrix* m, int dim) {
	if (dim == 1) {
		gsl_vector* sum_v = gsl_vector_calloc(m->size2);
		for (int i = 0; i < (int)m->size1; i++) {
			gsl_vector_const_view m_row_i = gsl_matrix_const_row(m, i);
			gsl_vector_add(sum_v, &m_row_i.vector);
		}
		return sum_v;
	} else if (dim == 2) {
		gsl_vector* sum_v = gsl_vector_calloc(m->size1);
		for (int i = 0; i < (int)m->size2; i++) {
			gsl_vector_const_view m_col_i = gsl_matrix_const_column(m, i);
			gsl_vector_add(sum_v, &m_col_i.vector);
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
	for (int i = 0; i < (int)m->size1; i++) {
		for (int j = i + k + 1; j < (int)m->size2; j++) {
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
	for (int i = 0; i < (int)m->size1; i++) {
		for (int j = i + k - 1; j >= 0; j--) {
			if (j < (int)m->size2) {
				gsl_matrix_set(triu_m, i, j, 0.0);
			}
		}
	}
	return triu_m;
}

gsl_vector* matlab::unique(const gsl_vector* v, const char* first_or_last, gsl_vector** i, gsl_vector** j) {
	if (std::strcmp(first_or_last, "first") != 0 && std::strcmp(first_or_last, "last") != 0) {
		return NULL;
	}
	gsl_vector* sort_v = sort(v);
	gsl_vector* unsized_v = gsl_vector_alloc(v->size);
	gsl_vector_set(unsized_v, 0, gsl_vector_get(sort_v, 0));
	int n = 1;
	for (int x = 1; x < (int)v->size; x++) {
		double prev_value = gsl_vector_get(sort_v, x - 1);
		double value = gsl_vector_get(sort_v, x);
		if (fp_compare(prev_value, value) != 0) {
			gsl_vector_set(unsized_v, n++, value);
		}
	}
	gsl_vector_free(sort_v);
	gsl_vector* unique_v = gsl_vector_alloc(n);
	gsl_vector_view unsized_subv = gsl_vector_subvector(unsized_v, 0, n);
	gsl_vector_memcpy(unique_v, &unsized_subv.vector);
	gsl_vector_free(unsized_v);
	if (i != NULL) {
		*i = gsl_vector_alloc(n);
		for (int x = 0; x < n; x++) {
			for (int y = 0; y < (int)v->size; y++) {
				if (fp_compare(gsl_vector_get(unique_v, x), gsl_vector_get(v, y)) == 0) {
					gsl_vector_set(*i, x, y);
					if (std::strcmp(first_or_last, "first") == 0) {
						break;
					}
				}
			}
		}
	}
	if (j != NULL) {
		*j = gsl_vector_alloc(v->size);
		for (int x = 0; x < (int)v->size; x++) {
			for (int y = 0; y < n; y++) {
				if (fp_compare(gsl_vector_get(v, x), gsl_vector_get(unique_v, y)) == 0) {
					gsl_vector_set(*j, x, y);
					break;
				}
			}
		}
	}
	return unique_v;
}

gsl_vector* matlab::unique(const gsl_matrix* m, const char* first_or_last, gsl_vector** i, gsl_vector** j) {
	gsl_vector* v = to_vector(m);
	gsl_vector* unique_v = unique(v, first_or_last, i, j);
	gsl_vector_free(v);
	return unique_v;
}

/*
 * Emulates (unique(m, "rows", first_or_last)).
 */
gsl_matrix* matlab::unique_rows(const gsl_matrix* m, const char* first_or_last, gsl_vector** i, gsl_vector** j) {
	if (std::strcmp(first_or_last, "first") != 0 && std::strcmp(first_or_last, "last") != 0) {
		return NULL;
	}
	gsl_matrix* sort_m = sortrows(m);
	gsl_matrix* unsized_m = gsl_matrix_alloc(m->size1, m->size2);
	gsl_vector_view first_row = gsl_matrix_row(sort_m, 0);
	gsl_matrix_set_row(unsized_m, 0, &first_row.vector);
	int n_unique = 1;
	for (int x = 1; x < (int)m->size1; x++) {
		gsl_vector_view sort_m_row_x_sub_1 = gsl_matrix_row(sort_m, x - 1);
		gsl_vector_view sort_m_row_x = gsl_matrix_row(sort_m, x);
		if (compare_vectors(&sort_m_row_x_sub_1.vector, &sort_m_row_x.vector) != 0) {
			gsl_matrix_set_row(unsized_m, n_unique++, &sort_m_row_x.vector);
		}
	}
	gsl_matrix_free(sort_m);
	gsl_matrix* unique_m = gsl_matrix_alloc(n_unique, m->size2);
	gsl_matrix_view unsized_subm = gsl_matrix_submatrix(unsized_m, 0, 0, n_unique, m->size2);
	gsl_matrix_memcpy(unique_m, &unsized_subm.matrix);
	gsl_matrix_free(unsized_m);
	if (i != NULL) {
		*i = gsl_vector_alloc(n_unique);
		for (int x = 0; x < n_unique; x++) {
			gsl_vector_view unique_m_row_x = gsl_matrix_row(unique_m, x);
			for (int y = 0; y < (int)m->size1; y++) {
				gsl_vector_const_view m_row_y = gsl_matrix_const_row(m, y);
				if (compare_vectors(&unique_m_row_x.vector, &m_row_y.vector) == 0) {
					gsl_vector_set(*i, x, y);
					if (std::strcmp(first_or_last, "first") == 0) {
						break;
					}
				}
			}
		}
	}
	if (j != NULL) {
		*j = gsl_vector_alloc(m->size1);
		for (int x = 0; x < (int)m->size1; x++) {
			gsl_vector_const_view m_row_x = gsl_matrix_const_row(m, x);
			for (int y = 0; y < n_unique; y++) {
				gsl_vector_view unique_m_row_y = gsl_matrix_row(unique_m, y);
				if (compare_vectors(&m_row_x.vector, &unique_m_row_y.vector) == 0) {
					gsl_vector_set(*j, x, y);
					break;
				}
			}
		}
	}
	return unique_m;
}

gsl_matrix* matlab::zeros(int size) {
	return gsl_matrix_calloc(size, size);
}

gsl_matrix* matlab::zeros(int size1, int size2) {
	return gsl_matrix_calloc(size1, size2);
}

/*
 * Emulates (zeros(size, 1)) or (zeros(1, size)).
 */
gsl_vector* matlab::zeros_vector(int size) {
	return gsl_vector_calloc(size);
}
