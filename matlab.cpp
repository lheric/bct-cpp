#include <cmath>
#include <cstdio>
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

/* 
 * Emulates 'scalar:scalar'
 * Ex: >>2:5 in matlab yields [2,3,4,5]
 */
gsl_vector* matlab::sequence(int start, int end) {
	if(end < start) {
		return NULL;
	}
	gsl_vector* v = gsl_vector_alloc(end-start+1);
	for(int i=0,val=start;val <= end;i++,val++) {
		gsl_vector_set(v, i, val);
	}
	return v;
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

// TODO: Single-argument zeros should return a vector?
gsl_matrix* matlab::zeros(int size) {
	return gsl_matrix_calloc(size, size);
}

gsl_matrix* matlab::zeros(int size1, int size2) {
	return gsl_matrix_calloc(size1, size2);
}

/*
 * Emulates ([v1 x])
 */
gsl_vector* matlab::concatenate(const gsl_vector* v, const double x) {
	if(v == NULL) {
		return NULL;
	}
	gsl_vector* cat_v = gsl_vector_alloc(v->size + 1);
	gsl_vector_view sub_v = gsl_vector_subvector(cat_v, 0, v->size);
	gsl_vector_memcpy(&sub_v.vector, v);
	gsl_vector_set(cat_v, v->size, x);
	return cat_v;
}
 
/*
 * Emulates ([v1 v2]) for row vectors or ([v1 ; v2]) for column vectors.
 */
gsl_vector* matlab::concatenate(const gsl_vector* v1, const gsl_vector* v2) {
	if (v1 == NULL && v2 == NULL) {
		return NULL;
	} else if (v1 == NULL) {
		return copy(v2);
	} else if (v2 == NULL) {
		return copy(v1);
	}
	gsl_vector* cat_v = gsl_vector_alloc(v1->size + v2->size);
	gsl_vector_view cat_v1 = gsl_vector_subvector(cat_v, 0, v1->size);
	gsl_vector_view cat_v2 = gsl_vector_subvector(cat_v, v1->size, v2->size);
	gsl_vector_memcpy(&cat_v1.vector, v1);
	gsl_vector_memcpy(&cat_v2.vector, v2);
	return cat_v;
}

/*
 * Emulates ([m1 ; m2]).
 */
gsl_matrix* matlab::concatenate_columns(const gsl_matrix* m1, const gsl_matrix* m2) {
	if (m1 == NULL && m2 == NULL) {
		return NULL;
	} else if (m1 == NULL) {
		return copy(m2);
	} else if (m2 == NULL) {
		return copy(m1);
	} else if (m1->size2 != m2->size2) {
		return NULL;
	}
	gsl_matrix* cat_m = gsl_matrix_alloc(m1->size1 + m2->size1, m1->size2);
	gsl_matrix_view cat_m1 = gsl_matrix_submatrix(cat_m, 0, 0, m1->size1, m1->size2);
	gsl_matrix_view cat_m2 = gsl_matrix_submatrix(cat_m, m1->size1, 0, m2->size1, m2->size2);
	gsl_matrix_memcpy(&cat_m1.matrix, m1);
	gsl_matrix_memcpy(&cat_m2.matrix, m2);
	return cat_m;
}

/*
 * Emulates ([m ; v]).
 */
gsl_matrix* matlab::concatenate_columns(const gsl_matrix* m, const gsl_vector* v) {
	if (m == NULL && v == NULL) {
		return NULL;
	} else if (m == NULL) {
		return to_row_matrix(v);
	} else if (v == NULL) {
		return copy(m);
	} else if (m->size2 != v->size) {
		return NULL;
	}
	gsl_matrix* cat_mv = gsl_matrix_alloc(m->size1 + 1, m->size2);
	gsl_matrix_view cat_m = gsl_matrix_submatrix(cat_mv, 0, 0, m->size1, m->size2);
	gsl_vector_view cat_v = gsl_matrix_row(cat_mv, cat_mv->size1-1);
	gsl_matrix_memcpy(&cat_m.matrix, m);
	gsl_vector_memcpy(&cat_v.vector, v);
	return cat_mv;
}

/*
 * Emulates ([v ; m]).
 */
gsl_matrix* matlab::concatenate_columns(const gsl_vector* v, const gsl_matrix* m) {
	if (m == NULL && v == NULL) {
		return NULL;
	} else if (m == NULL) {
		return to_row_matrix(v);
	} else if (v == NULL) {
		return copy(m);
	} else if (m->size2 != v->size) {
		return NULL;
	}
	gsl_matrix* cat_mv = gsl_matrix_alloc(m->size1 + 1, m->size2);
	gsl_vector_view cat_v = gsl_matrix_row(cat_mv, 0);
	gsl_matrix_view cat_m = gsl_matrix_submatrix(cat_mv, 1, 0, m->size1, m->size2);
	gsl_vector_memcpy(&cat_v.vector, v);
	gsl_matrix_memcpy(&cat_m.matrix, m);
	return cat_mv;
}

/*
 * Emulates ([v1 ; v2]) for row vectors
 */
gsl_matrix* matlab::concatenate_columns(const gsl_vector* v1, const gsl_vector* v2) {
	if (v1 == NULL && v2 == NULL) {
		return NULL;
	} else if (v1 == NULL) {
		return to_row_matrix(v2);
	} else if (v2 == NULL) {
		return to_row_matrix(v1);
	} else if (v1->size != v2->size) {
		return NULL;
	}
	gsl_matrix* cat_vv = gsl_matrix_alloc(2, v1->size);
	gsl_vector_view cat_v1 = gsl_matrix_row(cat_vv, 0);
	gsl_vector_view cat_v2 = gsl_matrix_row(cat_vv, 1);
	gsl_vector_memcpy(&cat_v1.vector, v1);
	gsl_vector_memcpy(&cat_v2.vector, v2);
	return cat_vv;
}

/*
 * Emulates ([m1 m2]).
 */
gsl_matrix* matlab::concatenate_rows(const gsl_matrix* m1, const gsl_matrix* m2) {
	if (m1 == NULL && m2 == NULL) {
		return NULL;
	} else if (m1 == NULL) {
		return copy(m2);
	} else if (m2 == NULL) {
		return copy(m1);
	} else if (m1->size1 != m2->size1) {
		return NULL;
	}
	gsl_matrix* cat_m = gsl_matrix_alloc(m1->size1, m1->size2 + m2->size2);
	gsl_matrix_view cat_m1 = gsl_matrix_submatrix(cat_m, 0, 0, m1->size1, m1->size2);
	gsl_matrix_view cat_m2 = gsl_matrix_submatrix(cat_m, 0, m1->size2, m2->size1, m2->size2);
	gsl_matrix_memcpy(&cat_m1.matrix, m1);
	gsl_matrix_memcpy(&cat_m2.matrix, m2);
	return cat_m;
}

/*
 * Emulates ([m v]).
 */
gsl_matrix* matlab::concatenate_rows(const gsl_matrix* m, const gsl_vector* v) {
	if (m == NULL && v == NULL) {
		return NULL;
	} else if (m == NULL) {
		return to_column_matrix(v);
	} else if (v == NULL) {
		return copy(m);
	} else if (m->size1 != v->size) {
		return NULL;
	}
	gsl_matrix* cat_mv = gsl_matrix_alloc(m->size1, m->size2 + 1);
	gsl_matrix_view cat_m = gsl_matrix_submatrix(cat_mv, 0, 0, m->size1, m->size2);
	gsl_vector_view cat_v = gsl_matrix_column(cat_mv, cat_mv->size2-1);
	gsl_matrix_memcpy(&cat_m.matrix, m);
	gsl_vector_memcpy(&cat_v.vector, v);
	return cat_mv;
}

/*
 * Emulates ([v m]).
 */
gsl_matrix* matlab::concatenate_rows(const gsl_vector* v, const gsl_matrix* m) {
	if (m == NULL && v == NULL) {
		return NULL;
	} else if (m == NULL) {
		return to_column_matrix(v);
	} else if (v == NULL) {
		return copy(m);
	} else if (m->size2 != v->size) {
		return NULL;
	}
	gsl_matrix* cat_mv = gsl_matrix_alloc(m->size1, m->size2 + 1);
	gsl_vector_view cat_v = gsl_matrix_column(cat_mv, 0);
	gsl_matrix_view cat_m = gsl_matrix_submatrix(cat_mv, 0, 1, m->size1, m->size2);
	gsl_vector_memcpy(&cat_v.vector, v);
	gsl_matrix_memcpy(&cat_m.matrix, m);
	return cat_mv;
}

/*
 * Emulates copy assignment.
 */
gsl_vector* matlab::copy(const gsl_vector* v) {
	gsl_vector* copy_v = gsl_vector_alloc(v->size);
	gsl_vector_memcpy(copy_v, v);
	return copy_v;
}

/*
 * Emulates copy assignment.
 */
gsl_matrix* matlab::copy(const gsl_matrix* m) {
	gsl_matrix* copy_m = gsl_matrix_alloc(m->size1, m->size2);
	gsl_matrix_memcpy(copy_m, m);
	return copy_m;
}

/*
 * Returns the truth value of a matrix: 1 if all the elements are non-zero,
 * 0 if at least 1 element is zero
 */
bool matlab::truth(const gsl_matrix* m) {
	for(int i = 0;i < m->size1;i++) {
		for(int j = 0;j < m->size2;j++) {  
			if(fp_zero(gsl_matrix_get(m, i,j ))) {
				return false;
			}
		}
	}
	return true;
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
 * Emulates (m1 & m2).
 */
gsl_matrix* matlab::logical_and(const gsl_matrix* m1, const gsl_matrix* m2) {
	if (m1->size1 != m2->size1 || m1->size2 != m2->size2) {
		return NULL;
	}
	gsl_matrix* and_m = gsl_matrix_alloc(m1->size1, m1->size2);
	for (int i = 0; i < m1->size1; i++) {
		for (int j = 0; j < m1->size2; j++) {
			bool nz1 = fp_nonzero(gsl_matrix_get(m1, i, j));
			bool nz2 = fp_nonzero(gsl_matrix_get(m2, i, j));
			gsl_matrix_set(and_m, i, j, (double)(nz1 && nz2));
		}
	}
	return and_m;
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
 * Emulates (~m)
 */
gsl_matrix* matlab::logical_not(const gsl_matrix* m) {
	gsl_matrix* not_m = gsl_matrix_alloc(m->size1, m->size2);
	for(int i = 0;i < m->size1;i++) {
		for(int j = 0;j < m->size2;j++) {
			bool z = fp_zero(gsl_matrix_get(m, i, j));
			gsl_matrix_set(not_m, i, j, z);
		}
	}
	return not_m;
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
 * Emulates (m1 | m2).
 */
gsl_matrix* matlab::logical_or(const gsl_matrix* m1, const gsl_matrix* m2) {
	if (m1->size1 != m2->size1 || m1->size2 != m2->size2) {
		return NULL;
	}
	gsl_matrix* or_m = gsl_matrix_alloc(m1->size1, m1->size2);
	for(int i = 0;i < m1->size1;i++) {
		for(int j = 0;j < m1->size2;j++) {
			bool nz1 = fp_nonzero(gsl_matrix_get(m1, i, j));
			bool nz2 = fp_nonzero(gsl_matrix_get(m2, i, j));
			gsl_matrix_set(or_m, i, j, (double)(nz1 || nz2));
		}
	}
	return or_m;
}

/*
 * Emulates (m1 * m2).
 */
gsl_matrix* matlab::mul(const gsl_matrix* m1, const gsl_matrix* m2) {
	gsl_matrix* product_m = gsl_matrix_alloc(m1->size1, m2->size2);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, m1, m2, 0.0, product_m);
	return product_m;
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
	gsl_matrix* pow_m = copy(m);
	for (int i = 2; i <= power; i++) {
		gsl_matrix* temp = mul(pow_m, m);
		gsl_matrix_free(pow_m);
		pow_m = temp;
	}
	return pow_m;
}

/*
 * Emulates (v op x), where op is a binary comparison operator.
 */
gsl_vector* matlab::compare_elements(const gsl_vector* v, compare_fn cmp_fn, double x) {
	gsl_vector* cmp_v = gsl_vector_alloc(v->size);
	for (int i = 0; i < v->size; i++) {
		double value = gsl_vector_get(v, i);
		gsl_vector_set(cmp_v, i, (double)cmp_fn(value, x));
	}
	return cmp_v;
}

/*
 * Emulates (v1 op v2), where op is a binary comparison operator.
 */
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

/*
 * Emulates (m op x), where op is a binary comparison operator.
 */
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

/*
 * Emulates (m1 op m2), where op is a binary comparison operator.
 */
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
 * Emulates vector indexing and assignment (v(indices) = x).
 */
void matlab::index_assign(gsl_vector* v, const gsl_vector* indices, double x) {
	for (int i = 0; i < indices->size; i++) {
		int index = (int)gsl_vector_get(indices, i);
		if (index < v->size) {
			gsl_vector_set(v, index, x);
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
 * Emulates vector logical indexing and assignment (v(lv) = x).
 */
void matlab::logical_index_assign(gsl_vector* v, const gsl_vector* lv, double x) {
	for (int i = 0; i < lv->size; i++) {
		if (fp_nonzero(gsl_vector_get(lv, i)) && i < v->size) {
			gsl_vector_set(v, i, x);
		}
	}
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

/*
 * Emulates matrix-to-vector conversion.  The vector is constructed by
 * consecutively appending columns.
 */
gsl_vector* matlab::to_vector(const gsl_matrix* m) {
	gsl_vector* v = gsl_vector_alloc(m->size1 * m->size2);
	for (int j = 0; j < m->size2; j++) {
		for (int i = 0; i < m->size1; i++) {
			double value = gsl_matrix_get(m, i, j);
			gsl_vector_set(v, j * m->size1 + i, value);
		}
	}
	return v;
}

/*
 * Converts a vector to a column matrix
 */
gsl_matrix* matlab::to_column_matrix(const gsl_vector* v) {
	gsl_matrix* m = gsl_matrix_alloc(v->size, 1);
	for(int i=0;i < v->size;i++) {
		gsl_matrix_set(m, i, 0, gsl_vector_get(v, i));
	}
	return m;
}

/*
 * Converts a vector to a row matrix
 */
gsl_matrix* matlab::to_row_matrix(const gsl_vector* v) {
	gsl_matrix* m = gsl_matrix_alloc(1, v->size);
	for(int i=0;i < v->size;i++) {
		gsl_matrix_set(m, 0, i, gsl_vector_get(v, i));
	}
	return m;
}
