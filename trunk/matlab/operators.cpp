#include <cmath>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "matlab.h"

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
 * Emulates ([v1 v2]) for column vectors.
 */
gsl_matrix* matlab::concatenate_rows(const gsl_vector* v1, const gsl_vector* v2) {
	if (v1 == NULL && v2 == NULL) {
		return NULL;
	} else if (v2 == NULL) {
		return to_column_matrix(v1);
	} else if (v1 == NULL) {
		return to_column_matrix(v2);
	} else if (v1->size != v2->size) {
		return NULL;
	}
	gsl_matrix* cat_vv = gsl_matrix_alloc(v1->size, 2);
	gsl_matrix_set_col(cat_vv, 0, v1);
	gsl_matrix_set_col(cat_vv, 1, v2);
	return cat_vv;
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
