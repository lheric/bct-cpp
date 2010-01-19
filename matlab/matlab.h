#ifndef MATLAB_H
#define MATLAB_H

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "quicksort.h"

namespace matlab {
	typedef bool (*fp_cmp_fn)(double, double);
	const double EPSILON = 1e-6;

	// Functions
	int all(const gsl_vector*);
	gsl_vector* all(const gsl_matrix*, int = 1);
	int any(const gsl_vector*);
	gsl_vector* any(const gsl_matrix*, int = 1);
	void dec2bin(int, char*);
	void dec2bin(int, int, char*);
	gsl_matrix* eye(int);
	gsl_matrix* eye(int, int);
	gsl_vector* find(const gsl_vector*, int = GSL_POSINF, const char* = "first");
	gsl_vector* find(const gsl_matrix*, int = GSL_POSINF, const char* = "first");
	gsl_matrix* find_ij(const gsl_matrix*, int = GSL_POSINF, const char* = "first");
	double max(const gsl_vector*);
	gsl_vector* max(const gsl_matrix*, int = 1);
	double min(const gsl_vector*);
	gsl_vector* min(const gsl_matrix*, int = 1);
	int nnz(const gsl_vector*);
	int nnz(const gsl_matrix*);
	gsl_vector* nonzeros(const gsl_matrix*);
	gsl_matrix* ones(int);
	gsl_matrix* ones(int, int);
	gsl_vector* reverse(gsl_vector*);
	gsl_matrix* sortrows(const gsl_matrix*, gsl_vector* = NULL);
	double sum(const gsl_vector*);
	gsl_vector* sum(const gsl_matrix*, int = 1);
	gsl_matrix* tril(const gsl_matrix*, int = 0);
	gsl_matrix* triu(const gsl_matrix*, int = 0);
	gsl_matrix* unique_rows(const gsl_matrix*, gsl_vector* = NULL, gsl_vector* = NULL);
	gsl_matrix* unique_rows(const gsl_matrix*, const char*, gsl_vector* = NULL, gsl_vector* = NULL);
	gsl_matrix* zeros(int);
	gsl_matrix* zeros(int, int);
	
	// Operators
	gsl_vector* concatenate(const gsl_vector*, double);
	gsl_vector* concatenate(double, const gsl_vector*);
	gsl_vector* concatenate(const gsl_vector*, const gsl_vector*);
	gsl_matrix* concatenate_columns(const gsl_vector*, const gsl_vector*);
	gsl_matrix* concatenate_columns(const gsl_matrix*, const gsl_vector*);
	gsl_matrix* concatenate_columns(const gsl_vector*, const gsl_matrix*);	
	gsl_matrix* concatenate_columns(const gsl_matrix*, const gsl_matrix*);
	gsl_matrix* concatenate_rows(const gsl_vector*, const gsl_vector*);
	gsl_matrix* concatenate_rows(const gsl_matrix*, const gsl_vector*);
	gsl_matrix* concatenate_rows(const gsl_vector*, const gsl_matrix*);
	gsl_matrix* concatenate_rows(const gsl_matrix*, const gsl_matrix*);
	gsl_vector* copy(const gsl_vector*);
	gsl_matrix* copy(const gsl_matrix*);
	gsl_vector* logical_and(const gsl_vector*, const gsl_vector*);
	gsl_matrix* logical_and(const gsl_matrix*, const gsl_matrix*);
	gsl_vector* logical_not(const gsl_vector*);
	gsl_matrix* logical_not(const gsl_matrix*);
	gsl_vector* logical_or(const gsl_vector*, const gsl_vector*);
	gsl_matrix* logical_or(const gsl_matrix*, const gsl_matrix*);
	gsl_matrix* mul(const gsl_matrix*, const gsl_matrix*);
	gsl_matrix* pow(const gsl_matrix*, int);
	gsl_vector* pow_elements(const gsl_vector*, double);
	gsl_vector* pow_elements(const gsl_vector*, const gsl_vector*);
	gsl_matrix* pow_elements(const gsl_matrix*, double);
	gsl_matrix* pow_elements(const gsl_matrix*, const gsl_matrix*);
	gsl_vector* sequence(int, int);
	gsl_vector* sequence(int, int, int);
	
	// Floating-point comparison
	int fp_compare(double, double);
	bool fp_zero(double);
	bool fp_nonzero(double);
	bool fp_positive(double);
	bool fp_negative(double);
	bool fp_equal(double, double);
	bool fp_not_equal(double, double);
	bool fp_less(double, double);
	bool fp_less_or_equal(double, double);
	bool fp_greater(double, double);
	bool fp_greater_or_equal(double, double);
	
	// Vector/matrix comparison
	int compare_vectors(const gsl_vector*, const gsl_vector*);
	int compare_vectorps(const gsl_vector**, const gsl_vector**);
	int compare_matrices(const gsl_matrix*, const gsl_matrix*);
	int compare_matrixps(const gsl_matrix**, const gsl_matrix**);
	gsl_vector* compare_elements(const gsl_vector*, fp_cmp_fn, double);
	gsl_vector* compare_elements(const gsl_vector*, fp_cmp_fn, const gsl_vector*);
	gsl_matrix* compare_elements(const gsl_matrix*, fp_cmp_fn, double);
	gsl_matrix* compare_elements(const gsl_matrix*, fp_cmp_fn, const gsl_matrix*);
	
	// Vector-by-vector indexing
	gsl_vector* ordinal_index(const gsl_vector*, const gsl_vector*);
	void ordinal_index_assign(gsl_vector*, const gsl_vector*, double);
	void ordinal_index_assign(gsl_vector*, const gsl_vector*, const gsl_vector*);
	gsl_vector* logical_index(const gsl_vector*, const gsl_vector*);
	void logical_index_assign(gsl_vector*, const gsl_vector*, double);
	void logical_index_assign(gsl_vector*, const gsl_vector*, const gsl_vector*);
	
	// Matrix-by-integer indexing
	double ordinal_index(const gsl_matrix*, int);
	void ordinal_index_assign(gsl_matrix*, int, double);
	
	// Matrix-by-vector indexing
	gsl_vector* ordinal_index(const gsl_matrix*, const gsl_vector*);
	void ordinal_index_assign(gsl_matrix*, const gsl_vector*, double);
	void ordinal_index_assign(gsl_matrix*, const gsl_vector*, const gsl_vector*);
	gsl_vector* logical_index(const gsl_matrix*, const gsl_vector*);
	void logical_index_assign(gsl_matrix*, const gsl_vector*, double);
	void logical_index_assign(gsl_matrix*, const gsl_vector*, const gsl_vector*);
	
	// Matrix-by-two-vectors indexing (non-mixed)
	gsl_matrix* ordinal_index(const gsl_matrix*, const gsl_vector*, const gsl_vector*);
	void ordinal_index_assign(gsl_matrix*, const gsl_vector*, const gsl_vector*, double);
	void ordinal_index_assign(gsl_matrix*, const gsl_vector*, const gsl_vector*, const gsl_matrix*);
	gsl_matrix* logical_index(const gsl_matrix*, const gsl_vector*, const gsl_vector*);
	void logical_index_assign(gsl_matrix*, const gsl_vector*, const gsl_vector*, double);
	void logical_index_assign(gsl_matrix*, const gsl_vector*, const gsl_vector*, const gsl_matrix*);
	
	// Matrix-by-two-vectors indexing (mixed)
	gsl_matrix* ord_log_index(const gsl_matrix*, const gsl_vector*, const gsl_vector*);
	void ord_log_index_assign(gsl_matrix*, const gsl_vector*, const gsl_vector*, double);
	void ord_log_index_assign(gsl_matrix*, const gsl_vector*, const gsl_vector*, const gsl_matrix*);
	gsl_matrix* log_ord_index(const gsl_matrix*, const gsl_vector*, const gsl_vector*);
	void log_ord_index_assign(gsl_matrix*, const gsl_vector*, const gsl_vector*, double);
	void log_ord_index_assign(gsl_matrix*, const gsl_vector*, const gsl_vector*, const gsl_matrix*);
	
	// Matrix-by-matrix indexing
	gsl_matrix* ordinal_index(const gsl_matrix*, const gsl_matrix*);
	void ordinal_index_assign(gsl_matrix*, const gsl_matrix*, double);
	void ordinal_index_assign(gsl_matrix*, const gsl_matrix*, const gsl_matrix*);
	gsl_vector* logical_index(const gsl_matrix*, const gsl_matrix*);
	void logical_index_assign(gsl_matrix*, const gsl_matrix*, double);
	void logical_index_assign(gsl_matrix*, const gsl_matrix*, const gsl_vector*);
	
	// Vector/matrix conversion
	bool to_bool(const gsl_vector*);
	bool to_bool(const gsl_matrix*);
	gsl_vector* to_vector(const gsl_matrix*);
	gsl_matrix* to_column_matrix(const gsl_vector*);
	gsl_matrix* to_row_matrix(const gsl_vector*);
};

#endif
