#ifndef MATLAB_H
#define MATLAB_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <limits>

namespace matlab {
	typedef bool (*compare_fn)(double, double);
	const double EPSILON = 1e-6;

	// Functions
	gsl_matrix* eye(int);
	gsl_matrix* eye(int, int);
	gsl_vector* find(const gsl_vector*, int = std::numeric_limits<int>::max(), const char* = NULL);
	gsl_vector* find(const gsl_matrix*, int = std::numeric_limits<int>::max(), const char* = NULL);
	gsl_matrix* find_ij(const gsl_matrix*, int = std::numeric_limits<int>::max(), const char* = NULL);
	int nnz(const gsl_vector*);
	int nnz(const gsl_matrix*);
	gsl_matrix* ones(int);
	gsl_matrix* ones(int, int);
	gsl_matrix* yens(int, double);
	gsl_matrix* yens(int, int, double);
	gsl_vector* sequence(int, int);
	double sum(const gsl_vector*);
	gsl_vector* sum(const gsl_matrix*, int = 1);
	gsl_matrix* tril(const gsl_matrix*, int = 0);
	gsl_matrix* triu(const gsl_matrix*, int = 0);
	gsl_matrix* zeros(int);
	gsl_matrix* zeros(int, int);
	
	// Operators
	gsl_vector* concatenate(const gsl_vector*, const double);
	gsl_vector* concatenate(const gsl_vector*, const gsl_vector*);
	gsl_matrix* concatenate_columns(const gsl_matrix*, const gsl_matrix*);
	gsl_matrix* concatenate_columns(const gsl_matrix*, const gsl_vector*);
	gsl_matrix* concatenate_columns(const gsl_vector*, const gsl_matrix*);	
	gsl_matrix* concatenate_rows(const gsl_matrix*, const gsl_matrix*);
	gsl_matrix* concatenate_rows(const gsl_matrix*, const gsl_vector*);
	gsl_matrix* concatenate_rows(const gsl_vector*, const gsl_matrix*);	
	gsl_vector* copy(const gsl_vector*);
	gsl_matrix* copy(const gsl_matrix*);
	gsl_vector* logical_and(const gsl_vector*, const gsl_vector*);
	gsl_vector* logical_not(const gsl_vector*);
	gsl_vector* logical_or(const gsl_vector*, const gsl_vector*);
	gsl_matrix* logical_or(const gsl_matrix*, const gsl_matrix*);
	gsl_matrix* mul(const gsl_matrix*, const gsl_matrix*);
	gsl_vector* pow_elements(const gsl_vector*, double);
	gsl_vector* pow_elements(const gsl_vector*, const gsl_vector*);
	gsl_matrix* pow_elements(const gsl_matrix*, double);
	gsl_matrix* pow_elements(const gsl_matrix*, const gsl_matrix*);
	gsl_matrix* pow(const gsl_matrix*, int);
	
	// Vector/matrix comparisons
	gsl_vector* compare_elements(const gsl_vector*, compare_fn, double);
	gsl_vector* compare_elements(const gsl_vector*, compare_fn, const gsl_vector*);
	gsl_matrix* compare_elements(const gsl_matrix*, compare_fn, double);
	gsl_matrix* compare_elements(const gsl_matrix*, compare_fn, const gsl_matrix*);
	bool cmp_equal(double, double);
	bool cmp_not_equal(double, double);
	bool cmp_greater(double, double);
	bool cmp_greater_or_equal(double, double);
	bool cmp_less(double, double);
	bool cmp_less_or_equal(double, double);
	
	// Floating-point comparisons
	int fp_compare(double, double);
	bool fp_equal(double, double);
	bool fp_not_equal(double, double);
	bool fp_zero(double);
	bool fp_nonzero(double);
	bool fp_positive(double);
	bool fp_negative(double);

	// Indexing
	gsl_vector* index(const gsl_vector*, const gsl_vector*);
	double index(const gsl_matrix*, int);
	gsl_matrix* index(const gsl_matrix*, const gsl_vector*, const gsl_vector*);
	gsl_matrix* index(const gsl_matrix*, const gsl_matrix*);
	void index_assign(gsl_vector*, const gsl_vector*, double);
	void index_assign(gsl_matrix*, int, double);
	void index_assign(gsl_matrix*, const gsl_vector*, double);
	void index_assign(gsl_matrix*, const gsl_matrix*, double);
	void index_assign(gsl_matrix*, const gsl_vector*, const gsl_vector*, double);
	gsl_vector* logical_index(const gsl_vector*, const gsl_vector*);
	gsl_vector* logical_index(const gsl_matrix*, const gsl_vector*);
	gsl_vector* logical_index(const gsl_matrix*, const gsl_matrix*);
	void logical_index_assign(gsl_vector*, const gsl_vector*, double);
	void logical_index_assign(gsl_matrix*, const gsl_vector*, double);
	void logical_index_assign(gsl_matrix*, const gsl_matrix*, double);
	
	// Conversions
	gsl_vector* to_vector(const gsl_matrix*);
	gsl_matrix* to_row_matrix(const gsl_vector*);
	gsl_matrix* to_column_matrix(const gsl_vector*);
};

#endif
