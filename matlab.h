#ifndef MATLAB_H
#define MATLAB_H

namespace matlab {	
	const double EPSILON = 1e-6;

	// Functions
	gsl_vector* find(const gsl_vector*, int = -1, const char* = NULL);
	gsl_vector* find(const gsl_matrix*, int = -1, const char* = NULL);
	int nnz(const gsl_vector*);
	int nnz(const gsl_matrix*);
	double sum(const gsl_vector*);
	gsl_vector* sum(const gsl_matrix*, int = 1);
	gsl_matrix* tril(const gsl_matrix*, int = 0);
	gsl_matrix* triu(const gsl_matrix*, int = 0);
	
	// Operators
	gsl_vector* logical_and(const gsl_vector*, const gsl_vector*);
	gsl_vector* logical_not(const gsl_vector*);
	gsl_vector* logical_or(const gsl_vector*, const gsl_vector*);
	
	// Comparisons
	int compare(double, double);
	bool is_equal(double, double);
	bool is_negative(double);
	bool is_nonzero(double);
	bool is_not_equal(double, double);
	bool is_positive(double);
	bool is_zero(double);

	// Indexing
	gsl_matrix* index(const gsl_matrix*, const gsl_vector*, const gsl_vector*);
	
	// Conversions
	gsl_vector* to_vector(const gsl_matrix*);
};

#endif
