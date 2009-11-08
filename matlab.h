#ifndef MATLAB_H
#define MATLAB_H

namespace matlab {	
	const double EPSILON = 1e-6;

	gsl_vector* find(const gsl_vector*);
	gsl_vector* logical_and(const gsl_vector*, const gsl_vector*);
	gsl_vector* logical_not(const gsl_vector*);
	int nnz(const gsl_vector*);
	int nnz(const gsl_matrix*);
	gsl_matrix* submatrix(const gsl_matrix*, const gsl_vector*, const gsl_vector*);
	double sum(const gsl_vector*);
	gsl_vector* sum(const gsl_matrix*, int = 1);
	gsl_matrix* tril(const gsl_matrix*, int = 0);
	gsl_matrix* triu(const gsl_matrix*, int = 0);
};

#endif
