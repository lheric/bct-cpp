#ifndef BCT_H
#define BCT_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

namespace bct {
	struct jdegree_outputs {
		gsl_matrix* J;
		int J_id, J_od, J_bl;
	};

	struct matching_ind_outputs {
		gsl_matrix* Min;
		gsl_matrix* Mout;
		gsl_matrix* Mall;
	};
	
	const double EPSILON = 1e-6;  // Used for floating-point equality comparisons

	// Density, degree, and assortativity
	double assortativity(const gsl_matrix*, int);
	gsl_vector* degrees_dir(const gsl_matrix*, gsl_vector* = NULL, gsl_vector* = NULL);
	gsl_vector* degrees_und(const gsl_matrix*);
	double density_dir(const gsl_matrix*);
	double density_und(const gsl_matrix*);
	struct jdegree_outputs jdegree(const gsl_matrix*);
	struct matching_ind_outputs matching_ind(const gsl_matrix*);
	gsl_vector* strengths_dir(const gsl_matrix*);
	gsl_vector* strengths_und(const gsl_matrix*);

	// Clustering
	gsl_vector* clustering_coef_bu(const gsl_matrix*);

	// Paths, distances, and cycles

	// Centrality

	// Motifs

	// Modularity and community structure

	// Data sets
	extern double cat_all[];
	extern double cat_ctx[];
	extern double fve30[];
	extern double fve32[];
	extern double macaque47[];
	extern double macaque71[];
	
	// MATLAB emulation
	gsl_vector* find(const gsl_vector*);
	int nnz(const gsl_vector*);
	int nnz(const gsl_matrix*);
	gsl_matrix* submatrix(const gsl_matrix*, const gsl_vector*, const gsl_vector*);
	double sum(const gsl_vector*);
	gsl_vector* sum(const gsl_matrix*, int);
	gsl_matrix* tril(const gsl_matrix*, int);
	gsl_matrix* triu(const gsl_matrix*, int);
	gsl_vector* vectorAnd(const gsl_vector*, const gsl_vector*);
	gsl_vector* vectorNot(const gsl_vector*);
	
	// Conversion
	gsl_vector* binary(const gsl_vector*);
	gsl_matrix* binary(const gsl_matrix*);
	gsl_matrix* zero_diagonal(const gsl_matrix*);
	
	// Utility
	gsl_matrix* find(const gsl_matrix*, int, double);	
	gsl_vector* pickCells(const gsl_vector*, const gsl_vector*);
	gsl_vector* splice(const gsl_vector*, const gsl_vector*);

	// Debugging
	void printf(const gsl_vector*, const char*);
	void printf(const gsl_matrix*, const char*);
};

#endif
