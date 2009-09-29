#ifndef BCT_H
#define BCT_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

namespace bct {
	const double EPSILON = 1e-6;  // Used for floating-point equality comparisons

	// Density, degree, and assortativity
	gsl_vector* degrees_und(const gsl_matrix*);
	gsl_vector* degrees_dir_in(const gsl_matrix*);
	gsl_vector* degrees_dir_out(const gsl_matrix*);
	gsl_vector* degrees_dir(const gsl_matrix*);

	// Clustering

	// Paths, distances, and cycles

	// Centrality

	// Motifs

	// Modularity and community structure

	// Utility
	gsl_vector* binary(const gsl_vector*);
	gsl_matrix* binary(const gsl_matrix*);

	// Debugging
	void printf(const gsl_vector*, const char*);
	void printf(const gsl_matrix*, const char*);
};

#endif
