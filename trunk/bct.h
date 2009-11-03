#ifndef BCT_H
#define BCT_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

namespace bct {
	class bct_exception { };
	class gsl_exception : public bct_exception { };
	class out_of_memory_exception : public bct_exception { };
	class size_exception : public bct_exception { };

	const double EPSILON = 1e-6;  // Used for floating-point equality comparisons
	extern bool safe_mode;  // Determines whether to check matrix status at the beginning of functions

	// Density, degree, and assortativity
	double assortativity_dir(const gsl_matrix*);
	double assortativity_und(const gsl_matrix*);
	double assortativity(const gsl_vector*, const gsl_vector*, int);
	gsl_vector* degrees_dir(const gsl_matrix*, gsl_vector* = NULL, gsl_vector* = NULL);
	gsl_vector* degrees_und(const gsl_matrix*);
	double density_dir(const gsl_matrix*);
	double density_und(const gsl_matrix*);
	gsl_matrix* jdegree(const gsl_matrix*); //returns the Joint distribution matrix
	double jdegree_id(gsl_matrix*); //returns other metrices from the Joint dstrbn matrix
	double jdegree_od(gsl_matrix*);
	double jdegree_bl(gsl_matrix*);
	gsl_matrix* matching_ind(const gsl_matrix*); //returns all indices
	gsl_matrix* matching_ind_in(const gsl_matrix*);
	gsl_matrix* matching_ind_out(const gsl_matrix*);
	gsl_vector* strengths_dir(const gsl_matrix*, gsl_vector* = NULL, gsl_vector* = NULL);
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
	gsl_vector* logical_and(const gsl_vector*, const gsl_vector*);
	gsl_vector* logical_not(const gsl_vector*);
	int nnz(const gsl_vector*);
	int nnz(const gsl_matrix*);
	gsl_matrix* submatrix(const gsl_matrix*, const gsl_vector*, const gsl_vector*);
	double sum(const gsl_vector*);
	gsl_vector* sum(const gsl_matrix*, int = 1);
	gsl_matrix* tril(const gsl_matrix*, int = 0);
	gsl_matrix* triu(const gsl_matrix*, int = 0);
	
	// Matrix status
	enum status {
		UNDIRECTED = 1, DIRECTED = 2,
		BINARY = 4, WEIGHTED = 8,
		POSITIVE = 16, SIGNED = 32,
		NO_LOOPS = 64, LOOPS = 128
	};
	bool check_status(const gsl_matrix*, int, const char* = NULL);
	bool has_loops(const gsl_matrix*);
	bool has_no_loops(const gsl_matrix*);
	bool is_binary(const gsl_matrix*);
	bool is_directed(const gsl_matrix*);
	bool is_positive(const gsl_matrix*);
	bool is_signed(const gsl_matrix*);
	bool is_undirected(const gsl_matrix*);
	bool is_weighted(const gsl_matrix*);
	
	// Matrix conversion
	gsl_matrix* binary(const gsl_matrix*);
	gsl_matrix* positive(const gsl_matrix*);
	gsl_matrix* remove_loops(const gsl_matrix*);
	gsl_matrix* undirected(const gsl_matrix*, bool = true);
	
	// Utility
	gsl_matrix* find(const gsl_matrix*, int, double);	
	gsl_vector* pick_cells(const gsl_vector*, const gsl_vector*);
	void set_safe_mode(bool);
	gsl_vector* splice(const gsl_vector*, const gsl_vector*);
	
	// Error handling
	void gsl_error_handler(const char*, const char*, int, int);

	// Debugging
	void printf(const gsl_vector*, const char*);
	void printf(const gsl_matrix*, const char*);
};

#endif
